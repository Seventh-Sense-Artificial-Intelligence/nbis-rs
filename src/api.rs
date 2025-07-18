use std::ffi::CStr;
use std::ptr::{null_mut, NonNull};
use std::sync::Arc;

use image::DynamicImage;
use std::os::raw::{c_int, c_uchar, c_void};

use image::Rgb;
use imageproc::drawing::{draw_filled_circle_mut, draw_filled_rect_mut};
use imageproc::rect::Rect;

use crate::consts::MM_PER_INCH;
use crate::encoding::decode_minutia;
use crate::errors::NbisError;
use crate::ffi::{
    comp_nfiq_featvctr, dflt_acfunc_hids, dflt_acfunc_outs, dflt_nHids, dflt_nInps, dflt_nOuts,
    dflt_wts, dflt_znorm_means, dflt_znorm_stds, free, free_minutiae, get_minutiae, runmlp2,
    sivv_ffi_from_bytes, znorm_fniq_featvctr, LFSPARMS, MINUTIAE, MIN_MINUTIAE, NFIQ_NUM_CLASSES,
    NFIQ_VCTRLEN,
};
use crate::imutils::{draw_arrow_with_head, png_bytes_from_rgb};
use crate::{Minutia, MinutiaKind, Minutiae};

/// Represents the result of the SIVV computation.
#[derive(Debug, Clone, uniffi::Record)]
pub struct SIVVResult {
    /// Index of the largest peak-valley pair (1-based)
    pub largest_pvp_index: i32,

    /// Total number of detected peak-valley pairs
    pub total_pvps: i32,

    /// Power difference between the peak and valley
    pub power_diff: f64,

    /// Frequency difference between the peak and valley
    pub freq_diff: f64,

    /// Slope between valley and peak (dy / dx)
    pub slope: f64,

    /// Frequency of the midpoint between valley and peak
    pub center_frequency: f64,

    /// Absolute frequency of the peak (undocumented in comments)
    pub peak_frequency: f64,
}

/// Represents the result of the NFIQ computation.
#[derive(Debug, Clone, uniffi::Record)]
pub struct NfiqResult {
    /// The NFIQ quality score.
    /// 1 = Excellent, 2 = Very Good, 3 = Good,
    /// 4 = Fair, 5 = Poor.
    /// See [`NfiqQuality`] for more details.
    pub nfiq: NfiqQuality,
    /// The confidence level of the NFIQ score.
    /// A value between 0.0 and 1.0, where 1.0 means very confident.
    /// This is a floating point value.
    pub confidence: f32,
}

/// Represents the quality of a fingerprint image as determined by NFIQ.
#[derive(Debug, Clone, PartialEq, PartialOrd, uniffi::Enum)]
pub enum NfiqQuality {
    /// Unknown quality fingerprint image. If the minutiae are loaded from a template,
    /// this means the quality is unknown.
    Unknown = 0,
    /// Excellent quality fingerprint image.
    /// This means the image is very clear and suitable for fingerprint recognition.
    Excellent = 1,
    /// Very good quality fingerprint image.
    /// This means the image is clear but may have some minor issues.
    VeryGood = 2,
    /// Good quality fingerprint image.
    /// This means the image is usable but has noticeable issues.
    Good = 3,
    /// Fair quality fingerprint image.
    /// This means the image is barely usable for fingerprint recognition.
    Fair = 4,
    /// Poor quality fingerprint image.
    /// This means the image is not usable for fingerprint recognition.
    Poor = 5,
}

impl NfiqQuality {
    pub fn from_i32(value: i32) -> Option<Self> {
        match value {
            1 => Some(NfiqQuality::Excellent),
            2 => Some(NfiqQuality::VeryGood),
            3 => Some(NfiqQuality::Good),
            4 => Some(NfiqQuality::Fair),
            5 => Some(NfiqQuality::Poor),
            _ => None,
        }
    }
}

#[uniffi::export]
pub fn sivv(image: &[u8]) -> Result<SIVVResult, NbisError> {
    // 0) Load the image ------------------------------------------------------
    let image = match image::load_from_memory(image) {
        Ok(img) => img,
        Err(_e) => return Err(NbisError::ImageLoadError),
    };

    // 1) Ensure 8‑bit grayscale ------------------------------------------------
    let gray = match image {
        DynamicImage::ImageLuma8(buf) => buf.clone(),
        _ => image.to_luma8(),
    };

    unsafe {
        let ptr = sivv_ffi_from_bytes(
            gray.as_ptr() as *mut c_uchar,
            gray.width() as i32,
            gray.height() as i32,
        );
        let str_result = CStr::from_ptr(ptr).to_string_lossy().into_owned();
        // Split the result into parts
        let parts: Vec<&str> = str_result.split(',').map(|s| s.trim()).collect();
        if parts.len() != 7 {
            return Err(NbisError::GenericError(
                "Invalid SIVV result format".to_string(),
            ));
        }

        // Parse the parts into the SIVVResult struct
        let result = SIVVResult {
            largest_pvp_index: parts[0].parse().unwrap_or_default(),
            total_pvps: parts[1].parse().unwrap_or_default(),
            power_diff: parts[2].parse().unwrap_or_default(),
            freq_diff: parts[3].parse().unwrap_or_default(),
            slope: parts[4].parse().unwrap_or_default(),
            center_frequency: parts[5].parse().unwrap_or_default(),
            peak_frequency: parts[6].parse().unwrap_or_default(),
        };

        Ok(result)
    }
}

/// Extracts minutiae from an 8‑bit grayscale fingerprint image using the
/// **NIST LFS v2** algorithm.
///
/// * `image` — the bytes of any image (PNG, JPEG, etc.) that can be converted
/// * `ppi`  — (Optional) scanner resolution in dpi. Default is 500 dpi.
///
/// Returns an [`Minutiae`] object containing the extracted minutiae.
///
/// If the image cannot be processed, returns an [`NbisError`].
///
/// This function is the main entry point for fingerprint minutiae extraction.
#[uniffi::export]
pub fn extract_minutiae(image: &[u8], ppi: Option<f64>) -> Result<Minutiae, NbisError> {
    let ppi = ppi.unwrap_or(500.0); // default to 500 dpi

    // 0) Load the image ------------------------------------------------------
    let image = match image::load_from_memory(image) {
        Ok(img) => img,
        Err(_e) => return Err(NbisError::ImageLoadError),
    };

    // 1) Ensure 8‑bit grayscale ------------------------------------------------
    let gray = match image {
        DynamicImage::ImageLuma8(buf) => buf.clone(),
        _ => image.to_luma8(),
    };
    let (iw, ih) = gray.dimensions();

    // Buffers and sizes returned by the C API -------------------------------
    let mut ominutiae: *mut MINUTIAE = null_mut();
    let mut oquality_map: *mut c_int = null_mut();
    let mut odirection_map: *mut c_int = null_mut();
    let mut olow_contrast_map: *mut c_int = null_mut();
    let mut olow_flow_map: *mut c_int = null_mut();
    let mut ohigh_curve_map: *mut c_int = null_mut();
    let mut map_w: c_int = 0;
    let mut map_h: c_int = 0;
    let mut obdata: *mut c_uchar = null_mut();
    let mut obw: c_int = 0;
    let mut obh: c_int = 0;
    let mut obd: c_int = 0;

    let ppmm = ppi / MM_PER_INCH; // convert to px/mm

    // 2) Call into C ---------------------------------------------------------
    let rc = unsafe {
        extern "C" {
            static lfsparms_V2: LFSPARMS;
        }
        get_minutiae(
            &mut ominutiae,
            &mut oquality_map,
            &mut odirection_map,
            &mut olow_contrast_map,
            &mut olow_flow_map,
            &mut ohigh_curve_map,
            &mut map_w,
            &mut map_h,
            &mut obdata,
            &mut obw,
            &mut obh,
            &mut obd,
            gray.as_ptr() as *mut c_uchar,
            iw as c_int,
            ih as c_int,
            8, // id = 8‑bit image
            ppmm,
            &lfsparms_V2 as *const _,
        )
    };

    if rc != 0 {
        return Err(NbisError::UnexpectedError(rc as i64));
    };

    // 3) Copy out the data we care about ------------------------------------
    // let quality = unsafe { std::slice::from_raw_parts(oquality_map, (map_w * map_h) as usize) };
    // let quality_vec = quality.to_vec();

    // let bin_slice = unsafe { std::slice::from_raw_parts(obdata, (obw * obh) as usize) };
    // let bin_img = ImageBuffer::<Luma<u8>, Vec<u8>>::from_raw(obw as u32, obh as u32, bin_slice.to_vec())
    //     .expect("size mismatch creating ImageBuffer");

    // Default to poor quality if we can't compute NFIQ
    let mut quality = NfiqResult {
        nfiq: NfiqQuality::Poor,
        confidence: 1.0,
    };
    // 3) compute the quality of the image -----------------------------
    // Only do quality assessment if there are enough minutiae
    if unsafe { (*ominutiae).num } > MIN_MINUTIAE as i32 {
        let mut featvctr = [0.0f32; NFIQ_VCTRLEN];
        let mut optflag = 0;

        let ret = unsafe {
            comp_nfiq_featvctr(
                featvctr.as_mut_ptr(),
                NFIQ_VCTRLEN as c_int,
                ominutiae,
                oquality_map,
                map_w,
                map_h,
                &mut optflag,
            )
        };

        if ret == 0 {
            // Z-normalize the feature vector
            unsafe {
                znorm_fniq_featvctr(
                    featvctr.as_mut_ptr(),
                    dflt_znorm_means.as_ptr(),
                    dflt_znorm_stds.as_ptr(),
                    NFIQ_VCTRLEN as c_int,
                )
            };

            // Call the MLP for NFIQ classification
            // Define the output arrays
            let mut outacs = [0.0f32; NFIQ_NUM_CLASSES];
            let mut class_idx: c_int = 0;
            let mut confidence: f32 = 0.0;
            let ret = unsafe {
                runmlp2(
                    dflt_nInps,
                    dflt_nHids,
                    dflt_nOuts,
                    dflt_acfunc_hids,
                    dflt_acfunc_outs,
                    dflt_wts.as_ptr() as *mut f32,
                    featvctr.as_mut_ptr(),
                    outacs.as_mut_ptr(),
                    &mut class_idx,
                    &mut confidence,
                )
            };

            if ret == 0 {
                // Map the class index to NFIQ quality
                quality = NfiqResult {
                    nfiq: NfiqQuality::from_i32(class_idx + 1).unwrap_or(NfiqQuality::Poor),
                    confidence,
                };
            }
        }
    }

    let minutiae = NonNull::new(ominutiae).expect("C returned null pointer");
    let minutiae_obj = unsafe {
        let mset = &*minutiae.as_ptr(); // &MINUTIAE
        let raw = std::slice::from_raw_parts(
            mset.list, // [*mut MINUTIA]
            mset.num as usize,
        );

        let minutiae_vec: Vec<Minutia> = raw
            .iter()
            .map(|ptr| {
                let m = &**ptr; // &MINUTIA
                Minutia {
                    x: m.x,
                    y: m.y,
                    direction: m.direction,
                    reliability: m.reliability,
                    // 0 = bifurcation, 1 = ridge ending
                    kind: if m.r#type == 0 {
                        MinutiaKind::Bifurcation
                    } else if m.r#type == 1 {
                        MinutiaKind::RidgeEnding
                    } else {
                        panic!("Unknown minutia type: {}", m.r#type);
                    },
                }
            })
            .collect();
        Minutiae::new(minutiae_vec, iw, ih, quality)
    };

    // 4) Free C allocations we no longer need -------------------------------
    unsafe {
        free(oquality_map as *mut c_void);
        free(odirection_map as *mut c_void);
        free(olow_contrast_map as *mut c_void);
        free(olow_flow_map as *mut c_void);
        free(ohigh_curve_map as *mut c_void);
        free(obdata as *mut c_void);
        free_minutiae(ominutiae);
    };

    Ok(minutiae_obj)
}

/// Loads an ISO/IEC 19794-2:2005 fingerprint template from bytes.
///
/// # Arguments
/// * `template_bytes` — the raw bytes of the ISO template.
///
/// Returns a [`Minutiae`] object containing the decoded minutiae.
///
/// If the template is invalid or cannot be parsed, returns an [`NbisError`].
#[uniffi::export]
pub fn load_iso_19794_2_2005(template_bytes: &[u8]) -> Result<Minutiae, NbisError> {
    if template_bytes.len() < 28 {
        return Err(NbisError::InvalidTemplate(
            "ISO template too short".to_string(),
        ));
    }

    // Check the header
    if &template_bytes[0..8] != b"FMR\0 20\0" {
        return Err(NbisError::InvalidTemplate("Invalid ISO header".to_string()));
    }

    let total_length = u32::from_be_bytes([
        template_bytes[8],
        template_bytes[9],
        template_bytes[10],
        template_bytes[11],
    ]) as usize;
    if total_length != template_bytes.len() {
        return Err(NbisError::InvalidTemplate(
            "Total length mismatch".to_string(),
        ));
    }

    let width = u16::from_be_bytes([template_bytes[14], template_bytes[15]]);
    let height = u16::from_be_bytes([template_bytes[16], template_bytes[17]]);

    let num_minutiae = i8::from_be_bytes([template_bytes[27]]) as usize;

    let mut minutiae = Vec::with_capacity(num_minutiae);
    for i in 0..num_minutiae {
        let start = 28 + 6 * i;
        let end = start + 6;
        if end > template_bytes.len() {
            return Err(NbisError::InvalidTemplate(
                "Minutia data overflow".to_string(),
            ));
        }
        let m_bytes: [u8; 6] = template_bytes[start..end].try_into().unwrap();
        minutiae.push(decode_minutia(&m_bytes));
    }

    Ok(Minutiae::new(
        minutiae,
        width as u32,
        height as u32,
        NfiqResult {
            nfiq: NfiqQuality::Unknown,
            confidence: 0.0,
        },
    ))
}

/// Extracts minutiae from an 8‑bit grayscale fingerprint image using the
/// **NIST LFS v2** algorithm.
///
/// * `image` — the bytes of any image (PNG, JPEG, etc.) that can be converted
/// * `ppi`  — (Optional) scanner resolution in dpi. Default is 500 dpi.
///
/// Returns an [`LfsResult`] on success or the negative error code propagated
/// from the C library.
#[uniffi::export]
pub fn extract_minutiae_from_image_file(
    path: &str,
    ppi: Option<f64>,
) -> Result<Minutiae, NbisError> {
    // Read the file bytes
    let image_bytes =
        std::fs::read(path).map_err(|_| NbisError::FileReadError(path.to_string()))?;

    // Call the main extraction function
    extract_minutiae(&image_bytes, ppi)
}

/// Annotates a fingerprint image with its extracted minutiae.
///
/// # Arguments
/// * `image` — the bytes of any image (PNG, JPEG, etc.) that can be converted
/// * `ppi`  — (Optional) scanner resolution in dpi. Default is 500 dpi.
/// * `min_quality` — (Optional) minimum reliability threshold for minutiae
///
/// Returns the bytes of a PNG image with the minutiae annotated.
#[uniffi::export]
pub fn annotate_minutiae(
    image: &[u8],
    ppi: Option<f64>,
    min_quality: Option<f64>,
) -> Result<Vec<u8>, NbisError> {
    // Try to load the image from bytes
    let mut image_rgb = match image::load_from_memory(image) {
        Ok(img) => match img {
            image::DynamicImage::ImageRgb8(rgb) => rgb,
            other => other.to_rgb8(),
        },
        Err(_) => return Err(NbisError::ImageLoadError),
    };

    let quality = min_quality.unwrap_or(0.5);
    if !(0.0..=1.0).contains(&quality) {
        return Err(NbisError::InvalidQuality(quality));
    }

    let minutiae = extract_minutiae(image, ppi)?;

    // Filter minutiae based on quality
    let minutiae: Vec<Arc<Minutia>> = minutiae
        .inner
        .into_iter()
        .filter(|m| m.reliability >= quality)
        .map(Arc::new)
        .collect();

    let img_w = image_rgb.width();
    let img_h = image_rgb.height();

    let square_radius = 2;
    let circle_radius = 2;

    // Draw each minutia on the annotated image
    for m in minutiae.iter() {
        let x = m.x;
        let y = m.y;
        if x >= 0 && y >= 0 && (x as u32) < img_w && (y as u32) < img_h {
            match m.kind {
                MinutiaKind::RidgeEnding => {
                    // Draw a red filled square
                    let rect = Rect::at(x - square_radius, y - square_radius).of_size(
                        (square_radius * 2 + 1) as u32,
                        (square_radius * 2 + 1) as u32,
                    );
                    draw_filled_rect_mut(&mut image_rgb, rect, Rgb([255, 0, 0]));
                }
                MinutiaKind::Bifurcation => {
                    // Draw a blue filled circle
                    draw_filled_circle_mut(&mut image_rgb, (x, y), circle_radius, Rgb([0, 0, 255]));
                }
            }

            let color = if m.kind == MinutiaKind::RidgeEnding {
                Rgb([255, 0, 0]) // Red for ridge ending
            } else {
                Rgb([0, 0, 255]) // Blue for bifurcation
            };

            draw_arrow_with_head(
                &mut image_rgb,
                (x as f32, y as f32),
                m.angle() as f32, // Convert to degrees
                15.0,             // Length of the arrow shaft
                6.0,              // Size of the arrowhead
                color,            // Green for direction
            );
        }
    }

    // Convert image_rgb to PNG bytes
    png_bytes_from_rgb(&image_rgb).map_err(|_| NbisError::ImageLoadError)
}

/// Annotates a fingerprint image with its extracted minutiae.
///
/// # Arguments
/// * `path` — the path to the image file.
/// * `ppi`  — (Optional) scanner resolution in dpi. Default is 500
/// * `min_quality` — (Optional) minimum reliability threshold for minutiae
///
/// Returns the bytes of a PNG image with the minutiae annotated.
///
/// If the image cannot be processed, returns an [`NbisError`].
#[uniffi::export]
pub fn annotate_minutiae_from_image_file(
    path: &str,
    ppi: Option<f64>,
    min_quality: Option<f64>,
) -> Result<Vec<u8>, NbisError> {
    // Read the image file bytes
    let image_bytes =
        std::fs::read(path).map_err(|_| NbisError::FileReadError(path.to_string()))?;

    annotate_minutiae(&image_bytes, ppi, min_quality)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_match() {
        // Load test images (these should be raw bytes of PNG/JPEG images)
        let p_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let p1_2 = fs::read("test_data/p1/p1_2.png").unwrap();
        let p1_3 = fs::read("test_data/p1/p1_3.png").unwrap();

        let res1 = extract_minutiae(&p_1, None).unwrap();
        let res2 = extract_minutiae(&p1_2, None).unwrap();
        let res3 = extract_minutiae(&p1_3, None).unwrap();
        let score1 = res1.compare(&res2);
        let score2 = res1.compare(&res3);
        let score3 = res2.compare(&res3);
        assert!(
            score1 > 50,
            "Match score between p1_1 and p1_2 should be greater than 50"
        );
        assert!(
            score2 > 50,
            "Match score between p1_1 and p1_3 should be greater than 50"
        );
        assert!(
            score3 > 50,
            "Match score between p1_2 and p1_3 should be greater than 50"
        );

        let p2_1 = fs::read("test_data/p2/p2_1.png").unwrap();
        let p2_2 = fs::read("test_data/p2/p2_2.png").unwrap();
        let p2_3 = fs::read("test_data/p2/p2_3.png").unwrap();

        let res4 = extract_minutiae(&p2_1, None).unwrap();
        let res5 = extract_minutiae(&p2_2, None).unwrap();
        let res6 = extract_minutiae(&p2_3, None).unwrap();
        let score4 = res4.compare(&res5);
        let score5 = res4.compare(&res6);
        let score6 = res5.compare(&res6);

        assert!(
            score4 > 50,
            "Match score between p2_1 and p2_2 should be greater than 50"
        );
        assert!(
            score5 > 50,
            "Match score between p2_1 and p2_3 should be greater than 50"
        );
        assert!(
            score6 > 50,
            "Match score between p2_2 and p2_3 should be greater than 50"
        );

        // Inter-fingerprint matching should yield lower scores
        let score7 = res1.compare(&res4);
        let score8 = res1.compare(&res5);
        let score9 = res1.compare(&res6);

        assert!(
            score7 < 50,
            "Match score between p1_1 and p2_1 should be less than 50"
        );
        assert!(
            score8 < 50,
            "Match score between p1_1 and p2_2 should be less than 50"
        );
        assert!(
            score9 < 50,
            "Match score between p1_1 and p2_3 should be less than 50"
        );
    }

    #[test]
    fn test_encode_to_iso() {
        let bryanc_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let res = extract_minutiae(&bryanc_1, None).unwrap();
        let encoded = res.to_iso_19794_2_2005();
        assert!(!encoded.is_empty(), "Encoded ISO data should not be empty");

        let minutiae = load_iso_19794_2_2005(&encoded).unwrap();
        assert_eq!(
            minutiae.inner.len(),
            res.inner.len(),
            "Decoded minutiae count should match original"
        );
        assert_eq!(
            minutiae.img_w, res.img_w,
            "Decoded image width should match original"
        );
        assert_eq!(
            minutiae.img_h, res.img_h,
            "Decoded image height should match original"
        );
        assert_eq!(
            minutiae.img_w, res.img_w,
            "Decoded image width should match original"
        );
        assert_eq!(
            minutiae.img_h, res.img_h,
            "Decoded image height should match original"
        );
        // All the minutiae should match
        for (m1, m2) in res.inner.iter().zip(minutiae.inner.iter()) {
            assert_eq!(m1.x, m2.x, "X coordinate should match");
            assert_eq!(m1.y, m2.y, "Y coordinate should match");
            assert_eq!(m1.direction, m2.direction, "Direction should match");
            // reliability is a float, so allow some tolerance
            assert!(
                (m1.reliability - m2.reliability).abs() < 1e-1,
                "Reliability should match within tolerance"
            );
            assert_eq!(m1.kind, m2.kind, "Kind should match");
        }
    }

    #[test]
    fn test_nfiq() {
        let p1_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let res = extract_minutiae(&p1_1, None).unwrap();
        assert!(
            (0.0..=1.0).contains(&res.quality().confidence),
            "Confidence should be between 0.0 and 1.0"
        );
        // Quality should be very good for this image
        assert!(
            res.quality().nfiq == NfiqQuality::Excellent,
            "NFIQ for p1_1 should be Excellent"
        );

        // Test a non-fingerprint image
        let random_image = fs::read("test_data/negative/landscape.jpg").unwrap();
        let res2 = extract_minutiae(&random_image, None).unwrap();
        // The quality should be poorest for non-fingerprint images
        assert!(
            res2.quality().nfiq == NfiqQuality::Poor,
            "NFIQ for non-fingerprint image should be Poor"
        );

        // Test a non-fingerprint image
        let random_image = fs::read("test_data/negative/face.jpeg").unwrap();
        let res2 = extract_minutiae(&random_image, None).unwrap();
        // The quality should be poorest for non-fingerprint images
        assert!(
            res2.quality().nfiq == NfiqQuality::Poor,
            "NFIQ for non-fingerprint image should be Poor"
        );
    }

    #[test]
    fn test_negative() {
        //Try to extract minutae from a file that is not an image
        let res1 = extract_minutiae_from_image_file("build.rs", None);

        // Check if the result is an error
        assert!(res1.is_err(), "Expected an error but got Ok");

        match res1 {
            Err(NbisError::ImageLoadError) => {
                // This is the expected variant — success!
            }
            Err(other) => panic!("Expected ImageLoadError but got: {:?}", other),
            Ok(_) => panic!("Expected error but got Ok"),
        }

        //Try to extract minutae from a file that does not exist
        let res2 = extract_minutiae_from_image_file("test_data/negative/x.png", None);

        // Check if the result is an error
        assert!(res2.is_err(), "Expected an error but got Ok");

        match res2 {
            Err(NbisError::FileReadError(_)) => {
                // This is the expected variant — success!
            }
            Err(other) => panic!("Expected FileReadError but got: {:?}", other),
            Ok(_) => panic!("Expected error but got Ok"),
        }

        // // Test with an image (neither face nor fingerprint)
        // let n_1 = fs::read("test_data/negative/no_face.jpeg").unwrap();

        // let res1 = extract_minutiae(&n_1, None).unwrap();
        // let res2 = extract_minutiae(&n_1, None).unwrap();
        // let score = res1.compare(&res2);
        // println!("{:?}", score);

        // Test with a face image
        let n_2 = fs::read("test_data/negative/varun_square.png").unwrap();

        let res1_n_2 = extract_minutiae(&n_2, None).unwrap();
        let res2_n_2 = extract_minutiae(&n_2, None).unwrap();
        let score_n_2 = res1_n_2.compare(&res2_n_2);
        println!("{:?}", score_n_2);
    }

    #[test]
    fn test_sivv() {
        let p1_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let sivv_result = sivv(&p1_1).unwrap();
        println!("SIVV result: {:?}", sivv_result);
    }
}
