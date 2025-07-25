use std::ffi::CStr;
use std::ptr::{null_mut, NonNull};
use std::sync::Arc;

use image::DynamicImage;
use std::os::raw::{c_int, c_uchar, c_void};

use image::Rgb;
use imageproc::drawing::{draw_filled_circle_mut, draw_filled_rect_mut};
use imageproc::rect::Rect;

use crate::consts::MM_PER_INCH;
use crate::errors::NbisError;
use crate::ffi::{
    comp_nfiq_featvctr, dflt_acfunc_hids, dflt_acfunc_outs, dflt_nHids, dflt_nInps, dflt_nOuts,
    dflt_wts, dflt_znorm_means, dflt_znorm_stds, free, free_minutiae, get_minutiae, runmlp2,
    sivv_ffi_free_bytes, sivv_ffi_from_bytes, znorm_fniq_featvctr, CPoint2i, LFSPARMS, MINUTIAE,
    MIN_MINUTIAE, NFIQ_NUM_CLASSES, NFIQ_VCTRLEN,
};
use crate::imutils::{draw_arrow_with_head, png_bytes_from_rgb};
use crate::{Minutia, MinutiaKind, Minutiae};

/// Represents the result of the SIVV computation.
pub(crate) struct SIVVResult {
    /// Index of the largest peak-valley pair (1-based)
    pub(crate) largest_pvp_index: i32,

    /// Total number of detected peak-valley pairs
    pub(crate) total_pvps: i32,

    /// Power difference between the peak and valley
    pub(crate) power_diff: f64,

    /// Frequency difference between the peak and valley
    pub(crate) freq_diff: f64,

    /// Slope between valley and peak (dy / dx)
    pub(crate) slope: f64,

    /// Frequency of the midpoint between valley and peak
    pub(crate) center_frequency: f64,

    /// Absolute frequency of the peak (undocumented in comments)
    pub(crate) peak_frequency: f64,
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

#[derive(Debug, Clone, uniffi::Record)]
pub struct Point {
    /// The x-coordinate of the point.
    pub x: i32,
    /// The y-coordinate of the point.
    pub y: i32,
}

#[derive(Debug, Clone, uniffi::Record)]
pub struct ROI {
    /// The x-coordinate of the top-left corner of the ROI.
    pub x1: i32,
    /// The y-coordinate of the top-left corner of the ROI.
    pub y1: i32,
    /// The x-coordinate of the bottom-right corner of the ROI.
    pub x2: i32,
    /// The y-coordinate of the bottom-right corner of the ROI.
    pub y2: i32,
    /// The center point of the ROI.
    pub center: Point,
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

    /// Encodes `NfiqQuality` as a 0–100 ISO quality byte.
    pub fn to_iso_quality(&self) -> u8 {
        match self {
            NfiqQuality::Excellent => 100,
            NfiqQuality::VeryGood => 80,
            NfiqQuality::Good => 60,
            NfiqQuality::Fair => 40,
            NfiqQuality::Poor => 20,
            NfiqQuality::Unknown => 0,
        }
    }

    /// Decodes a 0–100 ISO quality byte into `NfiqQuality`.
    pub fn from_iso_quality(value: u8) -> Self {
        match value {
            90..=100 => NfiqQuality::Excellent,
            70..=89 => NfiqQuality::VeryGood,
            50..=69 => NfiqQuality::Good,
            30..=49 => NfiqQuality::Fair,
            1..=29 => NfiqQuality::Poor,
            _ => NfiqQuality::Unknown,
        }
    }
}

fn is_fingerprint(result: &SIVVResult) -> bool {
    // The following values are from evaluation of the SIVV algorithm
    // on a mixed biometric dataset.
    let max_peak_freq = 0.15; // cycles/pixel
    let peak_height_threshold = 0.02;
    let _ = result.largest_pvp_index; // 1-based index, not used here
    let _ = result.total_pvps; // total number of peak-valley pairs, not used here
    let _ = result.freq_diff; // frequency difference, not used here
    let _ = result.slope; // slope, not used here
    let _ = result.center_frequency; // center frequency, not used here

    result.peak_frequency < max_peak_freq && result.power_diff > peak_height_threshold
}

// Safe Rust wrapper
#[allow(clippy::type_complexity)]
pub fn find_fingerprint_center(
    data: *const u8,
    width: c_int,
    height: c_int,
) -> Result<(CPoint2i, (i32, i32, i32, i32)), Box<dyn std::error::Error>> {
    let mut xbound_min: c_int = 0;
    let mut xbound_max: c_int = 0;
    let mut ybound_min: c_int = width;
    let mut ybound_max: c_int = height;

    // Call the C function
    let result = unsafe {
        crate::ffi::find_fingerprint_center_morph_c(
            data,
            width,
            height,
            &mut xbound_min,
            &mut xbound_max,
            &mut ybound_min,
            &mut ybound_max,
        )
    };

    // let point = opencv::core::Point2i::new(result.x, result.y);
    let bounds = (xbound_min, xbound_max, ybound_min, ybound_max);

    Ok((result, bounds))
}

fn sivv(image: *mut c_uchar, width: i32, height: i32) -> Result<SIVVResult, NbisError> {
    unsafe {
        let ptr = sivv_ffi_from_bytes(
            image,
            width as std::os::raw::c_int,
            height as std::os::raw::c_int,
        );
        let str_result = CStr::from_ptr(ptr).to_string_lossy().into_owned();
        sivv_ffi_free_bytes(ptr);

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

    // 2) Check SIVV result -----------------------------------
    let sivv_result = sivv(gray.as_ptr() as *mut c_uchar, iw as i32, ih as i32)?;
    if !is_fingerprint(&sivv_result) {
        // Early return if the image is not a fingerprint
        return Ok(Minutiae::new(
            Vec::new(),
            iw,
            ih,
            NfiqResult {
                nfiq: NfiqQuality::Unknown,
                confidence: 0.0,
            },
            None, // No ROI in this case
        ));
    }

    //get the finger print center
    let center = find_fingerprint_center(gray.as_ptr() as *mut c_uchar, iw as c_int, ih as c_int)
        .map_err(|e| NbisError::GenericError(e.to_string()))?;

    let roi = ROI {
        x1: center.1 .0,
        x2: center.1 .1,
        y1: center.1 .2,
        y2: center.1 .3,
        center: Point {
            x: center.0.x,
            y: center.0.y,
        },
    };

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
        Minutiae::new(minutiae_vec, iw, ih, quality, Some(roi))
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
    crate::encoding::load_iso_19794_2_2005(template_bytes)
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
        let encoded = res.to_iso_19794_2_2005(0.0);
        assert!(!encoded.is_empty(), "Encoded ISO data should not be empty");

        let minutiae = load_iso_19794_2_2005(&encoded).unwrap();

        // Qualiity should match the original
        assert_eq!(
            res.quality().nfiq,
            minutiae.quality().nfiq,
            "NFIQ quality should match original"
        );

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

        // Test encode with more than 255 minutiae
        let mut many_minutiae = res.inner.clone();
        // Add dummy minutiae to exceed 255
        for i in 0..300 {
            many_minutiae.push(Minutia {
                x: i as i32,
                y: i as i32,
                direction: 0,
                reliability: 0.0,
                kind: MinutiaKind::RidgeEnding,
            });
        }

        let many_res = Minutiae::new(many_minutiae, res.img_w, res.img_h, res.nfiq, None);
        let many_encoded = many_res.to_iso_19794_2_2005(0.0);
        assert!(
            !many_encoded.is_empty(),
            "Encoded ISO data should not be empty"
        );

        let many_minutiae_decoded = load_iso_19794_2_2005(&many_encoded).unwrap();
        assert_eq!(
            many_minutiae_decoded.inner.len(),
            255,
            "Decoded minutiae count should be capped at 255"
        );
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
            res2.quality().nfiq == NfiqQuality::Unknown,
            "NFIQ for non-fingerprint image should be Unknown"
        );

        // Test a non-fingerprint image
        let random_image = fs::read("test_data/negative/face.jpeg").unwrap();
        let res2 = extract_minutiae(&random_image, None).unwrap();
        // The quality should be poorest for non-fingerprint images
        assert!(
            res2.quality().nfiq == NfiqQuality::Unknown,
            "NFIQ for non-fingerprint image should be Unknown"
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
        assert_eq!(score_n_2, 0);
    }

    #[test]
    fn test_roi() {
        let p1_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let res = extract_minutiae(&p1_1, None).unwrap();
        assert!(res.roi().is_some(), "Expected ROI to be present");
        let roi = res.roi().unwrap();

        assert!(
            roi.x1 < roi.x2 && roi.y1 < roi.y2,
            "ROI coordinates should be valid"
        );
        assert!(roi.x1 == 0, "Expected ROI x1 to be 0");
        assert!(roi.y1 == 96, "Expected ROI y1 to be 96");
        assert!(roi.x2 == 382, "Expected ROI x2 to be 382");
        assert!(roi.y2 == 496, "Expected ROI y2 to be 496");
        assert_eq!(roi.center.x, 182, "Expected ROI center x to be 182");
        assert_eq!(roi.center.y, 296, "Expected ROI center y to be 296");

        // Uncomment the following lines to visualize the ROI on the image
        // // Load the original image to draw the ROI
        // let mut image_rgb = match image::load_from_memory(&p1_1) {
        //     Ok(img) => match img {
        //         image::DynamicImage::ImageRgb8(rgb) => rgb,
        //         other => other.to_rgb8(),
        //     },
        //     Err(_) => panic!("Failed to load image"),
        // };
        // let rect = Rect::at(roi.x1, roi.y1).of_size(
        //     (roi.x2 - roi.x1) as u32,
        //     (roi.y2 - roi.y1) as u32,
        // );
        // // Draw a red (hollow) rectangle around the ROI
        // imageproc::drawing::draw_hollow_rect_mut(
        //     &mut image_rgb,
        //     rect,
        //     Rgb([255, 0, 0]),
        // );

        // // Draw a cross at the center of the ROI
        // let center_x = roi.center.x;
        // let center_y = roi.center.y;

        // // Horizontal line
        // imageproc::drawing::draw_cross_mut(
        //     &mut image_rgb,
        //     Rgb([0, 255, 0]),
        //     center_x,
        //     center_y,
        // );

        // // Save the annotated image to verify the ROI visually
        // let annotated_path = "test_data/p1/p1_1_roi.png";
        // image_rgb
        //     .save(annotated_path)
        //     .expect("Failed to save annotated image with ROI");
    }
}
