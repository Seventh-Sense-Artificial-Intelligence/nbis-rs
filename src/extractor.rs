use std::{
    os::raw::{c_int, c_uchar, c_void},
    ptr::{null_mut, NonNull},
};

use image::{DynamicImage, Rgb};
use imageproc::{
    drawing::{draw_filled_circle_mut, draw_filled_rect_mut},
    rect::Rect,
};

use crate::{
    consts::MM_PER_INCH,
    ffi_nbis::{free, free_minutiae, get_minutiae, DEFAULT_BOZORTH_MINUTIAE, LFSPARMS, MINUTIAE},
    imutils::{draw_arrow_with_head, png_bytes_from_rgb},
    nfiq2_api::{new_nfiq2, Nfiq2},
    sivv::{find_fingerprint_center, is_fingerprint, sivv},
    structs::NbisExtractorSettings,
    Minutia, MinutiaKind, Minutiae, NbisError, Nfiq2Result, Point, ROI,
};

#[derive(Debug, Clone, uniffi::Object)]
pub struct NbisExtractor {
    settings: NbisExtractorSettings,
    nfiq2: Nfiq2,
}

#[uniffi::export]
pub fn new_nbis_extractor(settings: NbisExtractorSettings) -> Result<NbisExtractor, NbisError> {
    let nfiq2 = new_nfiq2()?;
    Ok(NbisExtractor { settings, nfiq2 })
}

impl NbisExtractor {
    pub fn new(settings: NbisExtractorSettings) -> Result<Self, NbisError> {
        Ok(NbisExtractor {
            settings,
            nfiq2: new_nfiq2()?,
        })
    }
}

#[uniffi::export]
impl NbisExtractor {
    pub fn settings(&self) -> NbisExtractorSettings {
        self.settings.clone()
    }

    pub fn load_iso_19794_2_2005(&self, template_bytes: &[u8]) -> Result<Minutiae, NbisError> {
        crate::encoding::load_iso_19794_2_2005(template_bytes)
    }

    pub fn annotate_minutiae_from_image_file(&self, path: &str) -> Result<Vec<u8>, NbisError> {
        // Load the image from the file
        let image = std::fs::read(path).map_err(|_| NbisError::FileReadError(path.to_string()))?;
        self.annotate_minutiae(&image)
    }

    pub fn annotate_minutiae(&self, image: &[u8]) -> Result<Vec<u8>, NbisError> {
        // Try to load the image from bytes
        let mut image_rgb = match image::load_from_memory(image) {
            Ok(img) => match img {
                image::DynamicImage::ImageRgb8(rgb) => rgb,
                other => other.to_rgb8(),
            },
            Err(_) => return Err(NbisError::ImageLoadError),
        };

        let minutiae = self.extract_minutiae(image)?;

        let img_w = image_rgb.width();
        let img_h = image_rgb.height();

        let square_radius = 2;
        let circle_radius = 2;

        // Draw each minutia on the annotated image
        for m in minutiae.inner.iter() {
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
                        draw_filled_circle_mut(
                            &mut image_rgb,
                            (x, y),
                            circle_radius,
                            Rgb([0, 0, 255]),
                        );
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

    pub fn extract_minutiae_from_image_file(&self, file_path: &str) -> Result<Minutiae, NbisError> {
        // Read the file bytes
        let image_bytes = std::fs::read(file_path)
            .map_err(|_| NbisError::FileReadError(file_path.to_string()))?;

        // Call the main extraction function
        self.extract_minutiae(&image_bytes)
    }

    pub fn extract_minutiae(&self, image_bytes: &[u8]) -> Result<Minutiae, NbisError> {
        let ppi = self.settings.ppi.unwrap_or(500.0); // default to 500 dpi

        // 1) Load the image ------------------------------------------------------
        let image = match image::load_from_memory(image_bytes) {
            Ok(img) => img,
            Err(_e) => return Err(NbisError::ImageLoadError),
        };

        // 2) Ensure 8‑bit grayscale ------------------------------------------------
        let gray = match image {
            DynamicImage::ImageLuma8(buf) => buf.clone(),
            _ => image.to_luma8(),
        };
        let (iw, ih) = gray.dimensions();

        // 3) Check SIVV result -----------------------------------
        if self.settings.check_fingerprint {
            let sivv_result = sivv(gray.as_ptr() as *mut c_uchar, iw as i32, ih as i32)?;
            if !is_fingerprint(&sivv_result) {
                // Early return if the image is not a fingerprint
                return Ok(Minutiae::new(
                    Vec::new(),
                    iw,
                    ih,
                    Nfiq2Result {
                        score: 0,
                        actionable: Vec::new(),
                        features: Vec::new(),
                    },
                    None, // No ROI in this case
                ));
            }
        }

        // 4) Get ROI -----------------------------------
        let roi = if self.settings.get_center {
            let center =
                find_fingerprint_center(gray.as_ptr() as *mut c_uchar, iw as c_int, ih as c_int)
                    .map_err(|e| NbisError::GenericError(e.to_string()))?;

            Some(ROI {
                x1: center.1 .0,
                x2: center.1 .1,
                y1: center.1 .2,
                y2: center.1 .3,
                center: Point {
                    x: center.0.x,
                    y: center.0.y,
                },
            })
        } else {
            None
        };

        // 5) Define buffers and sizes returned by the C API -------------------------------
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

        // 6) Call into C ---------------------------------------------------------
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

        // 7) Compute NFIQv2 quality assessment -----------------------------
        let quality = if self.settings.compute_nfiq2 {
            self.nfiq2.compute(image_bytes)?
        } else {
            Nfiq2Result {
                score: 0,
                actionable: Vec::new(),
                features: Vec::new(),
            }
        };

        // 8) Convert C results -----------------------------
        let minutiae = NonNull::new(ominutiae).expect("C returned null pointer");
        let mut minutiae_obj = unsafe {
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
            Minutiae::new(minutiae_vec, iw, ih, quality, roi)
        };

        // 7) filter, truncate and sort -----------------------------
        if self.settings.min_quality > 0.0 {
            minutiae_obj
                .inner
                .retain(|m| m.reliability >= self.settings.min_quality);
        }

        // Truncate to the default number of minutiae if necessary
        if minutiae_obj.inner.len() > DEFAULT_BOZORTH_MINUTIAE {
            minutiae_obj.inner.sort_by(|a, b| {
                b.reliability
                    .partial_cmp(&a.reliability)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            minutiae_obj.inner.truncate(DEFAULT_BOZORTH_MINUTIAE);
        }

        // Sort by x and then by y
        minutiae_obj
            .inner
            .sort_by(|a, b| a.x.cmp(&b.x).then(a.y.cmp(&b.y)));

        // 8) Free C allocations we no longer need -------------------------------
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
}

#[cfg(test)]
mod tests {
    use crate::ffi_nbis::DEFAULT_BOZORTH_MINUTIAE;

    use super::*;
    use std::fs;

    #[test]
    fn test_match() {
        // Load test images (these should be raw bytes of PNG/JPEG images)
        let p_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let p1_2 = fs::read("test_data/p1/p1_2.png").unwrap();
        let p1_3 = fs::read("test_data/p1/p1_3.png").unwrap();

        let extractor = new_nbis_extractor(NbisExtractorSettings::default()).unwrap();

        let res1 = extractor.extract_minutiae(&p_1).unwrap();
        let res2 = extractor.extract_minutiae(&p1_2).unwrap();
        let res3 = extractor.extract_minutiae(&p1_3).unwrap();
        let score1 = res1.compare(&res2);

        // Test that score is symmetric
        assert_eq!(score1, res2.compare(&res1), "Scores should be symmetric");

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

        let res4 = extractor.extract_minutiae(&p2_1).unwrap();
        let res5 = extractor.extract_minutiae(&p2_2).unwrap();
        let res6 = extractor.extract_minutiae(&p2_3).unwrap();
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
        let extractor = new_nbis_extractor(NbisExtractorSettings::default()).unwrap();
        let bryanc_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let res = extractor.extract_minutiae(&bryanc_1).unwrap();
        let encoded = res.to_iso_19794_2_2005();
        assert!(!encoded.is_empty(), "Encoded ISO data should not be empty");

        let minutiae = extractor.load_iso_19794_2_2005(&encoded).unwrap();

        // Qualiity should match the original
        assert_eq!(
            res.quality().score,
            minutiae.quality().score,
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
        let many_encoded = many_res.to_iso_19794_2_2005();
        assert!(
            !many_encoded.is_empty(),
            "Encoded ISO data should not be empty"
        );

        let many_minutiae_decoded = extractor.load_iso_19794_2_2005(&many_encoded).unwrap();
        assert_eq!(
            many_minutiae_decoded.inner.len(),
            DEFAULT_BOZORTH_MINUTIAE,
            "{}",
            &format!("Decoded minutiae count should be capped at {DEFAULT_BOZORTH_MINUTIAE}")
        );

        // Check if the match score before and after encoding is the same
        let bryanc_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let bryanc_2 = fs::read("test_data/p1/p1_2.png").unwrap();
        let r1 = extractor.extract_minutiae(&bryanc_1).unwrap();
        let r2 = extractor.extract_minutiae(&bryanc_2).unwrap();
        let e1 = r1.to_iso_19794_2_2005();
        let e2 = r2.to_iso_19794_2_2005();
        let reloaded_e1 = extractor.load_iso_19794_2_2005(&e1).unwrap();
        let reloaded_e2 = extractor.load_iso_19794_2_2005(&e2).unwrap();

        let s1 = r1.compare(&r2);
        let s2 = r1.compare(&reloaded_e2);
        let s3 = reloaded_e1.compare(&r2);
        let s4 = reloaded_e1.compare(&reloaded_e2);

        // Decode (Boolean, Boolean) as if the template was loaded from a file
        assert_eq!(
            r1.inner.len(),
            reloaded_e1.inner.len(),
            "Minutiae count should match after encoding"
        );
        // println!("r1 len: {}, reloaded_e1 len: {}", r1.inner.len(), reloaded_e1.inner.len());
        assert_eq!(
            s1, s2,
            "Match score should be the same for (False, False) vs (False True)"
        );
        assert_eq!(
            s1, s3,
            "Match score should be the same for (False, False) vs (True False)"
        );
        assert_eq!(
            s2, s4,
            "Match score should be the same for (False, True) vs (True True)"
        );
        assert_eq!(
            s3, s4,
            "Match score should be the same for (True, False) vs (True True)"
        );
        assert_eq!(
            s1, s4,
            "Match score should be the same for (False, False) vs (True True)"
        );
    }

    #[test]
    fn test_nfiq() {
        let extractor = new_nbis_extractor(NbisExtractorSettings::default()).unwrap();
        let p1_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let res = extractor.extract_minutiae(&p1_1).unwrap();

        // Quality should be very good for this image
        assert!(res.quality().score > 60, "NFIQ for p1_1 should > 60");

        // Test a non-fingerprint image
        let random_image = fs::read("test_data/negative/landscape.jpg").unwrap();
        let extractor = new_nbis_extractor(NbisExtractorSettings {
            min_quality: 0.0,
            get_center: false,
            check_fingerprint: true,
            compute_nfiq2: true,
            ppi: None,
        })
        .unwrap();
        let res2 = extractor.extract_minutiae(&random_image).unwrap();
        // The quality should be poorest for non-fingerprint images
        assert!(
            res2.quality().score == 0,
            "NFIQ for non-fingerprint image should be Unknown"
        );

        // Test a non-fingerprint image
        let random_image = fs::read("test_data/negative/face.jpeg").unwrap();
        let res2 = extractor.extract_minutiae(&random_image).unwrap();
        // The quality should be poorest for non-fingerprint images
        assert!(
            res2.quality().score == 0,
            "NFIQ for non-fingerprint image should be Unknown"
        );
    }

    #[test]
    fn test_negative() {
        let extractor = new_nbis_extractor(NbisExtractorSettings::default()).unwrap();
        //Try to extract minutae from a file that is not an image
        let res1 = extractor.extract_minutiae_from_image_file("build.rs");

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
        let res2 = extractor.extract_minutiae_from_image_file("test_data/negative/x.png");

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

        let extractor = new_nbis_extractor(NbisExtractorSettings {
            min_quality: 0.0,
            get_center: false,
            check_fingerprint: true,
            compute_nfiq2: false,
            ppi: None,
        })
        .unwrap();

        let res1_n_2 = extractor.extract_minutiae(&n_2).unwrap();
        let res2_n_2 = extractor.extract_minutiae(&n_2).unwrap();
        let score_n_2 = res1_n_2.compare(&res2_n_2);
        assert_eq!(score_n_2, 0);
    }

    #[test]
    fn test_roi() {
        let p1_1 = fs::read("test_data/p1/p1_1.png").unwrap();
        let extractor = new_nbis_extractor(NbisExtractorSettings {
            min_quality: 0.0,
            get_center: true,
            check_fingerprint: false,
            compute_nfiq2: false,
            ppi: None,
        })
        .unwrap();
        let res = extractor.extract_minutiae(&p1_1).unwrap();
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
