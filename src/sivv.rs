use std::ffi::CStr;

use std::os::raw::{c_int, c_uchar};

use crate::errors::NbisError;
use crate::ffi_nbis::{sivv_ffi_free_bytes, sivv_ffi_from_bytes, CPoint2i};
use crate::structs::SIVVResult;

pub(crate) fn is_fingerprint(result: &SIVVResult) -> bool {
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
pub(crate) fn find_fingerprint_center(
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
        crate::ffi_nbis::find_fingerprint_center_morph_c(
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

pub(crate) fn sivv(image: *mut c_uchar, width: i32, height: i32) -> Result<SIVVResult, NbisError> {
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
