use std::{
    ffi::CStr,
    os::raw::{c_char, c_uint, c_ushort},
    ptr,
};

use crate::{
    ffi_nfiq2::{
        nfiq2wrapper_compute, nfiq2wrapper_create, nfiq2wrapper_destroy, nfiq2wrapper_free_results,
        Nfiq2ResultsT, Nfiq2WrapperOpaque,
    },
    NbisError,
};

#[derive(Debug, Clone, uniffi::Record)]
pub struct Nfiq2Value {
    pub name: String,
    pub value: f64,
}

/// Safe Rust view of the results
#[derive(Debug, Clone, uniffi::Record)]
pub struct Nfiq2Result {
    pub score: u32,
    pub actionable: Vec<Nfiq2Value>,
    pub features: Vec<Nfiq2Value>,
}

/// The high‐level Rust handle
#[derive(Debug, Clone, uniffi::Object)]
pub struct Nfiq2 {
    ctx: *mut Nfiq2WrapperOpaque,
}

unsafe impl Send for Nfiq2 {}
unsafe impl Sync for Nfiq2 {}

/// Construct a new wrapper, or Err if allocation/initialization fails.
pub fn new_nfiq2() -> Result<Nfiq2, NbisError> {
    let ptr = unsafe { nfiq2wrapper_create() };
    if ptr.is_null() {
        Err(NbisError::Nfiq2CreateFailed)
    } else {
        Ok(Nfiq2 { ctx: ptr })
    }
}

impl Nfiq2 {
    /// Compute quality. Mirrors your C API.
    pub fn compute(&self, image_bytes: &[u8]) -> Result<Nfiq2Result, NbisError> {
        if self.ctx.is_null() {
            return Err(NbisError::Nfiq2NullContext);
        }

        // load the image from bytes
        let image =
            image::load_from_memory(image_bytes).map_err(|_| NbisError::Nfiq2ComputeFailed(-1))?;

        // convert to grayscale and get dimensions
        let image = image.to_luma8();

        let (cols, rows) = image.dimensions();
        let ppi = 500; // hardcoded PPI, can be adjusted as needed

        // zero the C struct
        let mut raw: Nfiq2ResultsT = unsafe { std::mem::zeroed() };

        let size = image.len() as c_uint;

        let rc = unsafe {
            nfiq2wrapper_compute(
                self.ctx,
                image.as_ptr(),
                size as c_uint,
                cols as c_uint,
                rows as c_uint,
                ppi as c_ushort,
                &mut raw,
            )
        };
        if rc != 0 {
            // free any partial allocations before returning
            unsafe { nfiq2wrapper_free_results(&mut raw) };
            return Err(NbisError::Nfiq2ComputeFailed(rc));
        }

        // helper to turn C arrays into Vec<(String,f64)>
        unsafe fn collect_pairs(
            ids_ptr: *const *const c_char,
            vals_ptr: *const f64,
            count: usize,
        ) -> Result<Vec<Nfiq2Value>, NbisError> {
            let mut out = Vec::with_capacity(count);
            let id_slice = std::slice::from_raw_parts(ids_ptr, count);
            let val_slice = std::slice::from_raw_parts(vals_ptr, count);
            for i in 0..count {
                let s = CStr::from_ptr(id_slice[i])
                    .to_str()
                    .map_err(|_| NbisError::Nfiq2ComputeFailed(-1))?
                    .to_string();

                out.push(Nfiq2Value {
                    name: s,
                    value: val_slice[i],
                });
            }
            Ok(out)
        }

        let actionable_count = raw.actionable_count as usize;
        let feature_count = raw.feature_count as usize;

        // collect actionable + features
        let actionable = unsafe {
            collect_pairs(
                raw.actionable_ids,
                raw.actionable_values as *const f64,
                actionable_count,
            )?
        };
        let features = unsafe {
            collect_pairs(
                raw.feature_ids,
                raw.feature_values as *const f64,
                feature_count,
            )?
        };

        let score = raw.score;

        // free C allocations
        unsafe { nfiq2wrapper_free_results(&mut raw) };

        Ok(Nfiq2Result {
            score,
            actionable,
            features,
        })
    }
}

impl Drop for Nfiq2 {
    fn drop(&mut self) {
        if !self.ctx.is_null() {
            unsafe { nfiq2wrapper_destroy(self.ctx) };
            self.ctx = ptr::null_mut();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nfiq2() {
        // construct wrapper
        let nfiq = new_nfiq2().expect("failed to create wrapper");

        let expected_scores = vec![54, 45, 53, 52, 57];
        let input_images = vec![
            "ext/NFIQ2-2.3.0/examples/images/SFinGe_Test01.pgm",
            "ext/NFIQ2-2.3.0/examples/images/SFinGe_Test02.pgm",
            "ext/NFIQ2-2.3.0/examples/images/SFinGe_Test03.pgm",
            "ext/NFIQ2-2.3.0/examples/images/SFinGe_Test04.pgm",
            "ext/NFIQ2-2.3.0/examples/images/SFinGe_Test05.pgm",
        ];

        for (i, img_path) in input_images.iter().enumerate() {
            // load test image bytes
            let img_bytes = std::fs::read(img_path).expect("failed to read test image");

            // call compute
            let res = nfiq.compute(&img_bytes).expect("compute failed");

            // check score
            assert_eq!(res.score, expected_scores[i]);
        }
    }
}
