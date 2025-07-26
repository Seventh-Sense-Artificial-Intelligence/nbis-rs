use once_cell::sync::Lazy;
use std::sync::{Arc, Mutex};

use crate::api::NfiqResult;
use crate::minutia::Minutia;
use crate::{bozorth::bz_match_score, encoding::to_nist_xyt_set};
/// A set of minutiae extracted from a fingerprint image.
#[derive(Debug, Clone, uniffi::Object)]
pub struct Minutiae {
    pub(crate) inner: Vec<Minutia>,
    pub(crate) img_w: u32, // original image width for XYT conversion
    pub(crate) img_h: u32, // original image height for XYT conversion
    pub(crate) nfiq: NfiqResult,
}

impl Minutiae {
    pub(crate) fn new(minutiae: Vec<Minutia>, img_w: u32, img_h: u32, nfiq: NfiqResult) -> Self {
        Minutiae {
            inner: minutiae.to_vec(),
            img_w,
            img_h,
            nfiq,
        }
    }
}

static BOZORTH_MUTEX: Lazy<Mutex<()>> = Lazy::new(|| Mutex::new(()));

#[uniffi::export]
impl Minutiae {
    /// Similarity via Bozorth‑3 (higher = more similar). A score > 50 is a likely match.
    ///
    /// # Arguments
    /// * `other` — another `Minutiae` object to compare against.
    ///
    /// Returns an `i32` score representing the similarity between the two sets of minutiae.
    /// A higher score indicates more similarity.
    pub fn compare(&self, other: &Minutiae) -> i32 {
        // Ensure we have a lock to prevent concurrent access issues as Bozorth is not thread-safe.
        let _lock = BOZORTH_MUTEX.lock().unwrap();
        let p = to_nist_xyt_set(self);
        let g = to_nist_xyt_set(other);
        //println!("Matching {} vs {}", p.xs.len(), g.xs.len());
        let score = bz_match_score(&p, &g);
        // #define QQ_SIZE 4000
        // #define QQ_OVERFLOW_SCORE QQ_SIZE
        // If the score is QQ_SIZE, it indicates an overflow condition.
        if score == 4000 {
            0 // Just return 0 for overflow
        } else {
            score // Return the actual score
        }
    }

    pub fn quality(&self) -> NfiqResult {
        self.nfiq.clone()
    }

    /// Returns a vector of `Minutia` objects representing the minutiae in this set.
    pub fn get(&self) -> Vec<Arc<Minutia>> {
        self.inner.iter().cloned().map(Arc::new).collect()
    }

    pub fn to_iso_19794_2_2005(&self, min_quality: f64) -> Vec<u8> {
        crate::encoding::to_iso_19794_2_2005(self, min_quality)
    }
}
