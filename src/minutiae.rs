use once_cell::sync::Lazy;
use std::sync::{Arc, Mutex};

use crate::api::NfiqResult;
use crate::minutia::Minutia;
use crate::{
    bozorth::bz_match_score,
    encoding::{encode_minutia, to_nist_xyt_set},
};
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

    /// Convert this set of minutiae into an ISO/IEC 19794-2:2005 template.
    /// The result is a `Vec<u8>` containing the full template data.
    pub fn to_iso_19794_2_2005(&self) -> Vec<u8> {
        let total_bytes = self.inner.len() * 6 + 28 + 2;

        let total_bytes_u32 =
            u32::try_from(total_bytes).expect("template larger than 4 294 967 295 bytes");

        let width = self.img_w as u16;
        let height = self.img_h as u16;

        let mut buf = Vec::with_capacity(28);

        // "FMR\0 20\0"  – 8 bytes
        buf.extend_from_slice(b"FMR\0 20\0");

        // 4-byte total length (big-endian)
        buf.extend_from_slice(&total_bytes_u32.to_be_bytes());

        // Two 0x00 padding bytes
        buf.extend_from_slice(&[0x00, 0x00]);

        // Width & height (big-endian, 2 bytes each)
        buf.extend_from_slice(&width.to_be_bytes());
        buf.extend_from_slice(&height.to_be_bytes());

        // Constant block 00 C5 00 C5 01 00 00 00 64
        buf.extend_from_slice(&[0x00, 0xC5, 0x00, 0xC5, 0x01, 0x00, 0x00, 0x00, 0x64]);

        // Number of minutiae
        let num_minutiae = self.inner.len() as i8;
        buf.extend_from_slice(&num_minutiae.to_be_bytes());

        for m in self.inner.iter() {
            // Encode each minutia into 6 bytes
            let encoded = encode_minutia(m);
            buf.extend_from_slice(&encoded);
        }

        let zero: u16 = 0; // same value as Python’s `zero = 0`
        buf.extend_from_slice(&zero.to_be_bytes());

        buf
    }
}
