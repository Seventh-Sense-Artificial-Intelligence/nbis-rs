use std::cmp::min;

use crate::ffi::{self, ensure_bozorth_inited};

use ffi::{
    xyt_struct,           // Bozorth input layout
    MAX_BOZORTH_MINUTIAE, // ‑‑”‑‑
};

/// Idiomatic Rust container for a minutiae set
#[derive(Clone)]
pub(crate) struct MinutiaeSet {
    pub xs: Vec<i32>,
    pub ys: Vec<i32>,
    pub theta: Vec<i32>,
}

impl MinutiaeSet {
    /// Converts to the C layout, truncating at MAX_BOZORTH_MINUTIAE
    fn to_c_struct(&self) -> xyt_struct {
        let mut xs = [0; MAX_BOZORTH_MINUTIAE];
        let mut ys = [0; MAX_BOZORTH_MINUTIAE];
        let mut theta = [0; MAX_BOZORTH_MINUTIAE];

        let n = min(self.xs.len(), MAX_BOZORTH_MINUTIAE);
        xs[..n].copy_from_slice(&self.xs[..n]);
        ys[..n].copy_from_slice(&self.ys[..n]);
        theta[..n].copy_from_slice(&self.theta[..n]);

        xyt_struct {
            nrows: n as i32,
            xcol: xs,
            ycol: ys,
            theta,
        }
    }
}

/// Safe wrapper around the C implementation.
pub(crate) fn bz_match_score(probe: &MinutiaeSet, gallery: &MinutiaeSet) -> i32 {
    ensure_bozorth_inited();
    let p_c = probe.to_c_struct();
    let g_c = gallery.to_c_struct();
    unsafe { ffi::bozorth_main(&p_c, &g_c) }
}
