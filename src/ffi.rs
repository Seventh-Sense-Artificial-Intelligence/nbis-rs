// src/ffi.rs
#![allow(
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    clippy::upper_case_acronyms
)]
use libc::FILE;
use std::os::raw::{c_double, c_int, c_uchar, c_void};
use std::sync::Once;

/* -----------------------------------------------------------------------
libc’s stderr — linker name depends on platform
-------------------------------------------------------------------- */
extern "C" {
    #[cfg_attr(any(target_os = "macos", target_os = "ios"), link_name = "__stderrp")]
    #[cfg_attr(not(any(target_os = "macos", target_os = "ios")), link_name = "stderr")]
    static mut libc_stderr: *mut FILE;
}

extern "C" {
    // declared in bozorth3.c
    static mut errorfp: *mut libc::FILE;
}

static INIT_BOZORTH: Once = Once::new();

#[inline]
pub(crate) fn ensure_bozorth_inited() {
    // SAFETY: we’re the only code touching `errorfp`, and we only write stderr
    unsafe {
        INIT_BOZORTH.call_once(|| {
            if errorfp.is_null() {
                errorfp = libc_stderr; // ← core fix
            }
        });
    }
}

pub const MAX_BOZORTH_MINUTIAE: usize = 200; // <-- match bozorth.h
pub const TOO_FEW_MINUTIAE: c_int = 2; // <-- match nfiq.h
pub const EMPTY_IMG: c_int = 1; // <-- match nfiq.h

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub(crate) struct xyt_struct {
    pub nrows: c_int,
    pub xcol: [c_int; MAX_BOZORTH_MINUTIAE],
    pub ycol: [c_int; MAX_BOZORTH_MINUTIAE],
    pub theta: [c_int; MAX_BOZORTH_MINUTIAE],
}

extern "C" {
    pub(crate) fn bozorth_main(probe: *const xyt_struct, gallery: *const xyt_struct) -> c_int;
}

////////////////////////////////////////////////////////////////////////////////
//  Structs copied verbatim from <lfs.h>
////////////////////////////////////////////////////////////////////////////////

#[repr(C)]
#[derive(Debug)]
pub(crate) struct MINUTIA {
    pub x: c_int,
    pub y: c_int,
    pub ex: c_int,
    pub ey: c_int,
    pub direction: c_int,
    pub reliability: c_double,
    pub r#type: c_int, // `type` is a Rust keyword – use a raw identifier
    pub appearing: c_int,
    pub feature_id: c_int,
    pub nbrs: *mut c_int,
    pub ridge_counts: *mut c_int,
    pub num_nbrs: c_int,
}

#[repr(C)]
#[derive(Debug)]
pub(crate) struct MINUTIAE {
    pub alloc: c_int,
    pub num: c_int,
    pub list: *mut *mut MINUTIA, // C: MINUTIA **list;
}

#[repr(C)]
#[derive(Debug)]
pub(crate) struct LFSPARMS {
    /* ------------ Image controls ------------------ */
    pub pad_value: c_int,
    pub join_line_radius: c_int,

    /* ------------ Map controls -------------------- */
    pub blocksize: c_int,
    pub windowsize: c_int,
    pub windowoffset: c_int,
    pub num_directions: c_int,
    pub start_dir_angle: c_double,
    pub rmv_valid_nbr_min: c_int,
    pub dir_strength_min: c_double,
    pub dir_distance_max: c_int,
    pub smth_valid_nbr_min: c_int,
    pub vort_valid_nbr_min: c_int,
    pub highcurv_vorticity_min: c_int,
    pub highcurv_curvature_min: c_int,
    pub min_interpolate_nbrs: c_int,
    pub percentile_min_max: c_int,
    pub min_contrast_delta: c_int,

    /* ------------ DFT controls -------------------- */
    pub num_dft_waves: c_int,
    pub powmax_min: c_double,
    pub pownorm_min: c_double,
    pub powmax_max: c_double,
    pub fork_interval: c_int,
    pub fork_pct_powmax: c_double,
    pub fork_pct_pownorm: c_double,

    /* ------------ Binarisation controls ----------- */
    pub dirbin_grid_w: c_int,
    pub dirbin_grid_h: c_int,
    pub isobin_grid_dim: c_int,
    pub num_fill_holes: c_int,

    /* ------------ Minutiae detection -------------- */
    pub max_minutia_delta: c_int,
    pub max_high_curve_theta: c_double,
    pub high_curve_half_contour: c_int,
    pub min_loop_len: c_int,
    pub min_loop_aspect_dist: c_double,
    pub min_loop_aspect_ratio: c_double,

    /* ------------ Linking ------------------------- */
    pub link_table_dim: c_int,
    pub max_link_dist: c_int,
    pub min_theta_dist: c_int,
    pub maxtrans: c_int,
    pub score_theta_norm: c_double,
    pub score_dist_norm: c_double,
    pub score_dist_weight: c_double,
    pub score_numerator: c_double,

    /* ------------ False-minutiae removal ---------- */
    pub max_rmtest_dist: c_int,
    pub max_hook_len: c_int,
    pub max_half_loop: c_int,
    pub trans_dir_pix: c_int,
    pub small_loop_len: c_int,
    pub side_half_contour: c_int,
    pub inv_block_margin: c_int,
    pub rm_valid_nbr_min: c_int,
    pub max_overlap_dist: c_int,
    pub max_overlap_join_dist: c_int,
    pub malformation_steps_1: c_int,
    pub malformation_steps_2: c_int,
    pub min_malformation_ratio: c_double,
    pub max_malformation_dist: c_int,
    pub pores_trans_r: c_int,
    pub pores_perp_steps: c_int,
    pub pores_steps_fwd: c_int,
    pub pores_steps_bwd: c_int,
    pub pores_min_dist2: c_double,
    pub pores_max_ratio: c_double,

    /* ------------ Ridge counting ------------------ */
    pub max_nbrs: c_int,
    pub max_ridge_steps: c_int,
}

////////////////////////////////////////////////////////////////////////////////
//  FFI function prototypes (unchanged)
////////////////////////////////////////////////////////////////////////////////

extern "C" {
    // --  Core entry point ---------------------------------------------------
    pub(crate) fn get_minutiae(
        ominutiae: *mut *mut MINUTIAE,
        oquality_map: *mut *mut c_int,
        odirection_map: *mut *mut c_int,
        olow_contrast_map: *mut *mut c_int,
        olow_flow_map: *mut *mut c_int,
        ohigh_curve_map: *mut *mut c_int,
        omap_w: *mut c_int,
        omap_h: *mut c_int,
        obdata: *mut *mut c_uchar,
        obw: *mut c_int,
        obh: *mut c_int,
        obd: *mut c_int,
        idata: *mut c_uchar,
        iw: c_int,
        ih: c_int,
        id: c_int,
        ppmm: c_double,
        lfsparms: *const LFSPARMS,
    ) -> c_int;

    // --  Library-provided destructors --------------------------------------
    pub(crate) fn free_minutiae(ptr: *mut MINUTIAE);
    pub(crate) fn free(ptr: *mut c_void); // libc::free
}

extern "C" {
    /// Computes NFIQ value and confidence from an 8-bit grayscale fingerprint image.
    /// - `idata`: pointer to grayscale image (8-bit)
    /// - `iw`, `ih`, `id`: width, height, depth (depth = 8)
    /// - `ippi`: resolution in pixels per inch (-1 for default 500)
    /// - `onfiq`: output pointer to NFIQ integer
    /// - `oconf`: output pointer to confidence (max MLP activation)
    /// - `optflag`: optional flags (can be null or &mut 0)
    ///
    /// Returns 0 on success, negative on error, or special codes for empty/low-quality images.
    pub(crate) fn comp_nfiq(
        onfiq: *mut c_int,
        oconf: *mut f32,
        idata: *mut c_uchar,
        iw: c_int,
        ih: c_int,
        id: c_int,
        ippi: c_int,
        optflag: *mut c_int,
    ) -> c_int;
}
