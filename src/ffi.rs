// src/ffi.rs
#![allow(
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    clippy::upper_case_acronyms
)]
use libc::FILE;
use std::os::raw::{c_char, c_double, c_int, c_uchar, c_void};
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
                                             //pub const TOO_FEW_MINUTIAE: c_int = 2; // <-- match nfiq.h
                                             //pub const EMPTY_IMG: c_int = 1; // <-- match nfiq.h
pub const MIN_MINUTIAE: usize = 5; // <-- match nfiq.h
pub const NFIQ_VCTRLEN: usize = 11; // <-- match nfiq.h
pub const NFIQ_NUM_CLASSES: usize = 5; // <-- match nfiq.h
pub const NFIQ_ML_WEIGHTS_LEN: usize = 379;

extern "C" {
    pub static dflt_znorm_means: [f32; NFIQ_VCTRLEN]; // 11 = NFIQ_VCTRLEN
    pub static dflt_znorm_stds: [f32; NFIQ_VCTRLEN]; // 11 = NFIQ_VCTRLEN
    pub static dflt_nInps: ::std::os::raw::c_int; // 11 = NFIQ_VCTRLEN
    pub static dflt_nHids: ::std::os::raw::c_int;
    pub static dflt_nOuts: ::std::os::raw::c_int;
    pub static dflt_acfunc_hids: ::std::os::raw::c_char;
    pub static dflt_acfunc_outs: ::std::os::raw::c_char;
    pub static dflt_wts: [f32; NFIQ_ML_WEIGHTS_LEN]; // 379 = NFIQ_ML_WEIGHTS_LEN
}

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
    /// Run the NFIQ MLP model on a feature vector
    ///
    /// # Arguments
    /// - `ninps`: number of input nodes
    /// - `nhids`: number of hidden nodes
    /// - `nouts`: number of output nodes
    /// - `acfunc_hids_code`: activation function code for hidden layer (e.g., 'L', 'S', 'N')
    /// - `acfunc_outs_code`: activation function code for output layer
    /// - `w`: pointer to MLP weights
    /// - `featvec`: pointer to input feature vector
    /// - `outacs`: pointer to output activation vector (length `nouts`)
    /// - `hypclass`: pointer to integer for output class (0 to nouts - 1)
    /// - `confidence`: pointer to float for confidence score (max activation)
    ///
    /// # Returns
    /// - 0 on success
    pub fn runmlp2(
        ninps: ::std::os::raw::c_int,
        nhids: ::std::os::raw::c_int,
        nouts: ::std::os::raw::c_int,
        acfunc_hids_code: ::std::os::raw::c_char,
        acfunc_outs_code: ::std::os::raw::c_char,
        w: *mut f32,
        featvec: *mut f32,
        outacs: *mut f32,
        hypclass: *mut ::std::os::raw::c_int,
        confidence: *mut f32,
    ) -> ::std::os::raw::c_int;
}

extern "C" {
    /// Z-normalize an NFIQ feature vector in-place
    ///
    /// # Arguments
    /// - `featvctr`: pointer to the feature vector (modified in-place)
    /// - `znorm_means`: pointer to the means for each coefficient
    /// - `znorm_stds`: pointer to the stddevs for each coefficient
    /// - `vctrlen`: length of the vectors
    pub fn znorm_fniq_featvctr(
        featvctr: *mut f32,
        znorm_means: *const f32,
        znorm_stds: *const f32,
        vctrlen: ::std::os::raw::c_int,
    );
}

extern "C" {
    /// Computes NFIQ feature vector from Mindtct output
    ///
    /// # Arguments
    /// - `featvctr`: pointer to output buffer (length = `vctrlen`)
    /// - `vctrlen`: length of the `featvctr` array
    /// - `minutiae`: pointer to C `MINUTIAE` struct
    /// - `quality_map`: pointer to quality map (int*)
    /// - `map_w`, `map_h`: dimensions of the quality map   
    /// - `optflag`: pointer to optional flags
    ///
    /// # Returns
    /// - 0 on success
    /// - EMPTY_IMG or error code otherwise
    pub fn comp_nfiq_featvctr(
        featvctr: *mut f32,
        vctrlen: ::std::os::raw::c_int,
        minutiae: *const MINUTIAE,
        quality_map: *mut ::std::os::raw::c_int,
        map_w: ::std::os::raw::c_int,
        map_h: ::std::os::raw::c_int,
        optflag: *mut ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}

extern "C" {
    pub(crate) fn sivv_ffi_from_bytes(data: *const u8, width: c_int, height: c_int) -> *mut c_char;
    pub(crate) fn sivv_ffi_free_bytes(ptr: *mut c_char);
}

#[repr(C)]
pub struct CPoint2i {
    pub x: c_int,
    pub y: c_int,
}

extern "C" {
    pub(crate) fn find_fingerprint_center_morph_c(
        data: *const u8,
        width: c_int,
        height: c_int,
        xbound_min: *mut c_int,
        xbound_max: *mut c_int,
        ybound_min: *mut c_int,
        ybound_max: *mut c_int,
    ) -> CPoint2i;
}
