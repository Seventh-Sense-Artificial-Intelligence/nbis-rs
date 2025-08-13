use std::os::raw::{c_char, c_int, c_uchar, c_uint, c_ushort};

#[repr(C)]
pub(crate) struct Nfiq2ResultsT {
    pub(crate) score: c_uint,
    pub(crate) actionable_count: c_uint,
    pub(crate) actionable_ids: *const *const c_char,
    pub(crate) actionable_values: *mut f64,
    pub(crate) feature_count: c_uint,
    pub(crate) feature_ids: *const *const c_char,
    pub(crate) feature_values: *mut f64,
}

/// Opaque C++ wrapper handle
#[repr(C)]
pub struct Nfiq2WrapperOpaque {
    _private: [u8; 0],
}

// FFI imports
extern "C" {
    pub(crate) fn nfiq2wrapper_create() -> *mut Nfiq2WrapperOpaque;
    pub(crate) fn nfiq2wrapper_destroy(ctx: *mut Nfiq2WrapperOpaque);

    pub(crate) fn nfiq2wrapper_compute(
        ctx: *mut Nfiq2WrapperOpaque,
        data: *const c_uchar,
        size: c_uint,
        cols: c_uint,
        rows: c_uint,
        ppi: c_ushort,
        out: *mut Nfiq2ResultsT,
    ) -> c_int;

    pub(crate) fn nfiq2wrapper_free_results(out: *mut Nfiq2ResultsT);
}
