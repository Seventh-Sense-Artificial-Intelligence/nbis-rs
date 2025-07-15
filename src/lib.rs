#![doc = include_str!("../README.md")]
// src/lib.rs
uniffi::setup_scaffolding!();

mod api;
mod bozorth;
mod consts;
mod encoding;
mod errors;
pub(crate) mod ffi;
mod imutils;
mod minutia;
mod minutiae;

pub use api::{extract_minutiae, extract_minutiae_from_image_file, load_iso_19794_2_2005};
pub use errors::NbisError;
pub use minutia::{Minutia, MinutiaKind, Position};
pub use minutiae::Minutiae;
