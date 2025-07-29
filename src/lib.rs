#![doc = include_str!("../README.md")]
// src/lib.rs
uniffi::setup_scaffolding!();

mod bozorth;
mod consts;
mod encoding;
mod errors;
mod extractor;
pub(crate) mod ffi;
mod imutils;
mod minutia;
mod minutiae;
mod sivv;
mod structs;

pub use structs::{NbisExtractorSettings, NfiqQuality, Point, ROI};

pub use errors::NbisError;
pub use extractor::new_nbis_extractor;
pub use extractor::NbisExtractor;
pub use minutia::{Minutia, MinutiaKind, Position};
pub use minutiae::Minutiae;
