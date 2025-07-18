use thiserror::Error;

#[derive(Error, Debug, uniffi::Enum)]
pub enum NbisError {
    #[error("Failed to interpret provided bytes as a PNG/JPEG image")]
    ImageLoadError,
    #[error("Failed to read file at the specified path: {0}")]
    FileReadError(String),
    #[error("quality must be between 0.0 and 1.0, got: {0}")]
    InvalidQuality(f64),
    #[error("An unexpected NBIS error occurred: {0}")]
    UnexpectedError(i64),
    #[error("Template could not be parsed: {0}")]
    InvalidTemplate(String),
    #[error("Generic error: {0}")]
    GenericError(String),
}
