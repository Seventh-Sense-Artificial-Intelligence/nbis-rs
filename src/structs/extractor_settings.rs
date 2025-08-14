#[derive(Debug, Clone, uniffi::Record)]
pub struct NbisExtractorSettings {
    /// The minimum quality of the minutiae to be extracted.
    pub min_quality: f64,

    /// Whether to extract the center point / ROI of the fingerprint.
    pub get_center: bool,

    /// Whether to use SIVV to check for a valid fingerprint.
    pub check_fingerprint: bool,

    /// Whether to compute the NFIQ2 quality of the fingerprint
    pub compute_nfiq2: bool,

    /// The PPI (pixels per inch) of the image. Default is 500.
    pub ppi: Option<f64>,
}

// Implementation of default settings for ExtractorSettings
impl Default for NbisExtractorSettings {
    fn default() -> Self {
        NbisExtractorSettings {
            min_quality: 0.0,
            get_center: false,
            check_fingerprint: false,
            compute_nfiq2: true,
            ppi: None,
        }
    }
}
