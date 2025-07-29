use crate::structs::nfiq_quality::NfiqQuality;

/// Represents the result of the NFIQ computation.
#[derive(Debug, Clone, uniffi::Record)]
pub struct NfiqResult {
    /// The NFIQ quality score.
    /// 1 = Excellent, 2 = Very Good, 3 = Good,
    /// 4 = Fair, 5 = Poor.
    /// See [`NfiqQuality`] for more details.
    pub nfiq: NfiqQuality,
    /// The confidence level of the NFIQ score.
    /// A value between 0.0 and 1.0, where 1.0 means very confident.
    /// This is a floating point value.
    pub confidence: f32,
}
