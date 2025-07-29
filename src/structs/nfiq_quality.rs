/// Represents the quality of a fingerprint image as determined by NFIQ.
#[derive(Debug, Clone, PartialEq, PartialOrd, uniffi::Enum)]
pub enum NfiqQuality {
    /// Unknown quality fingerprint image. If the minutiae are loaded from a template,
    /// this means the quality is unknown.
    Unknown = 0,
    /// Excellent quality fingerprint image.
    /// This means the image is very clear and suitable for fingerprint recognition.
    Excellent = 1,
    /// Very good quality fingerprint image.
    /// This means the image is clear but may have some minor issues.
    VeryGood = 2,
    /// Good quality fingerprint image.
    /// This means the image is usable but has noticeable issues.
    Good = 3,
    /// Fair quality fingerprint image.
    /// This means the image is barely usable for fingerprint recognition.
    Fair = 4,
    /// Poor quality fingerprint image.
    /// This means the image is not usable for fingerprint recognition.
    Poor = 5,
}

impl NfiqQuality {
    pub fn from_i32(value: i32) -> Option<Self> {
        match value {
            1 => Some(NfiqQuality::Excellent),
            2 => Some(NfiqQuality::VeryGood),
            3 => Some(NfiqQuality::Good),
            4 => Some(NfiqQuality::Fair),
            5 => Some(NfiqQuality::Poor),
            _ => None,
        }
    }

    /// Encodes `NfiqQuality` as a 0–100 ISO quality byte.
    pub fn to_iso_quality(&self) -> u8 {
        match self {
            NfiqQuality::Excellent => 100,
            NfiqQuality::VeryGood => 80,
            NfiqQuality::Good => 60,
            NfiqQuality::Fair => 40,
            NfiqQuality::Poor => 20,
            NfiqQuality::Unknown => 0,
        }
    }

    /// Decodes a 0–100 ISO quality byte into `NfiqQuality`.
    pub fn from_iso_quality(value: u8) -> Self {
        match value {
            90..=100 => NfiqQuality::Excellent,
            70..=89 => NfiqQuality::VeryGood,
            50..=69 => NfiqQuality::Good,
            30..=49 => NfiqQuality::Fair,
            1..=29 => NfiqQuality::Poor,
            _ => NfiqQuality::Unknown,
        }
    }
}
