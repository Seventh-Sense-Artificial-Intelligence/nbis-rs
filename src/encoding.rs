use std::cmp::Ordering;

use crate::{
    api::NfiqResult, bozorth::MinutiaeSet, consts::NUM_DIRECTIONS, ffi::MAX_BOZORTH_MINUTIAE,
    Minutia, MinutiaKind, Minutiae, NbisError, NfiqQuality,
};

/// Quantise an angle (degrees) into the 8-bit ISO/IEC 19794-2 orientation unit.
///
/// The ISO unit represents 360 ° / 256 ≈ 1.40625 ° per code value.
/// Any input (positive or negative) is first wrapped into the range [0, 360).
///
/// ```text
/// angle_deg  →  iso_unit (0‥=255)
/// 0.0        →  0
/// 1.40625    →  1
/// 360.0      →  0
/// -1.40625   →  255
/// ```
#[inline]
pub(crate) fn encode_iso_angle(angle_deg: f64) -> u8 {
    // Unquantize the NIST angle
    let angle_deg = 90.0 - angle_deg * 11.25;

    // Quantize to ISO/IEC 19794-2 orientation unit
    const ISO_STEP: f64 = 360.0 / 256.0; // 1.40625 °
    let norm = angle_deg.rem_euclid(360.0); // wrap into [0, 360)
    let quantised = (norm / ISO_STEP).round() as i32 & 0xFF;
    quantised as u8
}

/// Encode a single minutia into its 6-byte ISO/IEC 19794-2 representation.
///
/// * `x`, `y`          – image coordinates (0-based, ≤ 65 535)
/// * `angle_deg`       – orientation in ISO units (0-255, already quantised)
/// * `quality`         – 0-100 quality score (or vendor-defined)
/// * `min_type`        – 0 = ridge ending, 1 = bifurcation, 2 = other (uses top 2 bits)
///
/// Returns `[u8; 6]` matching the original Python `byte_list`.
#[inline]
pub(crate) fn encode_minutia(m: &Minutia) -> [u8; 6] {
    let min_type = if m.kind == MinutiaKind::Bifurcation {
        2
    } else {
        1
    }; // 0x02 = bifurcation, 0x01 = ridge ending
    let angle = encode_iso_angle(m.direction as f64);
    let quality = (m.reliability * 63.0).min(63.0) as u8; // 0-63 range
    let x = m.x;
    let y = m.y;
    let mut bytes = [0u8; 6];

    bytes[1] = (x & 0x00FF) as u8; // low-byte of X
    bytes[0] = ((x >> 8) as u8) + (min_type << 6); // high-byte of X + type*64
    bytes[2] = (y >> 8) as u8; // high-byte of Y
    bytes[3] = (y & 0x00FF) as u8; // low-byte of Y
    bytes[4] = angle; // orientation (ISO unit)
    bytes[5] = quality; // quality score

    bytes
}

fn nist_xyt(minutiae: &Minutiae, minutia: &Minutia) -> (i32, i32, i32) {
    // 1. Bottom‑left origin for Y.
    let x = minutia.x;
    let y = minutiae.img_h as i32 - minutia.y;

    // 2. Direction → degrees.
    let degrees_per_unit = 180.0 / NUM_DIRECTIONS; // 11.25 for 16 dirs
    let t = (270 - ((minutia.direction as f64 * degrees_per_unit).round() as i32)) % 360;
    let t = if t < 0 { t + 360 } else { t };

    (x, y, t)
}
/// Convert this result into a Bozorth‑ready [`MinutiaeSet`].
pub(crate) fn to_nist_xyt_set(minutiae: &Minutiae) -> MinutiaeSet {
    // 1. Copy, sort by reliability (desc), truncate
    let mut top = minutiae.inner.clone();
    top.sort_by(|a, b| b.reliability.partial_cmp(&a.reliability).unwrap());
    top.truncate(MAX_BOZORTH_MINUTIAE); // hard Bozorth limit

    // 2. Convert to NIST XYT
    let mut xs = Vec::with_capacity(top.len());
    let mut ys = Vec::with_capacity(top.len());
    let mut ts = Vec::with_capacity(top.len());

    for m in &top {
        let (ox, oy, ot) = nist_xyt(minutiae, m);
        xs.push(ox);
        ys.push(oy);
        ts.push(ot);
    }

    MinutiaeSet { xs, ys, theta: ts }
}

/// De-quantise an ISO/IEC 19794-2 orientation code back to the
/// *NIST orientation unit* (0‥15).
///
/// ```text
/// iso_code  →  nist_unit
/// 0         →  8
/// 64        →  0
/// 255       →  9           (because 255·1.40625° ≈ 358.59°, etc.)
/// ```
///
/// The result is returned as `u8`; change the cast if you prefer `f64`.
#[inline]
pub(crate) fn decode_iso_angle(iso_code: u8) -> i32 {
    const ISO_STEP: f64 = 360.0 / 256.0; // 1.40625°
    const NIST_STEP: f64 = 11.25; // one NIST unit in degrees

    // 1. Expand the 8-bit code to a real angle (centre of its bin).
    let iso_deg = iso_code as f64 * ISO_STEP;

    // 2. Undo the “90° − θ” swirl applied in the forward path.
    //    rem_euclid keeps the value positive before the division.
    let nist_unit = (90.0 - iso_deg).rem_euclid(360.0) / NIST_STEP;

    // 3. Re-quantise and wrap into 0-31 (5 bits)
    (nist_unit.round() as u32 & 0x1F) as i32
}

/// Decode the 6-byte minutia record used by `encode_minutia`.
#[inline]
pub(crate) fn decode_minutia(bytes: &[u8; 6]) -> Minutia {
    let x = (u16::from(bytes[0] & 0x3F) << 8) | u16::from(bytes[1]);
    let y = (u16::from(bytes[2]) << 8) | u16::from(bytes[3]);
    let direction = bytes[4];
    let reliability = f64::from(bytes[5]) / 63.0; // 0-63 → 0.0-1.0
    let kind = if bytes[0] & 0xC0 == 0x80 {
        MinutiaKind::Bifurcation
    } else {
        MinutiaKind::RidgeEnding
    };

    let direction = decode_iso_angle(direction);

    Minutia {
        x: x as i32,
        y: y as i32,
        direction,
        reliability,
        kind,
    }
}

/// Convert this set of minutiae into an ISO/IEC 19794-2:2005 template.
///
/// # Arguments
/// * `minutiae_obj` — the `Minutiae` object to convert.
/// * `min_quality` — minimum quality threshold for minutiae to be included (0.0 to 1.0).
///
/// Returns a `Vec<u8>` containing the ISO template bytes.
pub fn to_iso_19794_2_2005(minutiae_obj: &Minutiae, min_quality: f64) -> Vec<u8> {
    // The maximum number of minutiae is 255, so we can use u8 for the count.
    // therefore, first filter the top 255 minutiae by quality.
    let mut minutiae = minutiae_obj.inner.clone();

    // If a minimum quality is specified, filter out minutiae below that quality.
    if min_quality > 0.0 {
        minutiae.retain(|m| m.reliability >= min_quality);
    }

    if minutiae.len() > 255 {
        minutiae.sort_by(|a, b| {
            b.reliability
                .partial_cmp(&a.reliability)
                .unwrap_or(Ordering::Equal)
        });
        minutiae.truncate(255);
    }

    const ISO_HEADER_LENGTH: usize = 26;
    let total_bytes = minutiae.len() * 6 + ISO_HEADER_LENGTH;

    let total_bytes_u32 =
        u32::try_from(total_bytes).expect("template larger than 4 294 967 295 bytes");

    let width = minutiae_obj.img_w as u16;
    let height = minutiae_obj.img_h as u16;

    let mut buf = Vec::with_capacity(ISO_HEADER_LENGTH + minutiae.len() * 6);

    // "FMR\0 20\0"  – 8 bytes
    buf.extend_from_slice(b"FMR\0 20\0");

    // 4-byte total length (big-endian)
    buf.extend_from_slice(&total_bytes_u32.to_be_bytes());

    // Two 0x00 padding bytes
    buf.extend_from_slice(&[0x00, 0x00]);

    // Width & height (big-endian, 2 bytes each)
    buf.extend_from_slice(&width.to_be_bytes());
    buf.extend_from_slice(&height.to_be_bytes());

    let finger_position: u8 = 0; // 0 = unknown (or 1–10 = right fingers, 11–20 = left fingers)
    let view_number: u8 = 0; // Usually 0 for single view
    let impression_type: u8 = 0; // 0 = live-scan plain
    let view_and_impression = (view_number << 4) | (impression_type & 0x0F);

    // Convert NFIQ to ISO quality (0-100)
    let finger_quality: u8 = minutiae_obj.nfiq.nfiq.to_iso_quality();

    // Reserved 4 bytes
    let reserved = [0x00, 0x00, 0x00, 0x00];

    // Metadata bytes
    buf.extend_from_slice(&[finger_position]);
    buf.extend_from_slice(&[view_and_impression]);
    buf.extend_from_slice(&[finger_quality]);
    buf.extend_from_slice(&reserved);

    // Number of minutiae
    let num_minutiae = minutiae.len() as u8;
    buf.extend_from_slice(&num_minutiae.to_be_bytes());

    for m in minutiae.iter() {
        // Encode each minutia into 6 bytes
        let encoded = encode_minutia(m);
        buf.extend_from_slice(&encoded);
    }

    assert_eq!(buf.len(), total_bytes);

    buf
}

/// Loads an ISO/IEC 19794-2:2005 fingerprint template from bytes.
///
/// # Arguments
/// * `template_bytes` — the raw bytes of the ISO template.
///
/// Returns a [`Minutiae`] object containing the decoded minutiae.
///
/// If the template is invalid or cannot be parsed, returns an [`NbisError`].
pub fn load_iso_19794_2_2005(template_bytes: &[u8]) -> Result<Minutiae, NbisError> {
    if template_bytes.len() < 28 {
        return Err(NbisError::InvalidTemplate(
            "ISO template too short".to_string(),
        ));
    }

    // Check the header
    if &template_bytes[0..8] != b"FMR\0 20\0" {
        return Err(NbisError::InvalidTemplate("Invalid ISO header".to_string()));
    }

    let total_length = u32::from_be_bytes([
        template_bytes[8],
        template_bytes[9],
        template_bytes[10],
        template_bytes[11],
    ]) as usize;
    if total_length != template_bytes.len() {
        return Err(NbisError::InvalidTemplate(
            "Total length mismatch".to_string(),
        ));
    }

    let width = u16::from_be_bytes([template_bytes[14], template_bytes[15]]);
    let height = u16::from_be_bytes([template_bytes[16], template_bytes[17]]);

    let finger_quality = NfiqQuality::from_iso_quality(template_bytes[20]);

    let num_minutiae = template_bytes[25] as usize;
    let minutiae_start = 26;

    let mut minutiae = Vec::with_capacity(num_minutiae);
    for i in 0..num_minutiae {
        let start = minutiae_start + 6 * i;
        let end = start + 6;
        if end > template_bytes.len() {
            return Err(NbisError::InvalidTemplate(
                "Minutia data overflow".to_string(),
            ));
        }
        let m_bytes: [u8; 6] = template_bytes[start..end].try_into().unwrap();
        minutiae.push(decode_minutia(&m_bytes));
    }

    Ok(Minutiae::new(
        minutiae,
        width as u32,
        height as u32,
        NfiqResult {
            nfiq: finger_quality,
            confidence: 1.0,
        },
    ))
}
