use image::{codecs::png::PngEncoder, ExtendedColorType, ImageBuffer, ImageEncoder};
use image::{Rgb, RgbImage};
use imageproc::drawing::draw_line_segment_mut;
use std::io::Cursor;

/// Encode an RGB `ImageBuffer` to PNG and return the bytes.
pub(crate) fn png_bytes_from_rgb(
    img: &ImageBuffer<Rgb<u8>, Vec<u8>>,
) -> Result<Vec<u8>, image::ImageError> {
    let mut buf = Vec::with_capacity(img.width() as usize * img.height() as usize * 3);
    {
        // `Cursor` lets the encoder write into the Vec
        let mut writer = Cursor::new(&mut buf);
        let encoder = PngEncoder::new(&mut writer);
        encoder.write_image(
            img.as_raw(), // &[u8]
            img.width(),
            img.height(),
            ExtendedColorType::Rgb8, // 8-bit RGB
        )?;
    } // writer drops here; `buf` now contains PNG data

    Ok(buf)
}

/// Simulate a thick antialiased line by drawing multiple offset lines.
fn draw_thick_line_segment(
    img: &mut RgbImage,
    (x0, y0): (f32, f32),
    (x1, y1): (f32, f32),
    thickness: f32,
    color: Rgb<u8>,
) {
    let steps = thickness.ceil() as i32;
    let half = steps / 2;

    let dx = x1 - x0;
    let dy = y1 - y0;
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return;
    }

    // Perpendicular vector
    let (nx, ny) = (-dy / len, dx / len);

    for i in -half..=half {
        let offset = i as f32 * 0.5; // subpixel spacing
        let ox = nx * offset;
        let oy = ny * offset;

        draw_line_segment_mut(img, (x0 + ox, y0 + oy), (x1 + ox, y1 + oy), color);
    }
}

/// Draws a directional arrow from a center point, with arrowhead.
/// Angle is CCW from +X (mathematical convention), in degrees.
pub(crate) fn draw_arrow_with_head(
    img: &mut RgbImage,
    center: (f32, f32),
    angle_deg: f32,
    length: f32,
    head_size: f32,
    color: Rgb<u8>,
) {
    let theta = angle_deg.to_radians();

    // Compute arrow tip
    let x_tip = center.0 + length * theta.cos();
    let y_tip = center.1 - length * theta.sin(); // Y axis inverted in image

    // Draw the main shaft
    draw_thick_line_segment(img, center, (x_tip, y_tip), 2.0, color);

    // Arrowhead legs: symmetric around tip, ±160° from shaft
    let angle1 = theta + 150.0_f32.to_radians(); // adjust angle for arrowhead leg 1
    let angle2 = theta - 150.0_f32.to_radians(); // adjust angle for arrowhead leg 2

    let x1 = x_tip + head_size * angle1.cos();
    let y1 = y_tip - head_size * angle1.sin();

    let x2 = x_tip + head_size * angle2.cos();
    let y2 = y_tip - head_size * angle2.sin();

    draw_thick_line_segment(img, (x_tip, y_tip), (x1, y1), 2.0, color);
    draw_thick_line_segment(img, (x_tip, y_tip), (x2, y2), 2.0, color);
}
