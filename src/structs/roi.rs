use crate::structs::point::Point;

#[derive(Debug, Clone, uniffi::Record)]
pub struct ROI {
    /// The x-coordinate of the top-left corner of the ROI.
    pub x1: i32,
    /// The y-coordinate of the top-left corner of the ROI.
    pub y1: i32,
    /// The x-coordinate of the bottom-right corner of the ROI.
    pub x2: i32,
    /// The y-coordinate of the bottom-right corner of the ROI.
    pub y2: i32,
    /// The center point of the ROI.
    pub center: Point,
}
