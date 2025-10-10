#[derive(Debug, Clone, uniffi::Record)]
pub struct Point {
    /// The x-coordinate of the point.
    pub x: i32,
    /// The y-coordinate of the point.
    pub y: i32,
}
