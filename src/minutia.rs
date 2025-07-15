#[derive(Debug, Clone, PartialEq, uniffi::Enum)]
pub enum MinutiaKind {
    RidgeEnding = 1, // 1 = ridge ending
    Bifurcation = 2, // 0 = bifurcation
}

#[derive(uniffi::Record)]
pub struct Position {
    pub x: i32,
    pub y: i32,
}

/// A lightweight, fully-safe copy of one minutia.
/// (Add more fields if you need them.)
#[derive(Debug, Clone, uniffi::Object)]
pub struct Minutia {
    pub(crate) x: i32,
    pub(crate) y: i32,
    pub(crate) direction: i32,
    pub(crate) reliability: f64, // 0.0 to 1.0
    pub(crate) kind: MinutiaKind,
}

#[uniffi::export]
impl Minutia {
    /// Get the angle in degrees (0-360) for this minutia starting from the X axis (east).
    pub fn angle(&self) -> f64 {
        let angle = 90.0 - self.direction as f64 * 11.25;
        (angle % 360.0 + 360.0) % 360.0 // Ensures result is in [0, 360)
    }

    pub fn x(&self) -> i32 {
        self.x
    }

    pub fn y(&self) -> i32 {
        self.y
    }

    /// Get the position as a tuple (x, y).
    pub fn position(&self) -> Position {
        Position {
            x: self.x,
            y: self.y,
        }
    }

    /// Get the reliability score (0.0 to 1.0).
    pub fn reliability(&self) -> f64 {
        self.reliability
    }

    /// Get the kind of minutia (ridge ending or bifurcation).
    pub fn kind(&self) -> MinutiaKind {
        self.kind.clone()
    }
}
