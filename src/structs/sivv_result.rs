/// Represents the result of the SIVV computation.
pub(crate) struct SIVVResult {
    /// Index of the largest peak-valley pair (1-based)
    pub(crate) largest_pvp_index: i32,

    /// Total number of detected peak-valley pairs
    pub(crate) total_pvps: i32,

    /// Power difference between the peak and valley
    pub(crate) power_diff: f64,

    /// Frequency difference between the peak and valley
    pub(crate) freq_diff: f64,

    /// Slope between valley and peak (dy / dx)
    pub(crate) slope: f64,

    /// Frequency of the midpoint between valley and peak
    pub(crate) center_frequency: f64,

    /// Absolute frequency of the peak (undocumented in comments)
    pub(crate) peak_frequency: f64,
}
