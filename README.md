## NBIS-rs

This is a Rust binding to the NIST Biometric Image Software (NBIS) library, which is used for processing biometric images, particularly in the context of fingerprint recognition.

## Features

- Bindings to NBIS functions for minutia extraction, and matching
- Exports minutiae templates in ISO/IEC 19794-2:2005 format
- Matches minutiae templates against each other using the NBIS Bozorth3 algorithm

## Installation

To use NBIS-rs, add the following to your `Cargo.toml`:

```toml
[dependencies]
nbis = "0.1"
```

## Usage

Here's a simple example of how to use NBIS-rs in your project:

```rust
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use nbis;
    use nbis::Minutiae;

    // Read the bytes from a file (you could also use nbis::extract_minutiae_from_image_file)
    // but here we just load the image bytes as image paths on mobile platforms can be tricky.
    let image_bytes = std::fs::read("test_data/p1/p1_1.png")?;
    let minutiae_1 = nbis::extract_minutiae(&image_bytes, None)?;

    let image_bytes = std::fs::read("test_data/p1/p1_2.png")?;
    let minutiae_2 = nbis::extract_minutiae(&image_bytes, None)?;

    let image_bytes = std::fs::read("test_data/p1/p1_3.png")?;
    let minutiae_3 = nbis::extract_minutiae(&image_bytes, None)?;

    // Compare the two sets of minutiae
    let score = minutiae_1.compare(&minutiae_2);
    assert!(score > 50, "Expected a high similarity score between p1_1 and p1_2");
    let score = minutiae_1.compare(&minutiae_3);
    assert!(score > 50, "Expected a high similarity score between p1_1 and p1_3");
    let score = minutiae_2.compare(&minutiae_3);
    assert!(score > 50, "Expected a high similarity score between p1_2 and p1_3");

    // Next we will demonstrate conversion to ISO/IEC 19794-2:2005 format
    // and back to a `Minutiae` object.
    // First, convert the minutiae to ISO template bytes
    let iso_template: Vec<u8> = minutiae_1.to_iso_19794_2_2005();
    // And load it back
    let minutiae_from_iso = nbis::load_iso_19794_2_2005(&iso_template)?;
    // Compare the original minutiae with the one loaded from ISO template
    for (a, b) in minutiae_from_iso.get().iter().zip(minutiae_1.get().iter()) {
        assert_eq!(a.x(), b.x());
        assert_eq!(a.y(), b.y());
        assert_eq!(a.angle(), b.angle());
        assert_eq!(a.kind(), b.kind());
        // Reliability is quantized in the round-trip conversion,
        // so we allow a small margin of error.
        assert!((a.reliability() - b.reliability()).abs() < 1e-1);
    }

    // Finally we demonstrate loading from a file and comparing a negative match
    let minutiae_4 = nbis::extract_minutiae_from_image_file("test_data/p2/p2_1.png", None)?;
    let score = minutiae_1.compare(&minutiae_4);
    assert!(score < 50, "Expected a low similarity score between p1_1 and p2_1");

    Ok(())
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.