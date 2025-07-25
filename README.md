## NBIS-rs

[![CI](https://github.com/Seventh-Sense-Artificial-Intelligence/nbis-rs/actions/workflows/ci.yaml/badge.svg)](https://github.com/Seventh-Sense-Artificial-Intelligence/nbis-rs/actions/workflows/ci.yaml)

This is a Rust/Python binding to the [NIST Biometric Image Software](https://www.nist.gov/services-resources/software/nist-biometric-image-software-nbis) (NBIS) library, which is used for processing biometric images, particularly in the context of fingerprint recognition.

## Features

- Bindings to NBIS functions for minutia extraction, matching, and image quality assessment
- Exports minutiae templates in ISO/IEC 19794-2:2005 format
- Matches minutiae templates against each other using the NBIS Bozorth3 algorithm

## Installation (Rust)

To use NBIS-rs, add the following to your `Cargo.toml`:

```toml
[dependencies]
nbis = "0.1.2"
```

## Usage (Rust)

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
    let minimum_minutia_quality = 0.0; // Set minimum quality to 0.0 for no filtering
    let iso_template: Vec<u8> = minutiae_1.to_iso_19794_2_2005(minimum_minutia_quality);              
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

## Installation (Python)
To install the Python bindings, you can use pip:

```bash
pip install nbis-py
```

## Usage (Python)

Here's a simple example of how to use the NBIS Python bindings:

```python
import nbis

# Read the bytes from a file
image_bytes = open("test_data/p1/p1_1.png", "rb").read()
minutiae_1 = nbis.extract_minutiae(image=image_bytes, ppi=None)
image_bytes = open("test_data/p1/p1_2.png", "rb").read()
minutiae_2 = nbis.extract_minutiae(image=image_bytes, ppi=None)
image_bytes = open("test_data/p1/p1_3.png", "rb").read()
minutiae_3 = nbis.extract_minutiae(image=image_bytes, ppi=None)

# Compare the two sets of minutiae
score = minutiae_1.compare(minutiae_2)
assert score > 50, "Expected a high similarity score between p1_1 and p1_2"
score = minutiae_1.compare(minutiae_3)
assert score > 50, "Expected a high similarity score between p1_1 and p1_3"
score = minutiae_2.compare(minutiae_3)
assert score > 50, "Expected a high similarity score between p1_2 and p1_3"

# Convert minutiae to ISO/IEC 19794-2:2005 format
iso_template = minutiae_1.to_iso_19794_2_2005()
# Load it back
minutiae_from_iso = nbis.load_iso_19794_2_2005(iso_template)
# Compare the original minutiae with the one loaded from ISO template
for a, b in zip(minutiae_from_iso.get(), minutiae_1.get()):
    assert a.x() == b.x()
    assert a.y() == b.y()
    assert a.angle() == b.angle()
    assert a.kind() == b.kind()
    # Reliability is quantized in the round-trip conversion,
    # so we allow a small margin of error.
    assert abs(a.reliability() - b.reliability()) < 0.1

# Finally we demonstrate loading from a file and comparing a negative match
minutiae_4 = nbis.extract_minutiae_from_image_file("test_data/p2/p2_1.png", ppi=None)
score = minutiae_1.compare(minutiae_4)
assert score < 50, "Expected a low similarity score between p1_1 and p2_1"
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.