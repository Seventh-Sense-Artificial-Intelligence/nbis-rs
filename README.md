## NBIS-rs

This is a Rust/Python binding to the [NIST Biometric Image Software](https://www.nist.gov/services-resources/software/nist-biometric-image-software-nbis) (NBIS) library, which is used for processing biometric images, particularly in the context of fingerprint recognition.

For convenience, this library also binds to the [NIST Fingerprint Image Quality](https://www.nist.gov/services-resources/software/nfiq-2) (NFIQ) version 2. 

## Features

- Bindings to NBIS functions for minutia extraction, matching
- Exports minutiae templates in ISO/IEC 19794-2:2005 format
- Matches minutiae templates against each other using the NBIS Bozorth3 algorithm
- Provides support for NFIQ2 quality assessment

## Installation (Rust)

To use NBIS-rs, add the following to your `Cargo.toml`:

```toml
[dependencies]
nbis-rs = { git = "https://github.com/Seventh-Sense-Artificial-Intelligence/nbis-rs", branch = "main", version = "0.1.2" }
```
Or you can run the following command on the terminal of your new rust project:

cargo add nbis-rs --git https://github.com/Seventh-Sense-Artificial-Intelligence/nbis-rs --branch main

Running the above command will add the above dependency in your Cargo.toml

Now you can use the nbis-rs rust library in your project as mentioned in next section.

## Usage (Rust)

Here's a simple example of how to use NBIS-rs in your project:

```rust
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use nbis;
    use nbis::Minutiae;
    use nbis::NbisExtractorSettings;
    // Configuration for the NbisExtractor
    let settings = NbisExtractorSettings {
        // No filtering on minutiae quality (all minutiae will be included)
        min_quality: 0.0,
        // Do not compute ROI or center to save computing resources
        get_center: false,
        // Do not check if the image is a fingerprint using SIVV
        check_fingerprint: false,
        // compute the NFIQ score
        compute_nfiq2: true,
        // No specific PPI, use the default
        ppi: None,
    };

    let extractor = nbis::NbisExtractor::new(settings)?;

    // Read the bytes from a file (you could also use nbis::extract_minutiae_from_image_file)
    // but here we just load the image bytes as image paths on mobile platforms can be tricky.
    let image_bytes = std::fs::read("test_data/p1/p1_1.png")?;

    let minutiae_1 = extractor.extract_minutiae(&image_bytes)?;

    let image_bytes = std::fs::read("test_data/p1/p1_2.png")?;
    let minutiae_2 = extractor.extract_minutiae(&image_bytes)?;

    let image_bytes = std::fs::read("test_data/p1/p1_3.png")?;
    let minutiae_3 = extractor.extract_minutiae(&image_bytes)?;

    // Compare the two sets of minutiae
    let score = minutiae_1.compare(&minutiae_2);
    assert!(score > 35, "Expected a high similarity score between p1_1 and p1_2");
    let score = minutiae_1.compare(&minutiae_3);
    assert!(score > 35, "Expected a high similarity score between p1_1 and p1_3");
    let score = minutiae_2.compare(&minutiae_3);
    assert!(score > 35, "Expected a high similarity score between p1_2 and p1_3");

    // Next we will demonstrate conversion to ISO/IEC 19794-2:2005 format
    // and back to a `Minutiae` object.
    // First, convert the minutiae to ISO template bytes
    let iso_template: Vec<u8> = minutiae_1.to_iso_19794_2_2005();              
    // And load it back
    let minutiae_from_iso = extractor.load_iso_19794_2_2005(&iso_template)?;
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
    let minutiae_4 = extractor.extract_minutiae_from_image_file("test_data/p2/p2_1.png")?;
    let score = minutiae_1.compare(&minutiae_4);
    assert!(score < 35, "Expected a low similarity score between p1_1 and p2_1");

    // We can access the NFIQ2 quality via:
    let nfiq2_quality = minutiae_1.quality();
    assert!(nfiq2_quality.score > 50, "Expected a positive NFIQ2 quality score");

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
from nbis import NbisExtractor, NbisExtractorSettings

 #Configuration for the NbisExtractor
settings = NbisExtractorSettings(
    # Do not filter on minutiae quality (get all minutiae)
    min_quality=0.0,
    # Do not get the fingerprint center or ROI
    get_center=False,
    # Do not use SIVV to check if the image is a fingerprint
    check_fingerprint=False,
    # Compute the NFIQ2 quality score
    compute_nfiq2=True,
    # No specific PPI, use the default
    ppi=None,
)

extractor = nbis.new_nbis_extractor(settings)

# Read the bytes from a file
image_bytes = open("test_data/p1/p1_1.png", "rb").read()
minutiae_1 = extractor.extract_minutiae(image_bytes)
image_bytes = open("test_data/p1/p1_2.png", "rb").read()
minutiae_2 = extractor.extract_minutiae(image_bytes)
image_bytes = open("test_data/p1/p1_3.png", "rb").read()
minutiae_3 = extractor.extract_minutiae(image_bytes)

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
minutiae_from_iso = extractor.load_iso_19794_2_2005(iso_template)
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
minutiae_4 = extractor.extract_minutiae_from_image_file("test_data/p2/p2_1.png")
score = minutiae_1.compare(minutiae_4)
assert score < 50, "Expected a low similarity score between p1_1 and p2_1"

# We can access the NFIQ2 quality via:
nfiq2_quality = minutiae_1.quality()
assert nfiq2_quality.score > 50, "Expected a positive NFIQ2 quality score"
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.