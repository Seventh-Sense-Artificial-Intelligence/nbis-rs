## Python Bindings for NBIS

This package provides Python bindings for the NIST Biometric Image Software (NBIS) library, enabling minutiae extraction and matching functionalities on Linux and macOS platforms. The bindings allow you to work with biometric images, particularly for fingerprint recognition, using Python.

## Features
- Bindings to NBIS functions for minutiae extraction and matching
- Exports minutiae templates in ISO/IEC 19794-2:2005 format
- Matches minutiae templates using the NBIS Bozorth3 algorithm

## Installation
To install the Python bindings, you can use pip:

```bash
pip install nbis-py
```

## Usage
Here's a simple example of how to use the NBIS Python bindings:
```python
import nbis

# Read the bytes from a file
image_bytes = open("test_data/p1/p1_1.png", "rb").read()
minutiae_1 = nbis.extract_minutiae(image_bytes)
image_bytes = open("test_data/p1/p1_2.png", "rb").read()
minutiae_2 = nbis.extract_minutiae(image_bytes)
image_bytes = open("test_data/p1/p1_3.png", "rb").read()
minutiae_3 = nbis.extract_minutiae(image_bytes)

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
```

This example demonstrates how to extract minutiae from fingerprint images, compare them, and convert them to and from the ISO format.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you find a bug or have a feature request.

## License

This project is licensed under the MIT License - see the [LICENSE](../../../LICENSE) file for details.

## Additional Information
For more information on the NIST Biometric Image Software (NBIS) library, you can visit the [NBIS GitHub repository](https://github.com/lessandro/nbis).

