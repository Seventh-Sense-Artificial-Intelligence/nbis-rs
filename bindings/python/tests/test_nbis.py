import nbis
import numpy as np
from PIL import Image
import pytest
from io import BytesIO

# --- Fixtures ---

@pytest.fixture
def minutiae_p1_1():
    image_bytes = open("../../../test_data/p1/p1_1.png", "rb").read()
    return nbis.extract_minutiae(image=image_bytes, ppi=None)

@pytest.fixture
def minutiae_p1_2():
    image_bytes = open("../../../test_data/p1/p1_2.png", "rb").read()
    return nbis.extract_minutiae(image=image_bytes, ppi=None)

@pytest.fixture
def minutiae_p1_3():
    image_bytes = open("../../../test_data/p1/p1_3.png", "rb").read()
    return nbis.extract_minutiae(image=image_bytes, ppi=None)

@pytest.fixture
def minutiae_p2_1():
    return nbis.extract_minutiae_from_image_file("../../../test_data/p2/p2_1.png", ppi=None)

# --- Tests ---

def test_similarity_within_class(minutiae_p1_1, minutiae_p1_2, minutiae_p1_3):
    score_1_2 = minutiae_p1_1.compare(minutiae_p1_2)
    score_1_3 = minutiae_p1_1.compare(minutiae_p1_3)
    score_2_3 = minutiae_p1_2.compare(minutiae_p1_3)

    print("Score p1_1 vs p1_2:", score_1_2)
    print("Score p1_1 vs p1_3:", score_1_3)
    print("Score p1_2 vs p1_3:", score_2_3)

    assert score_1_2 > 50, "Expected high similarity between p1_1 and p1_2"
    assert score_1_3 > 50, "Expected high similarity between p1_1 and p1_3"
    assert score_2_3 > 50, "Expected high similarity between p1_2 and p1_3"

def test_iso_template_round_trip(minutiae_p1_1):
    # Set minimum quality to 0.0 for no filtering
    iso_template = minutiae_p1_1.to_iso_19794_2_2005(0.0)
    minutiae_loaded = nbis.load_iso_19794_2_2005(iso_template)

    for original, loaded in zip(minutiae_p1_1.get(), minutiae_loaded.get()):
        assert original.x() == loaded.x()
        assert original.y() == loaded.y()
        assert original.angle() == loaded.angle()
        assert original.kind() == loaded.kind()
        assert abs(original.reliability() - loaded.reliability()) < 0.1

def test_dissimilarity_between_classes(minutiae_p1_1, minutiae_p1_2, minutiae_p1_3, minutiae_p2_1):
    score = minutiae_p1_1.compare(minutiae_p2_1)
    print("Score p1_1 vs p2_1:", score)
    assert score < 50, "Expected low similarity between p1_1 and p2_1"
    score = minutiae_p1_2.compare(minutiae_p2_1)
    print("Score p1_2 vs p2_1:", score)
    assert score < 50, "Expected low similarity between p1_2 and p2_1"
    score = minutiae_p1_3.compare(minutiae_p2_1)
    print("Score p1_3 vs p2_1:", score)
    assert score < 50, "Expected low similarity between p1_3 and p2_1"


def images_equal(img1: Image.Image, img2: Image.Image) -> bool:
    """Compare two images pixel by pixel."""
    if img1.size != img2.size:
        return False
    arr1 = np.array(img1)
    arr2 = np.array(img2)
    return np.array_equal(arr1, arr2)

def test_annotated_image_matches_ground_truth():
    # This testcase will also cover `annotate_minutiae` which takes image bytes directly instead. 
    # Load the annotated image from the function (as bytes)
    annotated_bytes = nbis.annotate_minutiae_from_image_file(
        "../../../test_data/p1/p1_1.png",
        ppi=None,
        min_quality=None,
    )
    
    # Convert returned bytes into a PIL Image
    annotated_img = Image.open(BytesIO(annotated_bytes))

    # Load ground truth image
    ground_truth_img = Image.open("../../../test_data/p1/p1_1_annotated_ground_truth.png")

    # Optionally convert both to the same mode
    annotated_img = annotated_img.convert("RGB")
    ground_truth_img = ground_truth_img.convert("RGB")

    # Compare images
    assert images_equal(annotated_img, ground_truth_img), "Annotated image does not match ground truth"