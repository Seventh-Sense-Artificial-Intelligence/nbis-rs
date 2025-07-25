#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>      // brings in cv::morphologyEx, cv::Mat, etc.
#include <opencv2/opencv.hpp>
#include <cstring> // for memcpy
#include <string>
#include "fingerprint_wrapper.h"
#include "SIVVCore.h"

// original SIVV entry point
extern std::string sivv(const cv::Mat &src);
extern "C" {
    char* sivv_ffi_from_bytes(const unsigned char* data, int width, int height) {
        // Create cv::Mat from raw data (8-bit, 1 channel, grayscale)
        cv::Mat img(height, width, CV_8UC1, const_cast<unsigned char*>(data));
        // Call original function
        std::string result = sivv(img);
        // allocate and copy
        char* out = (char*)std::malloc(result.size() + 1);
        std::memcpy(out, result.c_str(), result.size() + 1);
        return out;
    }

    // free what sivv_ffi_from_bytes() allocated
    void sivv_ffi_free_bytes(char* ptr) {
        std::free(ptr);
    }
}

extern "C" {
CPoint2i find_fingerprint_center_morph_c(
        const unsigned char* data,
        int width,
        int height,
        int* xbound_min,
        int* xbound_max,
        int* ybound_min,
        int* ybound_max
) {
    // Create cv::Mat from raw data (8-bit, 1 channel, grayscale)
    cv::Mat img(height, width, CV_8UC1, const_cast<unsigned char*>(data));
    // Call the original C++ function
    cv::Point2i result = find_fingerprint_center_morph(
            img, xbound_min, xbound_max, ybound_min, ybound_max
    );

    // Convert result to C-compatible structure
    CPoint2i c_result;
    c_result.x = result.x;
    c_result.y = result.y;

    return c_result;
}
}