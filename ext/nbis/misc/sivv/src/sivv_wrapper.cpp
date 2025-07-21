#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>      // brings in cv::morphologyEx, cv::Mat, etc.
#include <opencv2/imgproc/imgproc_c.h>   // for cvGetHistValue_1D
#include <opencv2/core/core_c.h>         // for cvMinMaxLoc
#include <cstring> // for memcpy
#include <string>

// original SIVV entry point
extern std::string sivv(IplImage* src);

extern "C" {
    char* sivv_ffi_from_bytes(const unsigned char* data, int width, int height) {
        // Create IplImage header (8-bit, 1 channel)
        IplImage* img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);

        // Copy pixel data into image
        std::memcpy(img->imageData, data, width * height);

        // Call original function
        std::string result = sivv(img);

        // Cleanup
        cvReleaseImage(&img);

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