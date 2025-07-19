#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>      // brings in cv::morphologyEx, cv::Mat, etc.
#include <opencv2/imgproc/imgproc_c.h>   // for cvGetHistValue_1D
#include <opencv2/core/core_c.h>         // for cvMinMaxLoc
#include <cstring> // for memcpy
#include <string>

// original SIVV entry point
extern std::string sivv(IplImage* src);

extern "C" {
    const char* sivv_ffi_from_bytes(const unsigned char* data, int width, int height) {
        static std::string result;

        // Create IplImage header (8-bit, 1 channel)
        IplImage* img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);

        // Copy pixel data into image
        std::memcpy(img->imageData, data, width * height);

        // Call original function
        result = sivv(img);

        // Cleanup
        cvReleaseImage(&img);

        return result.c_str();
    }
}