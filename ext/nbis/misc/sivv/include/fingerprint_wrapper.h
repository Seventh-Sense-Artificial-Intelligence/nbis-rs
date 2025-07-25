//
// Created by jiong on 25/7/25.
//

#ifndef SIVV_FINGERPRINT_WRAPPER_H
#define SIVV_FINGERPRINT_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

// C-compatible point structure
typedef struct {
    int x;
    int y;
} CPoint2i;

// C-compatible matrix structure (simplified)
typedef struct {
    void* data;
    int rows;
    int cols;
    int type;
} CMat;

// C wrapper function
CPoint2i find_fingerprint_center_morph_c(
        const unsigned char* data,
        int width,
        int height,
        int* xbound_min,
        int* xbound_max,
        int* ybound_min,
        int* ybound_max
);

#ifdef __cplusplus
}
#endif

#endif //SIVV_FINGERPRINT_WRAPPER_H
