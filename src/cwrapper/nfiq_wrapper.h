// nfiq_wrapper.h
#pragma once

#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif

/// Opaque handle to our NFIQ2 wrapper object
typedef struct Nfiq2Wrapper Nfiq2Wrapper;

/// Quality‐score + feature arrays
typedef struct {
    uint32_t score;

    uint32_t actionable_count;
    const char** actionable_ids;
    double*      actionable_values;

    uint32_t feature_count;
    const char** feature_ids;
    double*      feature_values;
} nfiq2_results_t;

/// Create a new wrapper (allocates + initializes the embedded model)
Nfiq2Wrapper* nfiq2wrapper_create();

/// Destroy the wrapper (frees the model)
void nfiq2wrapper_destroy(Nfiq2Wrapper* ctx);

/// Compute quality on the given raw‐pixel buffer.
/// Returns 0 on success, 1 on invalid args, 2 on unexpected error.
int nfiq2wrapper_compute(Nfiq2Wrapper*    ctx,
                         const uint8_t*   data,
                         uint32_t         size,
                         uint32_t         cols,
                         uint32_t         rows,
                         uint16_t         ppi,
                         nfiq2_results_t* out);

/// Free any malloc’ed arrays inside results and zero it out.
void nfiq2wrapper_free_results(nfiq2_results_t* out);

#ifdef __cplusplus
}
#endif