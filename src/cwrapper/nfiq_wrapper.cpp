// nfiq_wrapper.cpp
#include "nfiq_wrapper.h"
#include <nfiq2.hpp>
#include <cstdlib>
#include <cstring>

struct Nfiq2Wrapper {
    NFIQ2::Algorithm model;
};

extern "C" {

Nfiq2Wrapper* nfiq2wrapper_create() {
    try {
        return new Nfiq2Wrapper{};
    } catch (...) {
        return nullptr;
    }
}

void nfiq2wrapper_destroy(Nfiq2Wrapper* ctx) {
    delete ctx;
}

int nfiq2wrapper_compute(Nfiq2Wrapper*    ctx,
                         const uint8_t*   data,
                         uint32_t         size,
                         uint32_t         cols,
                         uint32_t         rows,
                         uint16_t         ppi,
                         nfiq2_results_t* out)
{
    if (!ctx || !data || !out || size != cols * rows) {
        return 1;
    }

    try {
        // build the image data
        NFIQ2::FingerprintImageData img(data, size, cols, rows, 0 /*dpi units*/, ppi);

        // native measures
        auto algos = NFIQ2::QualityMeasures::computeNativeQualityMeasureAlgorithms(img);

        // unified score (reuse the same model each call!)
        out->score = ctx->model.computeUnifiedQualityScore(img);

        // actionable feedback
        auto act_ids = NFIQ2::QualityMeasures::getActionableQualityFeedbackIDs();
        auto act_map = NFIQ2::QualityMeasures::getActionableQualityFeedback(algos);

        out->actionable_count  = static_cast<uint32_t>(act_ids.size());
        out->actionable_ids    = (const char**)std::malloc(sizeof(char*) * act_ids.size());
        out->actionable_values = (double*)     std::malloc(sizeof(double) * act_ids.size());
        for (size_t i = 0; i < act_ids.size(); ++i) {
            const auto& id = act_ids[i];
            char* copy = (char*)std::malloc(id.size()+1);
            std::memcpy(copy, id.c_str(), id.size()+1);
            out->actionable_ids[i]    = copy;
            out->actionable_values[i] = act_map.at(id);
        }

        // native features
        auto feat_ids = NFIQ2::QualityMeasures::getNativeQualityMeasureIDs();
        auto feat_map = NFIQ2::QualityMeasures::getNativeQualityMeasures(algos);

        out->feature_count  = static_cast<uint32_t>(feat_ids.size());
        out->feature_ids    = (const char**)std::malloc(sizeof(char*) * feat_ids.size());
        out->feature_values = (double*)     std::malloc(sizeof(double) * feat_ids.size());
        for (size_t i = 0; i < feat_ids.size(); ++i) {
            const auto& id = feat_ids[i];
            char* copy = (char*)std::malloc(id.size()+1);
            std::memcpy(copy, id.c_str(), id.size()+1);
            out->feature_ids[i]    = copy;
            out->feature_values[i] = feat_map.at(id);
        }

        return 0;
    }
    catch (...) {
        return 2;
    }
}

void nfiq2wrapper_free_results(nfiq2_results_t* out) {
    if (!out) return;

    for (uint32_t i = 0; i < out->actionable_count; ++i) {
        std::free((void*)out->actionable_ids[i]);
    }
    std::free(out->actionable_ids);
    std::free(out->actionable_values);

    for (uint32_t i = 0; i < out->feature_count; ++i) {
        std::free((void*)out->feature_ids[i]);
    }
    std::free(out->feature_ids);
    std::free(out->feature_values);

    std::memset(out, 0, sizeof(*out));
}

} // extern "C"