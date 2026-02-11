#ifndef __MOD_BATCH_HPP
#define __MOD_BATCH_HPP

#include "eval_kmer_features.hpp"
#include "mod_options.hpp"

#include <corelib/5mc_motif_finder.hpp>
#include <openvino/openvino.hpp>

namespace ns_mods {

struct ModBatch {
    int kmer_size;
    int features_per_base;
    int features_per_sample;
    int samples_per_batch;

    int num_samples;
    float* batch_features;
    float* current_features;
    MolMethyCall* batch_samples;

    void* output_buffer;
    float* output_logits;

    size_t total_samples;

    std::vector<MolMethyCall>* samples;

    ov::InferRequest infer_request;

    ~ModBatch();

    ModBatch(int _kmer_size,
        int _features_per_base,
        int _samples_per_batch,
        std::vector<MolMethyCall>* _samples,
        ov::CompiledModel compiled_model);

    void call_current_batch();

    void call_mods_for_one_read(const int qid, EvalKmerFeaturesGenerator* ekfg);
};

} // ns_mods

#endif // __MOD_BATCH_HPP