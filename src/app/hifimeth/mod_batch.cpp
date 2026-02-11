#include "mod_batch.hpp"

namespace ns_mods {

ModBatch::ModBatch(int _kmer_size,
        int _features_per_base,
        int _samples_per_batch,
        std::vector<MolMethyCall>* _samples,
        ov::CompiledModel compiled_model)
{
        kmer_size = _kmer_size;
        features_per_base = _features_per_base;
        features_per_sample = kmer_size * features_per_base;
        samples_per_batch = _samples_per_batch;
        
        num_samples = 0;
        int features_per_batch = samples_per_batch * features_per_sample;
        batch_features = (float*)calloc(features_per_batch, sizeof(float));
        current_features = batch_features;

        batch_samples = (MolMethyCall*)calloc(samples_per_batch, sizeof(MolMethyCall));
        samples = _samples;

        ov::Shape input_shape = { (size_t)samples_per_batch,  (size_t)kmer_size,  (size_t)features_per_base };
        ov::Tensor input_tensor(ov::element::f32, input_shape, batch_features);
        infer_request = compiled_model.create_infer_request();
        infer_request.set_input_tensor(0, input_tensor);

        ov::Shape out_shape = compiled_model.output().get_shape();
        size_t out_shape_bytes = sizeof(float) * ov::shape_size(out_shape);
        output_buffer = malloc(out_shape_bytes);
        output_logits = (float*)(output_buffer);
        ov::Tensor out_tensor(ov::element::f32, out_shape, output_logits);
        infer_request.set_output_tensor(out_tensor);

        total_samples = 0;
}

ModBatch::~ModBatch()
{
    free(batch_features);
    free(batch_samples);
    free(output_buffer);
}

static void
s_logits_to_methy_probs(const float* data, size_t batch_size, MolMethyCall* calls)
{
    for (size_t i = 0; i < batch_size; ++i) {
        const float* logits = data + i * 2;
        float v0 = logits[0];
        float v1 = logits[1];
        
        float max_v = std::max(v0, v1);
        float exp0 = std::exp(v0 - max_v);
        float exp1 = std::exp(v1 - max_v);
        float sum = exp0 + exp1;
	    float p1 = exp1 / sum;
	    int v = 255 * p1;
	    if (v > 255) v = 255;

        calls[i].scaled_prob = v;
    }
}

void ModBatch::call_current_batch()
{
    infer_request.infer();
    s_logits_to_methy_probs(output_logits, num_samples, batch_samples);

    samples->insert(samples->end(), batch_samples, batch_samples + num_samples);
    total_samples += num_samples;
    num_samples = 0;
    current_features = batch_features;    
}

void ModBatch::call_mods_for_one_read(const int qid, EvalKmerFeaturesGenerator* ekfg) {
    int qoff;
    int strand;
    while (ekfg->get_next_sample_features(kmer_size, features_per_base,
        current_features, qoff, strand)) {
        MolMethyCall* p = batch_samples + num_samples;
        p->qid = qid;
        p->qoff = qoff;
        p->strand = strand;
        p->scaled_prob = 0;

        ++num_samples;
        current_features += features_per_sample;

        if (num_samples == samples_per_batch) call_current_batch();
    }
}

} // ns_mods
