#include "eval_kmer_features.hpp"

#include <vector>
#include <algorithm>
#include <cmath>

namespace ns_mods {

static void
s_extract_kmer_features(BamQuerySequence& query,
    BamKinetics& kinetics,
    MethylationContext& ctx,
    const int kmer_size,
    const int features_per_base,
    const int fwd_offset,
    int& strand,
    float* features)
{
    const int features_per_kmer = kmer_size * features_per_base;
    constexpr int NA_SIZE = EvalKmerFeaturesGenerator::NA_SIZE;
    auto ONE_HOT_ENCODING = EvalKmerFeaturesGenerator::ONE_HOT_ENCODING;

    const u8* seq;
    int offset;
    if (query.fwd_qs[fwd_offset] == ctx.FWD_MOD_BASE_CODE) {
        seq = query.fwd_qs;
        offset = fwd_offset;
        strand = FWD;
    } else {
        hbn_assert(query.fwd_qs[fwd_offset] == ctx.REV_MOD_BASE_CODE);
        seq = query.rev_qs;
        offset = query.size - 1 - fwd_offset;
        strand = REV;
        hbn_assert(seq[offset] == ctx.FWD_MOD_BASE_CODE);
    }
    const int HK = kmer_size / 2;
    const int qfrom = (offset >= HK) ? offset - HK : 0;
    const int qto = (offset + HK + 1 <= query.size) ? offset + HK + 1 : query.size;
    int fi = (HK > offset) ? (HK - offset) * features_per_base : 0;
    std::fill(features, features + features_per_kmer, 0.0f);

    for (int i = qfrom; i < qto; ++i) {
        u8 base = seq[i];
        for (int k = 0; k < NA_SIZE; ++k) features[fi++] = ONE_HOT_ENCODING[base][k];

        float v = kinetics.decoded_ipd(strand, i);
        v /= kinetics.MAX_KINETIC_VALUE;
        features[fi++] = v;

        v = kinetics.decoded_pw(strand, i);
        v /= kinetics.MAX_KINETIC_VALUE;
        features[fi++] = v;

        v = kinetics.decoded_ipd(REVERSED_STRAND(strand), query.size - 1 - i);
        v /= kinetics.MAX_KINETIC_VALUE;
        features[fi++] = v;

        v = kinetics.decoded_pw(REVERSED_STRAND(strand), query.size - 1 - i);
        v /= kinetics.MAX_KINETIC_VALUE;
        features[fi++] = v;

	//features[fi++] = fn;
	//features[fi++] = rn;
    }
}

void EvalKmerFeaturesGenerator::extract_chh_samples()
{
    _sample_offsets.clear();

    for (int i = 0; i <= M_query.size - M_ctx.CHH_MOTIF_SIZE; ++i) {
        uint8_t mi = M_ctx.get_fwd_chh_motif_idx(M_query.fwd_rqs + i);
        if (mi != M_ctx.CHH_INVALID_MOTIF_IDX) {
            _sample_offsets.push_back(i);
            continue;
        }
        mi = M_ctx.get_rev_chh_motif_idx(M_query.fwd_rqs + i);
        if (mi != M_ctx.CHH_INVALID_MOTIF_IDX) {
            _sample_offsets.push_back(i + M_ctx.CHH_MOTIF_SIZE - 1);
            continue;
        }
    }

    M_sample_offsets = _sample_offsets.data();
    M_num_samples = _sample_offsets.size();
    M_sample_idx = 0;
}

void EvalKmerFeaturesGenerator::extract_cpg_samples()
{
    _sample_offsets.clear();

    for (int i = 0; i <= M_query.size - M_ctx.CPG_MOTIF_SIZE; ++i) {
        if (M_query.fwd_rqs[i] == 'C' && M_query.fwd_rqs[i+1] == 'G') {
            _sample_offsets.push_back(i);
        }
    }

    M_sample_offsets = _sample_offsets.data();
    M_num_samples = _sample_offsets.size();
    M_sample_idx = 0;
}

void EvalKmerFeaturesGenerator::extract_chg_samples()
{
    _sample_offsets.clear();

    for (int i = 0; i <= M_query.size - M_ctx.CHG_MOTIG_SIZE; ++i) {
        if (strncasecmp(M_query.fwd_rqs + i, "CCG", 3) == 0) {
            _sample_offsets.push_back(i);
            continue;
        }
        if (strncasecmp(M_query.fwd_rqs + i, "CAG", 3) == 0) {
            _sample_offsets.push_back(i);
            continue;
        }
        if (strncasecmp(M_query.fwd_rqs + i, "CTG", 3) == 0) {
            _sample_offsets.push_back(i);
            continue;
        }
    }

    M_sample_offsets = _sample_offsets.data();
    M_num_samples = _sample_offsets.size();
    M_sample_idx = 0;
}

bool EvalKmerFeaturesGenerator::get_next_sample_features(
    const int kmer_size, const int features_per_base,
    float* features, int& offset, int& strand)
{
    if (M_sample_idx >= M_num_samples) return false;
    offset = M_sample_offsets[M_sample_idx++];
    s_extract_kmer_features(M_query, M_kinetics, M_ctx, kmer_size, features_per_base, offset, strand, features);
    return true;
}

} // ns_mods
