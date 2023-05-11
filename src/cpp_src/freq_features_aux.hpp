#ifndef __FREQ_FEATURES_AUX_HPP
#define __FREQ_FEATURES_AUX_HPP

#include "hbn_aux.h"

#include <cmath>

#include <algorithm>
#include <vector>

#define kMinCov 5
#define kNumSites 11
#define kMethProb 0.5
#define kUnmethProb 0.5

class ProbHist
{
public:
    ProbHist(int bins) {
        hbn_assert(bins > 0);
        m_bins = bins;
        m_step = 1.0 / bins;
        m_hist = new double[bins];
    }

    ~ProbHist() {
        delete m_hist;
    }

    void process_prob_list(double* prob_list, int list_size, std::vector<double>& norm_hist) {
        std::fill(m_hist, m_hist + m_bins, 0.0);
        for (int i = 0; i < list_size; ++i) {
            double p = prob_list[i];
            int b = p / m_step;
            if (b >= m_bins) b = m_bins - 1;
            m_hist[b] += 1.0;
        }
        double norm = 0.0;
        for (int i = 0; i < m_bins; ++i) norm += m_hist[i] * m_hist[i];
        norm = sqrt(norm);
        hbn_assert(norm > 0.001);
        for (int i = 0; i < m_bins; ++i) m_hist[i] /= norm;
        norm_hist.insert(norm_hist.end(), m_hist, m_hist + m_bins);
        //fprintf(stderr, "norm = %g\n", norm);
        //for (int i = 0; i < m_bins; ++i) fprintf(stderr, "%g\t", m_hist[i]);
        //fprintf(stderr, "\n");
    }

private:
    int m_bins;
    double m_step;
    double* m_hist;
};


#endif // __FREQ_FEATURES_AUX_HPP