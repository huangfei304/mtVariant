#ifndef _MTVARIANT_COVERAGECOUNTER_H
#define _MTVARIANT_COVERAGECOUNTER_H

namespace mtVariant {

static const int COVS_ALIGN = 6;

enum covs_idx {
    COVSIDX_DP   =   0,   // site reference depth
    COVSIDX_LQ   =   1,   // site reference depth with low quality
    COVSIDX_RS   =   2,  // reference positive strand depth
    COVSIDX_FT   =   3,  // site for depth with filtered
    COVSIDX_ALL  =   4,  // site for depth including variant
    COVSIDX_RP    =  5   // reference site in reads position
    // COVSIDX_INS  =   3  // indel insert
   // COVSIDX_SNP  =   4,  // snp variant depth
   // COVSIDX_INDEL =  5   // indel variant depth
};
typedef covs_idx covs_idx_e;

class CoverageCounter
{
  private:
    int (*_covs)[COVS_ALIGN];
    int _len;
    int _max_cov;

    CoverageCounter();

  public:
    CoverageCounter(int len, int max_cov);
    CoverageCounter(const CoverageCounter &);
    ~CoverageCounter();
    void reset();
    void add_coverage(int pos, covs_idx_e idx);
    void add_coverage(int pos, covs_idx_e idx, int add_cov);
    void add_coverage_range(int pos, covs_idx_e idx, int len);
    void add_coverage_range(int pos, covs_idx_e idx, int len, int add_cov);
    void set_coverage(int pos, covs_idx_e idx, int cov);
    int get_coverage(int pos, covs_idx_e idx);
    bool test_high_coverage(int pos, covs_idx_e idx, int l_qseq);
    int get_max_coverage();
    int get_length();
};

} // namepace mtVariant

#endif /* _MTVARIANT_COVERAGECOUNTER_H */
