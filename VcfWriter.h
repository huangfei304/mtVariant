#ifndef _MTVARIANT_VCFWRITER_H
#define _MTVARIANT_VCFWRITER_H

#include <cstdint>
#include <string>
#include <iomanip>
#include <sstream>

#include "htslib/vcf.h"
#include "htslib/kfunc.h"
#include "Bam.h"
#include "CoverageCounter.h"
#include "EventScanner.h"
#include "ReferenceSequence.h"
#include "SnpEvent.h"
#include "IndelEvent.h"
#include "mtVariant.h"


namespace mtVariant {
typedef int32_t qual_t;

enum var_filter_idx {
    _FILTER_PASS = 0,
    _FILTER_NO_DATA,
    _FILTER_NO_VAR,
    _FILTER_LOW_COVERAGE,
    _FILTER_LOW_QUAL,
    _FILTER_LOW_VARIANTREADS,
    _FILTER_LOW_VARIANTRATIO,
    _FILTER_PRIMER_SPECIFIC,
    _FILTER_READ_END_RATIO,
    _FILTER_VRFROMDELETION
};
typedef enum var_filter_idx var_filter_idx_e;

class VcfWriter
{
  private:
    CoverageCounter *_coverages;
    ReferenceSequence *_sequences;
    EventScanner *_events;
    std::string _region;
    int _rid;
    int _pass_int;
    bcf_hdr_t *_vcf_hdr;
    htsFile *_vcf_fp;
    const opts_s *_opts;
    bcf1_t *_vcf_rec;
    int _gt_missing;

    char _snp_alleles[4];
    const char *_indel_alleles[2];

    std::map< char, uint32_t > _snp_counts;
    std::map< char, double > _snp_prs;

    std::vector< int > _filter_idx;
    //std::vector< int > _indel_filter_idx;

    VcfWriter();
    void setup_filter_idx();
    void setup_vcf_header(const char *sn,int argc,char **argv,const char *version,const char *ref,bcf_hdr_t *hdr);

  public:
    VcfWriter(CoverageCounter &coverages,ReferenceSequence &sequences,EventScanner &events,bcf_hdr_t *vcf_hdr,htsFile *vcf_fp,const opts_s *opts);
    VcfWriter(const VcfWriter &);
    ~VcfWriter();
    void set_region(const char *reg);
    void setup_vcf(const char *sn,int argc,char **argv,const char *version,const char *ref);
    void set_pl(int32_t pl_arr[3], double r, double v);
    void print_variant_buffer(const std::string chr, const int32_t chr_len);
    bool print_snp_buffer(int32_t var_pos,float _read_pos);
    bool print_indel_buffer(int32_t var_pos,float _read_pos);
    void print_novar_buffer(int32_t var_pos);
    std::vector<int> get_unique_starts(std::vector<int> first, std::vector<int> second);
};

} // namepace xatlas

#endif /* _XATLAS_VCFWRITER_H */
