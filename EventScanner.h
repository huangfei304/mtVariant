#ifndef _MTVARIANT_EVENTHOLDER_H
#define _MTVARIANT_EVENTHOLDER_H

#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include "htslib/sam.h"

#include "IndelEvent.h"
#include "SnpEvent.h"
#include "Primer.h"

namespace mtVariant {
typedef std::map< int, std::map< char, SnpEvent > > SnpMap;
typedef std::map< int, std::map< std::string, IndelEvent > > IndelMap;
typedef std::vector< std::pair< int, SnpEvent > > AReadsSnpList;
typedef std::vector< IndelEvent > AReadsIndelList;
typedef std::vector<std::pair<int,int>> AReadsPrimerList;
typedef std::map< int, std::map< int, int> > PrimerMap;
class EventScanner
{
  private:
    int _near_end=5;
    float _max_sub;
    float _max_gap;
    int _max_chr_size=16569;
    int _read_min_len;
    int _min_rlen_count;
    int _base_score=20;
    int _distance = 10; // distance for primer
    int _filter_flag; // 0 for ok, 1 for mis or gap , 2 for primer

    EventScanner();
    double snp_nqs(int p, int qlen, const uint8_t *qquals);
    int32_t left_align(const ReferenceSequence &refseq, std::string& seq, int indel_pos, int read_start_pos);
    int32_t right_align(const ReferenceSequence &refseq, std::string& seq, int indel_pos, int read_last_pos);

  public:
    SnpMap _snps;
    IndelMap _indels;
    float* _read_pos = nullptr;
    bool _primer_filter = false;
    PrimerMap _primers;

    EventScanner(float max_mis_rate, float max_gap_rate,int min_base_score, int min_read, bool primer_filter);
    EventScanner(const EventScanner &);
    ~EventScanner();
    void reset();
    void _read_pos_init(int size);
    //void collect_indels(bam1_t *read, const ReferenceSequence &sequences, CoverageCounter &coverages);
    //void collect_snps(bam1_t *read, const ReferenceSequence &sequences, CoverageCounter &coverages);
    void collect_variant(bam1_t *read, const ReferenceSequence &sequences, CoverageCounter &coverages, Primer &primer);
    int primer_start(int start);
    void primer_add_count(int pos, int start);
    int get_primer_total(int pos);
    std::vector<int> get_all_starts(int pos);
    int get_primer_start_count(int pos, int start);
    bool _isHVR(int start, int end);
    int get_filter_flag(){return _filter_flag; }
    int get_short_read_count() {return _min_rlen_count;}
};

} // namepace mtVariant

#endif /* _MTVARIANT_EVENTHOLDER_H */
