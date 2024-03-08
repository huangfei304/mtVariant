#ifndef _MTVARIANT_H
#define _MTVARIANT_H

#include <cstdint>
#include <deque>
#include <map>
#include <string>
#include <vector>

#define MTVARIANT_NAME    "mtVariant"
#define MTVARIANT_VERSION "0.4"

typedef std::vector< std::pair< std::string, uint32_t >> regions_list_t;
typedef struct opts {
    uint8_t min_var_dp=4;   //min variant depth
    int min_site_dp=10;  //min site depth
    uint8_t min_base_score=20; //min base quality
    int min_low_depth = 0; // 0.2 * median_depth

    float snp_het_min=0.10; // min read ratio for snp het
    float snp_het_max=0.85; // max read ratio for snp het
    //double snp_strand_ratio_cutoff=0.1;

    float indel_het_min=0.10; //min read ratio for indel het
    float indel_het_max=0.80; //max read ratio for indel het
    float indel_min_var_ratio=0.1; // min variant read ratio for indel

    //int near_end = 5;  //near read end
    int max_cov =65535; //2^16=65536, 2^32 = 2147483647

    int extend = 65;  // extend the chrM last with the first 65 bp base

    //D-loop region
    //HVR-I 16024-16365, HVRII 73-340, HVRIII 340-576
    float read_max_mis_ratio = 0.05;   // reads max mismatch rate
    float read_max_gap_ratio = 0.02;   // reads max gap number
    float read_flank_mis = 4;  // start/end <25 bp reads mismatch and gap number
    uint8_t min_mapq = 1;    // min reads mapping quality

    //mitochondrial hotspots around 309 and 315 as well as 3107 according to the rCRS
    float multi_var_ratio = 0.90;   // the reads rate of the two variant
    bool all_variant = false; // only output variant site

    //primer parameter
    bool filter_primer = false;
    int distance = 5;  // primer site +/- 5 bp

    // read length
    int read_len_min = 40; // min read length
    int read_len_max = 100; // max read length


    // map statistic
    int total_reads = 0;
    int fail_unmap_reads=0;
    int map_reads=0;
    int paired_reads=0;
    int no_primary=0;
    int fwd_reads=0;
    int rev_reads=0;
    int mapqual_reads=0;
    int read_short = 0;
} opts_s;

#endif /* _mtVariant_H */
