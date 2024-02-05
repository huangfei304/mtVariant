#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>

#include "VcfWriter.h"

namespace mtVariant {

std::vector<int> VcfWriter::get_unique_starts(std::vector<int> first, std::vector<int> second){
    std::map<int, int> tmp_map;
    for(auto ele: first){
        tmp_map.emplace(ele,1);
    }
    for(auto ele: second){
        tmp_map.emplace(ele,1);
    }
    std::vector<int> uniq_list;
    for(auto ele: tmp_map){
        uniq_list.push_back(ele.first);
    }
    return uniq_list;
}

VcfWriter::VcfWriter()
{
}

VcfWriter::VcfWriter(
        CoverageCounter &coverages,
        ReferenceSequence &sequences,
        EventScanner &events,
        bcf_hdr_t *vcf_hdr,
        htsFile *vcf_fp,
        const opts_s *opts)
    : _coverages(&coverages),
      _sequences(&sequences),
      _events(&events),
      _region(""),
      _rid(-1),
      _vcf_hdr(vcf_hdr),
      _vcf_fp(vcf_fp),
      _opts(opts)
{
    if ((_vcf_rec = bcf_init1()) != nullptr)
    {
        //_pass_int = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, "PASS");
        _gt_missing = bcf_int32_missing;

        std::strncpy(_snp_alleles, ".,.", 4);
    }else{
        std::cerr<<" BCF init failed."<<std::endl;
    }
}

VcfWriter::VcfWriter(const VcfWriter &) = default;

VcfWriter::~VcfWriter()
{
    bcf_destroy1(_vcf_rec);
}

void VcfWriter::set_region(const char *reg)
{
    _region = reg;
    _rid = bcf_hdr_name2id(_vcf_hdr, reg);
}

void VcfWriter::setup_vcf_header(const char *sn, int argc, char **argv, const char *version, const char *ref, bcf_hdr_t *hdr) {
    time_t today = time(nullptr);
    tm now;
    localtime_r(&today, &now);

    bcf_hdr_printf(hdr, "##fileDate=%d%02d%02d\n", 1900 + now.tm_year, 1 + now.tm_mon, now.tm_mday);
    //bcf_hdr_printf(hdr, "##source=%s v%s\n", *argv, version);
    bcf_hdr_printf(hdr, "##reference=%s\n", ref);

    // FOMRAT
    bcf_hdr_append(hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled Genotype Likelihood Scores\">\n");

    // FILTER
    bcf_hdr_printf(hdr, "##FILTER=<ID=lowVR,Description=\"Variant read depth < %d\">\n", _opts->min_var_dp);
    bcf_hdr_append(hdr, "##FILTER=<ID=lowVRT,Description=\"Variant read ratio is less than PL cutoff\">\n");
    bcf_hdr_printf(hdr, "##FILTER=<ID=lowC,Description=\"Total coverage < %d (0.2*median_depth)\">\n", _opts->min_site_dp);
    bcf_hdr_append(hdr, "##FILTER=<ID=lowQ,Description=\"Variant average base quality is low (Q<20)\">\n");
    //bcf_hdr_append(hdr, "##FILTER=<ID=highRE,Description=\"Variant with high ratio of read near start or end (< 5 bp) > 0.96\">\n");
    //bcf_hdr_append(hdr, "##FILTER=<ID=,Description=\"Variant with high ratio of read near start or end (< 5 bp) > 0.96\">\n");

    bcf_hdr_append(hdr, "##FILTER=<ID=primerSpec,Description=\"Most variant reads are in one single primer, fisher test p-value <0.01\">\n");
    bcf_hdr_append(hdr, "##FILTER=<ID=NoDT,Description=\"No valid reads on this site\">\n");
    bcf_hdr_append(hdr, "##FILTER=<ID=NoVAR,Description=\"No valid variants reads on this site\">\n");
    //bcf_hdr_printf(hdr, "##FILTER=<ID=RER,Description=\"Ratio of variants reads within %d bp of read end is greater than %.2f\">\n", _opts->snp_near_end_bases);
    //bcf_hdr_append(hdr, "##FILTER=<ID=VRD,Description=\"Variant read coverage recorded in non-variant block originates from deletion called upstream\">\n");

    // INFO
    //bcf_hdr_append(hdr, "##INFO=<ID=P,Number=1,Type=Float,Description=\"p-value\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=EM,Number=0,Type=Flag,Description=\"The called SNP has an equal number of reads indicating another variant call and base was chosen by highest summed quality score\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=SRC,Number=4,Type=Integer,Description=\"Strand Reads Count, format: ref_fwd, ref_rev, var_fwd, var_rev\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=VBQ,Number=1,Type=Integer,Description=\"Variant Average base quaility\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=VFBQ,Number=1,Type=Integer,Description=\"Variant flank(+/- 5 bp) Average base quaility\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=NER,Number=2,Type=Float,Description=\"Near end reads (<5 bp) ratio for ref and allele\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=VARC, Number=1,Type=Float,Description=\"Average mismatch and gap rate by reads length\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=PRR, Number=4,Type=Float,Description=\"Read ratio in top two primers(>10 bp), format: ref_primer1,allele_primer1,ref_primer2,allele_primer2\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=PRC, Number=4,Type=Integer,Description=\"Read Count in top two Primers(>10 bp), format: ref_primer1,allele_primer1,ref_primer2,allele_primer2\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=LowQR, Number=2,Type=Float,Description=\"Low quality (<20) (average) read ratio, format: ref_lowqual_ratio,allele_lowqual_ratio\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=PSpval, Number=1,Type=Float,Description=\"Primer Specific fisher extract test p-value\">\n");
    bcf_hdr_append(hdr, "##INFO=<ID=BaseNum, Number=1, Type=String,Description=\"Each base depth\">\n");


    // contig
    //for (const auto &contig : regions) {
    //    if(contig.first == "chrM"){
    //        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%u>\n", contig.first.c_str(), (contig.second-_opts->shift));
    //    }else{
    //        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%u>\n", contig.first.c_str(), contig.second);
    //    }
    //}
    bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%u>\n", "chrM", 16569);

    // command
    const size_t BUFFER_SIZE = 0x1000;
    char *buffer = new char[BUFFER_SIZE];

    std::strncpy(buffer, "##command=", BUFFER_SIZE);
    for (int i = 0; i < argc; ++i) {
        std::strncat(buffer, argv[i], BUFFER_SIZE - 1);
        if (i != argc - 1) {
            std::strncat(buffer, " ", 2);
        }
    }
    std::strncat(buffer, "\n", 2);
    bcf_hdr_append(hdr, buffer);

    delete[] buffer;

    bcf_hdr_add_sample(hdr, sn);
    bcf_hdr_add_sample(hdr, nullptr);
}

void VcfWriter::setup_filter_idx() {
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "PASS"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "NoDT"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "NoVAR"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "lowC"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "lowQ"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "lowVR"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "lowVRT"));
    _filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "primerSpec"));
    //_filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "highRE"));
    //_filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "RER"));
    //_filter_idx.push_back(bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "VRD"));
}

void VcfWriter::setup_vcf(
        const char *sn,
        int argc,
        char **argv,
        const char *version,
        const char *ref)
{
    this->setup_vcf_header(sn, argc, argv, version, ref, _vcf_hdr);
    this->setup_filter_idx();
    _pass_int = bcf_hdr_id2int(_vcf_hdr, BCF_DT_ID, "PASS");
}

/**
 * Set PL values fro SNPs
 */
void VcfWriter::set_pl(int32_t pl_arr[3], double r, double v)
{
    const double error_rate = 0.001351248680180233;
    const double var_rate = 0.00333333333333;
    const double het_hom_ratio = 2.0;
    const double log_e = std::log10(error_rate);
    const double log_1_e = std::log10(1.0 - error_rate);

    double hom_prob = 1.0 / (1.0 + het_hom_ratio);
    double het_prob = 1.0 - hom_prob;
    double ref_prior = 1.0 - var_rate;
    double hom_prior = var_rate * hom_prob;
    double het_prior = var_rate * het_prob;
    double g[3];
    std::vector< int32_t > g_vec;

    g[0] = (int32_t)std::round(-10.0 * (std::log10(hom_prior) + r * log_e + v * log_1_e));   //  1/1
    g[1] = (int32_t)std::round(-10.0 * (std::log10(het_prior) + (r + v) * std::log10(0.5))); //  0/1
    g[2] = (int32_t)std::round(-10.0 * (std::log10(ref_prior) + r * log_1_e + v * log_e));   //  0/0

    g_vec.push_back(g[0]);
    g_vec.push_back(g[1]);
    g_vec.push_back(g[2]);

    std::sort(g_vec.begin(), g_vec.end());

    pl_arr[2] = g[0] - g_vec[0];
    pl_arr[1] = g[1] - g_vec[0];
    pl_arr[0] = g[2] - g_vec[0];
}

/**
 * Print SNPs
 */
bool VcfWriter::print_snp_buffer(int32_t var_pos,float _read_pos)
{
    //std::cerr<<"==========SNP start====================\n";
    bool equal_majority;
    char refbase;
    double p_value;
    float tmpf;
    int32_t tmpi, gt_arr[2], pl_arr[3],base_qual, flank_base_qual;
    int var_cov, dp_cov, rr_cov, ar_cov, total_coverage;
    std::string basecount="";

    // for each snp event up to next variant or end
    SnpMap::iterator sv_it = _events->_snps.find(var_pos);
    ar_cov = 0;
    _snp_counts.clear();
    _snp_prs.clear();
    for (auto &allele_it : sv_it->second) {
        size_t num_allele = allele_it.second._read_count;
        _snp_counts[allele_it.first] = num_allele;
        ar_cov += num_allele;
        basecount.push_back(allele_it.second._allele_base);
        //basecount.push_back(':');
        basecount = basecount + ":"+std::to_string(num_allele)+",";
        //basecount.push_back(',');
        //std::cerr<<"baseCount:"<<basecount<<std::endl;
    }

    refbase = _sequences->_seq[sv_it->first];
    rr_cov = _coverages->get_coverage(sv_it->first, COVSIDX_DP);
    total_coverage = rr_cov + ar_cov;
    std::memset(pl_arr, 0, sizeof pl_arr);

    // record allele with most supporting reads and its base quality score sum
    char high_base = '\0';
    int high_base_cov = 0;
    equal_majority = false;

    for (auto &allele_it : sv_it->second) {
        if (_snp_counts[allele_it.first] > high_base_cov) {
            high_base = allele_it.first;
            high_base_cov = _snp_counts[high_base];
            equal_majority = false;
        } else if (_snp_counts[allele_it.first] == high_base_cov) {
            equal_majority = true;
        }
    }

    // select equal majority snp with highest base quality sum
    if (equal_majority) {
        for (auto &allele_it : sv_it->second) {
          _snp_prs[allele_it.first] = allele_it.second.get_qual();
        }

        double high_pr = _snp_prs[high_base];
        for (auto &allele_it : sv_it->second) {
            if (_snp_counts[allele_it.first] == high_base_cov && _snp_prs[allele_it.first] > high_pr) {
                high_pr = _snp_prs[allele_it.first];
            }
        }
    }else{
        int cov_cutoff = int(high_base_cov * _opts->multi_var_ratio);
        for (auto &allele_it : sv_it->second) {
            if( _snp_counts[allele_it.first] == high_base_cov ) continue;
            if (_snp_counts[allele_it.first] >= cov_cutoff ) {
                equal_majority = true;
                break;
            }
        }
    }

    // coverage for allele with most supporting reads
    auto &called_snp = sv_it->second.at(high_base);

    var_cov = std::min(called_snp._read_count, _coverages->get_max_coverage());
    //_coverages->set_coverage(sv_it->first, COVSIDX_VR, var_cov);
    //_coverages->add_coverage(sv_it->first, COVSIDX_DP, var_cov);

    bcf_clear1(_vcf_rec);

    // set filter and genotype
    bool filter_flag = true;
    bool genotype = false;
    bool _het_flag = false;
    if ((rr_cov+var_cov) == 0) {
        bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_NO_DATA], 1);
        //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_NO_DATA]);
        gt_arr[0] = bcf_gt_missing;
        gt_arr[1] = bcf_gt_missing;
    } else {
        this->set_pl(pl_arr, (double)(rr_cov), (double)var_cov);
        /**
        if (pl_arr[2] < pl_arr[0]) {
            if (pl_arr[2] < pl_arr[1]) {
                gt_arr[0] = bcf_gt_unphased(1);
                gt_arr[1] = bcf_gt_unphased(1);
            } else {
                gt_arr[0] = bcf_gt_unphased(0);
                gt_arr[1] = bcf_gt_unphased(1);
            }
            genotype = true;
        } else {
            if (pl_arr[1] < pl_arr[0]) {
                gt_arr[0] = bcf_gt_unphased(0);
                gt_arr[1] = bcf_gt_unphased(1);
                genotype = true;
            } else {
                gt_arr[0] = bcf_gt_unphased(0);
                gt_arr[1] = bcf_gt_unphased(0);
                //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_LOW_VARIANTRATIO]);
            }
        }
        **/
        double var_rate = (double)var_cov / (double)(rr_cov+var_cov);
        if( var_rate >= _opts->snp_het_min ){
            if( var_rate >= _opts->snp_het_max){
                gt_arr[0] = bcf_gt_unphased(1);
            }else{
                gt_arr[0] = bcf_gt_unphased(0);
                _het_flag = true;
            }
            gt_arr[1] = bcf_gt_unphased(1);
            genotype = true;
        }else{
            gt_arr[0] = bcf_gt_unphased(0);
            gt_arr[1] = bcf_gt_unphased(0);
        }
        if( genotype ){
            int pos_cov = _coverages->get_coverage(sv_it->first, COVSIDX_RS);
            std::string strand_info = std::to_string(pos_cov)+","+std::to_string(rr_cov-pos_cov)+","+std::to_string(called_snp._pos_strand)+","+std::to_string(called_snp._read_count - called_snp._pos_strand);
            bcf_update_info_string(_vcf_hdr, _vcf_rec, "SRC", strand_info.c_str());

            base_qual =  int32_t(called_snp.get_qual());
            flank_base_qual=int32_t(called_snp.get_nqs());
            bcf_update_info_int32(_vcf_hdr, _vcf_rec, "VBQ", &base_qual, 1);
            bcf_update_info_int32(_vcf_hdr, _vcf_rec, "VFBQ", &flank_base_qual, 1);
            //float rel_pos = float(called_snp.get_rel_pos());
            float ref_rel_pos = (rr_cov ==0) ? 0.0 : (float)_read_pos/(float)rr_cov;
            float rel_pos[2] = { ref_rel_pos, float(called_snp.get_rel_pos())};
            bcf_update_info_float(_vcf_hdr, _vcf_rec, "NER", &rel_pos, 2);
            double fisher_left_p, fisher_right_p, fisher_twosided_p;
            if( _opts->filter_primer ){
                float vpr_info[4]={0.0,0.0,0.0,0.0}; //ref_p1,alt_p1, ref_p2, alt_p2
                int prc_info[4]={0,0,0,0}; //ref_pr1,alt_pr1, ref_pr2, alt_pr2
                int snp_primer_total = called_snp.get_total_start();
                int ref_primer_total = _events->get_primer_total(var_pos);
                std::vector<int> primer_list=get_unique_starts(called_snp.get_all_starts(),_events->get_all_starts(var_pos));
                if( primer_list.size()==1 ){
                    prc_info[0] = _events->get_primer_start_count(var_pos, primer_list[0]);
                    prc_info[1] = called_snp.get_start_count(primer_list[0]);
                }else {
                    std::vector<int> prc_total;
                    int first_idx,second_idx;
                    first_idx=0;
                    for(int i=0;i<primer_list.size();++i){
                        int prc_tmp = called_snp.get_start_count(primer_list[i])+_events->get_primer_start_count(var_pos, primer_list[i]);
                        prc_total.push_back(prc_tmp);
                        if( prc_tmp > prc_total[first_idx] ) first_idx = i;
                    }
                    second_idx = (first_idx == 0 ) ? 1 : 0;
                    for(int i=0;i<prc_total.size();++i){
                        if( i == first_idx ) continue;
                        if( prc_total[i] > prc_total[second_idx]) second_idx = i;
                    }
                    prc_info[0] = _events->get_primer_start_count(var_pos, primer_list[first_idx]);
                    prc_info[1] = called_snp.get_start_count(primer_list[first_idx]);
                    prc_info[2] = _events->get_primer_start_count(var_pos, primer_list[second_idx]);
                    prc_info[3] = called_snp.get_start_count(primer_list[second_idx]);
                }
                if( ref_primer_total > 0 ){
                    vpr_info[0] = (float)prc_info[0] /(float)ref_primer_total;
                    vpr_info[2] = (float)prc_info[2] /(float)ref_primer_total;
                }
                if( snp_primer_total > 0 ){
                    vpr_info[1] = (float)prc_info[1] /(float)snp_primer_total;
                    vpr_info[3] = (float)prc_info[3] /(float)snp_primer_total;
                }
                //int n11, int n12, int n21, int n22,
                kt_fisher_exact(prc_info[0], prc_info[2], prc_info[1], prc_info[3], &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
                //kt_fisher_exact(prc_info[0], prc_info[1], prc_info[2], prc_info[3], &fisher_left_p, &fisher_right_p, &fisher_twosided_p);

                bcf_update_info_float(_vcf_hdr, _vcf_rec, "PRR", &vpr_info, 4);
                bcf_update_info_int32(_vcf_hdr, _vcf_rec, "PRC", &prc_info, 4);
                if( _het_flag ){
                    float pval_tmp = (float) fisher_twosided_p;
                    bcf_update_info_float(_vcf_hdr, _vcf_rec, "PSpval", &pval_tmp, 1);
                }
            }

            int ref_lowq_count = _coverages->get_coverage(sv_it->first, COVSIDX_LQ);
            float ref_lowq_rate = (rr_cov == 0 ) ? 0.0 : (float)ref_lowq_count/(float)rr_cov;
            float low_qual[2] = {ref_lowq_rate, called_snp.get_lowqual_rate()};
            bcf_update_info_float(_vcf_hdr, _vcf_rec, "LowQR", &low_qual, 2);

            if( base_qual < 20 ) {
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_QUAL], 1);
            }else if( _opts->filter_primer && _het_flag && fisher_twosided_p < 0.01 ){
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_PRIMER_SPECIFIC], 1);
            }else if ( (rr_cov+var_cov) < _opts->min_site_dp ) {
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_COVERAGE],1);
                //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_LOW_COVERAGE]);
            }else if (var_cov <= _opts->min_var_dp ) {
                if( var_cov == 0 ) bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_NO_VAR],1);
                else  bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_VARIANTREADS],1);
            }else{
                filter_flag = false;
            }
        }
    }

    basecount.push_back(refbase);
    basecount = basecount + ":" + std::to_string(rr_cov);
    bcf_update_info_string(_vcf_hdr, _vcf_rec, "BaseNum", basecount.c_str());
    //std::cerr<<"final base count:"<<basecount<<std::endl;


    if( !filter_flag )  bcf_update_filter(_vcf_hdr, _vcf_rec, &_pass_int, 1);
    // write variant
    _vcf_rec->rid = _rid;
    _vcf_rec->pos = sv_it->first;
    _snp_alleles[0] = refbase;
    _snp_alleles[2] = called_snp._allele_base;
    bcf_update_alleles_str(_vcf_hdr, _vcf_rec, _snp_alleles);
    //_vcf_rec->qual = int(called_snp.get_qual());
    //_vcf_rec->qual = 60;
    bcf_update_genotypes(_vcf_hdr, _vcf_rec, gt_arr, 2);

    if (genotype && equal_majority) {
        bcf_update_info_flag(_vcf_hdr, _vcf_rec, "EM", nullptr, 1);
    }

    //tmpi = (int32_t)(rr_cov - var_cov);
    int32_t ads[2]={rr_cov, (int32_t)var_cov};
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "AD", &ads, 2);
    //bcf_update_format_int32(_vcf_hdr, _vcf_rec, "RR", &tmpi, 1);
    tmpi = (int32_t)total_coverage;
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "DP", &tmpi, 1);
    if( !filter_flag ){
        int tmp_i = std::max(pl_arr[1],pl_arr[2]);
        tmpi = std::min(pl_arr[0],tmp_i);
        tmpi = std::min(99, tmpi);
        bcf_update_format_int32(_vcf_hdr, _vcf_rec, "GQ", &tmpi, 1);
        bcf_update_format_int32(_vcf_hdr, _vcf_rec, "PL", &pl_arr, 3);
    }

    //if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
    //    std::cerr << "Failed to write snp VCF record" << std::endl;
    //    exit(EXIT_FAILURE);
    //}
    //std::cerr<<"==========SNP end====================\n";
    if( equal_majority && genotype ) genotype = false;
    //if( genotype && called_snp.get_rel_pos() >0.96 ){
    //    bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_READ_END_RATIO],1);
        // genotype = false;
    //}
    return genotype;
}

/**
 * Process indels
 */
bool VcfWriter::print_indel_buffer(int32_t var_pos,float _read_pos)
{
    //std::cerr<<"==========InDel start====================\n";
    char refbase;
    double p_value, var_ratio, ref_ratio;
    float tmpf, near_read_end_ratio;
    int32_t tmpi, var_start, gt_arr[2], pl_arr[3], base_qual, flank_base_qual;
    std::string alt, ref;
    int ar_cov, total_coverage, rr_cov;

    IndelMap::iterator iv_it = _events->_indels.find(var_pos);
    //for (iv_it = _events->_indels.begin(); iv_it != _events->_indels.end() && iv_it->first < next_var_pos; ++iv_it) {
    ar_cov = 0;
    auto ivg_it = iv_it->second.begin();
    IndelEvent *called_indel_ptr = &ivg_it->second;

    auto called_indel_read_count = std::min(called_indel_ptr->_read_count, _coverages->get_max_coverage());

    // find variant with maximum supporting reads
    while (ivg_it != iv_it->second.end()) {
        ar_cov += ivg_it->second._read_count;
        if (ivg_it->second._read_count > called_indel_read_count) {
            called_indel_ptr = &ivg_it->second;
            called_indel_read_count = std::min(called_indel_ptr->_read_count, _coverages->get_max_coverage());
        }
        ++ivg_it;
    }

    IndelEvent &called_indel = *called_indel_ptr;
    // "start position" for coverage
    var_start = called_indel._var_start;
    // what happens if var_start == 0, is that even valid?
    if (!called_indel._isdel && var_start > 0) {
        --var_start;
    }

    rr_cov = _coverages->get_coverage(var_start, COVSIDX_DP);
    if (!called_indel._isdel) rr_cov = (rr_cov > ar_cov ) ? (rr_cov - ar_cov) : 0;
    total_coverage = ar_cov + rr_cov;

    bcf_clear1(_vcf_rec);
    //bcf_update_filter(_vcf_hdr, _vcf_rec, &_pass_int, 1);

    //set filter
    bool filter_flag = true;
    bool _het_flag = false;
    bool genotype = false;
    if (total_coverage == 0) {
        bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_NO_DATA],1);
        gt_arr[0] = bcf_gt_missing;
        gt_arr[1] = bcf_gt_missing;
    }else{
        var_ratio = (double)called_indel_read_count / (double)total_coverage;
        ref_ratio = (double) rr_cov / (double)total_coverage;
        if ( var_ratio <= _opts->indel_het_min) {
            gt_arr[0] = bcf_gt_unphased(0);
            gt_arr[1] = bcf_gt_unphased(0);
        } else {
            genotype = true;
            if (var_ratio < _opts->indel_het_max) {
                gt_arr[0] = bcf_gt_unphased(0);
                gt_arr[1] = bcf_gt_unphased(1);
                _het_flag = true;
            } else {
                gt_arr[0] = bcf_gt_unphased(1);
                gt_arr[1] = bcf_gt_unphased(1);
            }
        }
        if( genotype ){
            flank_base_qual = called_indel.get_mean_avg_nqs();
            near_read_end_ratio = float(called_indel.get_near_read_end_ratio());
            bcf_update_info_int32(_vcf_hdr, _vcf_rec, "VFBQ", &flank_base_qual, 1);
            float read_rel_pos = (_coverages->get_coverage(var_start, COVSIDX_DP)>0) ? _read_pos/(float)_coverages->get_coverage(var_start, COVSIDX_DP) : 0.0;
            float rel_pos[2] = {read_rel_pos, near_read_end_ratio};
            bcf_update_info_float(_vcf_hdr, _vcf_rec, "NER", &rel_pos, 2);
            float mean_var_rate = float(called_indel.get_mean_var_rate());
            bcf_update_info_float(_vcf_hdr, _vcf_rec, "VARC", &mean_var_rate, 1);
            double fisher_left_p, fisher_right_p, fisher_twosided_p;
            if( _opts->filter_primer ){
                float vpr_info[4]={0.0,0.0,0.0,0.0}; //ref_p1,indel_p1, ref_p2, indel_p2
                int prc_info[4]={0,0,0,0}; //ref_pr1, indel_pr1, ref_pr2, indel_pr2
                int indel_primer_total = called_indel.get_total_start();
                int ref_primer_total = _events->get_primer_total(var_start);
                std::vector<int> primer_list=get_unique_starts(called_indel.get_all_starts(),_events->get_all_starts(var_start));
                if( primer_list.size()==1 ){
                    prc_info[0] = _events->get_primer_start_count(var_pos, primer_list[0]);
                    prc_info[1] = called_indel.get_start_count(primer_list[0]);
                }else {
                    std::vector<int> prc_total;
                    int first_idx,second_idx;
                    first_idx=0;
                    for(int i=0;i<primer_list.size();++i){
                        int prc_tmp = called_indel.get_start_count(primer_list[i])+_events->get_primer_start_count(var_pos, primer_list[i]);
                        prc_total.push_back(prc_tmp);
                        if( prc_tmp > prc_total[first_idx] ) first_idx = i;
                    }
                    second_idx = (first_idx == 0 ) ? 1 : 0;
                    for(int i=0;i<prc_total.size();++i){
                        if( i == first_idx ) continue;
                        if( prc_total[i] > prc_total[second_idx]) second_idx = i;
                    }
                    prc_info[0] = _events->get_primer_start_count(var_pos, primer_list[first_idx]);
                    prc_info[1] = called_indel.get_start_count(primer_list[first_idx]);
                    prc_info[2] = _events->get_primer_start_count(var_pos, primer_list[second_idx]);
                    prc_info[3] = called_indel.get_start_count(primer_list[second_idx]);
                }
                if( ref_primer_total > 0 ){
                    vpr_info[0] = (float)prc_info[0] /(float)ref_primer_total;
                    vpr_info[2] = (float)prc_info[2] /(float)ref_primer_total;
                }
                if( indel_primer_total > 0 ){
                    vpr_info[1] = (float)prc_info[1] /(float)indel_primer_total;
                    vpr_info[3] = (float)prc_info[3] /(float)indel_primer_total;
                }
                //int n11, int n12, int n21, int n22,
                kt_fisher_exact(prc_info[0], prc_info[2], prc_info[1], prc_info[3], &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
                //kt_fisher_exact(prc_info[0], prc_info[1], prc_info[2], prc_info[3], &fisher_left_p, &fisher_right_p, &fisher_twosided_p);

                bcf_update_info_float(_vcf_hdr, _vcf_rec, "PRR", &vpr_info, 4);
                bcf_update_info_int32(_vcf_hdr, _vcf_rec, "PRC", &prc_info, 4);
                if( _het_flag ){
                    float pval_tmp = (float) fisher_twosided_p;
                    bcf_update_info_float(_vcf_hdr, _vcf_rec, "PSpval", &pval_tmp, 1);
                }
            }

            int ref_lowq_count = _coverages->get_coverage(var_start, COVSIDX_LQ);
            float ref_lowq_rate = (rr_cov == 0 ) ? 0.0 : (float)ref_lowq_count/(float)rr_cov;
            float low_qual[2] = {ref_lowq_rate, called_indel.get_lowqual_rate()};
            bcf_update_info_float(_vcf_hdr, _vcf_rec, "LowQR", &low_qual, 2);

            if (called_indel_read_count < _opts->min_var_dp) {
                //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_LOW_VARIANTREADS]);
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_VARIANTREADS],1);
            }else if( _opts->filter_primer && _het_flag && fisher_twosided_p < 0.01 ){
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_PRIMER_SPECIFIC], 1);
            }else if((called_indel_read_count+rr_cov) < _opts->min_site_dp) {
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_COVERAGE],1);
                genotype = false;
            }else if(var_ratio < _opts->indel_min_var_ratio) {
                bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_VARIANTRATIO],1);
            }else{
                filter_flag=false;
            }
        }
    }

    if( !filter_flag ) bcf_update_filter(_vcf_hdr, _vcf_rec, &_pass_int,1);

    double score, read_end_score[3];
    double near_read_end_sum[3];
    uint32_t near_read_end_n[3];

    // ref and alt sequences
    ref.clear();
    alt.clear();
    if (called_indel._isdel) {
        ref.append(_sequences->_seq + called_indel._var_start - 1, (size_t)called_indel._var_len + 1);
        alt.push_back(_sequences->_seq[called_indel._var_start - 1]);
    } else {
        refbase = _sequences->_seq[called_indel._var_start - 1];
        ref.push_back(refbase);
        alt.push_back(refbase);
        alt.append(called_indel._seq);
    }

    //if ( genotype ) { // 1M1D74M or 74M1D1M
        // if(called_indel._isdel && (near_read_end_ratio>0.98 || near_read_end_ratio <0.02) && ref.length()==2 ){
        //if( genotype && near_read_end_ratio>0.96 ){
       //     bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_READ_END_RATIO],1);
        //}
    //}

    _vcf_rec->rid = _rid;
    _vcf_rec->pos = (called_indel._var_start - 1);
    _indel_alleles[0] = ref.c_str();
    _indel_alleles[1] = alt.c_str();
    bcf_update_alleles(_vcf_hdr, _vcf_rec, _indel_alleles, 2);
    bcf_update_genotypes(_vcf_hdr, _vcf_rec, gt_arr, 2);

    //if( ref =="CN" ){
    //    std::cerr<<"var_pos:"<<std::to_string(var_pos)<<std::endl;
    //}
    //tmpi = (int32_t)called_indel_read_count;
    int32_t ads[2]={(int32_t)rr_cov, (int32_t)called_indel_read_count};
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "AD", &ads, 2);
    tmpi = (int32_t)(total_coverage);
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "DP", &tmpi, 1);
    if( !filter_flag){
        this->set_pl(pl_arr, (double)rr_cov, (double)called_indel_read_count);
        int tmp_i = std::max(pl_arr[1],pl_arr[2]);
        tmpi = std::min(pl_arr[0],tmp_i);
        tmpi = std::min(99, tmpi);
        bcf_update_format_int32(_vcf_hdr, _vcf_rec, "GQ", &tmpi, 1);
        bcf_update_format_int32(_vcf_hdr, _vcf_rec, "PL", &pl_arr, 3);
    }
    //if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
    //    std::cerr << "Failed to write indel VCF record" << std::endl;
    //    exit(EXIT_FAILURE);
    //}
    //std::cerr<<"==========InDel end====================\n";
    return genotype;
}
void VcfWriter::print_novar_buffer(int32_t var_pos)
{
    //std::cerr<<"==========NoVar start====================\n";
    bool equal_majority;
    char refbase;
    double p_value;
    float tmpf;
    int32_t tmpi, gt_arr[2], pl_arr[3];
    int var_cov, total_coverage,ref_coverage;
    qual_t qual;
    std::string basecount="";

    var_cov = 0;
    // for each snp event up to next variant or end
    SnpMap::iterator sv_it = _events->_snps.find(var_pos);
    if( sv_it != _events->_snps.end() ){
        for (auto &allele_it : sv_it->second) {
            size_t num_allele = allele_it.second._read_count;
            var_cov += int(num_allele);
            basecount.push_back(allele_it.second._allele_base);
            //basecount.push_back(':');
            basecount = basecount + ":"+std::to_string(num_allele)+",";
        }
    }

    ref_coverage = _coverages->get_coverage(var_pos, COVSIDX_DP);
    refbase = _sequences->_seq[var_pos];
    //var_cov = _coverages->get_coverage(var_pos, COVSIDX_ALL) - _coverages->get_coverage(var_pos, COVSIDX_DP);
    total_coverage = var_cov + ref_coverage;

    bcf_clear1(_vcf_rec);
    // set filter and genotype
    //bcf_update_filter(_vcf_hdr, _vcf_rec, &_pass_int, 1);

    if (total_coverage == 0) {
        //bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_NO_DATA],1);
        //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_NO_DATA]);
        gt_arr[0] = bcf_gt_missing;
        gt_arr[1] = bcf_gt_missing;
    } else {
        //if (total_coverage < _opts->min_var_dp) {
            //bcf_update_filter(_vcf_hdr, _vcf_rec, &_filter_idx[_FILTER_LOW_COVERAGE],1);
            //bcf_add_filter(_vcf_hdr, _vcf_rec, _filter_idx[_FILTER_LOW_COVERAGE]);
        //}
        gt_arr[0] = bcf_gt_unphased(0);
        gt_arr[1] = bcf_gt_unphased(0);
        //bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_VARIANTRATIO]);
    }

    basecount.push_back(refbase);
    basecount = basecount + ":" + std::to_string(ref_coverage);
    bcf_update_info_string(_vcf_hdr, _vcf_rec, "BaseNum", basecount.c_str());

    _vcf_rec->rid = _rid;
    _vcf_rec->pos = var_pos;
    _snp_alleles[0] = refbase;
    _snp_alleles[2] = '.';
    bcf_update_alleles_str(_vcf_hdr, _vcf_rec, _snp_alleles);
    //_vcf_rec->qual = 60;
    bcf_update_genotypes(_vcf_hdr, _vcf_rec, gt_arr, 2);

    int32_t ads[2]={(int32_t)ref_coverage, (int32_t)var_cov};
    tmpi = (int32_t)total_coverage;
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "AD", &ads, 2);
    bcf_update_format_int32(_vcf_hdr, _vcf_rec, "DP", &tmpi, 1);

    //std::cerr<<"==========NoVar end====================\n";
    //if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
    //    std::cerr << "Failed to write non-variant VCF record" << std::endl;
    //    exit(EXIT_FAILURE);
    //}
}
void VcfWriter::print_variant_buffer(std::string chr, int32_t chr_len){
    set_region(chr.c_str());
    for(int32_t i=0; i<chr_len; ++i){
        //std::cerr<<"Pos: "<<std::to_string(i)<<std::endl;
        bool var_flag = false;
        if (_opts->all_variant ){
            if( _events->_indels.find(i+1) != _events->_indels.end()  ){
                var_flag = print_indel_buffer(i+1, _events->_read_pos[i]);
                if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
                    std::cerr << "Failed to write InDel record in site: " <<std::to_string(i)<<std::endl;
                    exit(EXIT_FAILURE);
                }
                var_flag = true;
            }
            if( _events->_snps.find(i) != _events->_snps.end() ){
                var_flag = print_snp_buffer(i,_events->_read_pos[i]);
                if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
                    std::cerr << "Failed to write SNP record in site: " <<std::to_string(i)<<std::endl;
                    exit(EXIT_FAILURE);
                }
                var_flag = true;
            }
            if( !var_flag ){
                print_novar_buffer(i);
                if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
                    std::cerr << "Failed to write SNP record in site: " <<std::to_string(i)<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }else{
            if (i == 3105 || i== 3106 ){
                print_novar_buffer(i);
            }else if( _events->_indels.find(i+1) != _events->_indels.end() ){
                var_flag = print_indel_buffer(i+1,_events->_read_pos[i]);
                //int ndp_arr = 0;
                //int total_dp=0;
                //int32_t *dp_arr = NULL;
                //if( bcf_get_info_int32(_vcf_hdr, _vcf_rec, "DP", &dp_arr, &ndp_arr) > 0 ){
                //    int total_dp = *((int32_t*)dp_arr);
                //}
                //if( total_dp < _opts->min_site_dp && _events->_snps.find(i) != _events->_snps.end()){
                //        var_flag = print_snp_buffer(i);
                //}else
                if( var_flag &&  _events->_snps.find(i) != _events->_snps.end() ){
                    bool var_flag_snp = print_snp_buffer(i,_events->_read_pos[i]);
                    if( !var_flag_snp ){
                        var_flag = print_indel_buffer(i+1,_events->_read_pos[i]);
                    }
                }else if( !var_flag && _events->_snps.find(i) != _events->_snps.end() ){
                    var_flag = print_snp_buffer(i,_events->_read_pos[i]);
                    if( !var_flag ){
                        print_novar_buffer(i);
                    }
                }
            }else if( _events->_snps.find(i) != _events->_snps.end() ){
                var_flag = print_snp_buffer(i,_events->_read_pos[i]);
                if( !var_flag ){
                    print_novar_buffer(i);
                }
            } else {
                print_novar_buffer(i);
            }
            if (bcf_write1(_vcf_fp, _vcf_hdr, _vcf_rec)) {
                std::cerr << "Failed to write VCF record in site: " <<std::to_string(i)<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}
} // namepace mtVariant
