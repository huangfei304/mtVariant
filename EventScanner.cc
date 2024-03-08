#include <algorithm>
#include <deque>

#include "EventScanner.h"

namespace mtVariant {

EventScanner::EventScanner()
{
}

EventScanner::EventScanner(float max_mis_rate, float max_gap_rate,int min_base_score, int min_read, bool primer_filter)
    : _max_sub(max_mis_rate)
    , _max_gap(max_gap_rate)
    , _base_score(min_base_score)
    , _read_min_len(min_read)
    , _primer_filter(primer_filter)
{
    _indels.clear();
    _snps.clear();
    _primers.clear();
    _max_chr_size = 16569; // chrM size
    _min_rlen_count = 0;
    _read_pos = new float[_max_chr_size];
    std::memset(_read_pos, 0, _max_chr_size * sizeof(float));
}

EventScanner::~EventScanner() = default;

void EventScanner::reset()
{
    _indels.clear();
    _snps.clear();
    _primers.clear();
    delete[] _read_pos;
}

int32_t EventScanner::left_align(const ReferenceSequence &refseq, std::string& seq, int32_t indel_pos, int32_t prev_ref_pos)
{
    int32_t align_pos = indel_pos;
    size_t len = seq.length();
    size_t len_1 = len - 1;
    uint32_t ref_idx = (uint32_t)indel_pos - 1;
    uint32_t seq_idx = 0;

    while ( align_pos > prev_ref_pos && refseq._seq[ref_idx] == seq[len_1 - seq_idx]) {
        --align_pos;
        --ref_idx;
        seq_idx = (seq_idx + 1) % len;
    }

    return align_pos;
}

int32_t EventScanner::right_align(const ReferenceSequence &refseq, std::string& seq, int indel_ref_start_pos, int read_last_pos, bool del ) {
    size_t len = seq.length();
    uint32_t ref_idx = ( del ) ? (uint32_t)(indel_ref_start_pos+len) : (uint32_t)(indel_ref_start_pos+1);
    int align_pos = indel_ref_start_pos;
    uint32_t seq_idx = 0;

    while ( ref_idx < read_last_pos && (refseq._seq[ref_idx] == seq[seq_idx] || refseq._seq[ref_idx] == 'N' ) ) {
        //if( indel_ref_start_pos == 3105 ) std::cerr<<"+++++++++++++"<<std::endl;
        ++align_pos;
        ++ref_idx;
        seq_idx = (seq_idx + 1) % len;
    }
    //if( indel_ref_start_pos == 3105 ) std::cerr<<"*********** pos:"<<std::to_string(align_pos)<<std::endl;
    return (del) ? align_pos : (align_pos+1);
    //return (ref_idx-1);
}

bool EventScanner::_isHVR( int start_pos, int end_pos ){
    //HVR-I 16024-16365, HVRII 73-340, HVRIII 340-576
    bool _isHVR = false;
    if( (start_pos >= 16024  && start_pos <= 16365) || (start_pos >= 73 && start_pos <= 576 ) ){
        _isHVR = true;
    }else if( (end_pos >= 16024 && end_pos <=16365) || (end_pos >= 73 && end_pos <= 576) ){
        _isHVR = true;
    }
    return _isHVR;
}
void EventScanner::collect_variant(bam1_t *read, const ReferenceSequence &refseq, CoverageCounter &coverages,Primer &primer) {
    bool strand;
    char refbase, allele, N;
    uint8_t *bseq, *bqual;
    uint32_t *cigar, c_op, c_type;
    int nm, gap, c_len, len, r_pos, q_pos, right_align_pos, var_count, last_match_ref, primer_idx;
    int mis_start, mis_end,gap_start, gap_end, q_pos_start, q_pos_end, start_idx, end_idx, q_start_idx, q_end_idx;
    std::string seq,read_name;
    AReadsSnpList this_reads_snps;
    AReadsIndelList this_reads_indels;
    AReadsPrimerList this_reads_primers;

    nm = 0; gap=0; q_pos = 0; var_count = 0; mis_start=0; mis_end=0;
    len = read->core.l_qseq;
    r_pos = read->core.pos;
    bseq = bam_get_seq(read);
    cigar = bam_get_cigar(read);
    q_pos_start = (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) ? bam_cigar_oplen(cigar[0])+25 : 25;
    q_pos_end = (bam_cigar_op(cigar[read->core.n_cigar-1]) == BAM_CSOFT_CLIP) ? len-26-bam_cigar_oplen(cigar[read->core.n_cigar-1]) : len-26;

    read_name = bam_get_qname(read);
    for (int i = 0; i < read->core.n_cigar; ++i) {
        c_op = bam_cigar_op(cigar[i]);
        c_len = bam_cigar_oplen(cigar[i]);
        switch (c_op) {
            case BAM_CMATCH: // M
                for(int j=0;j<c_len;++j){
                    if( seq_nt16_str[bam_seqi(bseq, (q_pos+j))] != refseq._seq[r_pos+j]){
                        if( (q_pos+j) < q_pos_start ) ++mis_start;
                        if( (q_pos+j) > q_pos_end ) ++mis_end;
                        ++nm;
                    }
                }
                break;
            case BAM_CEQUAL:
                break;
            case BAM_CDIFF: // X
                if( q_pos < q_pos_start ) ++mis_start;
                if( q_pos > q_pos_end ) ++mis_end;
                nm += c_len;
                break;
            case BAM_CDEL: // D
            case BAM_CINS: // I
                if( q_pos < q_pos_start ) ++mis_start;
                if( q_pos > q_pos_end ) ++mis_end;
                ++gap;
                break;
            case BAM_CSOFT_CLIP: // S
            default:
                break;
        }
        if ((bam_cigar_type(cigar[i]) & 0x01) != 0 ){
            q_pos += c_len;
        }
        if((bam_cigar_type(cigar[i]) & 0x02) !=0 ){
            r_pos += c_len;
        }
    }
    last_match_ref = r_pos;

    _filter_flag = 0;
    bool HVR_flag = _isHVR(read->core.pos+1, (1+last_match_ref));
    if( mis_start >=4 || mis_end >=4 ){
        _filter_flag = 100;
        return;
    }
    if( HVR_flag ){
        if( nm  > int(0.5+1.5*len*_max_sub) || gap > int(0.5+1.5*len*_max_gap) ){
            _filter_flag = 1;
            return;
        }
    } else {
        if( nm  > int(0.5+len*_max_sub) || gap > int(0.5+len*_max_gap) ){
            _filter_flag = 1;
            return;
        }
    }

    r_pos = read->core.pos;
    q_pos = 0;

    primer_idx = r_pos;
    PrimerType primer_type;
    start_idx = r_pos; end_idx = last_match_ref;
    if( _primer_filter ){
        primer_type = primer.check_primer((r_pos+1), (last_match_ref+1), q_pos);
        //if( read_name == "FT100023325L1C009R00100457498" ){
        //    std::cerr<<"mis_start:"<<std::to_string(mis_start)<<",mis_end:"<<std::to_string(mis_end)<<std::endl;
        //    std::cerr<<"-----read name:"<<read_name<<",start:"<<std::to_string(start_idx)<<",end:"<<std::to_string(end_idx)<<std::endl;
        //    std::cerr<<"r_pos:"<<std::to_string(r_pos)<<",last_pos:"<<std::to_string(last_match_ref)<<", type:"<<std::to_string(int(primer_type))<<std::endl;
        //}
        if (primer_type == PrimerType::NoPrim ){
            _filter_flag = 2;
            return;
        }
        primer_idx = primer.get_primer_idx((r_pos+1), (last_match_ref+1));
        if( primer_idx < 0 ){
            _filter_flag = 3;
            return;
        }
        // 1-base
        PrimerInfo tmp_pinfo = primer.get_prime_info(primer_idx);
        //if( primer_type == PrimerType::LeftPrim ){
        //    start_idx = tmp_pinfo.fwd_end; // convert to 0-base
        //}else if( primer_type == PrimerType::RightPrim ){
        //    end_idx = tmp_pinfo.rev_start - 1; // convert to 0-base
        //}else if( primer_type == PrimerType::BothPrim ){
            start_idx = tmp_pinfo.fwd_end;
            end_idx = tmp_pinfo.rev_start - 1;
        //}

        // reads end at the pos: 3107, trimmed the primer list
        if( end_idx == 3182 ) start_idx = 3107; //0-base
        //if( HVR_flag ){
            //std::cerr<<"read name:"<<read_name<<",start:"<<std::to_string(start_idx)<<",end:"<<std::to_string(end_idx)<<std::endl;
        //    if(  mis_start >=4 )  start_idx += 25;
        //    if( mis_end >=4 )  end_idx -= 25;
            //std::cerr<<"read name:"<<read_name<<",trimed start:"<<std::to_string(start_idx)<<",trimed end:"<<std::to_string(end_idx)<<std::endl;
        //}
        //if( read_name == "FT100038547L1C007R00500842346" ){
        //    std::cerr<<"********read name:"<<read_name<<",start:"<<std::to_string(start_idx)<<",end:"<<std::to_string(end_idx)<<std::endl;
        //    std::cerr<<"primer_idx:"<<std::to_string(primer_idx)<<",start:"<<std::to_string(start_idx)<<", end:"<<std::to_string(int(end_idx))<<std::endl;
        //}
        coverages.add_coverage(tmp_pinfo.index, COVSIDX_PRIMER);
    }

    if( (end_idx - start_idx+1) < _read_min_len ) {
        _filter_flag = 4;
        ++_min_rlen_count;
        return;
    }
    bseq = bam_get_seq(read);
    bqual = bam_get_qual(read);
    std::vector<std::pair<int32_t, int32_t>> coverage_rr_vec, coverage_del_vec;
    std::vector<int> read_vec;
    q_start_idx=0; q_end_idx=len;
    size_t found;

    strand = (bam_is_rev(read) == 0);
    r_pos = read->core.pos;
    q_pos = 0;
    for (int i = 0; i < read->core.n_cigar; ++i) {
        c_op = bam_cigar_op(cigar[i]);
        c_len = bam_cigar_oplen(cigar[i]);
        switch (c_op) {
            case BAM_CINS: // I
                if( r_pos < start_idx || r_pos > end_idx ) break;
                ++var_count;
                seq.clear();
                for(int j=0;j<c_len;++j){
                    seq.push_back(seq_nt16_str[bam_seqi(bseq, q_pos+j)]);
                }
                //found = seq.find(N);
                //if (found == std::string::npos ) {
                    right_align_pos = right_align(refseq, seq, r_pos, last_match_ref,false);
                    //if( read_name == "FT100023325L1C012R00402258842" ){
                    //     std::cerr<<"read_name:"<<read_name<<",ref_base:"<<refseq._seq[r_pos]<<",after base:"<<refseq._seq[right_align_pos]<<std::endl;
                    //    std::cerr<<"r_pos:"<<std::to_string(r_pos)<<",last_match_ref:"<<std::to_string(last_match_ref)<<",right_pos:"<<std::to_string(right_align_pos)<<",seq:"<<seq<<std::endl;
                    //}
                    //rotate seq if different align position
                    if(  right_align_pos > (r_pos+1) ){
                        std::rotate(seq.begin(), seq.begin() + c_len - (right_align_pos - r_pos-1) % c_len, seq.end());
                    }
                    //std::cerr<<"insert read name:"<<read_name<<",pos:"<<std::to_string(right_align_pos)<<std::endl;
                    if( right_align_pos >= _max_chr_size ) right_align_pos -=  _max_chr_size;
                    //if(right_align_pos == 1890 ){
                    //    std::cerr<<"insert read name:"<<read_name<<std::endl;
                    //}
                    this_reads_indels.push_back(
                        IndelEvent(right_align_pos, c_len, q_pos, seq, false, strand, read->core.qual, primer_idx)
                    );
                //}
                break;
            case BAM_CDEL: // D
                if(  r_pos < start_idx || r_pos > end_idx ) break;
                ++var_count;
                seq.clear();
                for(int j=0;j<c_len;++j){
                    seq.push_back(refseq._seq[r_pos+j]);
                }
                //found = seq.find(N);
                //if (found == std::string::npos) {
                    right_align_pos = right_align(refseq, seq, r_pos, last_match_ref, true);
                    //if( read_name == "FT100023315L1C007R00102147106" ){
                    //    std::cerr<<"read_name:"<<read_name<<",ref_base:"<<refseq._seq[r_pos]<<",after base:"<<refseq._seq[right_align_pos]<<std::endl;
                    //    std::cerr<<"r_pos:"<<std::to_string(r_pos)<<",last_match_ref:"<<std::to_string(last_match_ref)<<",right_pos:"<<std::to_string(right_align_pos)<<",seq:"<<seq<<std::endl;
                    //}
                    //seq.clear();
                    coverage_del_vec.push_back(std::make_pair(right_align_pos, c_len));
                    if( right_align_pos > r_pos ){
                        coverage_rr_vec.push_back(std::make_pair(r_pos, (right_align_pos-r_pos)));
                        read_vec.push_back(q_pos);
                        std::rotate(seq.begin(), seq.begin() + c_len - (right_align_pos - r_pos) % c_len, seq.end());
                    }
                    if( right_align_pos >= _max_chr_size ) right_align_pos -=  _max_chr_size;
                    //if(right_align_pos == 1890 ){
                    //    std::cerr<<"indel read name:"<<read_name<<std::endl;
                    //}
                    this_reads_indels.push_back(
                        IndelEvent(right_align_pos, c_len, q_pos, seq, true, strand, read->core.qual, primer_idx)
                    );
                //}
                break;
            case BAM_CMATCH: // M
            //case BAM_CREF_SKIP: // N
            case BAM_CEQUAL: // =
            case BAM_CDIFF: // X
                if( _primer_filter ){
                    int tmp_start = (r_pos > start_idx ) ? r_pos : start_idx;
                    int tmp_end = ( (r_pos+c_len) < end_idx) ? (r_pos + c_len) : end_idx;
                    if( tmp_start >= r_pos && tmp_end <= (r_pos + c_len)){
                        coverage_rr_vec.push_back(std::make_pair(tmp_start, (tmp_end-tmp_start)));
                        read_vec.push_back(q_pos+(tmp_start - r_pos));
                    }
                }else{
                    coverage_rr_vec.push_back(std::make_pair(r_pos, c_len));
                    read_vec.push_back(q_pos);
                }
                break;
            case BAM_CSOFT_CLIP: // S
            default:
                break;
        }

        if((bam_cigar_type(cigar[i]) & 0x02) !=0 ){
            if( start_idx >= r_pos && start_idx <= (r_pos+c_len) ){
                if ( (bam_cigar_type(cigar[i]) & 0x01) == 0 ) q_start_idx = q_pos;
                else q_start_idx = q_pos + (start_idx - r_pos);
            }
            if( end_idx >= r_pos && end_idx <= (r_pos+c_len) ){
                if ( (bam_cigar_type(cigar[i]) & 0x01) == 0 ) q_end_idx = q_pos;
                else q_end_idx = q_pos + (end_idx - r_pos);
            }
            r_pos += c_len;
            if( r_pos > end_idx ) break;
        }
        if ((bam_cigar_type(cigar[i]) & 0x01) != 0 ){
            q_pos += c_len;
        }
    }

    // add RR coverage, adjusting for deletions
    if ( !coverage_rr_vec.empty() ){
        for(int i=0;i<coverage_rr_vec.size();++i){
            int start = coverage_rr_vec[i].first;
            int end = start + coverage_rr_vec[i].second;

            bool flag = false;
            if( !coverage_del_vec.empty() ){
                for (int j=0; j<coverage_del_vec.size();j++) {
                    int del_start = coverage_del_vec[j].first;
                    int del_end = del_start + coverage_del_vec[j].second;
                    if ( del_end < start ) continue;
                    if( del_start > end ) break;
                    if( del_end >= start && del_end < end ){
                        flag = true;
                        int r_idx = read_vec[i] + (del_end - start);
                        for(int k=del_end; k<end; ++k){
                            if( r_idx >= len ) break;
                            refbase = refseq._seq[k];
                            allele = seq_nt16_str[bam_seqi(bseq, r_idx)];
                            int tmp_k = (k < _max_chr_size ) ? k : (k - _max_chr_size);
                            int _is_in_end = ( (r_idx - q_start_idx) < _near_end || (q_end_idx - r_idx ) < _near_end ? 1 : 0);
                            if( (int)bqual[r_idx] < _base_score ){
                                r_idx++;
                                continue;
                            }
                            if( refbase == allele ){
                                if( strand ) coverages.add_coverage(tmp_k, COVSIDX_RS);
                                if( (int)bqual[r_idx] < _base_score ) coverages.add_coverage(tmp_k, COVSIDX_LQ);
                                coverages.add_coverage(tmp_k, COVSIDX_DP);
                                coverages.add_coverage(tmp_k, COVSIDX_ALL);
                                //_read_pos[tmp_k] += (float)(r_idx)/(float)len;
                                _read_pos[tmp_k] += _is_in_end;
                                //primer_add_count(k, primer_idx);
                                this_reads_primers.push_back(std::make_pair(tmp_k, primer_idx));
                                //if( tmp_k == 5737 || tmp_k == 5759 || tmp_k == 10667 ){ // 0-base
                                if( tmp_k == 9307   ){ // 0-base
                                   std::cerr<<"pos: "<<std::to_string(tmp_k)<<",read name:"<<read_name<<",read_idx:"<<std::to_string(r_idx)<<",base:"<<allele<<", base quality:"<<std::to_string(bqual[r_idx])<<std::endl;
                                }
                            }else{
                                ++var_count;
                                coverages.add_coverage(tmp_k, COVSIDX_ALL);
                                //double relpos = (double)(r_idx)/(double)len;
                                double relpos = _is_in_end;
                                if( tmp_k == 9307 ){ // 0-base
                                    std::cerr<<"pos: "<<std::to_string(tmp_k)<<",read name:"<<read_name<<",read_idx:"<<std::to_string(r_idx)<<",base:"<<allele<<", base quality:"<<std::to_string(bqual[r_idx])<<std::endl;
                                }
                                this_reads_snps.push_back(std::make_pair(
                                    tmp_k,
                                    SnpEvent(allele, strand, (int)bqual[r_idx], relpos, this->snp_nqs(r_idx, len, bqual),primer_idx)
                                ));
                                ++var_count;
                            }
                            r_idx++;
                        }
                    }
                }
            }
            if( !flag ){
                int r_idx = read_vec[i];
                for(int k=start; k<end; ++k){
                    refbase = refseq._seq[k];
                    allele = seq_nt16_str[bam_seqi(bseq, r_idx)];
                    int tmp_k = (k < _max_chr_size ) ? k : (k - _max_chr_size);
                    int _is_in_end = ( (r_idx - q_start_idx) < _near_end || (q_end_idx - r_idx ) < _near_end ? 1 : 0);
                    if( (int)bqual[r_idx] < _base_score ){
                        r_idx++;
                        continue;
                    }
                    if( refbase == allele ){
                        if( strand ) coverages.add_coverage(tmp_k, COVSIDX_RS);
                        if( (int)bqual[r_idx] < _base_score ) coverages.add_coverage(tmp_k, COVSIDX_LQ);
                        coverages.add_coverage(tmp_k, COVSIDX_DP);
                        coverages.add_coverage(tmp_k, COVSIDX_ALL);
                        //_read_pos[tmp_k] += (float)(r_idx)/(float)len;
                        _read_pos[tmp_k] += _is_in_end;
                        //primer_add_count(k, primer_idx);
                        this_reads_primers.push_back(std::make_pair(tmp_k, primer_idx));
                        //if( tmp_k == 3106  ){ // 0-base
                        if( tmp_k == 9307  ){ // 0-base
                            std::cerr<<"pos: "<<std::to_string(tmp_k)<<",read name:"<<read_name<<",read_idx:"<<std::to_string(r_idx)<<",base:"<<allele<<", base quality:"<<std::to_string(bqual[r_idx])<<std::endl;
                        }
                    }else{
                        ++var_count;
                        coverages.add_coverage(tmp_k, COVSIDX_ALL);
                        //double relpos = (double)(r_idx)/(double)len;
                        double relpos = _is_in_end;
                        if( tmp_k == 9307  ){// 0-base
                            std::cerr<<"pos: "<<std::to_string(tmp_k)<<",read name:"<<read_name<<",read_idx:"<<std::to_string(r_idx)<<",base:"<<allele<<", base quality:"<<std::to_string(bqual[r_idx])<<std::endl;
                        }
                        this_reads_snps.push_back(std::make_pair(
                            tmp_k,
                            SnpEvent(allele, strand, (int)bqual[r_idx], relpos, this->snp_nqs(r_idx, len, bqual),primer_idx)
                        ));
                        ++var_count;
                    }
                    r_idx++;
                }
            }
        }
    }

    // add each snp
    if( !this_reads_snps.empty()){
        for (const auto &sv : this_reads_snps) {
            const auto &sv_it = _snps.find(sv.first);
            if (sv_it == _snps.end()) {
                _snps[sv.first].insert(std::make_pair(sv.second._allele_base, sv.second));
            } else {
                const auto &svg_it = sv_it->second.find(sv.second._allele_base);
                if (svg_it == sv_it->second.end()) {
                    sv_it->second.insert(std::make_pair(sv.second._allele_base, sv.second));
                } else {
                    svg_it->second.add_snp_event(sv.second);
                }
            }
        }
    }

    // add each indel
    if (!this_reads_indels.empty()) {
        double var_rate_gap_and_mismatch = (double)var_count / (q_end_idx - q_start_idx);
        for (auto &iv_it : this_reads_indels) {
            int32_t end_pos = iv_it._q_pos;
            if (!iv_it._seq.empty()) {
                end_pos += iv_it._var_len;
            }

            int32_t pos1 = iv_it._q_pos - _near_end;
            if (pos1 < 0) {
                pos1 = 0;
            }

            int32_t pos2 = end_pos + _near_end;
            if (pos2 >= q_end_idx) {
                pos2 = q_end_idx - 1;
            }

            int qual_sum=0;
            for (q_pos = pos1; q_pos <= pos2; ++q_pos) {
                qual_sum += (int)bqual[q_pos];
            }

            //iv_it._near_read_end_count = (iv_it._q_pos + 1 <= _near_end || len - end_pos < _near_end + 2 ? 1 : 0);
            iv_it._near_read_end_count = ( (iv_it._q_pos + 1 - q_start_idx) <= _near_end || (q_end_idx - end_pos) < _near_end ? 1 : 0);
            iv_it._avg_nbq = qual_sum / (pos2 - pos1 + 1);
            iv_it._var_rate_gap_and_mismatch = var_rate_gap_and_mismatch;

            if( iv_it._avg_nbq < _base_score ) continue;
            const auto &iv_search = _indels.find(iv_it._var_start);
            if (iv_search == _indels.end()) {
                _indels[iv_it._var_start].insert(std::make_pair(iv_it._id, iv_it));
            } else {
                const auto &ivg_it = iv_search->second.find(iv_it._id);
                if (ivg_it == iv_search->second.end()) {
                    iv_search->second.insert(std::make_pair(iv_it._id, iv_it));
                } else {
                    ivg_it->second.add_indel_event(iv_it);
                }
            }
        }
    }
    //
    if( _primer_filter && !this_reads_primers.empty() ){
        for (const auto &sv : this_reads_primers ) {
            const auto &sv_it = _primers.find(sv.first);
            if (sv_it == _primers.end()) {
                _primers[sv.first].insert(std::make_pair(sv.second, 1));
            } else {
                const auto &svg_it = sv_it->second.find(sv.second);
                if (svg_it == sv_it->second.end()) {
                    sv_it->second.insert(std::make_pair(sv.second, 1));
                } else {
                    ++svg_it->second;
                }
            }
        }
    }
}

/**
 * Sum of base quality score within range [p-5, p+4]
 */
inline double EventScanner::snp_nqs(int p, int qlen, const uint8_t *bqual)
{
    int start, end;
    int qual_sum;

    start = (p < 5) ? 0 : p - 5;
    end = (p > qlen - 4) ? qlen : p + 4;
    qual_sum = 0;

    for (int qpos = start; qpos <= end; ++qpos) {
        qual_sum += (int)bqual[qpos];
    }
    return (double)qual_sum /(double)(end-start+1); // NOTE adjust for near read ends
}

int EventScanner::get_primer_total(int pos){
    int total = 0;
    const auto &pm = _primers.find(pos);
    if (pm != _primers.end()) { // Found
        for( const auto pm_it: pm->second){
            total += pm_it.second;
        }
    }
    return total;
}
std::vector<int> EventScanner::get_all_starts(int pos){
    std::vector<int> _starts;
    const auto &pm = _primers.find(pos);
    if (pm != _primers.end()) { // Found
        for( const auto pm_it: pm->second){
            _starts.push_back(pm_it.first);
        }
    }
    return _starts;
}
int EventScanner::get_primer_start_count(int pos, int start){
    int number = 0;
    const auto &pm = _primers.find(pos);
    if (pm != _primers.end()) { // Found
        const auto &pm_it = pm->second.find(start);
        if (pm_it != pm->second.end()) {// Found
            number = pm_it->second;
        }
    }
    return number;
}

} // namepace mtVariant
