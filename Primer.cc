#include <cmath>
#include "Primer.h"

namespace mtVariant {

Primer::Primer()
{
    _primer_list.clear();
    _pfwd_map.clear();
    _prev_map.clear();
    //_tfwd_map.clear();
    //_trev_map.clear();
}
Primer::Primer(int distance,int read_mlen,std::string seq_type)
: _distance(distance)
,_read_mlen(read_mlen)
,_seq_type(seq_type)
{
    _primer_list.clear();
    _pfwd_map.clear();
    _prev_map.clear();
    //_tfwd_map.clear();
    //_trev_map.clear();
}
/**
void Primer::primer_cluster(int* site_depth,int size)
{
    std::vector<int> tmp_site;
    std::map<int,int> primer_list;
    int cut_off = _depth_cut / 2;
    int total = 0;
    for (int i=0; i< size; ++i ) {
        if( site_depth[i] > cut_off ){
            if( !tmp_site.empty() && (tmp_site.back() + _distance) > i ){ // overlap
               if( site_depth[tmp_site.back()] < site_depth[i] ) tmp_site.push_back(i);
            }else{
                tmp_site.push_back(i);
            }
        }
        total += site_depth[i];
    }

    for (int i=0; i < tmp_site.size(); ++i) {
        primer_list.emplace(tmp_site[i],0);
        int start = (tmp_site[i] < _distance) ? 0 : (tmp_site[i] - _distance);
        int end = ((tmp_site[i] + _distance) > size) ? size : (tmp_site[i]+_distance) ;
        for(int j=start;j<end;j++){
            primer_list[tmp_site[i]] += site_depth[j];
            _primer_map.emplace(j,tmp_site[i]);
        }
    }

    int left = 0;
    for(auto iter=primer_list.begin(); iter != primer_list.end(); ++iter ){
        if( iter->second < _depth_cut ){
            int start = (iter->first < _distance) ? 0 :  (iter->first - _distance);
            int end = ((iter->first + _distance)>size) ? size : (iter->first + _distance);
            for(int i=start;i<end;++i){
                auto it = _primer_map.find(i);
                if( it !=_primer_map.end() ) _primer_map.erase(it);
            }
        }else{
            left += iter->second;
        }
    }

    //finaly site
    fprintf(stderr,"Total reads:%d\n", total);
    fprintf(stderr, "Left reads:%d, %.4f\n", left, (float)left/(float)total);
}
**/

void Primer::read_chr_primer(std::string primer_file, std::string chr_name)
{
    //format: chr<tab>start<tab>primer_id
    const int LENS = 1024;
    int idx = 0;
    char buf[LENS];
    const char *delims = "\t \n";
    FILE *fp = fopen(primer_file.c_str(),"r");
    if( !fp ) std::cerr<<"Fail to open primer file:"<<primer_file<<std::endl;
    while( fgets(buf,LENS,fp) ){
        char* p_chr = strdup(strtok(buf,delims));
        std::string chr_s(p_chr, strlen(p_chr));
        int fwd_start =  atoi(strtok(NULL,delims));
        int fwd_end  =   atoi(strtok(NULL,delims));
        int rev_start =  atoi(strtok(NULL,delims));
        int rev_end   = atoi(strtok(NULL,delims));
        if( chr_s == chr_name ){
            if( fwd_start > fwd_end || rev_start > rev_end ) {
                std::cerr<<"The fwd/rev start > fwd/rev end, please check!"<<std::endl;
                exit(1);
            }
            if( fwd_start > rev_start ){
                int tmp_start = fwd_start;
                int tmp_end = fwd_end;
                fwd_start = rev_start;
                fwd_end = rev_end;
                rev_start = tmp_start;
                rev_end = tmp_end;
            }

            PrimerInfo tmp_pinfo(fwd_start,fwd_end, rev_start, rev_end);
            _primer_list.push_back(tmp_pinfo);

            // raw primer fwd start
            int start = (fwd_start < _distance) ? 0 : (fwd_start - _distance);
            int end = (fwd_start + _distance);
            if( _seq_type == "fwd" || _seq_type =="unknown" ){
                for(int i=start;i<=end;++i){
                    _pfwd_map.emplace(i, idx);
                }
                // trimed primer fwd end
                //start = (fwd_end < _distance) ? 0 : (fwd_end - _distance);
                //end = (fwd_end + _distance);
                //for(int i=start;i<end;++i){
                //    _tfwd_map.emplace(i, idx);
                //}

                // trimed primer rev start
                //start = (rev_start < _distance) ? 0 : (rev_start - _distance);
                //end = (rev_start + _distance);
                //for(int i=start;i<end;++i){
                //    _trev_map.emplace(i, idx);
                //}
            }
            if( _seq_type == "rev" || _seq_type =="unknown" ){
                // raw primer rev end
                start = (rev_end < _distance) ? 0 : (rev_end - _distance);
                end = (rev_end + _distance);
                for(int i=start;i<=end;++i){
                    _prev_map.emplace(i, idx);
                }
            }
            idx++;
        }
   }
   fclose(fp);

    if( idx == 0 ){
        std::cerr<<"No primer in chr:"<<chr_name<<", please check."<<std::endl;
        exit(EXIT_FAILURE);
    }
}
void Primer::reset(){
    _primer_list.clear();
    _pfwd_map.clear();
    _prev_map.clear();
    //_tfwd_map.clear();
    //_trev_map.clear();
}

int Primer::get_primer_idx(int start,int end)
{
    int idx = -1;
    //if( _tfwd_map.find(start) != _tfwd_map.end() ){
    //    idx = _tfwd_map.at(start);
    //}else if( _trev_map.find(end) != _trev_map.end() ){
    //    idx = _trev_map.at(end);
    //}else if( _pfwd_map.find(start) != _pfwd_map.end() ){
    if( _pfwd_map.find(start) != _pfwd_map.end() ){
        idx = _pfwd_map.at(start);
    }else if( _prev_map.find(end)  != _prev_map.end() ){
        idx = _prev_map.at(end);
    }
    return idx;
}

PrimerType Primer::check_primer(int start, int end, int read_len)
{
    PrimerType tmp_ptype = PrimerType::NoPrim;
    int map_len = (end - start + 1);
    if( _prev_map.find(end) != _prev_map.end() ){ // Found
        if(_pfwd_map.find(start) != _pfwd_map.end() ){ // Found
            //std::cerr<<"*****start:"<<std::to_string(start)<<",end:"<<std::to_string(end)<<"\n";
            //std::cerr<<"rev_map:"<<std::to_string(_prev_map.at(end))<<",start_map:"<<std::to_string(_pfwd_map.at(start))<<"\n";
            if( _prev_map.at(end) == _pfwd_map.at(start) ) { // Same index, read length == primer length
                tmp_ptype = PrimerType::BothPrim;
            }else{
                tmp_ptype = PrimerType::NoPrim;
            }
        }else{
            PrimerInfo tmp_pinfo = _primer_list.at(_prev_map[end]);
            //std::cerr<<"primer start:"<<std::to_string(tmp_pinfo.fwd_start)<<",end:"<<std::to_string(tmp_pinfo.rev_end)<<",insert:"<<std::to_string(tmp_pinfo.p_insert)<<std::endl;
            if( _read_mlen < tmp_pinfo.p_insert ){ // read length < primer length
                if( (map_len + _distance) < _read_mlen ){ // mapped length is too short
                    tmp_ptype = PrimerType::NoPrim;
                }else{
                    // due to sequence length ,fwd primer is not sequence.
                    tmp_ptype = PrimerType::RightPrim;
                }
            }else{// primer length < read length
                if( start < tmp_pinfo.fwd_start || map_len < tmp_pinfo.p_insert+ _distance ){
                    // read:    -------------------------------
                    // primer:  -------|----------------------|
                    tmp_ptype = PrimerType::BothPrim;
                }else{
                    tmp_ptype = PrimerType::NoPrim;
                }
            }
        }
    } else if( _pfwd_map.find(start) != _pfwd_map.end() ) { // Found
        if(_prev_map.find(end) != _prev_map.end() ){ // Found
            if( _prev_map.at(end) == _pfwd_map.at(start) ) { // Same index, read length == primer length
                tmp_ptype = PrimerType::BothPrim;
            }else{
                tmp_ptype = PrimerType::NoPrim;
            }
        }else{
            PrimerInfo tmp_pinfo = _primer_list.at(_pfwd_map[start]);
            if( _read_mlen < tmp_pinfo.p_insert ){ // read length < primer length
                if( (map_len + _distance) < _read_mlen ){ // mapped length is too short
                    tmp_ptype = PrimerType::NoPrim;
                }else{
                    // due to sequence length ,fwd primer is not sequence.
                    tmp_ptype = PrimerType::RightPrim;
                }
            }else{// primer length < read length
                if( end < tmp_pinfo.rev_end || map_len < tmp_pinfo.p_insert+ _distance ){
                    // read:    -------------------------------
                    // primer:  -------|----------------------|
                    tmp_ptype = PrimerType::BothPrim;
                }else{
                    tmp_ptype = PrimerType::NoPrim;
                }
            }
        }
    /**
    }else if( _tfwd_map.find(start) != _tfwd_map.end() ){
        if(_prev_map.find(end) != _prev_map.end() ){ // Found
            if( _tfwd_map.at(start) == _prev_map.at(end )){ // Same index
                tmp_ptype = PrimerType::RnTrimPrim;
            }else{// Trimed error
                //tmp_ptype = PrimerType::NoPrim;
                tmp_ptype = PrimerType::RnTrimPrim;
            }
        }else if( _trev_map.find(end) != _trev_map.end()){ // Found
            if( _tfwd_map.at(start) == _trev_map.at(end) ){ // Same Primer Index
                tmp_ptype = PrimerType::TrimPrimer;
            }else{// Trimed error
                //tmp_ptype = PrimerType::NoPrim;
                tmp_ptype = PrimerType::TrimPrimer;
            }
        }else{// No Found
            PrimerInfo tmp_pinfo = _primer_list.at(_tfwd_map[start]);
            if( _read_mlen <= tmp_pinfo.t_insert ){// read > insert
                tmp_ptype = PrimerType::TrimPrimer;
            }else{
                tmp_ptype = PrimerType::NoPrim;
            }
        }
    }else if( _trev_map.find(end) != _trev_map.end() ){
        PrimerInfo tmp_pinfo = _primer_list.at(_trev_map[end]);
        // two rev end
        
        if( _read_mlen <= tmp_pinfo.t_insert ){
            tmp_ptype = PrimerType::TrimPrimer;
        }else{
            tmp_ptype = PrimerType::NoPrim;
        }
    **/
    }else{
        tmp_ptype = PrimerType::NoPrim;
    }

    return tmp_ptype;
}

int Primer::get_primer_site(PrimerType primer_type, int site)
{
    int index=-1;
    if( primer_type == PrimerType::LeftPrim ){
        index = _pfwd_map[site];
    }else if( primer_type == PrimerType::RightPrim ){
        index = _prev_map[site];
    }
    return index;
}

} // namepace mtVariant
