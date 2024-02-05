#include "MapStat.h"

namespace mtVariant {

MapStat::MapStat(
        CoverageCounter &coverages,
        const opts_s *opts,
        std::string prefix)
    : _coverages(&coverages),
      _opts(opts),
      _prefix(prefix)
{
    //_cover_ptr = new int[_idx];
    //std::memset(_cover_ptr, 0, _idx * sizeof(int));
    for(int i=0;i<_idx;++i){
        _cover_ptr[i]=0;
    }
    _min_depth=600000;
    _max_depth=0;
    _total_site=0;
    _total_depth=0;
    _site_depth.clear();
    std::string statfile = _prefix+".map_stat.xls";
    fout = fopen(statfile.c_str(),"w");
    if( !fout ) std::cerr<<"Fail to open map stat file:"<<statfile<<std::endl;
    std::string depthfile = _prefix+".depth.xls.gz";
    gzout = gzopen(depthfile.c_str(),"wb");
    if( gzout == nullptr ) std::cerr<<"Fail to open depth file:"<<depthfile<<std::endl;
    gzprintf(gzout,"chr\tpos\tdepth\n");
}

MapStat::~MapStat()
{
    fclose(fout);
    gzclose(gzout);
}

void MapStat::calCoverage(std::string chr_name,int size){
    //int size = _coverages->get_length();
    int64_t total_depth = 0;
    //int tmp_size = (_site_depth.empty()) ? 0 : _site_depth.size();
    //>0x,>=10x,>=50x,>=100x,>=200x, >=500x, >=1000x,>=2000x, >=5000x
    for(int i=0;i<size;++i){
        int depth = _coverages->get_coverage(i, COVSIDX_ALL);
        _site_depth.push_back(depth);
        if( _min_depth > depth ) _min_depth = depth;
        if( _max_depth < depth ) _max_depth = depth;
        total_depth += depth;
        int end = 0;
        if( depth >= 5000 ){
            end = _idx;
        }else if( depth >= 2000){
            end = (_idx-1);
        }else if( depth >= 1000 ){
            end = (_idx - 2);
        }else if( depth >= 500 ){
            end = (_idx - 3);
        }else if( depth >= 200 ){
            end = (_idx - 4);
        }else if( depth >= 100 ){
            end = (_idx - 5);
        }else if( depth >= 50 ){
            end = (_idx - 6);
        }else if( depth >=10 ){
            end = (_idx - 7);
        }else{
            end = (depth > 0 ) ? 1 : 0;
        }
        for(int j=0;j<end;++j){
            _cover_ptr[j]++;
        }
    }
    _total_depth += total_depth;
    _total_site += size;

    if( size % 2 == 0 ){
        _median_depth = (_site_depth[size/2]+_site_depth[1+size/2])/2;
    }else{
        _median_depth = _site_depth[(size+1)/2];
    }
}
void MapStat::printStat() {
    //std::ifstream fin(outfile); read a file
    //FILE *fout = fopen(outfile.c_str(),"w");
    //if( fout ){
        fprintf(fout,"Term\tNumber\tRate(%)\n");
        fprintf(fout,"Raw Total Reads\t%d\t100.00\n", _opts->total_reads);
        fprintf(fout,"Fail and Unmapped reads\t%d\t%.2f\n\n", _opts->fail_unmap_reads, 100.0*_opts->fail_unmap_reads/(1.0*_opts->total_reads));

        fprintf(fout,"Mapped Reads\t%d\t%.2f\n",_opts->map_reads, 100.0*_opts->map_reads/(1.0*_opts->total_reads));
        if(_opts->paired_reads > 0 ){
            fprintf(fout,"Mapped and properly paired reads\t%d\t%.2f\n",_opts->paired_reads, 100.0*_opts->paired_reads/(1.0*_opts->total_reads));
        }
        fprintf(fout,"Non-primary mapped reads\t%d\t%.2f\n\n", _opts->no_primary, 100.0*_opts->no_primary/(1.0*_opts->total_reads));

        fprintf(fout,"Length < %d bp reads\t%d\t%.2f\n\n", _opts->read_len_min, _opts->read_short, 100.0*_opts->read_short/float(_opts->total_reads) );

        fprintf(fout,"Mapped quality cutoff value\t%d\n", _opts->min_mapq);
        fprintf(fout,"Mapping quality >=%d reads\t%d\t%.2f\n", _opts->min_mapq, _opts->mapqual_reads, 100.0*_opts->mapqual_reads/(1.0*_opts->total_reads));
        fprintf(fout,"Forward strand reads\t%d\t%.2f\n",_opts->fwd_reads,100.0*_opts->fwd_reads/(1.0*_opts->total_reads));
        fprintf(fout,"Reverse strand reads\t%d\t%.2f\n\n",_opts->rev_reads,100.0*_opts->rev_reads/(1.0*_opts->total_reads));

        fprintf(fout,"Length of region\t%d\n", _total_site);
        float aver_depth = 1.0*_total_depth/(1.0*_total_site);
        fprintf(fout,"Average depth\t%.2f\n", aver_depth);
        fprintf(fout,"Minimum depth\t%d\n",_min_depth);
        //int size = _site_depth.size();
        //if( size % 2 == 0 ){
        //    _median_depth = (_site_depth[size/2]+_site_depth[1+size/2])/2;
        //}else{
        //    _median_depth = _site_depth[(size+1)/2];
        //}
        fprintf(fout,"Median depth\t%d\n", _median_depth);
        fprintf(fout,"Maximum depth\t%d\n\n", _max_depth);

        int cutoff[]={1,10,50,100,200,500,1000,2000,5000};
        for(int i=0;i<_idx;++i){
                fprintf(fout,"Coverage (>=%dx)\t%d\t%.2f\n", cutoff[i],_cover_ptr[i],100.0*_cover_ptr[i]/(1.0*_total_site));
        }
    //}else{
    //    std::cerr<<"Fail to open file:"<<outfile<<std::endl;
    //}
    //fclose(fout);
}

void MapStat::printDepth(std::string chr_name, int size, bool reverse) {
    for(int i=0;i<size;++i){
        int depth = _coverages->get_coverage(i, COVSIDX_ALL);
        //int filter_depth = _coverages->get_coverage(i,COVSIDX_FT);
        //fprintf(fout, "%s\t%d\t%d\n", chr_name.c_str(), i, depth);
        int idx = ( reverse ) ? (size - i) : (i+1);
        gzprintf(gzout,"%s\t%d\t%d\n", chr_name.c_str(), idx, depth);
    }
}

void MapStat::printCumDepth() {
}

} // namepace mtVariant
