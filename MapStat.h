#ifndef _MTVARIANT_MAPSTAT_H
#define _MTVARIANT_MAPSTAT_H

#include "CoverageCounter.h"
#include "mtVariant.h"
#include <iostream>
#include <fstream>
#include <zlib.h>

namespace mtVariant {
const int _idx = 9;
class MapStat
{
    private:
        //>0x,>=10x,>=50x,>=100x,>=200x, >=500x, >=1000x,>=2000x, >=5000x
        //int cutoff[]={1,10,50,100,200,500,1000,2000,5000}
        //int *_cover_ptr;
        int _cover_ptr[_idx];
        int _min_depth;
        int _max_depth;
        int _median_depth;
        int _total_site;
        int64_t _total_depth;
        std::string _prefix;
        FILE *fout;
        //FILE *gzout;
        gzFile gzout = nullptr;

        std::vector<int> _site_depth;

        const opts_s *_opts;
        CoverageCounter *_coverages;

    public:
        void calCoverage(std::string chr_name, int chr_len);
        MapStat(CoverageCounter &coverages,const opts_s *opts, std::string prefix);
        ~MapStat();
        void printStat();
        void printDepth(std::string chr_name, int size,bool reverse=false);
        int get_median_depth() {return _median_depth;}
        void printCumDepth();
};

} // namepace mtVariant

#endif /* _MTVARIANT_MAPSTAT_H */
