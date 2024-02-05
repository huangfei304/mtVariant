#ifndef _MTVARIANT_PRIMER_H
#define _MTVARIANT_PRIMER_H

#include<map>
#include<vector>
#include<iostream>
#include<cstring>

namespace mtVariant {

enum class PrimerType {
    LeftPrim=0, // left primer need trim
    RightPrim, // right primer need trim
    BothPrim, // both primer need trim
    //PrimShort, // shorter than expected
    //PrimLong, // longer than expected
    NoPrim   // no in primer or primer error
};

struct PrimerInfo {
    PrimerInfo()
    {
    }
    PrimerInfo(int f_start, int f_end, int r_start, int r_end)
        : fwd_start(f_start)
        , fwd_end(f_end)
        , rev_start(r_start)
        , rev_end(r_end)
    {
        p_insert =  (rev_end - fwd_start + 1);
        t_insert = (rev_start - fwd_end - 1);
    }

    int fwd_start;
    int fwd_end;
    int rev_start;
    int rev_end;

    int p_insert;
    int t_insert;
};

typedef std::vector<PrimerInfo> _PrimerList;
class Primer
{
  private:
    Primer();
    std::map<int, int> _pfwd_map;
    std::map<int, int> _prev_map;
    //std::map<int, int> _tfwd_map;
    //std::map<int, int> _trev_map;

  public:
    int _distance=6;
    int _read_mlen;
    _PrimerList _primer_list;


    Primer(int distance,int read_mlen);
    void reset();
    //Primer(std::string primer_file,std::string chr_name);
    //void primer_cluster(int *_site_depth,int size);
    //Primer(std::map<int,int> site_depth, int depth_cut,int distance=10);
    void read_chr_primer(std::string primer_file,std::string chr_name);
    int get_primer_idx(int start, int end);
    PrimerType check_primer(int start, int end, int read_len);
    int get_primer_site(PrimerType primer_type, int site);
    PrimerInfo get_prime_info(int idx) { return _primer_list[idx]; }
};

} // namepace mtVariant

#endif /* _MTVARIANT_PRIMER_H */
