#ifndef _MTVARIANT_SNPEVENT_H
#define _MTVARIANT_SNPEVENT_H

#include<map>
#include<vector>

namespace mtVariant {

class SnpEvent
{
  private:
    SnpEvent();

  public:
    char _allele_base;
    int _read_count;
    int _lowq_count; // low quality count
    int _pos_strand;
    int _qual;
    double _rel_pos;
    double _nqs;
    std::map<int,int>_read_start;

    SnpEvent(char alt, bool rs, int qual, double rp, double nq, int read_start);
    ~SnpEvent();
    void add_snp_event(const SnpEvent &sv);
    double get_qual();
    double get_rel_pos();
    double get_nqs();
    int get_total_start();
    std::vector<int> get_all_starts();
    int get_start_count(int start);
    float get_lowqual_rate();
};

} // namepace mtVariant

#endif /* _MTVARIANT_SNPEVENT_H */
