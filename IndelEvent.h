#ifndef _MTVARIANT_INDELEVENT_H
#define _MTVARIANT_INDELEVENT_H

#include <cstdint>
#include<map>
#include<vector>

#include "CoverageCounter.h"
#include "ReferenceSequence.h"

namespace mtVariant {

class IndelEvent
{
  private:
    int _map_qual;

    IndelEvent();

  public:
    int _var_start;
    int _var_len;
    int _q_pos;
    std::string _seq;
    std::string _id;
    int _read_count;
    int _lowq_count;
    bool _isdel;
    int _near_read_end_count;
    int _avg_nbq;
    double _var_rate_gap_and_mismatch;
    std::map<int,int>_read_start;

    IndelEvent(
        int start,
        int len,
        int q_pos,
        std::string &seq,
        bool strand,
        int mqual,
        int read_start,
        int nearend = 0,
        int avgnbq = 0,
        double vrate = 0.0);
    ~IndelEvent();
    void add_indel_event(const IndelEvent &iv);
    double get_near_read_end_ratio();
    double get_mean_avg_nqs();
    double get_mean_var_rate();
    int get_total_start();
    std::vector<int> get_all_starts();
    int get_start_count(int start);
    float get_lowqual_rate(){ return (float)_lowq_count / (float)_read_count; }
};

} // namepace mtVariant

#endif /* _MTVARIANT_INDELEVENT_H */
