#include <cmath>
#include <cstring>
#include <iostream>
#include "IndelEvent.h"

namespace mtVariant {

IndelEvent::IndelEvent()
    :  _map_qual(0),
      _var_start(0),
      _var_len(0),
      _q_pos(0),
      _seq(""),
      _id(""),
      _read_count(0),
      _lowq_count(0),
      _isdel(false),
      _near_read_end_count(0),
      _avg_nbq(0),
      _var_rate_gap_and_mismatch(0)
{
    _read_start.clear();
}

IndelEvent::IndelEvent(int start, int len, int q_pos,  std::string &seq, bool isdel, bool strand, int mqual, int read_start, int nearend, int avgnbq, double vrate)
     : _map_qual(mqual),
      _var_start(start),
      _var_len(len),
      _q_pos(q_pos),
      _seq(seq),
      _id(isdel ? std::to_string(len) : seq),
      _read_count(1),
      _isdel(isdel),
      _near_read_end_count(nearend),
      _avg_nbq(avgnbq),
      _var_rate_gap_and_mismatch(vrate)
{
    if( _avg_nbq < 20 ) _lowq_count=1;
    else _lowq_count=0;
    _read_start.emplace(read_start, 1);
}

IndelEvent::~IndelEvent() = default;

void IndelEvent::add_indel_event(const IndelEvent &iv)
{
    ++_read_count;
    _near_read_end_count += iv._near_read_end_count;
    _map_qual += iv._map_qual;
    _avg_nbq += iv._avg_nbq;
    if(iv._avg_nbq < 20 ) ++_lowq_count;
    _var_rate_gap_and_mismatch += iv._var_rate_gap_and_mismatch;

    if (!_read_start.empty()) {
        for (auto i = iv._read_start.begin(); i != iv._read_start.end(); ++i) {
            const auto &this_it = _read_start.find(i->first);
            if( this_it == _read_start.end()){ // No Found
                _read_start.emplace(i->first, i->second);
            }else{ // Found
                //this_it->second += i->second;
                _read_start[i->first] += i->second;
            }
        }
    }else{
        _read_start.insert(iv._read_start.begin(), iv._read_start.end());
        //_read_start = iv._read_start;
    }
}

double IndelEvent::get_near_read_end_ratio()
{
    return (double)_near_read_end_count / _read_count;
}

double IndelEvent::get_mean_avg_nqs()
{
    return (double)_avg_nbq / _read_count;
}

double IndelEvent::get_mean_var_rate()
{
    return (double)_var_rate_gap_and_mismatch / _read_count;
}
int IndelEvent::get_total_start()
{
    int total=0;
    int idx=0;
    std::string line="";
    for (auto i = _read_start.begin(); i != _read_start.end(); ++i) {
        total += i->second;
    }
    //for(auto i = _read_start.begin(); i != _read_start.end(); ++i){
    //    line += std::to_string(float(i->second)/float(total));
    //    if( idx < (_read_start.size()-1))  line +=",";
    //    idx++;
    //}
    return total;
}
std::vector<int> IndelEvent::get_all_starts()
{
    std::vector<int> _starts;
    for (auto i = _read_start.begin(); i != _read_start.end(); ++i) {
        _starts.push_back(i->first);
    }
    return _starts;
}
int IndelEvent::get_start_count(int start){
    int number = 0;
    const auto &pm_it = _read_start.find(start);
    if (pm_it != _read_start.end()) {// Found
        number = pm_it->second;
    }
    return number;
}
} // namepace mtVariant
