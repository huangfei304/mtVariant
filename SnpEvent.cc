#include <cmath>
#include "SnpEvent.h"

namespace mtVariant {

SnpEvent::SnpEvent()
    : _allele_base(0),
      _read_count(0),
      _lowq_count(0),
      _pos_strand(0),
      _qual(0),
      _rel_pos(0.0),
      _nqs(0.0)
{
    _read_start.clear();
}

SnpEvent::SnpEvent(char alt,
                   bool rs,
                   int qual,
                   double rp,
                   double nq,
                   int read_start)
    : _allele_base(alt),
      _read_count(0),
      _lowq_count(0),
      _pos_strand(0),
      _qual(qual),
      _rel_pos(rp),
      _nqs(nq)
{
    ++_read_count;
    if( _qual < 20 ) ++_lowq_count;
    if (rs) {
        ++_pos_strand;
    }
    _read_start.emplace(read_start, 1);
}

SnpEvent::~SnpEvent() = default;

void SnpEvent::add_snp_event(const SnpEvent &sv)
{
    _read_count += sv._read_count;
    _lowq_count += sv._lowq_count;
    _pos_strand += sv._pos_strand;
    _qual += sv._qual;
    _rel_pos += sv._rel_pos;
    _nqs += sv._nqs;
    if (!_read_start.empty()) {
        for (auto i = sv._read_start.begin(); i != sv._read_start.end(); ++i) {
           const auto &this_it = _read_start.find(i->first);
            if( this_it == _read_start.end()){ // No Found
                _read_start.emplace(i->first, i->second);
            }else{ // Found
                this_it->second += i->second;
            }
        }
    }else{
        _read_start = sv._read_start;
    }
}

double SnpEvent::get_qual()
{
    return (double)_qual / (double)_read_count;
}

double SnpEvent::get_rel_pos()
{
    //return 1.0 - std::abs(0.5 - (double)_rel_pos / (double)_read_count);
    return (double)_rel_pos / (double)_read_count;
}

double SnpEvent::get_nqs()
{
    return (double)_nqs / (double)_read_count;
}

int SnpEvent::get_total_start()
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
std::vector<int> SnpEvent::get_all_starts()
{
    std::vector<int> _starts;
    for (auto i = _read_start.begin(); i != _read_start.end(); ++i) {
        _starts.push_back(i->first);
    }
    return _starts;
}
int SnpEvent::get_start_count(int start){
    int number = 0;
    const auto &pm_it = _read_start.find(start);
    if (pm_it != _read_start.end()) {// Found
        number = pm_it->second;
    }
    return number;
}

float SnpEvent::get_lowqual_rate(){
    return (float)_lowq_count / (float)_read_count;
}

} // namepace mtVariant
