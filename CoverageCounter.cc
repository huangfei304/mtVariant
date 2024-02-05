#include <cstring>

#include "CoverageCounter.h"

namespace mtVariant {

CoverageCounter::CoverageCounter()
    : _covs(nullptr),
      _len(0),
      _max_cov(0)
{
}

CoverageCounter::CoverageCounter(int len, int max_cov)
    : _covs(new int[len][COVS_ALIGN]),
      _len(len),
      _max_cov(max_cov)
{
    this->reset();
}

CoverageCounter::CoverageCounter(const CoverageCounter &) = default;

CoverageCounter::~CoverageCounter()
{
    delete[] _covs;
}

void CoverageCounter::reset()
{
    std::memset(_covs, 0, COVS_ALIGN * _len * sizeof(int));
}

void CoverageCounter::add_coverage(
        int pos,
        covs_idx_e idx)
{
    size_t p = (size_t)pos;

    if (p < _len && _covs[p][idx] < _max_cov) {
        ++_covs[p][idx];
    }
}

void CoverageCounter::add_coverage(
        int pos,
        covs_idx_e idx,
        int add_cov)
{
    size_t p = (size_t)pos;

    if (p < _len && _covs[p][idx] + add_cov < _max_cov) {
        _covs[p][idx] += add_cov;
    }
}

void CoverageCounter::add_coverage_range(
        int pos,
        covs_idx_e idx,
        int len)
{
    size_t p = (size_t)pos;
    size_t end = p + (size_t)len;
    if (end > _len) {
        end = _len;
    }

    while (p < end) {
        if (_covs[p][idx] < _max_cov) {
            ++_covs[p][idx];
        }
        ++p;
    }
}

void CoverageCounter::add_coverage_range(
        int pos,
        covs_idx_e idx,
        int len,
        int add_cov)
{
    size_t p = (size_t)pos;
    size_t end = p + (size_t)len;
    if (end > _len) {
        end = _len;
    }

    while (p < end) {
        if (_covs[p][idx] + add_cov <= _max_cov) {
            _covs[p][idx] += add_cov;
        }
        ++p;
    }
}

void CoverageCounter::set_coverage(
        int pos,
        covs_idx_e idx,
        int cov)
{
    size_t p = (size_t)pos;

    if (p < _len && cov <= _max_cov) {
        _covs[p][idx] = cov;
    }
}

int CoverageCounter::get_coverage(
        int pos,
        covs_idx_e idx)
{
    size_t p = (size_t)pos;

    return (p < _len) ? _covs[p][idx] : 0;
}

bool CoverageCounter::test_high_coverage(
        int pos,
        covs_idx_e idx,
        int l_qseq)
{
    bool rv = true;
    size_t p = (size_t)pos;

    if (p + l_qseq / 2 >= _len ||
        _covs[p][idx] < _max_cov ||
        _covs[p + l_qseq][idx] < _max_cov ||
        _covs[p + l_qseq / 2][idx] < _max_cov)
    {
        rv = false;
    }

    return rv;
}

int CoverageCounter::get_max_coverage()
{
    return _max_cov;
}
int CoverageCounter::get_length()
{
    return _len;
}
} // namepace mtVariant
