#include <cstring>
#include <iostream>

#include "ReferenceSequence.h"

namespace mtVariant {

ReferenceSequence::ReferenceSequence()
    : _seq(nullptr),
      _len(0)
{
}

ReferenceSequence::ReferenceSequence(const char *refseq)
{
    if ((_fai = fai_load(refseq)) == nullptr) {
        std::cerr<<"Loading reference file:"<<refseq<<" failed."<<std::endl;
    } else {
        // loaded reference index successfully
        _seq = nullptr;
        _len = 0;
    }
}

ReferenceSequence::~ReferenceSequence()
{
    fai_destroy(_fai);
    std::free(_seq);
}

void ReferenceSequence::set_region(const char *region)
{
    int region_len;
    // free previous region
    std::free(_seq);
    // set current region
    _seq = fai_fetch(_fai, region, &region_len);
    _len = region_len;

    // set all reference sequence to upper case
    for (size_t i = 0; i < (size_t)region_len; ++i) {
        switch (_seq[i]) {
        case 'A': case 'C': case 'G': case 'T': case 'N':
            break;
        case 'a': case 'c': case 'g': case 't': case 'n':
            _seq[i] = std::toupper(_seq[i]);
            break;
        default:
            _seq[i] = 'N';
            break;
        }
    }
}

} // namepace mtVariant
