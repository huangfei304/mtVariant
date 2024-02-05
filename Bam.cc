#include <cstring>
#include <cstdlib>
#include <iostream>

#include "Bam.h"
namespace mtVariant {

Bam::Bam()
    : _idx(nullptr),
      _sf(nullptr),
      _iter(nullptr),
      _hdr(nullptr)
{
}

Bam::Bam(const char *sf_fn, const char *ref_fn)
    : _idx(nullptr),
      _sf(nullptr),
      _iter(nullptr),
      _hdr(nullptr)
{
    size_t fn_strlen = std::strlen(sf_fn);
    size_t ref_fn_strlen = std::strlen(ref_fn);

    if (strcmp(sf_fn + fn_strlen - 4, ".bam") == 0) {
        _sf = sam_open(sf_fn, "rb");
    } else if (strcmp(sf_fn + fn_strlen - 5, ".cram") == 0) {
        _sf = sam_open(sf_fn, "rc");

        char *fai_fn = new char[ref_fn_strlen + 5];
        std::strncpy(fai_fn, ref_fn, ref_fn_strlen + 1);
        fai_fn[ref_fn_strlen] = '\0';
        std::strncat(fai_fn, ".fai", 5);

        if( !hts_set_fai_filename(_sf, fai_fn) ){
            std::cerr << "Failed to load reference index for \"" << ref_fn << "\"" << std::endl;
            exit(EXIT_FAILURE);
        }
        delete[] fai_fn;
    } else {
        // load as BAM format
        _sf = sam_open(sf_fn, "rb");
    }

    if (_sf == nullptr) {
        std::cerr << "Failed to load alignment file \"" << sf_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
    } else if ((_hdr = sam_hdr_read(_sf)) == nullptr) {
        std::cerr << "Failed to load header for \"" << sf_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
    } else if ((_idx = sam_index_load(_sf, sf_fn)) == nullptr) {
        std::cerr << "Failed to load index for \"" << sf_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }
}

Bam::~Bam()
{
    hts_itr_destroy(_iter);
    hts_idx_destroy(_idx);
    bam_hdr_destroy(_hdr);
    sam_close(_sf);
}

void Bam::set_iter(const char *region_str)
{
    hts_itr_destroy(_iter);
    _iter = sam_itr_querys(_idx, _hdr, region_str);
}

} // namepace mtVariant
