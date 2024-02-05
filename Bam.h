#ifndef _MTVARIANT_BAM_H
#define _MTVARIANT_BAM_H

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

namespace mtVariant {

class Bam
{
  private:
    hts_idx_t *_idx;

    Bam();

  public:
    samFile *_sf;
    hts_itr_t *_iter;
    bam_hdr_t *_hdr;

    explicit Bam(const char *sf_fn, const char *ref_fn);
    ~Bam();
    void set_iter(const char *where);
};

} // namepace mtVariant

#endif /* _MTVARIANT_BAM_H */
