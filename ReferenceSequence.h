#ifndef _MTVARIANT_SEQUENCE_H
#define _MTVARIANT_SEQUENCE_H

#include "htslib/faidx.h"

namespace mtVariant {

class ReferenceSequence
{
  private:
    ReferenceSequence();

  public:
    faidx_t *_fai;
    char *_seq;
    int _len;

    ReferenceSequence(const char *refseq);
    ~ReferenceSequence();
    void set_region(const char *region);
};

} // namepace mtVariant

#endif /* _MTVARIANT_SEQUENCE_H */
