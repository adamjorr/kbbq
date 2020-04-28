#ifndef __KBBQ_KSEQ_H
#define __KBBQ_KSEQ_H

#include <htslib/bgzf.h>
#include <htslib/kseq.h>

namespace kseq{
	KSEQ_INIT(BGZF*, bgzf_read)
}

#endif