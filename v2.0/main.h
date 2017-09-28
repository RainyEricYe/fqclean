/*
 *  yerui@genomics.cn
 *  2012-10-31
 *  2013-04-09
 *  2014-05-08
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif

#include "fq_read.h"
#include "gzstream.h"

using namespace __gnu_cxx;
using namespace std;

#define VERSION "v2.0"

typedef unsigned long ulong;
typedef unsigned int uint;
typedef map<int, ulong> mIntLong;

void usage(char *prog) {
    cout << "Program: fqclean   filtrate/trim pair-end reads with adapter or with too many low-quality-base\n"
            "Version: " << VERSION << "\n"
            "Contact: yerui <yerui@genomics.cn>\n\n"
            "Usage: " << prog << " in1.fq[.gz] in2.fq[.gz] out1.fq.gz out2.fq.gz\n\n"

            "       -n [f]    read in which rate of N > [f] will be discarded. [0.1]\n"
            "       -q [f]    read in which rate of low-quality-base > [f] will be discarded. [0.5]\n"
            "       -c [c]    phred quality < [c] is treated as low quality base. [ 'F' ]\n"
            "       -p [i]    [i] reads as a part to output. [0]\n"
            "       -m [i]    minimum length of reads after cleaning. [30]\n"
            "       -u [s]    UID types. ( NNNNNACT  means 5bp barcode + 3bp fix\n"
            "                              ACGNNNNN  meads 3bp fix + 5bp barcode\n"
            "                              NNNNNNNN  meads 8bp barcode. )\n"
            "       -t [i]    trim first N bp of reads (if has UID, trim N bp after UID) [0]\n"

            "\n" << endl;
}

template<typename M> string mtoa(M &i)
{
    ostringstream a;
    a << i;
    return a.str();
}

#endif // MAIN_H_
