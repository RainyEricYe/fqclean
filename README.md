# fqclean
Program: fqclean   filtrate/trim pair-end reads with adapter or with too many low-quality-base
Version: v2.0
Contact: yerui <yerui@genomics.cn>

Usage: fqclean in1.fq[.gz] in2.fq[.gz] out1.fq.gz out2.fq.gz

       -n [f]    read in which rate of N > [f] will be discarded. [0.1]
       -q [f]    read in which rate of low-quality-base > [f] will be discarded. [0.5]
       -c [c]    phred quality < [c] is treated as low quality base. [ 'F' ]
       -p [i]    [i] reads as a part to output. [0]
       -m [i]    minimum length of reads after cleaning. [30]
       -u [s]    UID types. ( NNNNNACT  means 5bp barcode + 3bp fix
                              ACGNNNNN  meads 3bp fix + 5bp barcode
                              NNNNNNNN  meads 8bp barcode. )
       -t [i]    trim first N bp of reads (if has UID, trim N bp after UID) [0]

