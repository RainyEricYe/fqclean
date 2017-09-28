/*
 *  yerui@genomics.cn
 *  2012-12-31
 *  2013-04-09
 *
 */

#include "main.h"

int main( int argc, char **argv )
{
    float n_cut(0.1);
    float q_cut(0.5);
    char  qc_cut('F');
    ulong part_size(0);
    int   min_len(30);
    string bar_str(""); // UID type
    UIDt uid_t;
    int trim_h(0);

    mIntLong LenMap; // read length -> freqency


    int c;
    while ( (c=getopt(argc,argv,"n:q:c:p:m:u:t:h")) != -1 ) {
        switch (c) {
            case 'n': n_cut  = atof(optarg);    break;
            case 'q': q_cut  = atof(optarg);    break;
            case 'c': qc_cut = optarg[0];       break;
            case 'p': part_size = atoi(optarg); break;
            case 'm': min_len = atoi(optarg);   break;
            case 'u': bar_str = optarg;         break;
            case 't': trim_h  = atoi(optarg);   break;
            case 'h':
            default:  usage(argv[0]);           exit(1);
        }
    }

    cout << "# n: " << n_cut
         << "\tq: " << q_cut
         << "\tc: " << qc_cut
         << "\tp: " << part_size
         << "\tm: " << min_len
         << "\tu: " << bar_str
         << "\tt: " << trim_h
         << endl;

    int bar_size = bar_str.size();

    if ( bar_size ) {
        uid_t.assign( bar_str );
        min_len += bar_size * 2 + trim_h; // conside with adapter contamination
    }

    if ( argc - optind != 4 )
        usage(argv[0]), exit(1);

    char *f[4];
    for ( int i(0); i != 4; ++i ) {
        f[i] = argv[optind + i];
    }

    if ( access(f[0], R_OK) != 0 or access(f[1], R_OK) != 0 )
        cerr << "infile access error!" << endl, exit(1);

    igzstream in1( f[0] );
    igzstream in2( f[1] );

    string suffix(".fq.gz");
    string ofq1( f[2] ), ofq2( f[3] );

    if (   ofq1.rfind( suffix ) != ofq1.size() - 6
        or ofq2.rfind( suffix ) != ofq2.size() - 6 )
        cerr << "out fq files should have suffix as .fq.gz " << endl, exit(0);

    int part(1);
    string outpre1 = ofq1.substr(0, ofq1.size()-6 );
    string outpre2 = ofq2.substr(0, ofq2.size()-6 );

    if ( part_size > 0 ) {
        ofq1 = outpre1 + '.' + mtoa(part) + suffix;
        ofq2 = outpre2 + '.' + mtoa(part) + suffix;
    }

    ogzstream out1( ofq1.c_str() );
    ogzstream out2( ofq2.c_str() );

    if ( access(ofq1.c_str(), W_OK) != 0 or access(ofq2.c_str(), W_OK) != 0 )
        cerr << "output error!" << endl, exit(1);

    uint clean_rp(0);
    uint total_rp(0);

    uint adapter_c(0); // num of reads pairs trimed because of adapter
    uint short_c(0);   // num of reads pairs discarded because of adapter
    uint barcode_c(0); // num of reads pairs trimed because of barcode contamination
    uint lowqua_c(0);  // num of reads pairs discarded because of low quality
    uint fixseq_c(0);  // num of reads pairs discarded because of mismatch on fix_seq


    // fq1_raw  fq2_raw    fq1_clean  fq2_clean
    // 0        1          2          3

    ulong base_n[4] = {0,0,0,0};
    ulong N_n[4] = {0,0,0,0};
    ulong lq_n[4] = {0,0,0,0};  // low quality

    for (;;) {
        if ( !in1.good() or !in2.good() ) break;

        fq_read  a, b;
        in1 >> a;
        in2 >> b;

        ++total_rp;

        if ( a.id.size() < 1 ) {
            continue;
        }

        if ( !is_pair(a, b) ) {
            cerr << "ERROR: not pair reads in line: " << total_rp * 4 - 3 << endl;
            cerr << '[' << a.id << "]\t[" << b.id << ']' << endl;
            break;
        }

        uint a_len, b_len, a_N, b_N, a_low, b_low;

        a_len = a.seq.size();
        b_len = b.seq.size();
        a_N = a.base_n('N');
        b_N = b.base_n('N');
        a_low = a.low_qua_n(qc_cut);
        b_low = b.low_qua_n(qc_cut);

        base_n[0] += a_len;
        base_n[1] += b_len;
        N_n[0]    += a_N;
        N_n[1]    += b_N;
        lq_n[0]   += a_low;
        lq_n[1]   += b_low;

        // adapter contamination
        int ov_len = has_adapter(a,b);

        if ( ov_len >= min_len ) {
            a.trim(ov_len);
            b.trim(ov_len);

            ++adapter_c;
        }
        else if ( ov_len > 0 ) { // ov_len < min_len
            continue;
            ++short_c;
        }
        else { // ov_len == 0, no adapter contamination
            if ( bar_size ) {
                if ( barcode_clean(a, b, bar_size) ) // clean barcode contamination at the end of reads
                    ++barcode_c;
            }
        }

        if ( bar_size ) {
            if ( ! fetch_UID(a, b, uid_t, ov_len) ) {
                ++fixseq_c;
                continue;
            }
        }

        if ( trim_h > 0 ) {
            a.trim_head( trim_h );
            b.trim_head( trim_h );
        }

        a_len = a.seq.size();
        b_len = b.seq.size();
        a_N = a.base_n('N');
        b_N = b.base_n('N');
        a_low = a.low_qua_n(qc_cut);
        b_low = b.low_qua_n(qc_cut);


        if ( (float)a_N > a_len * n_cut or (float)b_N > b_len * n_cut ) {
            ++lowqua_c;
            continue;
        }

        if ( (float)a_low > a_len * q_cut or (float)b_low > b_len * q_cut ) {
            ++lowqua_c;
            continue;
        }

        ++LenMap[a_len];
        ++clean_rp;

        out1 << a;
        out2 << b;

        if ( part_size && clean_rp % part_size == 0 ) {
            ++part;
            ofq1 = outpre1 + '.' + mtoa(part) + suffix;
            ofq2 = outpre2 + '.' + mtoa(part) + suffix;

            out1.close();
            out1.clear();
            out1.open(ofq1.c_str());

            out2.close();
            out2.clear();
            out2.open(ofq2.c_str());
        }

        base_n[2] += a_len;
        base_n[3] += b_len;
        N_n[2]    += a_N;
        N_n[3]    += b_N;
        lq_n[2]   += a_low;
        lq_n[3]   += b_low;
    }

    in1.close();
    in2.close();
    out1.close();
    out2.close();

    cout << "\ntype\traw_data\tclean_data\n"
         << "Read_pairs\t" << total_rp << '\t' << clean_rp << '\n'
         << "Bases\t" << base_n[0] + base_n[1] << '\t' << base_n[2] + base_n[3] << '\n'
         << "N_in_fq1\t"  << N_n[0] << '\t' << N_n[2] << '\n'
         << "N_in_fq2\t"  << N_n[1] << '\t' << N_n[3] << '\n'
         << "Low_quality_base_in_fq1\t" << lq_n[0] << '\t' << lq_n[2] << '\n'
         << "Low_quality_base_in_fq2\t" << lq_n[1] << '\t' << lq_n[3] << "\n\n";

    printf("Rate(clean_base/raw)\t%.2f%%\n", (base_n[2]+base_n[3]) * 100.0 / (base_n[0]+base_n[1]) );
    printf("Rate(clean_read_pair/raw)\t%.2f%%\n\n", clean_rp * 100.0 / total_rp );

    cout << "Reason\tOperate\tReads_pairs\tRate_over_raw(%)\n";
    printf("adapter\ttrim\t%d\t%.2f\n", adapter_c, adapter_c * 100.0 / total_rp );
    printf("adapter\tdiscard\t%d\t%.2f\n", short_c, short_c * 100.0 / total_rp );
    printf("barcode\ttrim\t%d\t%.2f\n", barcode_c, barcode_c * 100.0 / total_rp );
    printf("low-quality\tdiscard\t%d\t%.2f\n", lowqua_c, lowqua_c * 100.0 / total_rp );
    printf("fix-seq-error\tdiscard\t%d\t%.2f\n\n", fixseq_c, fixseq_c * 100.0 / total_rp );

    cout << "# read_length\tfrequency\n";
    for ( mIntLong::iterator it = LenMap.begin(); it != LenMap.end(); ++it ) {
        cout << it->first << '\t' << it->second << endl;
    }

    exit(0);
}
