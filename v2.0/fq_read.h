
// yerui@genomics.cn

// 2012-10-31  v1.0
// 2013-04-10  v1.1
// 2014-05-05  v2.0

#ifndef FQ_READ_H_
#define FQ_READ_H_

#include <iostream>
#include <string>
#include <set>

#include "UID_type.h"

using namespace std;

class fq_read {
public:
    string id;
    string seq;
    string qua;

    fq_read():id(""),seq(""),qua(""){}

    inline int base_n( char ) const;
    inline int low_qua_n( char ) const;
    inline void trim( int );
    inline void trim_head( int );

};

istream & operator >> ( istream &in,  fq_read &a )
{
    string plus;
    in >> a.id >> a.seq >> plus >> a.qua;
    return in;
}

ostream & operator << ( ostream &out, fq_read &a )
{
    out << a.id << '\n' << a.seq << "\n+\n" << a.qua << '\n';
    return out;
}

int fq_read::base_n( char c ) const
{
    int num(0);

    for ( int i(0); i != seq.size(); ++i )
        if ( seq[i] == c )  ++num;

    return num;
}

int fq_read::low_qua_n( char c ) const
{
    int num(0);

    for ( int i(0); i != qua.size(); ++i )
        if ( qua[i] < c )  ++num;

    return num;
}

void fq_read::trim( int ov_len )
{
    if ( seq.size() <= ov_len ) return;

    seq.erase( seq.begin() + ov_len, seq.end() );
    qua.erase( qua.begin() + ov_len, qua.end() );

    return;
}

void fq_read::trim_head( int len )
{
    seq.erase(0, len);
    qua.erase(0, len);
}


inline bool is_pair( const fq_read &a, const fq_read &b )
{
    int len = a.id.size();

    if ( len != b.id.size() )
        return false;
    else
        return a.id.substr(0, len-2) == b.id.substr(0, len-2);

}

string re_complement( const string &str )
{
    string a("");
    for ( int i = str.size()-1; i >= 0; --i ) {
        switch( str[i] ) {
            case 'A': a += "T"; break;
            case 'C': a += "G"; break;
            case 'G': a += "C"; break;
            case 'T': a += "A"; break;
            case 'N': a += "N"; break;
            default: cerr << "base error: [" << str << "]" << endl; break;
        }
    }

    return a;
}

void overlap(const string &seq1, const string &seq2, int p1, int p2, string &as, string &bs)
{
    if ( p1 < p2 ) {
        as = seq1.substr(0, seq1.length() - p2 + p1);
        bs = seq2.substr( p2-p1 );
    }
    else if ( p1 == p2 ) {
        as = seq1;
        bs = seq2;
    }
    else {
        as = seq1.substr( p1-p2 );
        bs = seq2.substr(0, seq1.length() - p1 + p2 );
    }
}

//  mismatch_num > n ? true : false

bool mismatch_num( string &x, string &y, int n )
{
    int len = x.size();
/*
    if ( len != y.size() )
        return false;
*/

    int mis_n(0);
    for ( int i(0); i != len; ++i ) {
        if ( x[i] != y[i] ) ++mis_n;
        if ( mis_n > n ) return true;
    }

/*
    if ( n == 2 ) {
        cout << "tail: " << x << "\ntail: " << y << endl;
    }
*/
    return false;
}


// return length of overlapped region

int has_adapter( fq_read &a, fq_read &b )
{

    int a_len = a.seq.size();

    int len(0); // length of overlapped region
    string rb = re_complement( b.seq );

    int s_len(10); // seed length
    set<int> df;   // record offsets between two reads seq
    set<int> okdf; // offsets that OK

    string as(""), bs("");

    // skip first base because its low quality
    for ( int i(1); i < 10; i += 2 ) {
        string seed = a.seq.substr(i, s_len);

        int start(i);
        for (;;) {
            size_t p = rb.find(seed, start);
            if ( p == string::npos) break;

            int offset = p - i;
            if ( df.find( offset ) != df.end() ) {
                start = p + 1;
                continue;
            }
            else
                df.insert( offset );

            overlap(a.seq, rb, i, p, as, bs);

            len = as.size();

            if ( mismatch_num(as, bs, 3) ) {
                start = p + 1;
                continue;
            }

            string tail_a = a.seq.substr(len, 10);
            string tail_b = b.seq.substr(len, 10);

            if ( mismatch_num(tail_a, tail_b, 2) ) {
                start = p + 1;
                continue;
            }

            okdf.insert( offset );
            if ( len > a_len / 2 ) return len;
            start = p + 1;
        }

        if ( len > a_len / 2 ) return len;
    }

    return ( okdf.size() > 0 ? a_len - *(okdf.begin()) : 0);
}

bool barcode_clean( fq_read &a, fq_read &b, int bar_size )
{
    int a_len = a.seq.size();
    string rb = re_complement( b.seq );

    int s_len(10); // seed length
    set<int> df;   // record offsets between two reads seq
    set<int> okdf; // offset that OK

    // skip first base because its low quality
    for ( int i(1); i < 10; i += 2 ) {
        string seed = rb.substr(i, s_len);

        int start(i);
        for (;;) {
            size_t p = a.seq.find(seed, start);
            if ( p == string::npos) break;

            int offset = p - i;
            if ( offset >= bar_size ) break;  // no contamination

            if ( df.find( offset ) != df.end() ) {
                start = p + 1;
                continue;
            }
            else
                df.insert( offset );

            string as(""), bs("");
            overlap(a.seq, rb, p, i, as, bs);

            int mismatch(0);
            if ( mismatch_num(as, bs, 3) ) {
                start = p + 1;
                continue;
            }

            okdf.insert( offset );
            break;
        }

        if ( okdf.size() )
            break;
    }

    if ( okdf.size() ) {
        int rd_len = a_len - bar_size + *(okdf.begin());

        a.trim(rd_len);
        b.trim(rd_len);
/*
        cout << "trim: " << a.seq << endl;
        cout << "trim: " << "           " << re_complement( b.seq ) << endl;
*/
        return true;
    }
    else
        return false;
}

bool fetch_UID( fq_read &a, fq_read &b, UIDt &v, int ov_len )
{
    int trim_len = v.bar_len + v.fix_len;

    // if fix_seq is not consistent, return false
    string a_fix = a.seq.substr(v.fix_pos, v.fix_len);
    string b_fix = b.seq.substr(v.fix_pos, v.fix_len);

    if ( mismatch_num(v.fix, a_fix, 1) or mismatch_num(v.fix, b_fix, 1) )
        return false;

    string uid_pr("#");

    uid_pr += a.seq.substr(v.bar_pos, v.bar_len);
    uid_pr += "#";
    uid_pr += b.seq.substr(v.bar_pos, v.bar_len);

    int insP = a.id.size() - 2;
    a.id.insert(insP, uid_pr);
    b.id.insert(insP, uid_pr);

    //cout << "test id\n" << a.id << '\n' << b.id << endl;

    if ( ov_len ) { // has adapter
        int rd_len = ov_len - trim_len * 2;

        a.seq = a.seq.substr(trim_len, rd_len);
        a.qua = a.qua.substr(trim_len, rd_len);

        b.seq = b.seq.substr(trim_len, rd_len);
        b.qua = b.qua.substr(trim_len, rd_len);
    }
    else {
        a.seq.erase(0, trim_len);
        a.qua.erase(0, trim_len);

        b.seq.erase(0, trim_len);
        b.qua.erase(0, trim_len);
    }

    return true;
}


#endif // FQ_READ_H_
