// yerui@genomics.cn
// 2013-04-09

#ifndef UID_TYPE_H_
#define UID_TYPE_H_

#include <string>

using namespace std;

class UIDt
{
public:
    int bar_pos;
    int bar_len;
    int fix_pos;
    int fix_len;
    string bar;
    string fix;

    UIDt():bar_pos(0),bar_len(0),fix_pos(0),fix_len(0),bar(""),fix("") {};

    void assign( string &s );

};

void UIDt::assign( string &s )
{
    if ( s[0] == 'N' ) {  // barcode first

        size_t p = s.rfind("N");

        bar_pos = 0;
        bar_len = p + 1;

        fix_pos = bar_len;
        fix_len = s.size() - bar_len;
    }
    else {
        size_t p = s.find("N");

        bar_pos = p;
        bar_len = s.size() - p;

        fix_pos = 0;
        fix_len = p;
    }

    bar = s.substr(bar_pos, bar_len);
    fix = s.substr(fix_pos, fix_len);
}

ostream & operator << ( ostream &out, UIDt &a )
{
    out << a.bar << '\t' << a.bar_pos << '\t' << a.bar_len << '\n'
        << a.fix << '\t' << a.fix_pos << '\t' << a.fix_len << '\n';

    return out;
}

#endif // UID_TYPE_H_
