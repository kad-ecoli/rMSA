const char* docstring=""
"addQuery2a3m input.fasta input.a3m output.a3m\n"
"    add the query sequence (input.fasta) to the template alignment\n"
"    (input.a3m) based on sliding over the first template hit.\n"
"    print the portion of identical residues\n"
"\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <map>
#include <climits>
#include <iomanip>

using namespace std;

void writea3m(const vector<string>&a3m_header_list,
              const vector<string>&a3m_sequence_list, const string &outfile)
{
    ofstream fp_out;
    fp_out.open(outfile.c_str(),ofstream::out);
    size_t i;
    for (i=0;i<a3m_sequence_list.size();i++)
        fp_out<<a3m_header_list[i]<<'\n'<<a3m_sequence_list[i]<<endl;
    fp_out.close();
    return;
}

float addQuery2a3m(const string &infile_fas, const string &infile_a3m, const string &outfile)
{
    float cscore=0;
    /* read query fasta */
    vector<string> a3m_sequence_list(1,"");
    vector<string> a3m_header_list(1,">query");
    string sequence;
    ifstream fp;
    if (infile_fas!="-") fp.open(infile_fas.c_str(),ios::in);
    while ((infile_fas!="-")?fp.good():cin.good())
    {
        if (infile_fas!="-") getline(fp,sequence);
        else getline(cin,sequence);
        if (sequence.size()==0) continue;
        else if (sequence[0]=='>') a3m_header_list[0]=sequence;
        else a3m_sequence_list[0]+=sequence;
    }
    fp.close();

    /* read a3m file */
    string aln;
    size_t r,i,j;
    size_t L=0;
    if (infile_a3m!="-") fp.open(infile_a3m.c_str(),ios::in);
    while ((infile_a3m!="-")?fp.good():cin.good())
    {
        if (infile_a3m!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.length()==0) continue;
        else if (sequence[0]=='>')
        {
            a3m_header_list.push_back(sequence);
            continue;
        }

        aln="";
        for (r=0;r<sequence.size();r++)
            if ((sequence[r]<'a' || sequence[r]>'z') && sequence[r]!='.')
                aln+=sequence[r];
        if (L==0) L=aln.size();
        else if (L && L!=aln.size())
        {
            cerr<<"ERROR! length mismatch for sequence "<<a3m_sequence_list.size()+1
                <<". "<<L<<"!="<<aln.size()<<endl;
            exit(0);
        }

        a3m_sequence_list.push_back(sequence);
        if (a3m_sequence_list.size()==INT_MAX)
        {
            cerr<<"WARNING! Cannot read beyond sequence number"<<INT_MAX<<endl;
            break;
        }
    }
    fp.close();

    ofstream fp_out;
    if (a3m_sequence_list.size()<=1 || a3m_sequence_list[0].size()==L)
    {
        /* if no need to add gap */
        fp_out.open(outfile.c_str(),ofstream::out);
        for (i=0;i<a3m_sequence_list.size();i++)
            fp_out<<a3m_header_list[i]<<'\n'<<a3m_sequence_list[i]<<endl;
        fp_out.close();
        cscore=1.;
    }
    else
    {
        /* find out the best offset by sliding */
        int offset=0;
        int best_offset=0;
        int best_cscore=-1;
        aln  ="";
        sequence=a3m_sequence_list[1];
        for (r=0;r<sequence.size();r++)
            if ((sequence[r]<'a' || sequence[r]>'z') && sequence[r]!='.')
                aln+=sequence[r];
        sequence=a3m_sequence_list[0];
        for (offset=0;offset<=sequence.size()-aln.size();offset++)
        {
            cscore=0;
            for (r=0;r<L;r++) cscore+=(aln[r]==sequence[r+offset]);
            if (cscore>best_cscore)
            {
                best_cscore=cscore;
                best_offset=offset;
            }
        }
        cscore=1.*best_cscore/aln.size();
        offset=best_offset;
        string lgap="";
        string rgap="";
        for (r=0;r<offset;r++) lgap+="-";
        for (r=0;r<sequence.size()-aln.size()-offset;r++) rgap+="-";

        /* write output */
        fp_out.open(outfile.c_str(),ofstream::out);
        fp_out<<a3m_header_list[0]<<'\n'<<a3m_sequence_list[0]<<'\n';
        for (i=1;i<a3m_sequence_list.size();i++)
            fp_out<<a3m_header_list[i]<<'\n'<<lgap<<a3m_sequence_list[i]<<rgap<<endl;
        fp_out.close();

        lgap.clear();
        rgap.clear();
    }

    /* clean up */
    sequence.clear();
    aln.clear();
    vector<string>().swap(a3m_sequence_list);
    vector<string>().swap(a3m_header_list);
    return cscore;
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc!=4)
    {
        cerr<<docstring;
        return 0;
    }
    string infile_fas=argv[1];
    string infile_a3m=argv[2];
    string outfile   =argv[3];
    if (outfile=="-")
    {
        cerr<<"ERROR! '-' is not an acceptable output file name"<<endl;
        return 1;
    }
    cout<<setiosflags(ios::fixed)<<setprecision(4)
        <<addQuery2a3m(infile_fas, infile_a3m, outfile)<<endl;
    /* clean up */
    infile_fas.clear();
    infile_a3m.clear();
    outfile.clear();
    return 0;
}
