const char* docstring=""
"insAln input.a3m input.afa output.a3m\n"
"    add back insertion residues from input.a3m to input.afa. output the\n"
"    result to output.a3m. input.afa but not input.a3m has the query.\n"
"    all input file must be single line fasta\n"
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

size_t insAln(const string infile_a3m, const string infile_afa, const string outfile)
{
    /* read a3m file */
    vector<string> a3m_sequence_list;
    vector<string> a3m_header_list;
    string sequence;
    string aln;
    size_t r=0;
    size_t L=0;
    ifstream fp;
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
        if (a3m_sequence_list.size()==0) L=aln.size();
        else if (a3m_sequence_list.size() && L!=aln.size())
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

    /* read afa file */
    vector<string> afa_aln_list;
    vector<string> afa_header_list;
    if (infile_afa!="-") fp.open(infile_afa.c_str(),ios::in);
    while ((infile_afa!="-")?fp.good():cin.good())
    {
        if (infile_afa!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.length()==0) continue;
        else if (sequence[0]=='>')
        {
            afa_header_list.push_back(sequence);
            continue;
        }

        if (afa_aln_list.size()==0) L=sequence.size();
        else if (afa_aln_list.size() && L!=sequence.size())
        {
            cerr<<"ERROR! length mismatch for sequence "<<afa_aln_list.size()+1
                <<". "<<L<<"!="<<sequence.size()<<endl;
            exit(0);
        }

        afa_aln_list.push_back(sequence);
        if (afa_aln_list.size()==INT_MAX)
        {
            cerr<<"WARNING! Cannot read beyond sequence number"<<INT_MAX<<endl;
            break;
        }
    }

    /* check sequence number */
    if (afa_aln_list.size()!=a3m_sequence_list.size()+1)
    {
        cerr<<"ERROR! a3m input should have one less sequence compared to afa input\n"
            <<infile_a3m<<":"<<a3m_sequence_list.size()<<"\t"
            <<infile_afa<<":"<<afa_aln_list.size()<<endl;
        if (afa_aln_list.size()>1+a3m_sequence_list.size()) exit(1);
    }

    /* find out afa all gap columns */
    vector<bool> template_gap_vec(L,true);
    int i,j;
    for (j=0;j<L;j++)
    {
        for (i=1;i<afa_aln_list.size();i++)
        {
            if (afa_aln_list[i][j]=='-') continue;
            template_gap_vec[j]=false;
            break;
        }
    }
    vector<bool> query_gap_vec(L,true);
    char aa;
    sequence="";
    for (j=0;j<L;j++)
    {
        aa=afa_aln_list[0][j];
        if (aa=='-') continue;
        query_gap_vec[j]=false;
        sequence+=aa;
    }
    afa_aln_list[0]=sequence;

    /* add insertion states */
    L=aln.size();
    vector<string> insert_vec(L+1,"");
    L=afa_aln_list[1].size();
    for (i=1;i<afa_aln_list.size();i++)
    {
        sequence=a3m_sequence_list[i-1];
        for (r=0;r<insert_vec.size();r++) insert_vec[r].clear();
        r=0;
        for (j=0;j<sequence.size();j++)
        {
            aa=sequence[j];
            if (aa=='.') continue;
            else if ('a'<=aa && aa<='z') insert_vec[r]+=aa;
            else r++;
        }
        sequence=insert_vec[0];
        r=0;
        for (j=0;j<L;j++)
        {
            aa=afa_aln_list[i][j];
            if (query_gap_vec[j])
            {
                if (aa=='-') aa=0;
                else aa=tolower(aa);
            }
            if (aa) sequence+=aa;
            if (template_gap_vec[j]==false)
            {
                r++;
                sequence+=insert_vec[r];
            }
        }
        afa_aln_list[i]=sequence;
    }
        
    /* output */
    string txt;
    for (i=0;i<afa_aln_list.size();i++)
    {
        if (afa_header_list.size()>i) txt+=afa_header_list[i]+'\n';
        txt+=afa_aln_list[i]+'\n';
    }
    if (outfile=="-") cout<<txt<<flush;
    else
    {
        ofstream fp_out;
        fp_out.open(outfile.c_str(),ofstream::out);
        fp_out<<txt;
        fp_out.close();
    }

    /* clean up */
    sequence.clear();
    aln.clear();
    txt.clear();
    vector<string> ().swap(a3m_sequence_list);
    vector<string> ().swap(a3m_header_list);
    vector<string> ().swap(afa_aln_list);
    vector<string> ().swap(afa_header_list);
    vector<bool>   ().swap(template_gap_vec);
    vector<bool>   ().swap(query_gap_vec);
    return afa_aln_list.size();
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile_a3m=argv[1];
    string infile_afa=argv[2];
    string outfile=(argc<=3)?"-":argv[3];
    insAln(infile_a3m, infile_afa, outfile);
    return 0;
}
