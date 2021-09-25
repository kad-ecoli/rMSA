const char* docstring=""
"filla3m input.a3m output.a3m\n"
"    fill in gaps using insertion residues if the residue type is identical to\n"
"    the most frequent residue type in that position\n"
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
#include <cctype>

using namespace std;

size_t filla3m(const string infile, const string outfile)
{
    /* read aln file */
    vector<string> aln_list;
    vector<string> sequence_list;
    vector<string> header_list;
    string sequence;
    string aln;
    size_t r=0;
    size_t L=0;
    char aa;
    ifstream fp;
    if (infile!="-") fp.open(infile.c_str(),ios::in);
    while ((infile!="-")?fp.good():cin.good())
    {
        if (infile!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.size()==0) continue;
        else if (sequence[0]=='>')
        {
            header_list.push_back(sequence);
            continue;
        }

        aln.clear();
        for (r=0;r<sequence.size();r++)
            if ((sequence[r]<'a' || sequence[r]>'z') && sequence[r]!='.')
                aln+=sequence[r];
        if (aln_list.size()==0) L=aln.size();
        else if (aln_list.size() && L!=aln.size())
        {
            cerr<<"ERROR! length mismatch for sequence "<<aln_list.size()+1
                <<". "<<L<<"!="<<aln.size()<<endl;
            exit(0);
        }

        aln_list.push_back(aln);
        sequence_list.push_back(sequence);
        if (aln_list.size()==INT_MAX)
        {
            cerr<<"WARNING! Cannot read beyond sequence number"<<INT_MAX<<endl;
            break;
        }
    }
    fp.close();

    /* get residue type */
    string aa_list="";
    vector<vector<size_t> >aa_count_mat;
    vector<size_t> aa_count_vec(L,0);
    size_t seq_num=aln_list.size();
    size_t a,i,j;
    bool new_aa;
    for (i=0;i<seq_num;i++)
    {
        for (j=0;j<L;j++)
        {
            aa=aln_list[i][j];
            if (aa=='-') continue;
            new_aa=true;
            for (a=0;a<aa_list.size();a++)
            {
                if (aa==aa_list[a])
                {
                    new_aa=false;
                    break;
                }
            }
            if (new_aa)
            {
                aa_list+=aa;
                aa_count_mat.push_back(aa_count_vec);
            }
            aa_count_mat[a][j]++;
        }
    }

    /* get the most frequent aa type per column */
    for (j=0;j<L;j++)
    {
        for (a=0;a<aa_list.size();a++)
        {
            if (aln_list[0][j]==aa_list[a])
            {
                aa_count_vec[j]=a;
                break;
            }
        }
        for (a=0;a<aa_list.size();a++)
        {
            if (aa_count_mat[a][j]>aa_count_mat[aa_count_vec[j]][j])
                aa_count_vec[j]=a;
        }
    }
    vector<char> max_aa_vec(L,'-');
    for (j=0;j<L;j++) max_aa_vec[j]=aa_list[aa_count_vec[j]];
    
    vector<vector<size_t> >().swap(aa_count_mat);
    vector<string> ().swap(aln_list);
    vector<size_t> ().swap(aa_count_vec);
    aa_list.clear();
    
    /* fill gap using insertion */
    vector<char> match_vec(L,'-');
    vector<string> insert_vec(L+1,"");
    char max_aa; // most frequent aa at a position
    bool succeed=true;
    for (i=0;i<seq_num;i++)
    {
        sequence=sequence_list[i];
        for (r=0;r<L+1;r++) insert_vec[r].clear();
        r=0;
        for (j=0;j<sequence.size();j++)
        {
            aa=sequence[j];
            if ('a'<=aa && aa<='z') insert_vec[r]+=aa;
            else
            {
                match_vec[r]=aa;
                r++;
            }
        }

        succeed=true;
        while (succeed)
        {
            succeed=false;
            for (r=1;r<L;r++)
            {
                if (insert_vec[r].size()==0) continue;
                for (j=r;j>=0;j--)
                {
                    if (match_vec[j]!='-')
                    {
                        j++;
                        break;
                    }
                }
                if (0<=j && j<=r)
                {
                    max_aa=max_aa_vec[j];
                    if (toupper(insert_vec[r][0])==max_aa)
                    {
                        match_vec[j]=max_aa;
                        insert_vec[r]=insert_vec[r].substr(1);
                        succeed=true;
                        continue;
                    }
                }
                for (j=r;j<L;j++)
                {
                    if (match_vec[j]!='-')
                    {
                        j--;
                        break;
                    }
                }
                if (r<=j && j<L)
                {
                    max_aa=max_aa_vec[j];
                    if (toupper(insert_vec[r].back())==max_aa)
                    {
                        match_vec[j]=max_aa;
                        insert_vec[r]=insert_vec[r].substr(0,insert_vec[r].size()-1);
                        succeed=true;
                        continue;
                    }
                }
                if (match_vec[r-1]=='-')
                {
                    max_aa=max_aa_vec[r-1];
                    if (toupper(insert_vec[r].back())==max_aa)
                    {
                        match_vec[r-1]=max_aa;
                        insert_vec[r-1]+=insert_vec[r].substr(0,insert_vec[r].size()-1);
                        insert_vec[r].clear();
                        succeed=true;
                        continue;
                    }
                }
                if (match_vec[r]=='-')
                {
                    max_aa=max_aa_vec[r];
                    if (toupper(insert_vec[r][0])==max_aa)
                    {
                        match_vec[r]=max_aa;
                        insert_vec[r+1]=insert_vec[r].substr(1)+insert_vec[r+1];
                        succeed=true;
                        continue;
                    }
                }
            }
            /* parse initial insertion */
            if (insert_vec[0].size())
            {
                for (j=0;j<L;j++)
                {
                    if (match_vec[j]!='-')
                    {
                        j--;
                        break;
                    }
                }
                if (j>=0)
                {
                    max_aa=max_aa_vec[j];
                    if (toupper(insert_vec[0].back())==max_aa)
                    {
                        match_vec[j]=max_aa;
                        insert_vec[0]=insert_vec[0].substr(0,insert_vec[0].size()-1);
                        succeed=true;
                    }
                }
            }
            /* parse ending insertion */
            if (insert_vec[L].size())
            {
                for (j=L;j>0;j--)
                {
                    if (match_vec[j]!='-')
                    {
                        j++;
                        break;
                    }
                }
                if (j<L)
                {
                    max_aa=max_aa_vec[j];
                    if (toupper(insert_vec[L][0])==max_aa)
                    {
                        match_vec[j]=max_aa;
                        insert_vec[L]=insert_vec[L].substr(1);
                        succeed=true;
                    }
                }
            }
        }

        sequence.clear();
        for (r=0;r<L;r++) sequence+=insert_vec[r]+match_vec[r];
        sequence_list[i]=sequence;
    }


    /* output */
    string txt;
    for (i=0;i<seq_num;i++)
    {
        if (header_list.size()>i) txt+=header_list[i]+'\n';
        txt+=sequence_list[i]+'\n';
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
    vector<char>   ().swap(max_aa_vec);
    vector<char>   ().swap(match_vec);
    vector<string> ().swap(insert_vec);
    vector<string> ().swap(header_list);
    vector<string> ().swap(sequence_list);
    sequence.clear();
    aln.clear();
    txt.clear();
    return seq_num;
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile =argv[1];
    string outfile=(argc<=2)?"-":argv[2];
    filla3m(infile,outfile);
    return 0;
}
