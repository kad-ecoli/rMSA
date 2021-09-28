const char* docstring=""
"getRNAcentralTaxonID accession.list rnacentral.tsv output.tsv\n"
"    read RNAcentral accessions list and taxonID mapping file rnacentral.tsv,\n"
"    write output file to output.tsv\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <map>

using namespace std;

int getRNAcentralTaxonID(const string infile="-",
    const string mapfile="rnacentral.tsv", const string outfile="-")
{
    /* read accession list */
    ifstream fp_in;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    string accession,line;
    size_t r;
    map<string,string> taxon_map;
    vector<string> accession_list;
    vector<string> header_list;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        accession=line;
        if (line.size()>1 && line[0]=='>') accession=line.substr(1);   
        if (accession.size()<4) continue;
        if (accession.substr(0,3)!="URS") continue;
        for (r=3;r<accession.size();r++)
        {
            if (('A'<=accession[r] && accession[r]<='Z')||
                ('0'<=accession[r] && accession[r]<='9')) continue;
            accession=accession.substr(0,r);
            break;
        }
        taxon_map[accession]="";
        accession_list.push_back(accession);
        header_list.push_back(line);
    }
    fp_in.close();

    /* read mapping */
    if (mapfile!="-") fp_in.open(mapfile.c_str(),ios::in);
    string Lch,taxon;
    while ((mapfile!="-")?fp_in.good():cin.good())
    {
        if (mapfile!="-") fp_in>>accession>>Lch>>taxon;
        else cin>>accession>>Lch>>taxon;
        if (taxon_map.count(accession))
            taxon_map[accession]=taxon;
    }
    fp_in.close();

    /* output */
    size_t nseqs=0;
    size_t i;
    ofstream fp_out;
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    for (i=0;i<accession_list.size();i++)
    {
        accession=accession_list[i];
        if (taxon_map[accession].size()==0) continue;
        nseqs++;
        if (outfile=="-") cout<<header_list[i]+'\t'+taxon_map[accession]<<endl;
        else            fp_out<<header_list[i]+'\t'+taxon_map[accession]<<endl;
    }
    fp_out.close();

    /* clean up */
    accession.clear();
    line.clear();
    Lch.clear();
    taxon.clear();
    map<string,string> ().swap(taxon_map);
    vector<string> ().swap(accession_list);
    vector<string> ().swap(header_list);
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    string mapfile=argv[2];
    string outfile=(argc<=3)?"-":argv[3];
    getRNAcentralTaxonID(infile, mapfile, outfile);
    return 0;
}
