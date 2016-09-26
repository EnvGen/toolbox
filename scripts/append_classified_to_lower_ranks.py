#!/usr/bin/env python

import pandas as pd, sys


def append_classified(df):
    ## Append classified name to each Unclassified rank in turn
    p_unc = df[df.phylum=="Unclassified"]
    df.loc[p_unc.index,"phylum"] += "."+df.loc[p_unc.index,"superkingdom"]
    c_unc = df[df["class"]=="Unclassified"]
    df.loc[c_unc.index,"class"] += "."+df.loc[c_unc.index,"phylum"]
    o_unc = df[df["order"]=="Unclassified"]
    df.loc[o_unc.index,"order"] += "."+df.loc[o_unc.index,"class"]
    f_unc = df[df["family"]=="Unclassified"]
    df.loc[f_unc.index,"family"] += "."+df.loc[f_unc.index,"order"]
    g_unc = df[df["genus"]=="Unclassified"]
    df.loc[g_unc.index,"genus"] += "."+df.loc[g_unc.index,"family"]
    s_unc = df[df["species"]=="Unclassified"]
    df.loc[s_unc.index,"species"] += "."+df.loc[s_unc.index,"genus"]

    ## Search and replace multiple "Unclassified."
    df.replace(to_replace="Unclassified."*6,value="Unclassified.", inplace=True, regex=True)
    df.replace(to_replace="Unclassified."*5,value="Unclassified.", inplace=True, regex=True)
    df.replace(to_replace="Unclassified."*4,value="Unclassified.", inplace=True, regex=True)
    df.replace(to_replace="Unclassified."*3,value="Unclassified.", inplace=True, regex=True)
    df.replace(to_replace="Unclassified."*2,value="Unclassified.", inplace=True, regex=True)
    df.replace(to_replace="Unclassified.Unclassified",value="Unclassified", inplace=True, regex=True)
    
    return df

def main():
    infile = sys.argv[1]
    df = pd.read_csv(infile, header=0, sep=",", index_col=0)
    df.fillna("Unclassified",inplace=True)
    df = append_classified(df)
    df.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__': 
    main()
