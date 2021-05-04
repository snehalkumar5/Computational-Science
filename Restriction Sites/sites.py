import json
import sys
import pandas as pd

sites={}
MIN_SITE_LENGTH=4
MAX_SITE_LENGTH=8
    
def getcomplement(DNA):
    strand = DNA
    strand = strand.replace("G","%temp%").replace("C","G").replace("%temp%","C")
    strand = strand.replace("A","%var%").replace("T","A").replace("%var%","T")
    return strand

def get_restrictionsites(DNA):
    # reversedDNA = DNA[::-1]
    complementaryDNA = getcomplement(DNA)
    length=len(DNA)
    # print(complementaryDNA[::-1])
    for i in range(MIN_SITE_LENGTH,MAX_SITE_LENGTH+1):
        for j in range(length-i+1):
            k = j+i
            seq = DNA[j:k]
            # if DNA[j:k] == reversedDNA[length-k:length-j]:
            # rev = reversedDNA[j:k]
            rev = complementaryDNA[j:k]
            if seq == rev[::-1]:
                if seq in sites:
                    sites[seq]["Frequency"] = int(sites[seq]["Frequency"])+1
                    sites[seq]["Location"].append((j+1,k))
                else:
                    sites[seq]={}
                    sites[seq]["Frequency"] = 1
                    sites[seq]["Location"] = [(j+1,k)]
    return sites

if __name__=="__main__":
    # dna_file = "puc19.txt"
    # with open(dna_file,"r") as f:
    #     DNA = f.read()
    # DNA=str(input())
    DNA="""GAGAGAGAGAGAGA"""

    DNA = DNA.replace("\n","")
    DNA = DNA.replace("\r","")  
    DNA = DNA.replace(" ","")
    length = len(DNA)
    DNA = DNA.upper()
    # print(DNA)
    sites = get_restrictionsites(DNA)
    print("Restriction sites:\n")
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_colwidth', None)

    df = pd.DataFrame.from_dict(sites,orient='index')
    print(df)
    # print(sites)
    # print(json.dumps(sites,indent=1))
    # cnt = df["Frequency"].sum()
    print("\nTotal number of distinct restriction recognition sites:",len(sites))
    # print("Total number of restriction recognition sites:",cnt,)