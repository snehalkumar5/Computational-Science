def transcribe(DNA):
    rna = DNA.replace("T","U")
    return rna

def translate(DNA):
    flag=-1
    codon_table = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W'}
    protein = ""
    for i in range(0,len(DNA)):
        k=i+3
        codon = DNA[i:k]
        if (len(codon)!=3):
            break
        if codon_table[codon]=='M':
            flag=i
            break
    if flag==-1:
        return protein
    for i in range(flag,len(DNA),3):
        k=i+3
        codon = DNA[i:k]
        if (len(codon)!=3):
            break
        if codon_table[codon]=='STOP':
            break
        protein += codon_table[codon]
    return protein

if __name__=="__main__":
    # dna_file = "puc19.txt"
    # with open(dna_file,"r") as f:
    #     DNA = f.read()
    # DNA=str(input())
    DNA="""ATGATGGGGGCCCGACGTACGACGTAA"""
    DNA = DNA.replace("\n","").replace(" ","").replace("\r","").replace("\\","")
    DNA = DNA.upper()
    print("DNA strand:", DNA)

    RNA = transcribe(DNA)
    Protein = translate(DNA)
    print("\nRNA strand:", RNA)
    if Protein=="":
        print("\nError: No start codon found!")
        quit()
    print("\nProtein strand:",Protein)

