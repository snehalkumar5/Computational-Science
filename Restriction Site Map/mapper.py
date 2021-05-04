import pandas as pd
RE_DICT={
    "BamH1":"GGATCC","Bal2":"AGATCT","Bal1":"TGGCCA","EcoR1":"GAATTC","BstUI":"CGCG"
}

def restrictionmap(DNA,RE):
    re_map={}
    # re_map[RE]={}
    print("RE: ",RE, "\nRecognition sequence: ",RE_DICT[RE])
    re_map["Cut site"] = []
    re_map["Location"] = []
    ck=RE_DICT[RE]
    length = len(ck)
    for i in range(len(DNA)):
        if ck[0] == DNA[i]:
            if DNA[i:i+length] == ck:
                re_map["Location"].append((i,i+length))
                re_map["Cut site"].append(i+2)
                # re_map["Site"].append(DNA[-len(DNA)+i-8:-len(DNA)+i]+"-"+DNA[i:i+length]+"-"+DNA[i+length:i+length+8])
    return re_map

if __name__=="__main__":
    # dna_file = "puc19.txt"
    # with open(dna_file,"r") as f:
    #     DNA = f.read()
    DNA="""tcgcgcgttt cggtgatgac ggtgaaaacc tctgacacat gcagctcccg gagacggtca
cagcttgtct gtaagcggat gccgggagca gacaagcccg tcagggcgcg tcagcgggtg
ttggcgggtg tcggggctgg cttaactatg cggcatcaga gcagattgta ctgagagtgc
accatatgcg gtgtgaaata ccgcacagat gcgtaaggag aaaataccgc atcaggcgcc
attcgccatt caggctgcgc aactgttggg aagggcgatc ggtgcgggcc tcttcgctat
tacgccagct ggcgaaaggg ggatgtgctg caaggcgatt aagttgggta acgccagggt
tttcccagtc acgacgttgt aaaacgacgg ccagtgaatt cgagctcggt acccggggat
cctctagagt cgacctgcag gcatgcaagc ttggcgtaat catggtcata gctgtttcct
gtgtgaaatt gttatccgct cacaattcca cacaacatac gagccggaag cataaagtgt
aaagcctggg gtgcctaatg agtgagctaa ctcacattaa ttgcgttgcg ctcactgccc
gctttccagt cgggaaacct gtcgtgccag ctgcattaat gaatcggcca acgcgcgggg
agaggcggtt tgcgtattgg gcgctcttcc gcttcctcgc tcactgactc gctgcgctcg
gtcgttcggc tgcggcgagc ggtatcagct cactcaaagg cggtaatacg gttatccaca
gaatcagggg ataacgcagg aaagaacatg tgagcaaaag gccagcaaaa ggccaggaac
cgtaaaaagg ccgcgttgct ggcgtttttc cataggctcc gcccccctga cgagcatcac
aaaaatcgac gctcaagtca gaggtggcga aacccgacag gactataaag ataccaggcg
tttccccctg gaagctccct cgtgcgctct cctgttccga ccctgccgct taccggatac
ctgtccgcct ttctcccttc gggaagcgtg gcgctttctc aatgctcacg ctgtaggtat
ctcagttcgg tgtaggtcgt tcgctccaag ctgggctgtg tgcacgaacc ccccgttcag
cccgaccgct gcgccttatc cggtaactat cgtcttgagt ccaacccggt aagacacgac
ttatcgccac tggcagcagc cactggtaac aggattagca gagcgaggta tgtaggcggt
gctacagagt tcttgaagtg gtggcctaac tacggctaca ctagaaggac agtatttggt
atctgcgctc tgctgaagcc agttaccttc ggaaaaagag ttggtagctc ttgatccggc
aaacaaacca ccgctggtag cggtggtttt tttgtttgca agcagcagat tacgcgcaga
aaaaaaggat ctcaagaaga tcctttgatc ttttctacgg ggtctgacgc tcagtggaac
gaaaactcac gttaagggat tttggtcatg agattatcaa aaaggatctt cacctagatc
cttttaaatt aaaaatgaag ttttaaatca atctaaagta tatatgagta aacttggtct
gacagttacc aatgcttaat cagtgaggca cctatctcag cgatctgtct atttcgttca
tccatagttg cctgactccc cgtcgtgtag ataactacga tacgggaggg cttaccatct
ggccccagtg ctgcaatgat accgcgagac ccacgctcac cggctccaga tttatcagca
ataaaccagc cagccggaag ggccgagcgc agaagtggtc ctgcaacttt atccgcctcc
atccagtcta ttaattgttg ccgggaagct agagtaagta gttcgccagt taatagtttg
cgcaacgttg ttgccattgc tacaggcatc gtggtgtcac gctcgtcgtt tggtatggct
tcattcagct ccggttccca acgatcaagg cgagttacat gatcccccat gttgtgcaaa
aaagcggtta gctccttcgg tcctccgatc gttgtcagaa gtaagttggc cgcagtgtta
tcactcatgg ttatggcagc actgcataat tctcttactg tcatgccatc cgtaagatgc
ttttctgtga ctggtgagta ctcaaccaag tcattctgag aatagtgtat gcggcgaccg
agttgctctt gcccggcgtc aatacgggat aataccgcgc cacatagcag aactttaaaa
gtgctcatca ttggaaaacg ttcttcgggg cgaaaactct caaggatctt accgctgttg
agatccagtt cgatgtaacc cactcgtgca cccaactgat cttcagcatc ttttactttc
accagcgttt ctgggtgagc aaaaacagga aggcaaaatg ccgcaaaaaa gggaataagg
gcgacacgga aatgttgaat actcatactc ttcctttttc aatattattg aagcatttat
cagggttatt gtctcatgag cggatacata tttgaatgta tttagaaaaa taaacaaata
ggggttccgc gcacatttcc ccgaaaagtg ccacctgacg tctaagaaac cattattatc
atgacattaa cctataaaaa taggcgtatc acgaggccct ttcgtc"""
    DNA = DNA.replace("\n","").replace("\r","").replace(" ","")
    DNA = DNA.upper()
    print("DNA Sequence:",DNA)
    print("Number of bases:",len(DNA))
    RE="BstUI"
    sites = restrictionmap(DNA,RE)
    sites["Fragment"]=[]
    sites["Fragment size"]=[]
    pos=0
    ls=[]
    ps=[]
    check = sites["Cut site"]
    for i,val in enumerate(check):
        if i == len(check)-1:
            ps.append((val+1,check[0]))
            ll = len(DNA)-check[i]+check[0]
            ls.append(ll)
        else:
            ps.append((val+1,check[i+1]))
            ls.append(check[i+1]-val)
        
        pos = val
    sites["Fragment size"] = ls
    sites["Fragment"] = ps
    # df = pd.DataFrame.from_dict(sites,orient='index')
    df = pd.DataFrame(data=sites)
    df.sort_values(by=["Cut site"],ascending=True,inplace=True)
    df.sort_values(by=["Fragment size"],ascending=False,inplace=True)
    print("\nRestriction sites Map:")
    print(df.to_string())
    print("\nTotal number of restriction sites: ",len(df))