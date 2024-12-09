
def loadSiteVariantData(fileName):
    """Load our variant data in its by site form."""
    f=open(fileName,"r")

    header=f.readline().rstrip()
    sampleNamesL=header.split()[5:]

    siteDataL=[]

    while True:
        s=f.readline()
        if s=="":
            break

        s=s.rstrip()
        strL=s.split("\t")
        Chr=strL[0]
        pos=strL[1]
        ID=strL[2]
        ref=strL[3]
        chimp=strL[4]
        genDataL=strL[5:]
        
        genotypesL=[]
        for dipgeno in genDataL:
            a,b=dipgeno.split(",")
            genotypesL.append((a,b))

        siteDataL.append((Chr,int(pos),ID,ref,chimp,tuple(genotypesL)))

    f.close()

    return sampleNamesL,siteDataL


def getDataByHaplotype(sampleNamesL,siteDataL):
    """Take data organized by site (sites in rows) and convert to a
    format organized by haplotype (haplotypes in rows)."""

    hapNamesL=[]
    for sampleName in sampleNamesL:
        hapNamesL.append(sampleName+".A")
        hapNamesL.append(sampleName+".B")

    hapLL=[[] for i in range(len(hapNamesL))] # where we'll put the data

    for (Chr,pos,ID,ref,chimp,genotypesT) in siteDataL:
        i=0
        for a,b in genotypesT:
            hapLL[i].append(a)
            i+=1
            hapLL[i].append(b)
            i+=1

    # tuple-ize and return
    hapDataL=[]
    for hapL in hapLL:
        hapDataL.append(tuple(hapL))

    return hapNamesL,hapDataL

def getChimp(siteDataL):
    """Get chimp alleles at each site, and return as list."""
    chimpL=[]
    for (Chr,pos,ID,ref,chimp,genotypesT) in siteDataL:
        chimpL.append(chimp)
    return chimpL

finSampleNamesL,finSiteDataL=loadSiteVariantData("fin-2.136445439-136692143.tsv")
yorSampleNamesL,yorSiteDataL=loadSiteVariantData("yor-2.136445439-136692143.tsv")

chimpL = getChimp(finSiteDataL) # will be same for fin and yor

# convert to hap based
finHapNamesL,finHapDataL=getDataByHaplotype(finSampleNamesL,finSiteDataL)
yorHapNamesL,yorHapDataL=getDataByHaplotype(yorSampleNamesL,yorSiteDataL)

del finSampleNamesL
del finSiteDataL
del yorSampleNamesL
del yorSiteDataL

