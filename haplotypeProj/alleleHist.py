import matplotlib.pyplot as pyplot
from coalSim import *

## Plotting derived allele freq hist of real lac data

def alleleFreqs(hapDataL,chimpL):
    """Calculate the derived allele freq for every site in
    hapDataL. Return as list."""
    freqL=[]
    numSites=len(hapDataL[0])
    numHaps=len(hapDataL)
    for site in range(numSites):
        allelesL=[]
        chimp=chimpL[site]

        for hap in range(numHaps):
            allelesL.append(hapDataL[hap][site])

        allelesPresentL=list(set(allelesL))

        # we're only interested in sites with 2 alleles. skip
        # 1's. There are no 3's in fin yor (and we'd skip them if
        # there were).
        if len(allelesPresentL)==2:
            if allelesPresentL[0] == chimp:
                derivedAllele = allelesPresentL[1]
            else:
                derivedAllele = allelesPresentL[0]
            freq = float(allelesL.count(derivedAllele))/len(allelesL)
            if freq != 0 and freq != 1: # skip the 0's and 1's
                freqL.append(freq)

    return freqL

def lacPlots(hapDataL,chimpL):
    """Wrapper for plotting allele count histograms."""

    freqL=alleleFreqs(hapDataL,chimpL)
    pyplot.hist(freqL)
    pyplot.title("Derived allele frequency histogram")
    pyplot.ylabel("Number of occurances in each bin")
    pyplot.ylim(0,350)
    pyplot.show()

## Plotting derived allele freq hist of coalescent simulations

def coalAlleleFreqs(hapL,numMuts):
    """Allele freq histogram of coalSim output."""
    derAlleleFreqL=[]
    for i in range(numMuts):
        count=0
        for seq in hapL:
            if i in seq:
                count+=1
        freq=float(count)/len(hapL)
        if freq!=0 and freq!=1:
            derAlleleFreqL.append(freq)
    return derAlleleFreqL

def coalHist(derAlleleFreqLL, numBins):
    """Takes a list of derived allele frequency lists. Calculates
    histograms for each, then plots a histogram with the average and
    standard error of the bin heights for these. numBins specifies the
    number of bins which we divide the frequencies into."""

    numReps=len(derAlleleFreqLL)

    # determine bins
    spacer=1.0/numBins
    totBins=[spacer*x for x in range(numBins+1)]

    # put the hists for each entry in a numpy array
    coalHistAr=numpy.zeros((numReps,numBins))
    for i in range(numReps):
        coalHistAr[i,]=numpy.histogram(derAlleleFreqLL[i],totBins)[0]

    meanL=[]
    stderrL=[]
    for i in range(numBins):
        meanL.append(numpy.mean(coalHistAr[:,i]))
        stderrL.append(stats.tsem(coalHistAr[:,i]))

    wid=(totBins[1]-totBins[0])
    pyplot.bar(totBins[:-1],meanL,width=wid,yerr=stderrL,color="gray")

def coalPlots(popSize,numAlleles,numMuts,numReps):
    """Wrapper for coalescent simulations."""

    numBins=10

    ## coalescent analysis
    derAlleleFreqLL=[]

    for i in range(numReps):
        hapL=coalSim(numAlleles,popSize,numMuts) # numbers to match fin
        derAlleleFreqL=coalAlleleFreqs(hapL,numMuts)
        derAlleleFreqLL.append(derAlleleFreqL)

    # plot mean of hists (with errorbars)
    pyplot.figure()
    coalHist(derAlleleFreqLL,numBins)

    #pyplot.ylim(0,700)
    pyplot.title("Coalescent derived allele frequency histogram")
    pyplot.xlabel("Derived allele frequency bins")

    pyplot.show()
