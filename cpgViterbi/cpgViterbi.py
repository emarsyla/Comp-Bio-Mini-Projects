import math
import fasta
from context import *


def islandPedictionSummaryStats(predict,correct):
    """Print some statistics summarizing how good the prediction is
    relative to the correct solution. Our definitions of sensitivity
    and specificity assume the 'lower case' model represents the
    feature of interest."""
    trueLower=0
    trueUpper=0
    predictLowerAndTrueLower=0
    predictUpperAndTrueUpper=0
    for i in range(len(predict)):
        if correct[i].islower(): trueLower+=1
        else: trueUpper+=1

        if predict[i].islower() and correct[i].islower(): predictLowerAndTrueLower+=1
        if predict[i].isupper() and correct[i].isupper(): predictUpperAndTrueUpper+=1

    sensitivity = 1.0 * predictLowerAndTrueLower / trueLower
    specificity = 1.0 * predictUpperAndTrueUpper / trueUpper
    
    print("Sensitivity (TPR):", format(sensitivity,".3f"))
    print("Specificity (TNR):", format(specificity,".3f"))
    print("Youden's J       :", format(sensitivity + specificity - 1,".3f"))
    

def runExample():
    """Wrapper to run the provided cpg example."""
    
    states = ["A","C","G","T","a","c","g","t"]
    
    # calculate transition matrix from using A solution for training
    Asolution=fasta.load("A-cpgSolution.fa")[0][1]
    cpgTransD,cpgCountD=context([Asolution],1,states,1)
    cpgTransD = logDictValues(cpgTransD)
    
    # now use viterbi to get sol for B
    B=fasta.load("B-region.fa")[0][1]
    scM,btM=viterbi(B,cpgTransD)
    solB=bt(B,scM,btM)

    f=open("Bpredict.fa","w")
    print("> prediction",file=f)
    print(solB,file=f)
    f.close()

    # test out our results on B vs. the provided solution
    Bsolution=fasta.load("B-cpgSolution.fa")[0][1]
    islandPedictionSummaryStats(solB,Bsolution)

def viterbi(seq, transD):
    '''Runs viterbi algorithm to find most likely path through CPG island and non CPG island states for the
    given data'''
    # Initialize dp table and backtracking table
    cpgIsland = [_ for _ in range(len(seq))]
    nonCpgIsland = [_ for _ in range(len(seq))]
    btCpgIsland = [_ for _ in range(len(seq))]
    btNonCpgIsland = [_ for _ in range(len(seq))]
    cpgIsland[0] = math.log(0.5)
    nonCpgIsland[0] = math.log(0.5)
    btCpgIsland[0] = '-'
    btNonCpgIsland[0] = '-'

    # Run the viterbi algorithm
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        first = pair[0]
        second = pair[1]
        #calculate log likelihoods for each scenario
        nontonon = nonCpgIsland[i] + (transD[first.upper() + second.upper()])
        nontois = nonCpgIsland[i] + (transD[first.upper() + second.lower()])
        istonon = cpgIsland[i] + (transD[first.lower() + second.upper()])
        istois = cpgIsland[i] + (transD[first.lower() + second.lower()])
        # Add to cpg island dp table and backtracking table
        if istois > nontois:
            cpgIsland[i + 1] = istois
            btCpgIsland[i + 1] = "I"
        else: 
            cpgIsland[i + 1] = nontois
            btCpgIsland[i + 1] = "N"
        # Add to non cpg island dp and backtracking tables
        if istonon > nontonon:
            nonCpgIsland[i + 1] = istonon
            btNonCpgIsland[i + 1] = "I"
        else:
            nonCpgIsland[i + 1] = nontonon
            btNonCpgIsland[i + 1] = "N"

    scoreTable = [nonCpgIsland, cpgIsland]
    backtrackTable = [btNonCpgIsland, btCpgIsland]
    return scoreTable, backtrackTable

def bt(seq, scM, btM):
    '''Calculate the most probable path through seq, representing cpg islands as
    lowercase, and non cpg islands as uppercase'''
    cpgSeq = ""
    length = len(seq)
    # Find the final probability for cpg or non cpg islands
    finalProbNon = scM[0][length - 1]
    finalProbIs = scM[1][length - 1]
    # Set current state based on the maximum
    if finalProbIs > finalProbNon:
        currentState = "island"
        bt = btM[1]
    else:
        currentState = "nonIsland"
        bt = btM[0]
    # Trace backward through sequence
    for i in range(length):
        location = length - i - 1
        prevCPG = bt[location]
        char = seq[location]
        # add to cpg sequence based on current character and state
        if currentState == "island":
            cpgSeq = char.lower() + cpgSeq
        else:
            cpgSeq = char.upper() + cpgSeq
        # Change state based on backtrace table
        if prevCPG == "I":
            currentState = "island"
            bt = btM[1]
        else:
            currentState = "nonIsland"
            bt = btM[0]
    return cpgSeq


