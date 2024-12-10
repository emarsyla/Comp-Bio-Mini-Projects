from energyD import *
import sys

# Based on Zuker and Stiegler 1981: "Optimal computer folding of large
# RNA sequences using thermodynamics and auxiliary information".
# Charge 0 for bifurcation loops and all unpaired bases which aren't
# part of loops.

MINPAIRDIST=5

#########################
### YOUR FUNCTIONS HERE #
#########################

def efold(RNA,efD,fpD):
    # In dictionary
    if RNA in efD:
        return efD[RNA]
    # Base case
    if len(RNA) < MINPAIRDIST + 1:
        return (0, [])

    #Lose it
    loseE = efold(RNA[1:], efD, fpD)
    bestSoFar = (loseE[0], adjust(loseE[1], 1))

    #Use it
    char = RNA[0]
    if char == 'A':
        charpair = 'U'
    elif char == 'U':
        charpair = 'A'
    elif char == 'C':
        charpair = 'G'
    else:
        charpair = 'C'
    for i in range(MINPAIRDIST, len(RNA)):
        if RNA[i] == charpair:
            insideLoop = fpfold(RNA[:i+1], fpD)
            rest = efold(RNA[i+1:], efD, fpD)
            energy = insideLoop[0] + rest[0]
            if energy <= bestSoFar[0]:
                pairings = [(0, i)]
                pairings.extend(adjust(insideLoop[1], 1))
                pairings.extend(adjust(rest[1], i+1))
                bestSoFar = (energy, pairings)
    efD[RNA] = bestSoFar
    return bestSoFar

def fpfold(RNA,fpD):
    # In dictionary
    if RNA in fpD:
        return fpD[RNA]
    # Find min of hairpin, sbiloop, biloop
    hairpinEn = hairpin(RNA)
    hairpinE = hairpinEn[0]
    minEsbiLoop = sbiLoop(RNA, fpD)
    minEsbiLoopE = minEsbiLoop[0]
    minEbiLoop = biLoop(RNA, fpD)
    minEbiLoopE = minEbiLoop[0]
    energies = [hairpinE, minEsbiLoopE, minEbiLoopE]
    minE = min(energies)
    if minE == hairpinE:
        energy = hairpinEn
    elif minE == minEsbiLoopE:
        energy = minEsbiLoop
    else:
        energy = minEbiLoop
    fpD[RNA] = energy
    return energy

def hairpin(RNA):
    char = RNA[0]
    if char == 'A' or char == 'U':
        pairing = 'A'
    elif char == 'C' or char == 'G':
        pairing = 'G'

    newRNA = RNA[1:-1]
    if len(newRNA) < MINPAIRDIST - 1:
        return float('inf')
    length = str(len(newRNA))
    if len(newRNA) >= 30:
        length = '30+'
    hairpinE = energyD['hairpin' + pairing + length]
    return (hairpinE, [])


def sbiLoop(RNA,fpD):
    if len(RNA) < MINPAIRDIST:
        return (float('inf'), [])
    char = RNA[0]
    firstchar = RNA[0]
    lastchar = RNA[-1]
    firstpair = firstchar + lastchar
    bestSoFar = (float('inf'), [])
    newRNA = RNA[1:-1]
    for i in range(len(newRNA) - MINPAIRDIST):
        for j in range(i + MINPAIRDIST, len(newRNA)):
            # Stacking case
            if i == 0 and j == len(newRNA) - 1:
                stackE = energyD[firstpair + newRNA[i] + newRNA[j]]
                rest = fpfold(newRNA[i:j+1], fpD)
                energy = stackE + rest[0]
                if energy <= bestSoFar[0]:
                    pairings = [(i, j)]
                    pairings.extend(adjust(rest[1], 1))
                    bestSoFar = (energy, pairings)
            # Bulge cases
            elif i == 0:
                length = len(newRNA) - j - 1
                if length >= 30:
                    length = '30+'
                bulgeE = energyD['bulge' + str(length)]
                stackE = energyD[firstpair + newRNA[i] + newRNA[j]]
                rest = fpfold(newRNA[i:j+1], fpD)
                energy = bulgeE + stackE + rest[0]
                if energy <= bestSoFar[0]:
                    pairings = [(i, j)]
                    pairings.extend(adjust(rest[1], 1))
                    bestSoFar = (energy, pairings)
            elif j == len(newRNA) - 1:
                length = i
                if length >= 30:
                    length = '30+'
                bulgeE = energyD['bulge' + str(length)]
                stackE = energyD[firstpair + newRNA[i] + newRNA[j]]
                rest = fpfold(newRNA[i:j+1], fpD)
                energy = bulgeE + stackE + rest[0]
                if energy <= bestSoFar[0]:
                    pairings = [(i, j)]
                    pairings.extend(adjust(rest[1], i+1))
                    bestSoFar = (energy, pairings)
            # Interior face
            else:
                length = len(newRNA) - j + i - 1
                if length >= 30:
                    length = '30+'
                if (char == "G" or char == "C") and (newRNA[i] == 'G' or newRNA[i] == 'C'):
                    interior = 'interiorGG' + str(length)
                elif (char == "A" or char == "U") and (newRNA[i] == 'A' or newRNA[i] == 'U'):
                    interior = 'interiorAA' + str(length)
                else:
                    interior = "interiorGA" + str(length)
                interiorE = energyD[interior]
                rest = fpfold(newRNA[i:j+1], fpD)
                energy = interiorE + rest[0]
                if energy <= bestSoFar[0]:
                    pairings = [(i, j)]
                    pairings.extend(adjust(rest[1], i+1))
                    bestSoFar = (energy, pairings)
    return bestSoFar


def biLoop(RNA,fpD):
    newRNA = RNA[1:-1]
    if len(newRNA) < 2 * MINPAIRDIST + 2:
        return (float('inf'), [])
    bestSoFar = (float('inf'), [])
    for i in range(len(newRNA) - 2*MINPAIRDIST - 2):
        for j in range(i + MINPAIRDIST + 1, len(newRNA) - MINPAIRDIST - 2):
            for k in range(j + 1, len(newRNA) - MINPAIRDIST):
                for l in range(k + MINPAIRDIST, len(newRNA)):
                    firstRNA = newRNA[i:j+1]
                    secondRNA = newRNA[k:l+1]
                    firstRNAE = fpfold(firstRNA, fpD)
                    secondRNAE = fpfold(secondRNA, fpD)
                    energy = firstRNAE[0] + secondRNAE[0]
                    if energy <= bestSoFar[0]:
                        pairings = [(i, j), [k, l]]
                        pairings.extend(adjust(firstRNAE[1], i+1))
                        pairings.extend(adjust(secondRNAE[1], k+1))
                        bestSoFar = (energy, pairings)
    return bestSoFar


def toVienna(RNA,sol):
    '''Convert efold output to the string format used by Vienna
    RNAplot.'''
    outL=["."]*len(RNA)
    for pair in sol[1]:
        outL[pair[0]]='('
        outL[pair[1]]=')'
    return "".join(outL)

def adjust(pairsL, k):
    '''Add k to the coordinates in pairsL, returning a new list.'''
    newPairsL=[]
    for l,r in pairsL:
        newPairsL.append((l+k,r+k))
    return newPairsL

def energyFold(RNA):
    '''Wrapper to run energy based RNA folding.'''
    sol = efold(RNA,{},{})
    print(sol)
    print(RNA)
    print(toVienna(RNA,sol))
