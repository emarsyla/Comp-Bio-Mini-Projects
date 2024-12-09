import numpy, random
from scipy import stats
from lacData import *
from alleleHist import *

def timeToNextCoalescence(numAlleles,popSize):

    """Returns a geometric random variable, the time to the next

    coalescence event (in generations)."""

    p = (float(numAlleles) * (numAlleles-1)) / (4*popSize)

    return int(stats.geom.rvs(p))

def randomCoalescentTree(numAlleles,popSize):

    """Returns a random coalescent tree."""

    # Initialize the nodes of the tree
    allelesL = [(i, (), (), 0) for i in range(numAlleles)]
    nodeNum = numAlleles - 1
    timetocoalesce = 0

    # Choose two random nodes to combine until there is a fully combined tree.
    while len(allelesL) > 1:
        nodeNum += 1
        length = len(allelesL)
        randInts = random.sample(list(range(length)), 2)
        node1 = allelesL[randInts[0]]
        node2 = allelesL[randInts[1]]
        timetocoalesce = timetocoalesce + timeToNextCoalescence(length, popSize)
        newNode = [(nodeNum, node1, node2, timetocoalesce)]
        allelesL.remove(node1)
        allelesL.remove(node2)
        allelesL = allelesL + newNode
    return allelesL[0]
        
def convertToGensPerBranch(Tree, parentGensBeforePresent):

    """Converts variant 1 tree to variant 2 tree"""

    # Base case: tree is empty
    newTree = ()
    if Tree == ():
        return newTree
    
    # otherwise, find branchGens and then recurse through the tree
    else:
        nodeNum, leftTree, rightTree, totGensBeforePresent = Tree
        branchGens = parentGensBeforePresent - totGensBeforePresent
        newTree = (nodeNum, convertToGensPerBranch(leftTree, totGensBeforePresent), convertToGensPerBranch(rightTree, totGensBeforePresent), branchGens)
        return newTree

def getBranchProportionList(Tree):

    """Returns a list of the node numbers proportional to the length of the branches"""

    propList = []
    # Base case tree is nan empty tree
    if Tree == ():
        return propList
    
    # Otherwise add nto propList the nodeNum branchGens times, then rescurse
    else:
        nodeNum, leftTree, rightTree, branchGens = Tree
        for i in range(branchGens):
            propList.append(nodeNum)
        propList.extend(getBranchProportionList(leftTree))
        propList.extend(getBranchProportionList(rightTree))  
        return propList
    
def findNumNodes(Tree):
    """Helper function for assignMutsToBranch. Finds the number of nodes in the tree"""

    if not Tree or Tree == ():
        return 0
    
    else:
        nodeNum, leftTree, rightTree, branchGens = Tree
        left_count = findNumNodes(leftTree)
        right_count = findNumNodes(rightTree)
        return 1 + left_count + right_count
    
def assignMutsToBranch(Tree, numMuts):

    """Assigns mutations randomly to branches of the tree"""

    size = findNumNodes(Tree)
    mutList = [[] for _ in range(size)]
    propList = getBranchProportionList(Tree)
    for i in range(numMuts):
        index = random.choice(propList)
        mutList[index].append(i)
    return mutList
    
def createSeqs(Tree, branchMutL, seqT):

    """Creates sequences or haplotypes based on the mutations along the branches"""

    nodeNum, leftTree, rightTree, branchGens = Tree
    tupList = []
    muts = branchMutL[nodeNum]
    for mut in muts:
        seqT += (mut,)

    # Base case: Tree is a leaf node
    if leftTree == () and rightTree == ():
        tupList.append(seqT)
        return tupList

    # Otherwise, extend the tupList
    else:
        tupList.extend(createSeqs(leftTree, branchMutL, seqT))
        tupList.extend(createSeqs(rightTree, branchMutL, seqT))
        return tupList

def coalSim(numAlleles,popSize,numMuts):

    """Call the prvious functions for one coalescent simulation"""

    Tree = randomCoalescentTree(numAlleles,popSize)
    timeOfFinalCoalescence=Tree[-1]
    Tree = convertToGensPerBranch(Tree, timeOfFinalCoalescence)
    branchMutL = assignMutsToBranch(Tree, numMuts)
    return createSeqs(Tree, branchMutL, ())

def runCoalSim(popSize,numAlleles,numMuts,numReps,maxHapCountRealData):
    
    """wrapper function to run mutliple simulations"""

    simList = []
    for i in range(numReps):
        simulation = coalSim(numAlleles,popSize,numMuts)
        unique = list(set(simulation))
        simHapCounts = []
        for elem in unique:
            simHapCounts.append(len([tup for tup in simulation if tup == elem]))
        simList.append(max(simHapCounts))
    return len([count for count in simList if count > maxHapCountRealData])

