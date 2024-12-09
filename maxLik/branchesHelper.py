# Support functions for maximum likelihood branch length calculation
# Eliot Bush, Aug 2016
import math
from UnrootedTree import *
from branches import *

## load

def middleComma(s):
    '''Find the 'middle' comma. This is the one which has had the same
number of ( as ) before it. Return the index.'''
    parenBalance=0
    for i in range(len(s)):
        if s[i]=='(':
            parenBalance+=1
        elif s[i]==')':
            parenBalance-=1
        elif s[i]==',' and parenBalance==0:
            return i
    return None

    
def topologicalNewickToTree3(s):
    '''Given a Newick type string representation of a tree that lacks
branch lengths, convert to tuple tree3 format and return.
    '''
    if "(" not in s and ")" not in s:
        # its a tip
        return(s,(),())
    else:
        s=s[1:-1] # strip of outer perens
        mci=middleComma(s)
        leftTree=topologicalNewickToTree3(s[:mci].strip()) # strip any leading or
        rightTree=topologicalNewickToTree3(s[mci+1:].strip()) # trailing white space
        return('anc',leftTree,rightTree)

def loadTree(treeFN):
    '''Load a tree in newick form from file. Assumes no branch lengths given.'''
    f = open(treeFN,'r')
    s = f.read().rstrip()[:-1] # :-1 so we get rid of semicolon
    f.close()
    tree = topologicalNewickToTree3(s)
    return tree

def loadData(dataFN):
    '''Load sequence data from file, return as dictionary. Assumes one
species per line. Name, white space, sequence.'''
    f = open(dataFN,'r')
    D = {}
    while True:
        s = f.readline()
        if s == '':
            break
        species,seq = s.split()
        D[species] = seq
    f.close()
    return D


## Preliminary branch lengths

def jukesProb2Branch(propDiff):
    """Correction for multiple hits using Jukes-Cantor model. Goes from
proportion of changes to a branch length."""
    if propDiff < 0:
        raise ValueError("jukesProb2Branch was passed a negative number, which it doesn't know how to handle.")
    elif propDiff == 0:
        return 0
    elif propDiff < 0.75:
        return -3.0*math.log(1-(4.0*propDiff/3))/4
    else:
        raise ValueError("jukesProb2Branch was passed a number >= 0.75, which is too big.")


def propDifferent(seqA,seqB):
    """Calculate the proporton of sites which differ between aligned seqs
A and B. Ignore sites with a gap in one sequence or the other."""
    numNonGap=0
    numMisMatch=0
    for i in range(len(seqA)):
        if seqA[i]!="-" and seqB[i]!="-":
            numNonGap+=1
            if seqA[i]!=seqB[i]:
                numMisMatch+=1
    return(float(numMisMatch)/numNonGap)


def makeMatrix(dataD):
    '''Given a data matrix keyed by species name and with sequences as
values, calculate the Jukes-Cantor distance between all
sequences. Store in dictionary, keyed by tuples of two species
names.'''
    # make list of nodes, each being a tip of tree4
    nodeL=[]
    for sp in dataD:
        nodeL.append((sp,(),(),0))
    # distances
    distD={}
    for sp1 in dataD:
        sp1TupTip = (sp1,(),(),0)
        for sp2 in dataD:
            sp2TupTip = (sp2,(),(),0)
            dist=jukesProb2Branch(propDifferent(dataD[sp1],dataD[sp2]))
            distD[(sp1TupTip,sp2TupTip)] = dist
            distD[(sp2TupTip,sp1TupTip)] = dist
    return nodeL,distD

def createMergeOrderL(tree3):
    '''Given a three tuple type tree (without branch lengths) create a
list which shows the order in which one shuld merge nodes to produce
this tree, starting from the tips.'''
    if tree3[1] == ():
        return []
    else:
        leftL = createMergeOrderL(tree3[1])
        rightL = createMergeOrderL(tree3[2])

        thisNode = (tree3[1],tree3[2])

        return leftL + rightL + [thisNode]

def isomorphic(tree1,tree2):
    '''Takes two rooted tuple type trees, and returns True if they are
isomorphic in their topology. Accepts both 3 and 4 tuple trees.
    '''
    if tree1[0:3] == tree2[0:3]: return True # slice in case 4 tuple tree
    elif tree1[1] == () or tree2[1] == (): return False
    else:
        if isomorphic(tree1[1],tree2[1]) and isomorphic(tree1[2],tree2[2]):
            return True
        elif isomorphic(tree1[1],tree2[2]) and isomorphic(tree1[2],tree2[1]):
            return True
        else: return False

def findMatchingTree(node3,nodeL):
    '''Given a list of nodes, nodeL, containing trees in tree4 format,
find the one which mateches node3, which is in three tuple
format. Return the matching node.'''

    for node4 in nodeL:
        if isomorphic(node3,node4):
            return node4
    return None # this shouldn't happen
    
    
# NJ based functions

def nodeSep(node,nodeL,distD):
    '''A measure of the separation between a node of interest and
    the other nodes.'''
    distSum=0
    for iternode in nodeL:
        if node != iternode:
            distSum+=distD[(node,iternode)]
    return(float(distSum)/(len(nodeL)-2))

def branchLength(nodeA,nodeB,nodeL,distD):
    '''Takes two nodes we're planning to merge, nodeA and nodeB. Calculates the
branch lengths from their common ancestor to each.'''
    dist=distD[(nodeA,nodeB)]
    sepA=nodeSep(nodeA,nodeL,distD)
    sepB=nodeSep(nodeB,nodeL,distD)
    branchA=0.5*(dist+(sepA-sepB))
    branchB=0.5*(dist+(sepB-sepA))
    return(branchA,branchB)

def mergeNodes(nodeA,nodeB,branchLenA,branchLenB):
    '''Given two nodes to be merged, update their branch lengths, then merge.'''
    updatedNdA=(nodeA[0],nodeA[1],nodeA[2],branchLenA)
    updatedNdB=(nodeB[0],nodeB[1],nodeB[2],branchLenB)
    return ('anc',updatedNdA,updatedNdB,0)

def updateDistances(nodeA,nodeB,newNode,nodeL,distD):
    '''Updates the distance dictionary distD, given a newly merged
node. Modifies distD in place and returns None.'''
    for node in nodeL:
        if node!=nodeA and node!=nodeB:
            newdist=(distD[nodeA,node]+distD[nodeB,node]-distD[nodeA,nodeB])*0.5
            distD[(newNode,node)]=newdist
            distD[(node,newNode)]=newdist

def preliminaryBrLen(dataD,guideTree3):
    '''Determined a set of preliminary branch lengths. Uses topology from
guide tree and a neighbor joining style greedy approach. Returns tree in
four tuple format.'''


    nodeL,distD = makeMatrix(dataD)
    mergeOrderL = createMergeOrderL(guideTree3)

        
    for leftNode3,rightNode3 in mergeOrderL[:-1]: # do last merge outside loop

        # get the nodes in nodeL to merge by finding the ones
        # that match the next pair in mergeOrderL
        leftNode4 = findMatchingTree(leftNode3,nodeL)
        rightNode4 = findMatchingTree(rightNode3,nodeL)        


        lenA,lenB=branchLength(leftNode4,rightNode4,nodeL,distD)
        mergedNode=mergeNodes(leftNode4,rightNode4,lenA,lenB)
        
        updateDistances(leftNode4,rightNode4,mergedNode,nodeL,distD)

        nodeL.remove(leftNode4) # remove old nodes and add new
        nodeL.remove(rightNode4)
        nodeL.append(mergedNode)

    # merge last two nodes
    leftNode4,rightNode4=nodeL
    finalDist=distD[nodeL[0],nodeL[1]]
    finalTree = mergeNodes(leftNode4,rightNode4,finalDist/2.0,finalDist/2.0)
    return(finalTree)


## Making an unrooted tree

def getUnrootedFromNode(tree4,parentNode,branchCounter,nodeCounter,nodeNameD,nodeBranchD,branchNodeD,branchLenD):
    '''Given a four tuple tree, update the dictionaries
nodeNameD,nodeBranchD,branchNodeD,branchLenD so as to specify an
unrooted tree. parentNode is the number of the parent node for the
tree.'''

    thisBranch = branchCounter
    thisNode = nodeCounter
    
    # add for this branch and node. same stuff if tip or internal
    branchNodeD[thisBranch] = [thisNode,parentNode]
    branchLenD[thisBranch] = tree4[3]
    nodeBranchD[thisNode] = [thisBranch]
    nodeBranchD[parentNode].append(thisBranch) # add this branch to parent node
    nodeNameD[thisNode] = tree4[0]
    branchCounter += 1
    nodeCounter += 1

    #  only recurse if its not a tip
    if not tree4[1] == (): # not a tip
        # left
        branchCounter,nodeCounter = getUnrootedFromNode(tree4[1],thisNode,branchCounter,nodeCounter,nodeNameD,nodeBranchD,branchNodeD,branchLenD)
        # right
        branchCounter,nodeCounter = getUnrootedFromNode(tree4[2],thisNode,branchCounter,nodeCounter,nodeNameD,nodeBranchD,branchNodeD,branchLenD)

        
    return branchCounter,nodeCounter

def tupleIzeDict(D):
    '''For a dictionary with values consisting of lists, convert to tuples.'''
    newD = {}
    for key in D:
        newD[key] = tuple(D[key])

    return newD
    
def getUnrooted(tree4):
    '''Takes a four tuple tree and constructs a dictionary based
(unrooted) representation of it.
    '''
    nodeNameD = {}
    nodeBranchD = {}
    branchNodeD = {}
    branchLenD = {}

    rootNode = 0
    tempRootBranch = 0 # an branch attached to this root node we'll remove later
    
    branchCounter = 0
    nodeCounter = 1

    nodeBranchD[rootNode] = []
    
    branchCounter,nodeCounter = getUnrootedFromNode(tree4[1],rootNode,branchCounter,nodeCounter,nodeNameD,nodeBranchD,branchNodeD,branchLenD)

    branchCounter,nodeCounter = getUnrootedFromNode(tree4[2],rootNode,branchCounter,nodeCounter,nodeNameD,nodeBranchD,branchNodeD,branchLenD)


    # now edit these dictionaries to remove the root node and the
    # tempRootBranch branch attached to that root node. And update the
    # connections accordingly

    branch1,branch2 = nodeBranchD[rootNode]
    branchToKeep = branch1 if branch2==tempRootBranch else branch2

    # must update branchToKeep entry in branchNodeD
    node1,node2 = branchNodeD[tempRootBranch]
    nodeToReplaceWith = node1 if node2==rootNode else node2

    node1,node2 = branchNodeD[branchToKeep]
    if node1 == rootNode:
        branchNodeD[branchToKeep] = [nodeToReplaceWith,node2]
    else:
        branchNodeD[branchToKeep] = [node1,nodeToReplaceWith]

    # make branchToKeep entry in branchLenD longer by the length of tempRootBranch
    branchLenD[branchToKeep] += branchLenD[tempRootBranch]
        
    # update the nodeToReplaceWith entry in nodeBranchD. Take out
    # tempRootBranch and put in branchToKeep
    ind = nodeBranchD[nodeToReplaceWith].index(tempRootBranch)
    nodeBranchD[nodeToReplaceWith][ind] = branchToKeep

    # now get rid of rootNode and tempRootBranch entries in nodeBranchD,
    # branchNodeD, and branchLenD 
    del nodeBranchD[rootNode]
    del branchNodeD[tempRootBranch]
    del branchLenD[tempRootBranch]

    branchNodeD = tupleIzeDict(branchNodeD)
    nodeBranchD = tupleIzeDict(nodeBranchD)
    
    return UnrootedTree(nodeNameD,nodeBranchD,branchNodeD,branchLenD)

## ML

def jukesBranch2Prob(branchLength):
    '''Use jukes cantor model to go from a branch length to a proportion
of changes (or a probability).'''
    return 3.0/4.0*(1 - math.exp(-4.0/3.0 * branchLength))

def singlePosition(tree, rootChar, position, dataD, memo):
    '''Calculate the likelihood for a single position, given data in
dataD and a particular base at the root (rootChar). tree is a
rooted four tuple type tree. By RLH, with tweaks by EB.
    '''
    rootName = tree[0]
    leftTree = tree[1]
    rightTree = tree[2]
    branchLength = tree[3]
    if len(leftTree) == 0:  # Tip
        if rootChar == dataD[rootName][position]:
            return 1.0
        else:
            return 0.0
    elif (tree,rootChar,position) in memo:
        return memo[(tree,rootChar,position)]
    else:
        score = 0.0
        leftBranchLength = leftTree[3]
        rightBranchLength = rightTree[3]
        
        for leftChar in ["A", "T", "C", "G"]:
            for rightChar in ["A", "T", "C", "G"]:
                leftScore = singlePosition(leftTree, leftChar, position, dataD, memo)
                probObservedChange = jukesBranch2Prob(leftBranchLength)
                if rootChar == leftChar:
                    leftScore *= (1.0-probObservedChange)
                else:
                    leftScore *= (1.0/3.0) * probObservedChange
            

                rightScore = singlePosition(rightTree, rightChar, position, dataD, memo)
                probObservedChange = jukesBranch2Prob(rightBranchLength)
                if rootChar == rightChar:
                    rightScore *= (1.0 - probObservedChange)
                else:
                    rightScore *= (1.0/3.0) * probObservedChange
                score += leftScore * rightScore

        memo[(tree,rootChar,position)] = score
                    
        return score

## Output
    
def tree4ToNewick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string.'''
    if tree[1]==():
        return tree[0]+":"+str(tree[3])
    else:
        leftStr=tree4ToNewick(tree[1])
        rightStr=tree4ToNewick(tree[2])
        return "("+leftStr+","+rightStr+"):"+str(tree[3])

def utree2tree4(utree):
    '''Given a utree object, convert to tree4 format. Makes use of
(student written) getTree4. We arbitrarily root it.'''

    branch = utree.branches[0] # arbitrarily pick a branch to root it on
    nodesT = utree.getNodes(branch)

    ltree4 = getTree4(utree,nodesT[0],branch)
    ltree4 = (ltree4[0],ltree4[1],ltree4[2],utree.getBrlen(branch)/2.0)
    

    rtree4 = getTree4(utree,nodesT[1],branch)
    rtree4 = (rtree4[0],rtree4[1],rtree4[2],utree.getBrlen(branch)/2.0)
    
    outTree = ('anc',ltree4,rtree4,0)
    return outTree
    
def writeUTree(utree,fileName):
    '''Write utree to fileName (in newick format).'''
    tree4 = utree2tree4(utree)

    # create left and right tree strings separately, then
    # comebine. this avoids having a zero length branch at the end,
    # which our load functions don't like.

    outStr = "("+tree4ToNewick(tree4[1])+","+tree4ToNewick(tree4[2])+");"
    
    f=open(fileName,"w")
    f.write(outStr)
    f.close()


