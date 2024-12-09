# Maximum likelihood branch length calculation
import sys
from UnrootedTree import *
from branchesHelper import *




######


def branches(treeFN,dataFN,hillClimbDelta,outTreeFN):
    '''Wrapper function for maximum likelihood branch length optimization.'''

    tree3 = loadTree(treeFN) # topological tree
    dataD = loadData(dataFN)

    tree4 = preliminaryBrLen(dataD,tree3) # now with starter branch lengths
    utree = getUnrooted(tree4)

    finalLik,utree = optimizeBranches(utree,dataD,hillClimbDelta)

    print("Final branch lengths:")
    print(utree)

    print("Final Likelihood:", format(finalLik,".5e"))

    writeUTree(utree,outTreeFN)
    
    print("Tree (arbitrarily rooted) with branch lengths written to:",outTreeFN)

def singPosMluTree(tree4_1, tree4_2, position, branchLen, dataD):
    '''Helper function for mlutree. Finds the likelihood at a single position'''
    cases = []
    # Iterate over all possible combinations of nucleotides
    for X in ["A", "T", "G", "C"]:
        for Y in ["A", "T", "G", "C"]:
            # Find likelihoods for left tree, right tree, and jukes cantor on the branch
            tree1lik = singlePosition(tree4_1, X, position, dataD, {})
            tree2lik = singlePosition(tree4_2, Y, position, dataD, {})
            jukesLik = jukesBranch2Prob(branchLen)
            if X == Y:
                # If the nucleotides are the same, likelihood is calculated this way
                cases.append(0.25 * tree1lik * tree2lik * (1-jukesLik))
            else:
                cases.append(0.25 * tree1lik * tree2lik * 1/3 * jukesLik)
    return sum(cases)


def mlutree(tree4_1,tree4_2,dataD,branchLen):
    '''Finds the likelihood of the tree given the data across all positions'''
    sequencelist = list(dataD.values())
    likelihood = 1

    for position in range(len(sequencelist[0])):
        likelihood = likelihood * singPosMluTree(tree4_1, tree4_2, position, branchLen, dataD)

    return likelihood

def getTree4(utree, node, branch):
    '''Gets a rooted tree from an unrooted utree'''
    name = utree.getName(node)
    branchlen = utree.getBrlen(branch)
    if name != 'anc':
        # we're on a tip
        return (name, (), (), branchlen)
    else:
        # Get branches
        branch1, branch2, branch3 = utree.getBranches(node)
        branchList = [branch1, branch2, branch3]
        branchList.remove(branch)
        node1, node2 = utree.getNodes(branchList[0])
        node3, node4 = utree.getNodes(branchList[1])
        nodeList = [node1, node2, node3, node4]
        nodeList.remove(node)
        tree4_1 = getTree4(utree, nodeList[0], branchList[0])
        tree4_2 = getTree4(utree, nodeList[1], branchList[1])
        return (name, tree4_1, tree4_2, branchlen)
    
def singleBranchOptimize(utree, dataD, delta, branch):
    '''Helper func for optimizeBranches. Optimizes a single branch of the tree'''
    optimized = False
    branchLen = utree.getBrlen(branch)
    node1, node2 = utree.getNodes(branch)
    tree4_1 = getTree4(utree, node1, branch)
    tree4_2 = getTree4(utree, node2, branch)
    likelihood = 0
    while not optimized:
        branchLen1 = branchLen / delta
        branchLen2 = branchLen * delta
        likelihoodDown = mlutree(tree4_1,tree4_2,dataD,branchLen1)
        likelihoodUp = mlutree(tree4_1,tree4_2,dataD,branchLen2)
        # Longer branch better than shorter branch
        if (likelihoodUp > likelihoodDown):
            newlikelihood = likelihoodUp
            newbranchLen = branchLen2
        # Shorter branch better than longer branch
        else:
            newlikelihood = likelihoodDown
            newbranchLen = branchLen1
        # Only update branch length if likelihood improves
        if newlikelihood > likelihood:
            branchLen = newbranchLen
            likelihood = newlikelihood
            # End loop
        else:
            optimized = True
    return likelihood, branchLen



    
def optimizeBranches(utree,dataD,delta):
    # Fails test cases but is off by 0.00001 :(
    '''Given a tree, iterate over the branches, optimizing each branch until we findt he optimal tree'''
    optimized = False
    branches = utree.branches
    likelihood = 0

    while not optimized:
        # Run through all branches before testing new likelihood
        for branch in branches:
            newlikelihood, branchLen = singleBranchOptimize(utree, dataD, delta, branch)
            utree.setBrlen(branch, branchLen)
        # Only update likelihood if the likelihood improves
        if newlikelihood > likelihood:
            likelihood = newlikelihood
        # End loop
        else:
            optimized = True
    
    return (likelihood, utree)






#if __name__ == "__main__":



    tree3 = loadTree('example1.tre') # topological tree

    dataD = loadData('example1.data')



    tree4 = preliminaryBrLen(dataD,tree3) # now with starter branch lengths

    utree = getUnrooted(tree4)

