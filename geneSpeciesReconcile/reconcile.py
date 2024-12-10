# Parsimony algorithm for tree reconciliation

## test trees

# cospeciation with loss
gTreeA = ('i',('h',('f',('a',(),()),('b',(),())),('d',(),())),('e',(),()))
spTreeA =  ('I',('H',('G',('F',('A',(),()),('B',(),())),('C',(),())),('D',(),())),('E',(),()))
tipMapA = set([('a','A'),('b','B'),('d','D'),('e','E')])

# cospeciation, duplication
gTreeB = ('g',('e',('a',(),()),('b',(),())),('f',('c1',(),()),('c2',(),())))
spTreeB =  ('G',('E',('A',(),()),('B',(),())),('C',(),()))
tipMapB = set([('a','A'),('b','B'),('c1','C'),('c2','C')])

# cospeciation, transfer, loss
gTreeC = ('k',('h',('a',(),()),('c',(),())),('j',('m',('b',(),()),('d',(),())),('i',('e',(),()),('f',(),()))))
spTreeC = ('K',('H',('G',('A',(),()),('B',(),())),('C',(),())),('J',('I',('E',(),()),('F',(),())),('D',(),())))
tipMapC = set([('a','A'),('b','B'),('c','C'),('d','D'),('e','E'),('f','F')]) # host trans

# deep duplication with loss
gTreeD = ('j',('i',('f',('a',(),()),('b',(),())),('h',('g1',('d1',(),()),('e1',(),())),('c',(),()))),('g2',('d2',(),()),('e2',(),())))
spTreeD = ('I',('F',('A',(),()),('B',(),())),('H',('G',('D',(),()),('E',(),())),('C',(),())))
tipMapD = set([('a','A'),('b','B'),('c','C'),('d1','D'),('e1','E'),('d2','D'),('e2','E')])

# cospeciation with loss
gTreeE = ('h',('f',('a',(),()),('b',(),())),('g',('c',(),()),('e',(),())))
spTreeE =  ('Z',('W',('A',(),()),('B',(),())),('Y',('X',('D',(),()),('E',(),())),('C',(),())))
tipMapE = set([('a','A'),('b','B'),('c','C'),('e','E')])


## functions

def getTips(tree):
    '''Return list of tips of tree.'''
    if tree[1]==():
        return [tree[0]]
    else:
        return getTips(tree[1]) + getTips(tree[2])

def getInternalNodes(tree):
    '''Return list of internal nodes in tree.'''
    if tree[1]==():
        return []
    else:
        return getInternalNodes(tree[1])+getInternalNodes(tree[2]) + [tree[0]] 

def getAllNodes(tree):
    '''Return all nodes in tree'''
    if tree[1]==():
        return [tree[0]]
    else:
        return getAllNodes(tree[1])+getAllNodes(tree[2]) + [tree[0]] 

def subtree(node,tree):
    '''Return subtree of tree, rooted at node.'''
    if tree[0] == node:
        return tree
    elif tree[1] == ():
        return None # node isn't present
    else:
        lt=subtree(node,tree[1])
        rt=subtree(node,tree[2])
        
        if lt != None: return lt
        else: return rt

def initializeDpTables(gTree,spTree,tipMap):
    '''Given gene and species trees, create the dp tables in the form of
dicts. We create one which has gene tree nodes vs. species tree
nodes. And another which has gene tree nodes vs. species tree
branches. Populate with scores for the gene tree tips. Return as dict
of dicts, where the key 'n' has the node dict as value, and the key
'b' has the branch dict as value.
    '''
    nodeDpD = {}
    branchDpD = {}
    # do gTree tips
    for gNode in getTips(gTree):
        for spNB in getAllNodes(spTree):
            # regard spNB as a node
            if (gNode,spNB) in tipMap:
                nodeDpD[(gNode,spNB)] = (0,('t',))
            else:
                nodeDpD[(gNode,spNB)] = (float('inf'),())

            # regard spNB as a branch (gTree tips map to branches w/ inf cost)
            branchDpD[(gNode,spNB)] = (float('inf'),())

    dpD = {'n':nodeDpD,'b':branchDpD}
    
    return dpD

def printDp(gTree,spTree,dpD):
    '''Print out the dp table nicely.'''
    L=extractDp(gTree,spTree,dpD)
    printDpHelper(L)
    
def extractDp(gTree,spTree,dpD):
    '''Return dp table scores as list of lists..'''

    colNames=['']+getAllNodes(spTree) + getAllNodes(spTree)

    L=[colNames]
    
    for gNode in getTips(gTree) + getInternalNodes(gTree):
        rowL=[gNode]
        for spNB in getAllNodes(spTree):
            if (gNode,spNB) in dpD['n']:
                rowL.append(str(dpD['n'][(gNode,spNB)][0])) # regard as node
            else:
                rowL.append('')
        for spNB in getAllNodes(spTree):
            if (gNode,spNB) in dpD['n']:
                rowL.append(str(dpD['b'][(gNode,spNB)][0])) # regard as branch
            else:
                rowL.append('')
        L.append(rowL)
    return L
        
def printDpHelper(L,indent=0):
    '''Given tabular data in a list of lists (where sublists are rows)
print nicely so columns line up. Indent is an optional number of blank spaces to put in front of each row.'''
    # get max width for each column
    colMax=[]
    for col in range(len(L[0])):
        mx=0
        for row in L:
            if len(row[col]) > mx:
                mx = len(row[col])
        colMax.append(mx)
    
    # print
    for row in L:
        for col in range(len(row)):
            row[col]=row[col]+' ' * (colMax[col]-len(row[col]))
            
    for row in L:
        printStr = " "*indent + " | ".join(row)
        print(printStr.rstrip())

def getAncestors(spNode,spTree):
    '''Get list of ancestors nodes. Assume spNode is in spTree.'''
    
    if spNode == spTree[0]:
        return []
    elif spTree[1] == ():
        # hit a tip and spNode doesn't match
        return None
    else:
        lt=getAncestors(spNode,spTree[1])
        rt=getAncestors(spNode,spTree[2])

        # one or both must be None
        if lt!=None:
            return [spTree[0]] + lt
        elif rt!=None:
            return [spTree[0]] + rt
        else:
            return None

def getDescendants(node,tree):
    '''Return a list of all nodes in tree descending from node. Includes
node itself.'''
    st=subtree(node,tree)
    return getAllNodes(st)
    
def getPossibleTranferBranches(spNode,spTree):
    '''Given a node on the species tree, determine all the legal branches
one of its subnodes could horizontal transfer too. The restriction is,
it must not go to a branch in this same lineage, either forward
or backward.'''

    # these functions get nodes, but we can think of the branch
    # leading to each of those nodes
    allBrL=getAllNodes(spTree)
    parBrL=getAncestors(spNode,spTree)
    decBrL=getDescendants(spNode,spTree)

    branchL=[]
    for br in allBrL:
        if br not in parBrL and br not in decBrL:
            branchL.append(br)

    return branchL

def isTip(node,tree):
    '''Returns boolean indicating whether node is a tip.'''
    if tree[1] == ():
        if node==tree[0]: return True
        else: return False
    else:
        return isTip(node,tree[1]) or isTip(node,tree[2])

def lossCountHelper(tree,node):
    '''Count the number of nodes we pass to get to node starting at
root. Return tuple with that count, and tuple of loss nodes.'''
    if tree[0] == node:
        return (0,()) # we're not passing a node here
    elif tree[1]==():
        return (-float('inf'),())
    else:
        lt=lossCountHelper(tree[1],node)
        rt=lossCountHelper(tree[2],node)

        if lt[0] > rt[0]:
            lossT= lt[1] + (tree[0],) # add this one
            return (lt[0]+1,lossT)
        else:
            lossT= rt[1] + (tree[0],)
            return (rt[0]+1,lossT)

def lossCount(tree,parent,child,nb):
    '''Given tree, and a node placement of the parent at parent and a
child at child, figure out how many loss events are implied. nb is 'n'
or 'b' and tells us whether parent represents a node or the branch
leading to that node. Return tuple with count, and tuple of loss
nodes.
    '''
    if parent == child:
        return (0,())
    else:
        subtr=subtree(parent,tree)
        if nb=='n':
            # parent is a node, we don't want to count the root of
            # subtr itself, since that's parent
            lt = lossCountHelper(subtr[1],child)
            rt = lossCountHelper(subtr[2],child)
            if lt[0] > rt[0]: return lt
            else: return rt
        else:
            # parent is branch leading to node, so we should count it
            return lossCountHelper(subtr,child)
            
def findMinCostRootPlacement(root,spTree,dpD):
    '''Given filled out dp tables, return the minimum cost placement for
root root.
    '''
    rootPlacement = ()
    minCost=float('inf')
    for spNB in getAllNodes(spTree):
        # regard spNB as node
        if dpD['n'][(root,spNB)][0] < minCost:
            rootPlacement = (root,spNB,'n')
            minCost = dpD['n'][(root,spNB)][0]
        # now regard spNB as branch
        if dpD['b'][(root,spNB)][0] < minCost:
            rootPlacement = (root,spNB,'b')
            minCost = dpD['b'][(root,spNB)][0]

    return rootPlacement,minCost

def reconcile(gTree,spTree,tipMap,dup,trans,loss):
    '''Reconcile a gene tree with a species tree.'''

    dpD=initializeDpTables(gTree,spTree,tipMap)

    # fill dp table
    for gNode in getInternalNodes(gTree):
        for spNB in getAllNodes(spTree):
            # regard spNB as node
            dpD['n'][(gNode,spNB)]=scoreInternalGeneNodeToSpNode(gTree,spTree,dpD,loss,gNode,spNB)
            # now regard spNB as branch
            dpD['b'][(gNode,spNB)]=scoreInternalGeneNodeToSpBranch(gTree,spTree,dpD,dup,trans,loss,gNode,spNB)
            
    rootPlacement,minCost = findMinCostRootPlacement(gTree[0],spTree,dpD)
    
    return dpD,rootPlacement,minCost


def scoreInternalGeneNodeToSpNode(gTree,spTree,dpD,loss,gNode,spNode):
    '''Intakes a gene tree (gTree), a species tree (spTree), a dp table, loss cost, and the gNode being placed into spNode
    Returns the minimum cost and the solution tuple for where the subnodes should be placed'''
    #solution tuple format: (left node placement, right node placement, event, losses)
    solutionTuple = ()
    # Get the gene tree at the current node
    geneTree = subtree(gNode, gTree)
    gLeftNode = geneTree[1][0] 
    gRightNode = geneTree[2][0]
    # Get the species tree at the current node
    speciesTree = subtree(spNode, spTree)
    spLeftTree = speciesTree[1]
    spRightTree = speciesTree[2]
    # We are on a tip
    if spLeftTree == () or spRightTree == ():
        return (float('inf'), solutionTuple)
    # Get all nodes in left and right species trees
    spLeftTreeNodes = getAllNodes(spLeftTree)
    spRightTreeNodes = getAllNodes(spRightTree)
    costs = {}
    # Iterate through all possible placements of the gene subnodes
    for spRightNode in spRightTreeNodes:
        for spLeftNode in spLeftTreeNodes:
            # Find the costs and losses associated with placements
            rightLossCount, rightLossNodes = lossCount(spTree, spNode, spRightNode, 'n')
            leftLossCount, leftLossNodes = lossCount(spTree, spNode, spLeftNode, 'n')
            lossCost = loss * (rightLossCount + leftLossCount)
            lossNodes = rightLossNodes + leftLossNodes
            # There are 8 cases for the placements of the left and right nodes onto the sp tree nodes
            costs[(('n', gLeftNode, 'n', spLeftNode), ('n', gRightNode, 'n', spRightNode), "Cospeciation", lossNodes)] = dpD["n"][(gRightNode, spRightNode)][0] + dpD["n"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spRightNode), ('n', gRightNode, 'n', spLeftNode), "Cospeciation", lossNodes)] = dpD["n"][(gRightNode, spLeftNode)][0] + dpD["n"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spLeftNode), ('n', gRightNode, 'b', spRightNode), "Cospeciation", lossNodes)] = dpD["b"][(gRightNode, spRightNode)][0] + dpD["n"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spRightNode), ('n', gRightNode, 'b', spLeftNode),"Cospeciation", lossNodes)] = dpD["b"][(gRightNode, spLeftNode)][0] + dpD["n"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spRightNode), ('n', gRightNode, 'n', spLeftNode), "Cospeciation", lossNodes)] = dpD["n"][(gRightNode, spLeftNode)][0] + dpD["b"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spLeftNode), ('n', gRightNode, 'n', spRightNode), "Cospeciation", lossNodes)] = dpD["n"][(gRightNode, spRightNode)][0] + dpD["b"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spRightNode), ('n', gRightNode, 'b', spLeftNode), "Cospeciation", lossNodes)] = dpD["b"][(gRightNode, spLeftNode)][0] + dpD["b"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spLeftNode), ('n', gRightNode, 'b', spRightNode), "Cospeciation", lossNodes)] = dpD["b"][(gRightNode, spRightNode)][0] + dpD["b"][(gLeftNode, spLeftNode)][0] + lossCost
    # Iterate throught the dictionary to find the minimum cost
    mincost = float('inf')
    for solntuple, cost in costs.items():
        if cost <= mincost:
            mincost = cost
            solutionTuple = solntuple
    if mincost == float('inf'):
        solutionTuple = ()
    cost = mincost
    return (cost, solutionTuple)

def scoreInternalGeneNodeToSpBranch(gTree,spTree,dpD,dup,trans,loss,gNode,spBranch):
    '''Intakes a gene tree (gTree), a species tree (spTree), a dp table, loss cost, a duplication cost, a transfer cost, and the gNode being placed into sBranch.
    Since horizontal transfer and duplication are the only outcomes possible, it only considers these.
    Returns the minimum cost and the solution tuple for where the subnodes should be placed'''
    #solution tuple format: (left node placement, right node placement, event, losses)
    solutionTuple = ()
    # Get gene tree and species tree subtrees at the branch we are looking at
    geneTree = subtree(gNode, gTree)
    speciesTree = subtree(spBranch, spTree)
    # Get left and right gene tree nodes
    gLeftNode = geneTree[1][0]
    gRightNode = geneTree[2][0]
    # Get all descen
    descendantNodes = getDescendants(spBranch, spTree)
    transferBranches = getPossibleTranferBranches(spBranch, spTree)
    if descendantNodes == [] or transferBranches == []:
        transferCost = float('inf')
        transSolutionTuple = ()
    else:
        costs = {}
        # Horizontal transfer
        for descendantNode in descendantNodes:
            for transferBranch in transferBranches:
                lossCost, lossNodes = lossCount(spTree, spBranch, descendantNode, 'b')
                lossCost = lossCost * loss
                costs[(('n', gLeftNode, 'n', descendantNode, 'd'), ('n', gRightNode, 'n', transferBranch, 't'), "Horizontal Transfer", lossNodes)] = dpD['n'][(gLeftNode, descendantNode)][0] + dpD['n'][(gRightNode, transferBranch)][0] + lossCost
                costs[(('n', gLeftNode, 'n', transferBranch, 't'), ('n', gRightNode, 'n', descendantNode, 'd'), "Horizontal Transfer", lossNodes)] = dpD['n'][(gLeftNode, transferBranch)][0] + dpD['n'][(gRightNode, descendantNode)][0] + lossCost
                costs[(('n', gLeftNode, 'n', descendantNode, 'd'), ('n', gRightNode, 'b', transferBranch, 't'), "Horizontal Transfer", lossNodes)] = dpD['n'][(gLeftNode, descendantNode)][0] + dpD['b'][(gRightNode, transferBranch)][0] + lossCost
                costs[(('n', gLeftNode, 'n', transferBranch, 't'), ('n', gRightNode, 'b', descendantNode, 'd'), "Horizontal Transfer", lossNodes)] = dpD['n'][(gLeftNode, transferBranch)][0] + dpD['b'][(gRightNode, descendantNode)][0] + lossCost
                costs[(('n', gLeftNode, 'b', descendantNode, 'd'), ('n', gRightNode, 'n', transferBranch, 't'), "Horizontal Transfer", lossNodes)] = dpD['b'][(gLeftNode, descendantNode)][0] + dpD['n'][(gRightNode, transferBranch)][0] + lossCost
                costs[(('n', gLeftNode, 'b', transferBranch, 't'), ('n', gRightNode, 'n', descendantNode, 'd'), "Horizontal Transfer", lossNodes)] = dpD['b'][(gLeftNode, transferBranch)][0] + dpD['n'][(gRightNode, descendantNode)][0] + lossCost
                costs[(('n', gLeftNode, 'b', descendantNode, 'd'), ('n', gRightNode, 'b', transferBranch, 't'), "Horizontal Transfer", lossNodes)] = dpD['b'][(gLeftNode, descendantNode)][0] + dpD['b'][(gRightNode, transferBranch)][0] + lossCost
                costs[(('n', gLeftNode, 'b', transferBranch, 't'), ('n', gRightNode, 'b', descendantNode, 'd'), "Horizontal Transfer", lossNodes)] = dpD['b'][(gLeftNode, transferBranch)][0] + dpD['b'][(gRightNode, descendantNode)][0] + lossCost
        mincost = float('inf')
        for solntuple, cost in costs.items():
            if cost <= mincost:
                mincost = cost
                transSolutionTuple = solntuple
        if mincost == float('inf'):
            transSolutionTuple = ()
        transferCost = mincost + trans

    # Duplication
    spLeftTree = speciesTree[1]
    spRightTree = speciesTree[2]
    # Find the nodes included in the left and right species trees
    if spLeftTree == () and spRightTree == ():
        spLeftTreeNodes = [spBranch]
        spRightTreeNodes = [spBranch]
    elif spRightTree == (): 
        spLeftTreeNodes = getAllNodes(spLeftTree) + [spBranch]   
        spRightTreeNodes = [spBranch]
    elif spLeftTree == ():
        spLeftTreeNodes = [spBranch]
        spRightTreeNodes = getAllNodes(spRightTree) + [spBranch]
    else:
        spLeftTreeNodes = getAllNodes(spLeftTree) + [spBranch]
        spRightTreeNodes = getAllNodes(spRightTree) + [spBranch]
    costs = {}
    # Find cost for placing the left and right gene nodes in all possible sp tree nodes
    for spRightNode in spRightTreeNodes:
        for spLeftNode in spLeftTreeNodes:
            # Find loss cost and loss nodes
            leftLossCost, leftLossNodes = lossCount(spTree, spBranch, spLeftNode, 'b')
            rightLossCost, rightLossNodes = lossCount(spTree, spBranch, spRightNode, 'b')
            lossCost = loss * (leftLossCost + rightLossCost)
            lossNodes = leftLossNodes + rightLossNodes
            # 8 cases
            costs[(('n', gLeftNode, 'n', spLeftNode), ('n', gRightNode, 'n', spRightNode), "Duplication", lossNodes)] = dpD["n"][(gRightNode, spRightNode)][0] + dpD["n"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spRightNode), ('n', gRightNode, 'n', spLeftNode), "Duplication", lossNodes)] = dpD["n"][(gRightNode, spLeftNode)][0] + dpD["n"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spLeftNode), ('n', gRightNode, 'b', spRightNode), "Duplication", lossNodes)] = dpD["b"][(gRightNode, spRightNode)][0] + dpD["n"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'n', spRightNode), ('n', gRightNode, 'b', spLeftNode),"Duplication", lossNodes)] = dpD["b"][(gRightNode, spLeftNode)][0] + dpD["n"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spRightNode), ('n', gRightNode, 'n', spLeftNode), "Duplication", lossNodes)] = dpD["n"][(gRightNode, spLeftNode)][0] + dpD["b"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spLeftNode), ('n', gRightNode, 'n', spRightNode), "Duplication", lossNodes)] = dpD["n"][(gRightNode, spRightNode)][0] + dpD["b"][(gLeftNode, spLeftNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spRightNode), ('n', gRightNode, 'b', spLeftNode), "Duplication", lossNodes)] = dpD["b"][(gRightNode, spLeftNode)][0] + dpD["b"][(gLeftNode, spRightNode)][0] + lossCost
            costs[(('n', gLeftNode, 'b', spLeftNode), ('n', gRightNode, 'b', spRightNode), "Duplication", lossNodes)] = dpD["b"][(gRightNode, spRightNode)][0] + dpD["b"][(gLeftNode, spLeftNode)][0] + lossCost
    # Find minimum duplication cost
    mincost = float('inf')
    for solntuple, cost in costs.items():
        if cost <= mincost:
            mincost = cost
            dupSolutionTuple = solntuple
    if mincost == float('inf'):
        dupSolutionTuple = ()
    dupCost = mincost + dup

    # Compare transfer to duplication
    if transferCost < dupCost:
        cost = transferCost
        solutionTuple = transSolutionTuple
    else:
        cost = dupCost
        solutionTuple = dupSolutionTuple
    return (cost, solutionTuple)

def bt(placement, dpD, list):
    '''Backtrack through solution tuples'''
    # placement = (geneNode, speciesPosition, nb)
    if placement[2] == 'b':
        cost, solutionTuple = dpD['b'][placement[0], placement[1]]
    else:
        cost, solutionTuple = dpD['n'][placement[0], placement[1]]

    if solutionTuple == ('t',):
        return []

    # Get information from solution tuple    
    leftGeneNode = solutionTuple[0][1]
    leftSpeciesPos = solutionTuple[0][2]
    leftSpeciesPosName = solutionTuple[0][3]
    rightGeneNode = solutionTuple[1][1]
    rightSpeciesPos = solutionTuple[1][2]
    rightSpeciesPosName = solutionTuple[1][3]
    event = solutionTuple[2]
    losses = solutionTuple[3]

    # Fill out list to be passed back into bt
    if event == "Horizontal Transfer":
        list.append([placement, (leftGeneNode, leftSpeciesPosName, leftSpeciesPos, solutionTuple[0][4]), (rightGeneNode, rightSpeciesPosName, rightSpeciesPos, solutionTuple[1][4]), event, losses])
    else:
        list.append([placement, (leftGeneNode, leftSpeciesPosName, leftSpeciesPos), (rightGeneNode, rightSpeciesPosName, rightSpeciesPos), event, losses])

    return bt((leftGeneNode, leftSpeciesPosName, leftSpeciesPos), dpD, list),  bt((rightGeneNode, rightSpeciesPosName, rightSpeciesPos), dpD, list)

def printReconciliation(dpD,rootPlacement,minCost):
    list = []
    bt(rootPlacement, dpD, list)
    print("Minimum cost:", minCost)
    for i in range(len(list)):
        x = list[i]
        initNode = x[0]
        if initNode[2] == 'n':
            loc = "Node"
        else:
            loc = "Branch"
        leftNode = x[1]
        if leftNode[2] == 'n':
            leftloc = "Node"
        else:
            leftloc = "Branch"
        rightNode = x[2]
        if rightNode[2] == 'n':
            rightloc = "Node"
        else:
            rightloc = "Branch"
        event = x[3]
        losses = x[4]
        if event == "Horizontal Transfer":
            if leftNode[3] == 'd':
                leftEvent = 'normal descent to'
                rightEvent = 'transfer to'
            else:
                rightEvent = 'normal descent to'
                leftEvent = 'transfer to'
        else:
            leftEvent = ''
            rightEvent = ''
        print(event, ": (Node ", initNode[0], ", ", loc, initNode[1], ") -->", leftEvent, "(node ", leftNode[0], ", ", leftloc, leftNode[1], ")", rightEvent, "(node ", rightNode[0], ", ", rightloc, rightNode[1], ")")
        if losses != ():
            for loss in losses:
                print("gene loss at node", loss)

    
            


