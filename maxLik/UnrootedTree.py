# UnrootedTree class
# Bio 188
# Eliot Bush


class UnrootedTree:

    def __init__(self,nodeNameD,nodeBranchD,branchNodeD,branchLenD = None):
        '''Create an UnrootedTree object.'''

        # nodeBranchD and nodeNameD should have the same keys
        if sorted(nodeBranchD.keys()) != sorted(nodeNameD.keys()):
            raise ValueError("nodeBranchD and nodeNameD don't have the same keys.")

        # tip nodes should have one branch, internal nodes three
        for key in nodeBranchD:
            if nodeNameD[key] == 'anc':
                if len(nodeBranchD[key]) != 3:
                    raise ValueError("Internal nodes in nodeBranchD should be attached to 3 branches.")
            else:
                if len(nodeBranchD[key]) != 1:
                    raise ValueError("Tips in nodeBranchD should be attached to 1 branch.")

        self.nodeBranchD = nodeBranchD
        self.nodeNameD = nodeNameD

        # every branch should have two nodes
        if not all(len(branchNodeD[branch])==2 for branch in branchNodeD):
            raise ValueError("Error in branchNodeD, each branch should be attached to two nodes.")
        else:
            self.branchNodeD = branchNodeD


        # branch lengths
        if branchLenD == None:
            branchLenD = {}
            for key in branchNodeD:
                branchLenD[key] = None
            
                
        self.branchLenD = branchLenD

        # now make nodeList and branchList attributes
        branchList = list(branchLenD.keys())
        branchList.sort()
        self.branches = tuple(branchList)

        nodeList = list(nodeBranchD.keys())
        nodeList.sort()
        self.nodes = tuple(nodeList)

        
            
    def __repr__(self):
        rowsL=[]
        rowsL.append(['Sp.','Nnum)', 'branch','Length','(Nnum','Sp.'])
        rowsL.append(['---','-----', '------','------','-----','---'])
        for branch in self.branchLenD:
            node1,node2 = self.branchNodeD[branch]

            # get strings for nodes
            
            nm1Str = self.getName(node1) if self.getName(node1) != 'anc' else ''
            nm2Str = self.getName(node2) if self.getName(node2) != 'anc' else ''

            # get branch len str
            if self.branchLenD[branch] == None:
                brLenStr = str(None)
            else:
                brLenStr = format(self.branchLenD[branch],".4f")
            

            rowsL.append([nm1Str,str(node1)+")",str(branch),brLenStr,"("+str(node2),nm2Str])
        
        return tableString(rowsL)

    def getNodes(self,branch):
        '''Given a branch number, return the nodes it's attached to.'''
        return self.branchNodeD[branch]

    def getBranches(self,node):
        '''Given a node number, return the branches it's attached to.'''
        return self.nodeBranchD[node]

    def getName(self,node):
        '''Given a node number, return the species name it corresponds to.'''
        return self.nodeNameD[node]

    def getBrlen(self,branch):
        '''Given a branch number, return its length.'''
        return self.branchLenD[branch]
    
    def setBrlen(self,branch,branchLen):
        '''Set the length of branch to be branchLen.'''
        self.branchLenD[branch]=branchLen
        
    def printNodeNameSeq(self,dataD):
        '''Given a dict of sequences, print node number, name, seq for each node.'''
        print("Node number, name and sequence for tips:")
        rowL=[]
        for node in self.nodes:
            if self.getName(node) != 'anc':
                rowL.append([str(node),self.getName(node),dataD[self.getName(node)]])

        print(tableString(rowL))

def tableString(L,indent=2):
    '''Given tabular data in a list of lists (where sublists are rows)
format nicely so columns line up, and return as string. Indent is an
optional number of blank spaces to put in front of each row.
    '''
    # get max width for each column
    outStr=''
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
        outStr += " "*indent + " ".join(row).rstrip()+'\n'

    return outStr
