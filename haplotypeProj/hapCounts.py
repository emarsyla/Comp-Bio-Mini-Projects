def getHapData(filename):
    f = open(filename)
    num_tuples = 2 * (len(f.readline().split()) - 5)
    list_of_tuples = [() for _ in range(num_tuples)]
    iteration = 0
    while True:
        # Read each line
        s=f.readline()
        if s=="":
            break
        L=s.split()
        ID = L[2]
        L = L[5:]

        # Make a list of the nucleotides
        list = []
        for string in L:
            for char in string:
                if char != ",":
                    list.append(char)
        # Put into the appropriate tuple
        for enum, elem in enumerate(list):
            list_of_tuples[enum] = list_of_tuples[enum] + (elem,)

        if ID == "rs4988235":
            lactase_id_pos = iteration
        
        iteration += 1
    return list_of_tuples, lactase_id_pos



def hapCounts(hapDataL):
    unique = list(set(hapDataL))
    dict = {}

    # make dictionary
    for elem in unique:
        dict[elem] = len([tup for tup in hapDataL if tup == elem])

    return dict

def answers(filename1, filename2):
    data1, idpos1 = getHapData(filename1)
    numHaps1 = len(data1)
    data2, idpos2 = getHapData(filename2)
    numHaps2 = len(data2)
    print("There were", numHaps1, "haplotypes sequenced int he finninsh population, and", numHaps2, "Haplotypes sequenced in the yoruba population")
    dict1 = hapCounts(data1)
    dict2 = hapCounts(data2)

    # Find top 3 counts in dictionary
    counts1 = [0, 0, 0]
    haps1 = ["", "", ""]
    for key, value in dict1.items():
        if value > counts1[0]:
            counts1[2] = counts1[1]
            counts1[1] = counts1[0]
            counts1[0] = value
            haps1[2] = haps1[1]
            haps1[1] = haps1[0]
            haps1[0] = key[idpos1]
        elif value > counts1[1]:
            counts1[2] = counts1[1]
            counts1[1] = value
            haps1[2] = haps1[1]
            haps1[1] = key[idpos1]
        elif value > counts1[2]:
            counts1[2] = value
            haps1[2] = key[idpos1]

    counts2 = [0, 0, 0]
    haps2 = ["", "", ""]
    for key, value in dict2.items():
        if value > counts2[0]:
            counts2[2] = counts2[1]
            counts2[1] = counts2[0]
            counts2[0] = value
            haps2[2] = haps2[1]
            haps2[1] = haps2[0]
            haps2[0] = key[idpos2]
        elif value > counts2[1]:
            counts2[2] = counts2[1]
            counts2[1] = value
            haps2[2] = haps2[1]
            haps2[1] = key[idpos2]
        elif value > counts2[2]:
            counts2[2] = value
            haps2[2] = key[idpos2]
        
    print("The counts for the top three most common haplotypes in the finnish population are", counts1, "While the counts of the top most common haplotypes in the yoruba population are", counts2)
    print("The lactase persitance allele for finnish are", haps1, "and", haps2, "for yoruba")