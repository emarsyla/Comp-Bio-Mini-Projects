import numpy, random

def loadNeandReadD(neandFileName):
    """Loads the Neanderthal data file and then returns a dictionary of
    the alleles present for each chromosome and position."""
    f = open(neandFileName)
    dict = {}
    while True:
        s = f.readline()
        if s=="":
            break
        L = s.split()
        chromosome = L[0]
        position = int(L[1])
        nucleotides = L[2]
        tuple = ()
        for char in nucleotides:
            if char != ",":
                tuple += (char, )
        # Store alleles for each chromosome and position
        dict[chromosome, position] = tuple
    return dict

def loadModernHumanData(filename, neandReadD):
    """Loads the modern human data and adds information to a list of tuples """
    list_of_tuples = []
    f = open(filename)
    s = f.readline()
    while True:
        allelesL = []
        s = f.readline()
        if s == "":
            break
        L = s.split()
        Chr = L[0]
        pos = int(L[1])
        chimp = L[4]
        nucs = L[5:]
        # Only add information if it is included in the neanderthal dictionary
        if (Chr, pos) in neandReadD:
            for string in nucs:
                for char in string:
                    if char != ",":
                        allelesL.append(char)
            # Store information for each chromosome and position
            tuple = (Chr,pos,chimp,allelesL)
            list_of_tuples.append(tuple)
    return list_of_tuples

def derAlleleCount(neandReadD,h1SiteDataL,h2SiteDataL):
    """Loads the neanderthal dictionary and the list of tuples for two human populations. 
    This function compares two modern human populations"""
    human1count = 0
    human2count = 0
    #Iterates through each element in the list
    for i in range(len(h1SiteDataL)):
        human1 = h1SiteDataL[i]
        human2 = h2SiteDataL[i]
        chr = human1[0]
        pos = human1[1]
        chimp = human1[2]
        neand = neandReadD[chr, pos]
        # Choose random alleles
        hum1allele = random.choice(human1[3])
        hum2allele = random.choice(human2[3])
        neandRead = random.choice(neand)
        alleles = [chimp, neandRead, hum1allele, hum2allele]
        # Criteria 1) Are there a total of two types of alleles present in the 4 samples
        if len(set(alleles)) == 2:
            # Criteria 2) Do the two human alleles differ
            if alleles[2] != alleles[3]:
                # Criteria 3) Does neandertha; have the derived allele
                derived = [allele for allele in set(alleles) if allele != chimp][0]
                if neandRead == derived:
                    # Add to counts if the humans have the derived allele
                    if hum1allele == derived:
                        human1count += 1
                    if hum2allele == derived:
                        human2count += 1
    return human1count, human2count

def wrapper(finnishFilename, yorubaFilename, chineseFilename, neandFilename):
    # Wrapper function runs all 3 comparisons
    neandReadD = loadNeandReadD(neandFilename)
    finnish = loadModernHumanData(finnishFilename,neandReadD)
    yoruba = loadModernHumanData(yorubaFilename,neandReadD)
    chinese = loadModernHumanData(chineseFilename,neandReadD)
    # Finnish - Yoruba comparison
    difList = []
    for _ in range(100):
        finCount,yorCount=derAlleleCount(neandReadD,finnish,yoruba)
        difList.append(finCount - yorCount)
    mean = numpy.mean(difList)
    stddev = numpy.std(difList)
    FinYorComp = (mean, stddev)

    # Finnish - Chinese comparison
    difList = []
    for _ in range(100):
        finCount,chiCount=derAlleleCount(neandReadD,finnish,chinese)
        difList.append(finCount - chiCount)
    mean = numpy.mean(difList)
    stddev = numpy.std(difList)
    FinChiComp = (mean, stddev)
    
    # Chinese - Yoruba comparison
    difList = []
    for _ in range(100):
        chiCount,yorCount=derAlleleCount(neandReadD,chinese,yoruba)
        difList.append(chiCount - yorCount)
    mean = numpy.mean(difList)
    stddev = numpy.std(difList)
    ChiYorComp = (mean, stddev)

    print("The mean and standard deviations for the Finnish - Yoruba comparison are", FinYorComp)
    print("The mean and standard deviations for the Finnish - Chinese comparison are", FinChiComp)
    print("The mean and standard deviations for the Chinese - Yoruba comparison are", ChiYorComp)
    