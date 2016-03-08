from tabulate import tabulate

# --------------------------------------------------------------
# These functions are used with the "*.pdb.text" file.
# --------------------------------------------------------------


def is_atom(line):
    if line[0] == 'A' and line[1] == 'T' and line[2] == 'O':
        return True
    else:
        return False


def is_alpha_carbon(line):
    if line[13] == "C" and line[14] == "A":
        return True
    else:
        return False


def get_atom_number(line):
    atom_number = ''
    start_check = 10
    while line[start_check] != ' ':
        atom_number = line[start_check] + atom_number
        start_check -= 1
    atom_number = int(atom_number)
    return atom_number


def get_alpha_carbon_atom_numbers(proteinDBFilename):
    assert type(proteinDBFilename) is str
    alphaCarbonAtomNumbers = []
    proteinDB = open(proteinDBFilename, 'r')
    for line in proteinDB:
        if is_atom(line) and is_alpha_carbon(line):
            alphaCarbonAtomNumbers.append(get_atom_number(line))

    alphaCarbonAtomNumbers = tuple(alphaCarbonAtomNumbers)

    return alphaCarbonAtomNumbers


def get_residue_type(line):
    residue_type = ''
    for i in range(17,20):
        residue_type += line[i]
    return residue_type


def return_each_residue_type(filename):
    residueTypes = []
    proteinDB = open(filename, 'r')
    for line in proteinDB:
        if is_atom(line) and is_alpha_carbon(line):
            residueTypes.append(get_residue_type(line))
    proteinDB.close()
    return residueTypes

# --------------------------------------------------------------
# The next two functions are used with the "*.SurfReport" file.
# --------------------------------------------------------------

# This simply parses a "SurfReport" file and returns the size of each cavity.
# Output is a list whose index is cavity identity and value is cavity size.


def get_cavity_size(filename):
    assert type(filename) is str
    openFile = open(filename, 'r')
    resultCavitySize = []
    for line in openFile:
        rawSize = line[line.find('(')+1:line.find(')')-1]
        size = rawSize.split(':')
        size = float(size[1])
        resultCavitySize.append(size)
    openFile.close()
    return resultCavitySize

# This function returns a list of lists whose elements are residue numbers. For instance,
# when calling list[i[n]], "i" represents the Cavity number, and "n" represents a residue number
# in cavity "i"


def get_residues_by_cavity(filename):
    assert type(filename) is str
    cavityInfo = open(filename, 'r')
    duesByCavity = []
    for line in cavityInfo.readlines():
        residues = line[line.find('{')+1:line.find('}')-1]
        listResidues = [int(i) for i in residues.split(',')]
        duesByCavity.append(listResidues)
    cavityInfo.close()
    return duesByCavity


# --------------------------------------------------------------
# The next two functions are used with the "*.xml" file containing
# rigid structure data.
# --------------------------------------------------------------


def get_atoms_in_rigids(rigidFilename):

    from xml.etree import ElementTree

    rigidData = ElementTree.parse(rigidFilename)
    root = rigidData.getroot()
    atomsInRigids = []
    for child in root[0]:
        for pointSet in child:
            # print('This is a new rigid structure')
            atomsInPointSet = []
            for point in pointSet:
                atomsInPointSet.append(int(point.attrib['id']))
            atomsInRigids.append(atomsInPointSet)

    return atomsInRigids


def get_alpha_carbons_in_rigids(rawAtomsInRigids, alphaCarbonAtomNumbers):
    for rigid in rawAtomsInRigids:
        for atom in rigid[:]:
            if atom not in alphaCarbonAtomNumbers:
                rigid.remove(atom)
    return rawAtomsInRigids

def get_all_aminos(listByType):
    aminos = []
    while len(aminos) != 20:
        for list in listByType:
            for elem in list:
                if elem not in aminos:
                    aminos.append(elem)
    return aminos


proteinName = input('What is the alphanumeric designation of this protein? ')

pdbFile = proteinName + '.pdb.txt'
surfFile = str.lower(proteinName) + 'A.SurfReport'
xmlFile = str.lower(proteinName) + '_Processed_postPG_BBH.xml'


alphaNums = get_alpha_carbon_atom_numbers(pdbFile)
duesByCavity = get_residues_by_cavity(surfFile)
AtomsInRigids = get_atoms_in_rigids(xmlFile)

rigidSize = []
for rigid in AtomsInRigids:
    rigidSize.append(len(rigid))

rawAtomsInRigids = get_alpha_carbons_in_rigids(AtomsInRigids, alphaNums)
alphaCarbonAtomNumbers = get_alpha_carbon_atom_numbers(pdbFile)
cavSize = get_cavity_size(surfFile)
ResidueTypes = return_each_residue_type(pdbFile)

rigidsByCavity = []
ResidueTypeByCavity = []
ResiduesPerCavity = []
inCavity = False


for cavity in duesByCavity:
    rigidsThisCavity = 0
    for individualRigid in rawAtomsInRigids:
        stopLoop = int(len(cavity))
        residueCounter = 0
        while not inCavity and residueCounter != stopLoop:
            for residueNumber in cavity:
                if alphaCarbonAtomNumbers[residueNumber] in individualRigid:
                    inCavity = True
                    rigidsThisCavity += 1
                    break
                residueCounter += 1
        inCavity = False
    residuesThisCavity = []
    for residue in cavity:
        residuesThisCavity.append(ResidueTypes[residue-1])
    ResidueTypeByCavity.append(residuesThisCavity)
    ResiduesPerCavity.append(len(residuesThisCavity))
    rigidsByCavity.append(rigidsThisCavity)

outputFile = proteinName + "Output.txt"
outputFile = open(outputFile, 'w')

cavNums = []
for i in range(0, len(cavSize)):
    cavNums.append(i)

# --------------------------------------------------------------
# For the sake of not printing a table with a ton of columns,
# I completely arbitrarily chose to count the number
# of Alanines, Glycines and Valines in each cavity. However, the
# 'ResidueTypesByCavity' object contains counts for every
# residue type, so the code could very easily be modified to
# include print every residue type for every cavity.
# --------------------------------------------------------------

table = [['Cav #', 'Cav Size', '# Residues', "ALA's", "GLY's", "VAL's", '# Rigids', 'Total Size of Rigids']]
for i in range(0, len(cavSize)):
    table.append([cavNums[i], cavSize[i], ResiduesPerCavity[i], ResidueTypeByCavity[i].count('ALA'), ResidueTypeByCavity[i].count('GLY'), ResidueTypeByCavity[i].count('VAL'), rigidsByCavity[i], rigidSize[i]])

outputFile.write(tabulate(table))
outputFile.close()

