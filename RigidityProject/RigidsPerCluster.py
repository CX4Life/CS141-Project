# Here we open the file to then save the alpha-carbom
# Atom's numbers into a list so that we can then simplify
# our rigid cluster data by EXCLUDING anything that is no
# an alpha-carbon atom in the rigid cluster data.

proteinDB = open("1HE4.pdb.txt", 'r')


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

alphaCarbonAtomNumbers = []

for line in proteinDB:
    if is_atom(line):
        if is_alpha_carbon(line):
            alphaCarbonAtomNumbers.append(get_atom_number(line))

proteinDB.close()
# At this point, we have a list containing all the alpha carbon atom numbers. This is good.
# We can now use this list to simplify our rigid cluster data. The result that I want is a
# dictionary that has key values of rigid cluster ID's which are paired to lists containing
# all alpha carbons in those rigid clusters.

rigidFilename = '1he4_Processed_postPG_BBH.xml'

def get_atoms_in_rigids(rigidFilename):

    from xml.etree import ElementTree

    rigidData = ElementTree.parse(rigidFilename)
    root = rigidData.getroot()
    atomsInRigids = []
    i = 0
    for child in root[0]:
        for pointSet in child:
            # print('This is a new rigid structure')
            atomsInPointSet = []
            for point in pointSet:
                atomsInPointSet.append(int(point.attrib['id']))
            atomsInRigids.append(atomsInPointSet)
    
    return atomsInRigids

rawAtomsInRigids = get_atoms_in_rigids(rigidFilename)


# Our next step is to remove from that big list of all atom numbers in the the rigid structures any atom
# that is not an alpha carbon. This should hopefully speed up our final evaluation, but honestly, I don't
# know if it really matters. This is worth checking in with Filip from a Big O standpoint.

for rigid in rawAtomsInRigids:
    for atom in rigid[:]:
        if atom not in alphaCarbonAtomNumbers:
            rigid.remove(atom)

# Now we have a list organized by index:value where index = residue number and value = alpha carbon atom number
# and a list organized by index:value where index is rigid cluster index and value is a list of alpha carbons IN
# that rigid cluster.

cavityInfo = open('1he4A.SurfReport', 'r')
duesByCavity = []
for line in cavityInfo.readlines():
    residues = line[line.find('{')+1:line.find('}')-1]
    listResidues = [int(i) for i in residues.split(',')]
    duesByCavity.append(listResidues)

cavityInfo.close()
# Now we have all three of our important data sets created. duesByCavity gives lists of residues in each cavity,
# rawAtomsInRigids has been converted to only hold alpha carbon atom numbers, and alphaCarbonAtom... holds carbon atom
# numbers who's index is their residue number. We've got everything we need, we just need to evaluate cavity by cavity.

rigidsByCavity = []
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
    rigidsByCavity.append(rigidsThisCavity)

# The code below was my original attempt. It had a major logical error, in that
# it would check to see if an individual residue was contained in any rigid, then
# iterate a counter for that residue. The error was, if two residues were contained
# by the same rigid structure, it would count that rigid twice in the result. The
# above code solves that problem by checking each rigid structure to see if contains
# residues that are in the cavity, this way each rigid structure is only counted once
# per cavity.

    # for residueNumber in cavity:
    #     #The residueNuber element is an integer which represents one residue
    #     for individualRigid in rawAtomsInRigids:
    #         #IndividualRigid is a list of of integers
    #         if alphaCarbonAtomNumbers[residueNumber] in individualRigid:
    #             rigidsThisCavity += 1
    # rigidsInCavity.append(rigidsThisCavity)
    # rigidsThisCavity = 0

print(rigidsByCavity)

# print(len(rigidsInCavity) - len(duesByCavity))
# print(len(rawAtomsInRigids))