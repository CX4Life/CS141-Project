# For this guy, we can just steal the residue type generator from
# ResidueTypesByCavity...
# It's becoming clear at this point that I'm duplicating a TON of code
# and that really the only thing that differs substantially between each of these
# little evaluators is how the various data structures are combined or cross-referenced.

# If I knew more about how to effectively use classes, then this would be a much more flexible,
# extensible single program, rather than an amalgam of little, highly duplicitous modules.

# These first functions and the opening of the PDB are lifted from the "RigidsPerCluster" project.

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

# This next function is simply a modification of the "get_atom_number" function from
# the "RigidsPerCluster" project

def get_atom_number(line):
    atom_number = ''
    for i in range(7, 11):
        atom_number += line[i]
    atom_number = int(atom_number)
    return atom_number

def get_residue_type(line):
    residue_type = ''
    for i in range(17, 20):
        residue_type += line[i]
    return residue_type

# This section creates our list of residue types. The index of the list elements is the
# residue number minus 1.

ResidueTypes = []
atomsInThisResidue = []
atomsByResidue = []
firstOrLastAtom = True
dueChecker = ''

for line in proteinDB:
    # One thing that would be nice is to only have to iterate through
    # All the lines in proteinDB once, and get my list of atom numbers
    # by residue AND my list of residue types at the same time...

    if is_atom(line):
        if is_alpha_carbon(line):
            ResidueTypes.append(get_residue_type(line))
        if firstOrLastAtom:
            dueChecker = line[25]
            firstOrLastAtom = False
        if line[25] == dueChecker:
            atomsInThisResidue.append(get_atom_number(line))
        else:
            atomsByResidue.append(atomsInThisResidue)
            atomsInThisResidue = []
            dueChecker = line[25]
            atomsInThisResidue.append(get_atom_number(line))
    else:
        if not firstOrLastAtom and len(atomsByResidue) >= 1:
            firstOrLastAtom = True
            atomsByResidue.append(atomsInThisResidue)

# NICE. That felt good.

proteinDB.close()

# This time, we aren't concerned with cavities, just rigid clusters. We need to steal a
# big chunk from RigidsPerCluster to get a list of residues in the rigid clusters... by atom number.
# This could be a huge, slow, pain, especially if we forego the original condition that we're only
# concerned about the alpha carbon atoms of each residue. It can be done, though.

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

# A nice thing would be a list of residues in each rigid BY residue number. Let's make that happen...


def get_residues_in_rigids(rawAtomsInRigids, atomsByResidue):

    residuesInRigids = []
    for individualRigid in rawAtomsInRigids:
        residuesInThisRigid = []
        for residue in atomsByResidue:
            for atomNumber in individualRigid:
                if atomNumber > residue[-1]:
                    break
                if atomNumber in residue and (atomsByResidue.index(residue) not in residuesInThisRigid):
                    residuesInThisRigid.append(atomsByResidue.index(residue))
        residuesInRigids.append(residuesInThisRigid)
    return residuesInRigids

# Cool. Now just changing those integer values to residue types should be simple.

numericalResiduesInRigids = get_residues_in_rigids(rawAtomsInRigids, atomsByResidue)

bugCheck1 = len(numericalResiduesInRigids)

residueTypesByRigid = []

for rigid in numericalResiduesInRigids:
    residueTypesThisRigid = []
    for residue in rigid:
        residueTypesThisRigid.append(ResidueTypes[residue])
    residueTypesByRigid.append(residueTypesThisRigid)

bugCheck2 = len(residueTypesByRigid)

assert bugCheck1 == bugCheck2
assert len(rawAtomsInRigids) == len(residueTypesByRigid)

print(residueTypesByRigid)

