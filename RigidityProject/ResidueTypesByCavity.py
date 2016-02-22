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

def get_residue_type(line):
    residue_type = ''
    for i in range(17,20):
        residue_type += line[i]
    return residue_type

# This section creates our list of residue types. The index of the list elements is the
# residue number minus 1.

ResidueTypes = []

for line in proteinDB:
    if is_atom(line) and is_alpha_carbon(line):
        ResidueTypes.append(get_residue_type(line))

proteinDB.close()

# Once again, we lift the residue # in cavity section from the "RigidsPerCluster" project.

cavityInfo = open('1he4A.SurfReport', 'r')
duesByCavity = []
for line in cavityInfo.readlines():
    residues = line[line.find('{')+1:line.find('}')-1]
    listResidues = [int(i) for i in residues.split(',')]
    duesByCavity.append(listResidues)

cavityInfo.close()

# Boom. Let the cross referencing begin.

ResidueTypeByCavity = []


for cavity in duesByCavity:
    residuesThisCavity = []
    for residue in cavity:
        residuesThisCavity.append(ResidueTypes[residue-1])
    ResidueTypeByCavity.append(residuesThisCavity)

# Cool. Now we can do fun stuff with this list of lists.

for cavity in ResidueTypeByCavity:
    print(cavity.count('MET'))
    # We can even check to see if I screwed up...
    assert len(cavity) == len(duesByCavity[ResidueTypeByCavity.index(cavity)])

# Or just print the thing.

print(ResidueTypeByCavity)

### CAV number - size (angstroms) residue count, residue TYPES, # of rigid clusters, total size of rigid clusters