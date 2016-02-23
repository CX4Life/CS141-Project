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


# --------------------------------------------------------------
# The next two functions are used with the "*.SurfReport" file.
# --------------------------------------------------------------

# This simply parses a "SurfReport" file and returns the size of each cavity.
# Output is a list whose index is cavity identity and value is cavity size.

def getCavitySize(filename):
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

def getResiduesByCavity(filename):
    assert type(filename) is str
    cavityInfo = open(filename, 'r')
    duesByCavity = []
    for line in cavityInfo.readlines():
        residues = line[line.find('{')+1:line.find('}')-1]
        listResidues = [int(i) for i in residues.split(',')]
        duesByCavity.append(listResidues)
    cavityInfo.close()
    return duesByCavity
