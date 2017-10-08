"""Accept three files concerning one protein as input, and output text file
containing information about that specific protein.

The *.pdb.text file contains a list of atoms and their classification within
each 'residue' in the protein, only one of which is that alpha carbon atom.
This is relevant as the alpha carbon atom numbers are used to compare features
between the rigidity and cavity data.


"""
import sys

# --------------------------------------------------------------
# These functions are used with the "*.pdb.text" file.
# --------------------------------------------------------------


def is_atom(_line):
    """Return true if line in file concerns an atom, else false."""
    return _line[:4] == 'ATOM'

def is_alpha_carbon(_line):
    """Return true if line in file concerns an alpha carbon, else false."""
    return _line[13:15] == "CA"

def get_atom_number(_line):
    """Return the atom number given a line of text."""
    atom_number = ''
    start_check = 10
    while _line[start_check] != ' ':
        atom_number = _line[start_check] + atom_number
        start_check -= 1
    return int(atom_number)


def get_alpha_carbon_atom_numbers(_pdb_filename):
    """Return the all the alpha carbons by number from the
    protein database text file."""
    assert isinstance(_pdb_filename, str)
    ac_numbers = []
    all_res = []
    _pdb_file = open(_pdb_filename, 'r')
    for _line in _pdb_file:
        if is_atom(_line) and is_alpha_carbon(_line):
            atom_number = int(_line[23:27])
            if atom_number not in all_res:
                ac_numbers.append(get_atom_number(_line))
                all_res.append(atom_number)

    _pdb_file.close()
    return [ac_numbers, all_res]


def get_residue_type(_line):
    """Return the residue type (what amino acid) for a given line of atom data."""
    return _line[17:21]


def return_each_residue_type(_pdb_filename):
    """Return of list of all residues in a protein."""
    residue_types = []
    _pdb_file = open(_pdb_filename, 'r')
    for _line in _pdb_file:
        if is_atom(_line) and is_alpha_carbon(_line):
            residue_types.append(get_residue_type(_line))
    _pdb_file.close()
    return residue_types

# --------------------------------------------------------------
# The next two functions are used with the "*.SurfReport" file.
# --------------------------------------------------------------

def get_cavity_size(_surf_filename):
    """Return the size of each cavity in a SurfReport file."""
    assert isinstance(_surf_filename, str)
    __surf_file = open(_surf_filename, 'r')
    cavity_sizes = []
    for _line in __surf_file:
        _raw_size = _line[_line.find('(')+1:_line.find(')')-1]
        _size = float(_raw_size.split(':')[1])
        cavity_sizes.append(_size)
    __surf_file.close()
    return cavity_sizes


def get_residues_by_cavity(_surf_filename):
    """Return a list of lists whose elements are residue numbers.

    For instance, when calling list[i[n]], "i" represents the Cavity number,
    and "n" represents a residue number in cavity "i"
    """
    assert isinstance(_surf_filename, str)
    __surf_file = open(_surf_filename, 'r')
    residues_in_each_cavity = []
    for _line in __surf_file.readlines():
        _residues = _line[_line.find('{')+1:_line.find('}')-1]
        _list_of_residues = [int(i) for i in _residues.split(',')]
        residues_in_each_cavity.append(_list_of_residues)
    __surf_file.close()
    return residues_in_each_cavity


# --------------------------------------------------------------
# The next two functions are used with the "*.xml" file containing
# rigid structure data.
# --------------------------------------------------------------


def get_atoms_in_rigids(rigid_filename):
    """Return the atoms by number in each rigid structure as a
    list of lists. The outer list is each rigid structure, the
    inner list is each atom number."""
    from xml.etree import ElementTree

    _rigid_data = ElementTree.parse(rigid_filename)
    _root = _rigid_data.getroot()
    atoms_in_rigid = []
    for _child in _root[0]:
        for _pointset in _child:
            _atoms_in_pointset = []
            for _point in _pointset:
                _atoms_in_pointset.append(int(_point.attrib['id']))
            atoms_in_rigid.append(_atoms_in_pointset)

    return atoms_in_rigid


def get_alpha_carbons_in_rigids(rigids_to_mutate, ac_numbers):
    """Strip every atom which is not an alpha carbon from each rigid
    structure, and return only those atoms which are alpha carbons in
    and in a rigid structure."""
    for _rigid in rigids_to_mutate:
        for atom in _rigid[:]:
            if atom not in ac_numbers:
                _rigid.remove(atom)
    return rigids_to_mutate


def get_all_aminos(list_by_type):
    """Return all amino acids present in this protein."""
    aminos = []
    while len(aminos) != 20:
        for list_ in list_by_type:
            for _elem in list_:
                if _elem not in aminos:
                    aminos.append(_elem)
    return aminos

PDB_FILE = sys.argv[1]
SURF_FILE = sys.argv[2]
XML_FILE = sys.argv[3]
PROTEIN_NAME = PDB_FILE[6:10]

RESIDUES_BY_CAVITY = get_residues_by_cavity(SURF_FILE)
ATOMS_BY_RIGID = get_atoms_in_rigids(XML_FILE)

RIGID_SIZE = []
for rigid in ATOMS_BY_RIGID:
    RIGID_SIZE.append(len(rigid))

AC_ATOM_NUMBERS = get_alpha_carbon_atom_numbers(PDB_FILE)[0]
RAW_ATOMS_IN_RIGIDS = get_alpha_carbons_in_rigids(ATOMS_BY_RIGID, AC_ATOM_NUMBERS)
CAVITY_SIZE = get_cavity_size(SURF_FILE)
RESIDUE_TYPES = return_each_residue_type(PDB_FILE)

RIGIDS_BY_CAVITY = []
RESIDUES_TYPES_BY_CAVITY = []
RESIDUES_PER_CAVITY = []
RIGID_SIZE_BY_CAVITY = []
IN_CAVITY = False


for cavity in RESIDUES_BY_CAVITY:
    rigidsThisCavity = 0
    rigidSizeThisCavity = 0
    for individualRigid in RAW_ATOMS_IN_RIGIDS:
        stopLoop = int(len(cavity))
        residueCounter = 0
        while not IN_CAVITY and residueCounter != stopLoop:
            for residueNumber in cavity:
                if AC_ATOM_NUMBERS[cavity.index(residueNumber)] in individualRigid:
                    rigidSizeThisCavity += RIGID_SIZE[RAW_ATOMS_IN_RIGIDS.index(individualRigid)]
                    in_cavity = True
                    rigidsThisCavity += 1
                    break
                residueCounter += 1
        in_cavity = False
    residuesThisCavity = []
    for residue in cavity:
        residuesThisCavity.append(RESIDUE_TYPES[cavity.index(residue)])
    RESIDUES_TYPES_BY_CAVITY.append(residuesThisCavity)
    RESIDUES_PER_CAVITY.append(len(residuesThisCavity))
    RIGIDS_BY_CAVITY.append(rigidsThisCavity)
    RIGID_SIZE_BY_CAVITY.append(rigidSizeThisCavity)

OUTPUT_FILENAME = PROTEIN_NAME + "Output.txt"
OUTPUT_FILE = open(OUTPUT_FILENAME, 'w')

# --------------------------------------------------------------
# For the sake of not printing a table with a ton of columns,
# I completely arbitrarily chose to count the number
# of Alanines, Glycines and Valines in each cavity. However, the
# 'ResidueTypesByCavity' object contains counts for every
# residue type, so the code could very easily be modified to
# include print every residue type for every cavity.
# --------------------------------------------------------------

TABLE = []
OUTPUT_FILE.write("Output is as such:\nC# : Cavity Number\n")
OUTPUT_FILE.write("CS : Cavity Size\nRsC : Residues in this cavity\n")
OUTPUT_FILE.write("RgC : Rigid clusters containing >= 1 residue in this cavity\n")
OUTPUT_FILE.write("SumRgSiz : Total atoms in all Rigid clusters interacting with this cavity\n\n")

OUTPUT_FILE.write("C#      CS      RsC     RgC     SumRgSiz\n")
for i in range(0, len(CAVITY_SIZE)):
    TABLE.append(
        [i, CAVITY_SIZE[i],
         RESIDUES_PER_CAVITY[i],
         RIGIDS_BY_CAVITY[i],
         RIGID_SIZE_BY_CAVITY[i]]
        )
for line in TABLE:
    outString = ""
    for elem in line:
        outString += str(elem) + "\t"
    OUTPUT_FILE.write(outString + "\n")
OUTPUT_FILE.close()
