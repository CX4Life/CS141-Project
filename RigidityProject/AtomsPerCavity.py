proteinDB = open("1HE4.pdb.txt", 'r')


def is_atom(line):
    if line[0] == 'A' and line[1] == 'T' and line[2] == 'O':
        return True
    else:
        return False

# This next section is the portion of the code where we count atoms per
# residue, and append a list with that number of atoms. The means of determining
# if a particular line is part of a NEW residue or the same residue it just checked
# seem pretty clunky to me, because for every line it iterate 4 times just to
# check the residue number. Nothing crazy, but it could be better.

atomCounter = 0
residueCounter = 0
firstOrLastResidue = True
atomsPerResidue = []

for line in proteinDB:
    residueNumber = ''
    if is_atom(line):
        for i in range(22, 26):
            residueNumber += line[i]
        residueNumber = int(residueNumber)
        if firstOrLastResidue:
            residueCounter = residueNumber
            firstOrLastResidue = False
        if residueNumber == residueCounter:
            atomCounter += 1
        else:
            atomsPerResidue.append(atomCounter)
            residueCounter += 1
            atomCounter = 1
    else:
        if not firstOrLastResidue and len(atomsPerResidue) >= 1:
            atomsPerResidue.append(atomCounter)
            firstOrLastResidue = True

proteinDB.close()

# print(atomsPerResidue)
# print(len(atomsPerResidue))
# assert False

# atomsPerResidue now holds a list of integers, the index of which
# represents residue number - 1, and the value of which is the number
# of atoms in the residue.

# This code is directly lifted from the first... an example of
# something that could be generalized for different metrics, for sure.

cavityInfo = open('1he4A.SurfReport', 'r')
duesByCavity = []
for line in cavityInfo.readlines():
    residues = line[line.find('{')+1:line.find('}')-1]
    listResidues = [int(i) for i in residues.split(',')]
    duesByCavity.append(listResidues)

cavityInfo.close()

atomCountPerCavity = []

for cavity in duesByCavity:
    atomsInCavity = 0
    for residue in cavity:
        atomsInCavity += atomsPerResidue[residue - 1]
    atomCountPerCavity.append(atomsInCavity)

print(atomCountPerCavity)

## print(atomCountPerCavity)




