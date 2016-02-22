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
    return resultCavitySize

