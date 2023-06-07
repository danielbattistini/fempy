def Chunkyfy(array, nElemPerChunk=50):
    if len(array) == 0:
        return [[]]
    nChunks = (len(array) -1) // nElemPerChunk +1
    return [array[iChunk*nElemPerChunk:nElemPerChunk * (iChunk +1)] for iChunk in range(nChunks)]



print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 1))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 2))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 3))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 4))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 5))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 6))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 7))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 8))
print(Chunkyfy([1, 2, 3, 4, 5, 6, 7], 9))

