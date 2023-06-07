import fempy
from .analysis import Pair
from .io import *
from .format import *

def Chunkyfy(array, nElemPerChunk=200):
    if len(array) == 0:
        return [[]]
    nChunks = (len(array) -1) // nElemPerChunk +1
    return [array[iChunk*nElemPerChunk:nElemPerChunk * (iChunk +1)] for iChunk in range(nChunks)]


def MergeInChunks(mergeList, oFile, nElemPerChunk=200):
    # perform the mergin in 2 steps. Temporary merge files are not targets of the task
    tmpFiles = []
    for iMerge, mList in enumerate(Chunkyfy(mergeList, nElemPerChunk)):
        print(mList)
        tmpFiles.append(os.path.splitext(oFile)[0] + f'_{iMerge}.root ')
        command = 'hadd -f ' + tmpFiles[-1] + ' '.join(mList)
        print(command)
        os.system(command)

    # Final merging    
    command = 'hadd -f ' + oFile + ' ' + ' '.join(tmpFiles)
    os.system(command)


def DivideCanvas(canvas, n):
    if n <= 1:
        return canvas
    elif n<=3:
        canvas.Divide(n, 1)
        return canvas
    elif n == 4:
        canvas.Divide(2, 2)
        return canvas
    elif n <= 8:
        canvas.Divide((n+1)//2, 2)
        return canvas
    elif n <= 12:
        canvas.Divide((n+1)//3, 3)
        return canvas
    elif n <= 20:
        canvas.Divide((n+1)//5, 4)
        return canvas
    else:
        fempy.error("too many pads")
        return None


def GetNPanels(n):
    if n<=3:
        return (n, 1)
    elif n == 4:
        return (2, 2)
    elif n <= 8:
        return ((n+1)//2, 2)
    elif n <= 12:
        return ((n+1)//3, 3)
    elif n <= 20:
        return ((n+1)//5, 4)
    else:
        fempy.error("too many pads")
    