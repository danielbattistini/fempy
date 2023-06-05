import os
from ROOT import TFile, TH1, TH2, TH3, TDirectoryFile, TList

def GetObjectFromFile(inFile, pathToObj):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TFile
     - TList

    Parameters
    -----------
    inFile: TFile or name of the input file
    pathToObj: path of the object inside the root file

    Returns:
    -----------
    outObj: target root object
    '''

    pathToObj = os.path.normpath(pathToObj)
    pathElements = pathToObj.split(os.sep)

    outObj = inFile.Get(pathElements.pop(0))
    if isinstance(outObj, TH1) or isinstance(outObj, TH2) or isinstance(outObj, TH3):
        outObj.SetDirectory(0)
    for iContainer, containerName in enumerate(pathElements):
        if isinstance(outObj, TFile) or isinstance(outObj, TDirectoryFile):
            outObj = outObj.Get(containerName)
        elif isinstance(outObj, TList):
            outObj = outObj.FindObject(containerName)
        elif outObj == None:
            print(f'\033[31mError\033[0m: {containerName} is None. Exit!')
            sys.exit()

        else:
            print(f'\033[31mError\033[0m: instance of {type(outObj)} not implemented. Exit!')
            sys.exit()
    
    return outObj