import os
def newDirect(path):
    """
    This script opens the specified save path and finds the highest folder number
    It then returns the next highest number as a name for the currently running simulation.
    """
    i=1
    test=2
    files =os.listdir(path)
    n=len(files)
    numFiles=[]
    if n>0:
        for i in range(n):
            try:
                numFiles.append(float(files[i]))
            except ValueError:
                pass
        if len(numFiles)>0:
            k=max(numFiles)
        else:
            k=0
    else:
        k=0
    return k+1.0

