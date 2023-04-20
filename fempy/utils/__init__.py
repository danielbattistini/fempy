import fempy
from .analysis import Pair
from .io import *
from .format import *

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
    