import numpy as np
import matplotlib.pyplot as plt
from numba import jit, njit

def koch_generator(start:np.array, end:np.array)->np.array:
    start = np.array(start)
    end = np.array(end)
    
    dx, dy = (end - start) / 4  # steps

    p1 = start
    p2 = p1 + np.array([dx, dy]) 
    p3 = p2 + np.array([-dy, dx]) # raise
    p4 = p3 + np.array([dx, dy])  
    p5 = p4 + 2*np.array([dy, -dx])   # lower
    p6 = p5 + np.array([dx, dy])
    p7 = p6 + np.array([-dy, dx]) # raise
    p8 = end

    return np.array([p1, p2, p3, p4, p5, p6, p7, p8])


def generate_square(lwl:int , L:int) ->np.array:
    pass


print(koch_generator((np.array([0,0])), np.array([0,2])))