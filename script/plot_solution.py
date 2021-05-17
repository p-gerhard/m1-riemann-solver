#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":   
    try:
        filename = sys.argv[1]
    except IndexError:
        raise IndexError("the solution data file must be passed as argument")
    
    data = np.genfromtxt(filename, delimiter=' ')
    
    plt.plot(data[:,0], data[:,1], label=r'rho')
    plt.plot(data[:,0], data[:,2], label=r'I')
    plt.ylabel('w')
    plt.xlabel('x')
    
    plt.legend(loc=0)
    
    plt.show()