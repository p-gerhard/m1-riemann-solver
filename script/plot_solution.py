#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":   
    filename = "u_ex_s2.out"
    data = np.genfromtxt(filename, delimiter=' ')
    
    plt.plot(data[:,0], data[:,1], label=r'rho')
    plt.plot(data[:,0], data[:,2], label=r'I')
    plt.ylabel('w')
    plt.xlabel('x')
    
    plt.legend(loc=0)
    
    plt.show()