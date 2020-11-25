# When a degrading system starts to approach its doomsday, it acts "funny."
# The goal of this code is to detect this "funny" signal, if present, in
# biological population (drosophila).

# Some technical stuff:
# "funny" signal --> Early Warning Signal
# reduction in larvae food --> degrading environment

# Notes to self:
# 1. Remember that you are no longer passing x5 as an arguement for Egg2Fecund
#    fuction.
# 2.

# Properties of the system:
LarFood=1.6
LarFoodminus=0.001
AdNut=1
AdNutminus=0
hatchability=0.98
Mc=1.1  # Critical mass --> minimum mass for larvae to become an adult
sex_ratio=0.5  # male:female
SenDen=0.17
SenSize=1.7
NoG=100  #previously 100
NoR=10 #previously 10
N_Eggs_init=30

import random as ran
import numpy as np
import statistics

# Function to generate flies per generation
def Egg2Fecund(N_Eggs, LarFood, AdNut, hatchability, Mc, sex_ratio, SenDen,
               SenSize):
    # Some scaling "constants":
    x1 = 4.2
    x2 = 1
    x3 = 0.009
    x41 = 1
    x42 = 1
    x5 = 85
    x6 = 2

    N_larvae = int(round(hatchability * N_Eggs))
    SD = 0.31  # Standard deviation for larvae size
    fecundity = []
    if N_larvae != 0:
        MeanLarvaeSize = x1 *(1 -(1.0 /(x2 +np.exp(-x3 *N_larvae /LarFood))))
    else:
        MeanLarvaeSize = 0
    larvae_size = np.random.normal(MeanLarvaeSize, SD, N_larvae)
    adult_size = [x42 * j for j in larvae_size if j > Mc]
    N_Adult = int(round(x41 * len(adult_size)))
    adult_size = ran.sample(adult_size, N_Adult)
    if adult_size == []:
        MeanAdultSize = 0
    else:
        MeanAdultSize = statistics.mean(adult_size)
    SizeDepFecundity = AdNut * x5 * np.log(x6 + SenSize * np.array(adult_size))
    DenDepFecundity = 1 / (1 + SenDen * N_Adult)
    fecundity = SizeDepFecundity * DenDepFecundity
    for k in range(N_Adult):
        if ran.random() < sex_ratio:
            fecundity[k] = -1
    NMales = np.count_nonzero(fecundity == -1)
    fecundity = list(fecundity)
    return (fecundity, N_larvae, MeanLarvaeSize, MeanAdultSize)
