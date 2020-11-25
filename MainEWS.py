# When a degrading system starts to approach its doomsday, it acts "funny."
# The goal of this code is to detect this "funny" signal, if present, in
# biological population (drosophila).

# Some technical stuff:
# "funny" signal --> Early Warning Signal
# reduction in larvae food --> degrading environment

# Notes to self:
# 1. Remember that you are no longer passing x5 as an arguement for Egg2Fecund,
#    timeseries and replicates fuctions.
# 2. Remember that you changed position of ext variable in return statement of
#    timeseries function.
# 3. Remember that you changed position of ext_count variable in return
#    statement of replicates function.

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
    N_males = np.count_nonzero(fecundity == -1)
    fecundity = list(fecundity)
    return (fecundity, N_larvae, MeanLarvaeSize, MeanAdultSize)

# Generating timeseries
def timeseries(N_Eggs, LarFood, AdNut, NoG, hatchability, Mc, sex_ratio,
               SenDen, SenSize, LarFoodminus, AdNutminus):
    adult_tseries = []
    egg_tseries = []
    eggsperfemale_tseries = []
    larvae_tseries = []
    MeanLarvaeSize_tseries = []
    MeanAdultSize_tseries = []
    ext=0  # variable to count number of extinctions
    for t in range(NoG):
        egg_tseries += [N_Eggs]
        fecj, N_larvae, MeanLarvaeSize, MeanAdultSize = Egg2Fecund(N_Eggs,
                                                        LarFood, AdNut,
                                                        hatchability, Mc,
                                                        sex_ratio, SenDen,
                                                        SenSize)
        N_females = np.count_nonzero(fecundity != -1)
        adult_tseries += [len(fecj)]
        larvae_tseries += [N_larvae]
        MeanLarvaeSize_tseries += [MeanLarvaeSize]
        MeanAdultSize_tseries += [MeanAdultSize]
        N_Eggs = sum([int(round(x)) for x in fecj if x != -1])
        if N_females !=0:
            eggsperfemale_tseries += [N_Eggs/N_females]
        else:
            eggsperfemale_tseries += [0]
        if len(fecj) == 0 & ext == 0:
            ext = 1
        LarFood, AdNut = food_decrease(LarFood, LarFoodminus, AdNut,
                         AdNutminus)
    return [egg_tseries, adult_tseries, larvae_tseries,
            MeanLarvaeSize_tseries, MeanAdultSize_tseries,
            eggsperfemale_tseries, ext]

# Decreasing Larval and Adult food
def food_decrease(LarFood, LarFoodminus, AdNut, AdNutminus):
    if LarFood>LarFoodminus:
        LarFood = LarFood -LarFoodminus
    if AdNut>AdNutminus:
        AdNut = AdNut -AdNutminus
    return (LarFood, AdNut)

# Replicating timeseries
def replicates(LarFood, AdNut, hatchability, Mc, sex_ratio, SenDen, SenSize,
               NoG, NoR, N_Eggs_init, LarFoodminus, AdNutminus):
    all_egg = []
    all_adult = []
    all_larvae = []
    all_larvae_size = []
    all_adult_size = []
    all_eggperfemale = []
    N_Eggs = N_Eggs_init
    ext_count=0
    for rep in range(NoR):
        [N_egg, N_adult, N_larvae,
        MeanLarvaeSize_tseries,
        MeanAdultSize_tseries,
        eggsperfemale_tseries, ext]= timeseries(N_Eggs, LarFood, AdNut,
                                                NoG, hatchability, Mc,
                                                sex_ratio, SenDen, SenSize,
                                                LarFoodminus, AdNutminus)
        all_egg += [N_egg]
        all_adult += [N_adult]
        all_larvae += [N_larvae]
        all_larvae_size += [MeanLarvaeSize_tseries]
        all_adult_size += [MeanAdultSize_tseries]
        all_eggperfemale += [eggsperfemale_tseries]
        ext_count += ext
    return (all_egg, all_adult, all_larvae, all_larvae_size, all_adult_size,
            all_eggperfemale, ext_count)

all_egg, all_adult, all_larvae, all_larvae_size, all_adult_size, all_eggperfemale, ext_count = replicates(LarFood, AdNut, hatchability,
                                                                                               Mc, sex_ratio, SenDen, SenSize, NoG, NoR,
                                                                                               N_Eggs_init, LarFoodminus, AdNutminus)
