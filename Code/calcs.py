import pandas as pd
import numpy as np
from scipy import signal, constants

# Read CSVs
DF = [pd.read_csv(f'../Data/{i}.csv') for i in range(1,19)]

# Calculos generales

def period(df: pd.DataFrame):
    ''' T = t[θmax2] - t[θmax1] '''

    maxs, _ = signal.find_peaks(df["θ"]) # Encuentra picos del ángulo encuentro 2 θmax
    periods = []
    for i in range(len(maxs)-2):
        t1 = df["t"].iloc[maxs[i+1]]
        t0 = df["t"].iloc[maxs[i]]
        periods.append(t1-t0)  
    return(np.mean(periods))

def angularFreq(df: pd.DataFrame):
    ''' ω = 2π / T '''
    T = period(df)
    return((2*np.pi) / T)

def oscilationFreq(df: pd.DataFrame):
    ''' f = 1 / T '''
    T = period(df)
    return(1/T)

#Cálculos para pequeñas oscilaciones (Small Oscilations SO)

def angularFreqSO(df: pd.DataFrame):
    ''' ω = sqrt( g / l )'''
    l = np.mean(df["r"])
    return np.sqrt(constants.g / l)

def periodSO(df: pd.DataFrame):
    ''' T = ω / 2π'''
    w = angularFreqSO(df)
    return w / (2*np.pi)

def oscilationFreqSO(df: pd.DataFrame):
    ''' f = 2π / ω'''
    w = angularFreqSO(df)
    return (2*np.pi) / w

####################################################################################################################

def main():
    print(DF[0])
    print(period(DF[0]))
    print(angularFreq(DF[0]))
    print(oscilationFreq(DF[0]))


if __name__ == "__main__":
    main()