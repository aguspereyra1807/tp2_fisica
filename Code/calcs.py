import pandas as pd
import numpy as np
from scipy import signal, constants

# Read CSVs
DF = []
for i in range(1,19):
    with open(f'../Data/{i}.csv', encoding='utf-8') as f:
        m = float(f.readline().strip().split('=')[1])
    df = pd.read_csv(f'../Data/{i}.csv', skiprows=1)
    df.m = m
    DF.append(df)

############################## Calculos generales #############################

def period(df: pd.DataFrame):
    ''' T = t[θmax2] - t[θmax1] '''

    maxs, _ = signal.find_peaks(df["θ"]) # Encuentra picos del ángulo encuentro 2 θmax
    periods = []
    for i in range(len(maxs)-1):
        t1 = df["t"].iloc[maxs[i+1]]
        t0 = df["t"].iloc[maxs[i]]
        periods.append(t1-t0)  
    return np.mean(periods)

def angularFreq(df: pd.DataFrame, t: float):
    ''' ω = 2π / T '''
    return (2*np.pi) / t

def oscilationFreq(df: pd.DataFrame, t: float):
    ''' f = 1 / T '''
    return 1/t

########## Cálculos para pequeñas oscilaciones (Small Oscilations SO) ##########

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

############### Cargo valores T, ω, f

for df in DF:
        df.period = period(df)
        df.w = angularFreq(df, df.period)
        df.f = oscilationFreq(df, df.period)

if __name__ == '__main__':
    pass