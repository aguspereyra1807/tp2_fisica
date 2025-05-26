import pandas as pd
import numpy as np
from scipy import signal, constants

# Cargar DataFrames
DF = []
for i in range(1,19):
    with open(f'../Data/{i}.csv', encoding='utf-8') as f:
        m = float(f.readline().strip().split('=')[1])
    df = pd.read_csv(f'../Data/{i}.csv', skiprows=1)
    df.m = m
    DF.append(df)

DF2 = []
for i in range(19,22):
    with open(f'../Data/{i}.csv', encoding='utf-8') as f:
        m = float(f.readline().strip().split('=')[1])
    df = pd.read_csv(f'../Data/{i}.csv', skiprows=1)
    df.m = m
    DF2.append(df)

############################## Calculos generales #############################

def period(df: pd.DataFrame):
    ''' T = t[θmax2] - t[θmax1] '''

    maxs, _ = signal.find_peaks(df["θ"]) # Encuentra picos del ángulo encuentro 2 θmax
    periods = []
    for i in range(len(maxs)-1):
        t1 = df["t"].iloc[maxs[i+1]]
        t0 = df["t"].iloc[maxs[i]]
        periods.append(t1-t0)  
    return round(np.mean(periods), 2)

def angularFreq(df: pd.DataFrame, t: float):
    ''' ω = 2π / T '''
    
    return round((2*np.pi) / t, 2)

def oscilationFreq(df: pd.DataFrame, t: float):
    ''' f = 1 / T '''

    return round(1/t, 2) 

def radius(df: pd.DataFrame):
    return round(np.mean(df["r"]), 2)

########## Cálculos para pequeñas oscilaciones (Small Oscilations SO) ##########

def angularFreqSO(df: pd.DataFrame):
    ''' ω = sqrt( g / l )'''

    l = np.mean(df["r"])
    return round(np.sqrt(constants.g / l), 2)

def periodSO(df: pd.DataFrame):
    ''' T = ω / 2π'''

    w = angularFreqSO(df)
    return round(w / (2*np.pi), 2)

def oscilationFreqSO(df: pd.DataFrame):
    ''' f = 2π / ω'''

    w = angularFreqSO(df)
    return round((2*np.pi) / w, 2)

############### Cargo valores T, ω, f

for df in DF:
        df.period = period(df)
        df.w = angularFreq(df, df.period)
        df.f = oscilationFreq(df, df.period)

for df2 in DF2:
        df2.period = period(df2)
        df2.w = angularFreq(df2, df2.period)
        df2.f = oscilationFreq(df2, df2.period)

############### Cálculos varios

def estimateGravity():
    samples = []
    for df in DF:
        r = np.mean(df['r']) / 100 # cm -> m
        samples.append( (4 * r * np.pi ** 2) / (df.period ** 2) )
    return round(np.mean(samples), 2), samples

if __name__ == '__main__':
    print(estimateGravity()[0])

