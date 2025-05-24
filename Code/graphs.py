import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calcs import DF, angularFreq, oscilationFreq, radius

# DF por masas
m1 = [df for df in DF if df.m == 6.07]
m2 = [df for df in DF if df.m == 22.15]
m3 = [df for df in DF if df.m == 72.36]
dfM = [m1,m2,m3]

#### [1] (Información general)
##### ω en función de θ inicial

def AngleVsAngular():
    plt.figure(figsize=(12,8))

    col = { 
        6.07 : 'cyan',
        22.15 : 'violet',
        72.36 : 'yellow'
    }

    for m0 in dfM:
        x1 = [max(df['θ']) for df in m0 if df['r'][0] >= 30]
        y1 = [angularFreq(df, df.period) for df in m0 if df['r'][0] >= 30]

        x2 = [max(df['θ']) for df in m0 if df['r'][0] < 30]
        y2 = [angularFreq(df, df.period) for df in m0 if df['r'][0] < 30]

        mass = m0[0].m
        label1 = rf'$m = {mass} \pm 0.01$g, $L \approx 32.8cm$'
        label2 = rf'$m = {mass} \pm 0.01$g, $L \approx 21.2cm$'

        plt.plot(x1,y1, 'X', label=label1, markersize=8, color=col[mass], markeredgecolor='black')
        plt.plot(x2,y2, 'o', label=label2, markersize=8, color=col[mass], markeredgecolor='black')

    plt.xlabel(r"Ángulo inicial $\theta_0$ [ºC]")
    plt.ylabel(r"Frecuencia angular $\omega$ [rad/s]")
    plt.grid(True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.title(r"$\theta_0$ vs $\omega$")

    plt.savefig("../Graphs/AngleVsAngularFreq.png", bbox_inches='tight')

##### [2]