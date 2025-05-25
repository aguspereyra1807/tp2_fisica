import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calcs import DF, DF2, angularFreq, oscilationFreq, radius

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
##### m, θ constantes: L1, L2, L3 

def differentLenght():
    w_exp = [df.w for df in DF2]
    l_values = [radius(df) for df in DF2]  # l es el promedio de r en cada experimento

    plt.figure(figsize=(12,8))
    plt.plot(l_values, w_exp, '-o', label=r'$\omega$ experimental', markersize=8, color='red', markeredgecolor='black')

    plt.xlabel(r'Longitud del péndulo $L$ [cm]')
    plt.ylabel(r'Frecuencia angular $\omega$ [rad/s]')
    plt.title(r'$\omega$ vs $L$')
    plt.grid(True)
    plt.legend()
    plt.savefig("../Graphs/differentLength.png", bbox_inches='tight')

def plot_trajectories_grid():
    # Usar solo los primeros 18 CSVs
    DF_sub = DF[:18]
    
    masses = [6.07, 22.15, 72.36]
    styles = ['-', '--', ':']
    markers = ['o', 's', '^']
    
    fig, axs = plt.subplots(3, 2, figsize=(14, 12), sharex=True, sharey=True)
    fig.suptitle("Trayectorias angulares (θ vs t)", fontsize=16)

    for i, m in enumerate(masses):
        df_mass = [df for df in DF_sub if df.m == m]

        # L largo (r >= 30)
        long_l = [df for df in df_mass if df['r'][0] >= 30]
        # L corto (r < 30)
        long_c = [df for df in df_mass if df['r'][0] < 30]

        for j, long_group in enumerate([long_l, long_c]):
            ax = axs[i, j]

            for k, df in enumerate(long_group):
                style = styles[k % len(styles)]
                marker = markers[k % len(markers)]
                label = fr"$\theta_0 \approx {max(df['θ']):.1f}$ rad"

                ax.plot(df['t'], df['θ'], linestyle=style, marker=marker,
                        label=label, markevery=15)

            title = "Larga" if j == 0 else "Corta"
            ax.set_title(f"m = {m} g, L {title}")
            ax.grid(True)
            if j == 0:
                ax.set_ylabel(r"$\theta(t)$ [rad]")
            if i == 2:
                ax.set_xlabel("Tiempo [s]")
            ax.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("../Graphs/Trajectories_Grid.png", bbox_inches='tight')
    plt.close()


plot_trajectories_grid()