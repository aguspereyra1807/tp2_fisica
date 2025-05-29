import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calcs import DF, DF2, angularFreq, oscilationFreq, radius
from scipy.signal import find_peaks
from scipy.stats import linregress
from collections import defaultdict
import matplotlib.cm as cm


# DF por masas
m1 = [df for df in DF if df.m == 6.07]
m2 = [df for df in DF if df.m == 22.15]
m3 = [df for df in DF if df.m == 72.36]
dfM = [m1,m2,m3]

#### [1] ω vs θ inicial con propagación de errores
##### ω en función de θ inicial

def angleVsAngular():
    plt.figure(figsize=(12,8))

    col = { 
        6.07 : "#b00000",
        22.15 : "#17bc96",
        72.36 : "#cfcc0c"
    }

    ecol = {
        6.07 : "#B52121",
        22.15 : "#19c29a",
        72.36 : "#d9d626"
    }

    for m0 in dfM:
        x1 = [max(df['θ']) for df in m0 if df['r'][0] >= 30]
        y1 = [df.w for df in m0 if df['r'][0] >= 30]
        err1 = [df.angularError for df in m0 if df['r'][0] >= 30]

        x2 = [max(df['θ']) for df in m0 if df['r'][0] < 30]
        y2 = [df.w for df in m0 if df['r'][0] < 30]
        err2 = [df.angularError for df in m0 if df['r'][0] < 30]

        mass = m0[0].m
        label1 = rf'$m = {mass} \pm 0.01$g, $L \approx 32.8cm$'
        label2 = rf'$m = {mass} \pm 0.01$g, $L \approx 21.2cm$'

        plt.errorbar(x1,y1, err1, fmt='^', label=label1, markersize=6, color=col[mass], markeredgecolor='black', ecolor=ecol[mass], elinewidth=1, capsize=3)
        plt.errorbar(x2,y2, err2, fmt='v', label=label2, markersize=6, color=col[mass], markeredgecolor='black', ecolor=ecol[mass], elinewidth=1, capsize=3)

    plt.xlabel(r"Ángulo inicial $\theta_0$ [ºC]")
    plt.ylabel(r"Frecuencia angular $\omega$ [rad/s]")
    plt.grid(True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.title(r"$\theta_0$ vs $\omega$ con propagación de errores")

    plt.savefig("../Graphs/AngleVsAngularFreq2.png", bbox_inches='tight')


##### [2] Variación de ω con el largo 
##### m, θ constantes: L1, L2, L3 

def differentLenght():
    w_exp = [df.w for df in DF2]
    l_values = [radius(df) for df in DF2]  # l es el promedio de r en cada experimento
    initialMean = np.mean([max(DF2[0]['θ']), max(DF2[1]['θ']), max(DF2[2]['θ'])])

    plt.figure(figsize=(12,8))

    slope, intercept, _, _, _ = linregress(l_values, w_exp)
    t_fit = np.linspace(min(l_values), max(l_values), 100)
    theta_fit = slope * t_fit + intercept
    plt.plot(t_fit, theta_fit, linestyle='--', color='#ff484d', label=rf'Ajuste lineal')

    plt.plot(l_values, w_exp, 'o', markersize=8, color='red', markeredgecolor='black', label=rf'$\theta_0 \approx {initialMean}°, m = $')

    plt.xlabel(r'Longitud del péndulo $L$ [cm]')
    plt.ylabel(r'Frecuencia angular $\omega$ [rad/s]')
    plt.title(r'$\omega$ vs $L$')
    plt.grid(True)
    plt.legend()
    plt.savefig("../Graphs/differentLength.png", bbox_inches='tight')

##### [3] Todas las trayectorias
##### θ(t)

def trajectoriesGrid():
    colors = [ "#b0009e", "#1772bc", "#0ccf43"]
    
    DF_sub = DF[:18]
    
    masses = [6.07, 22.15, 72.36]
    styles = ['-', '--', ':']
    markers = ['o', 's', '^']
    
    fig, axs = plt.subplots(6, 1, figsize=(8, 24), sharex=True)
    fig.suptitle(r"Trayectorias angulares ($\theta$ vs $t$)", fontsize=16)

    for i, m in enumerate(masses):
        df_mass = [df for df in DF_sub if df.m == m]

        # L largo (r >= 30)
        long_l = [df for df in df_mass if df['r'][0] >= 30]
        # L corto (r < 30)
        long_c = [df for df in df_mass if df['r'][0] < 30]

    
        idx = 0
        for i, m in enumerate(masses):
            df_mass = [df for df in DF_sub if df.m == m]
            long_l = [df for df in df_mass if df['r'][0] >= 30]
            long_c = [df for df in df_mass if df['r'][0] < 30]
            for j, long_group in enumerate([long_l, long_c]):
                ax = axs[idx]
                for k, df in enumerate(long_group):
                    style = styles[k % len(styles)]
                    marker = markers[k % len(markers)]
                    col = colors[k % len(colors)]
                    label = fr"$\theta_0 \approx {max(df['θ']):.1f}$°"
                    ax.plot(df['t'], df['θ'], linestyle=style, marker=marker,
                            label=label, markevery=15, color=col)
                title = "Larga" if j == 0 else "Corta"
                ax.set_title(rf"$m$ = {m} g, $L$ {title}")
                ax.grid(True)
                ax.set_ylabel(r"$\theta(t)$ [$°$]")
                if idx == 5:
                    ax.set_xlabel(r"Tiempo [$s$]")
                ax.legend(fontsize=8)
                idx += 1

    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig("../Graphs/TrajectoriesGrid.png", bbox_inches='tight')

###### [3] ω 
###### L vs ω y m vs ω

def angularVsMassLarge():
    w = [angularFreq(df.assign(θ=np.deg2rad(df['θ'])), df.period) for df in DF]
    m = [df.m for df in DF]  # masa en gramos
    l = [radius(df) for df in DF]  # longitud en cm

    _, axs = plt.subplots(1, 2, figsize=(14, 6))

    # ω vs masa
    axs[0].plot(m, w, 'o', color='blue', markeredgecolor='black')
    axs[0].set_xlabel("Masa [g]")
    axs[0].set_ylabel("Frecuencia angular ω [rad/s]")
    axs[0].set_title("ω vs masa")
    axs[0].grid(True)

    # ω vs longitud
    axs[1].plot(l, w, 'o', color='green', markeredgecolor='black')
    axs[1].set_xlabel("Longitud L [cm]")
    axs[1].set_ylabel("Frecuencia angular ω [rad/s]")
    axs[1].set_title("ω vs longitud")
    axs[1].grid(True)

    plt.tight_layout()
    plt.savefig("../Graphs/OmegaVsMasaYLongitud.png")

###### [4] Máximos y ajuste lineal
###### θ(t)

def gridMaxLinealRegression():
    masas = [6.07, 22.15, 72.36]
    largos_labels = ['Long. Mayor', 'Long. Menor']

    fig, axs = plt.subplots(3, 2, figsize=(16, 14), sharex=True, sharey=True)
    fig.suptitle(r"Amplitud de oscilacion de $\theta(t)$° con ajuste lineal", fontsize=19)

    for i, m in enumerate(masas):
        dfs_masa = [df for df in DF if df.m == m]

        # Separar por largo
        dfs_larga = [df for df in dfs_masa if df['r'].iloc[0] >= 30]
        dfs_corta = [df for df in dfs_masa if df['r'].iloc[0] < 30]

        colores = ["#0C1CAB", "#750567", "#2c94a0"] 
        
        for j, dfs_largo in enumerate([dfs_larga, dfs_corta]):
            ax = axs[i, j]

            for k, df in enumerate(dfs_largo):
                t = df['t'].values
                theta = df['θ'].values

                # Encontrar máximos
                peaks, _ = find_peaks(theta)
                t_peaks = t[peaks]
                theta_peaks = theta[peaks]

                # Color y estilo para esta curva
                color = colores[k % len(colores)]
                
                # Plot curva original sin marker, solo línea
                ax.plot(t, theta, '-', color=color, label=rf'$\theta_0 \approx$ {max(theta):.1f}°')

                # Plot máximos como puntos X (sin marker en curva, solo aquí)
                ax.plot(t_peaks, theta_peaks, color=color, marker="o", markersize=8, linestyle='None',
                        label=rf'Máximos $\theta_0 \approx$ {max(theta):.1f}°')

                # Ajuste lineal a los máximos
                if len(t_peaks) > 1:
                    slope, intercept, _, _, _ = linregress(t_peaks, theta_peaks)
                    t_fit = np.linspace(t_peaks.min(), t_peaks.max(), 100)
                    theta_fit = slope * t_fit + intercept
                    ax.plot(t_fit, theta_fit, linestyle='--', color=color,
                            label=rf'Ajuste max $\theta_0 \approx${max(theta):.1f}°')

            ax.set_title(rf'$m$ = {m} g, {largos_labels[j]}')
            ax.grid(True)
            if j == 0:
                ax.set_ylabel(r'Ángulo $\theta$ [°]')
            if i == 2:
                ax.set_xlabel(r'Tiempo [$s$]')
            ax.legend(fontsize=7)

    plt.tight_layout()
    plt.savefig("../Graphs/MaximosConFit.png")