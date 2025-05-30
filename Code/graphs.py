import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calcs import DF, DF2, angularFreq, oscilationFreq, radius
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.interpolate import make_interp_spline
from scipy import constants
from collections import defaultdict
import matplotlib.cm as cm


# DF por masas
m1 = [df for df in DF if df.m == 6.07]
m2 = [df for df in DF if df.m == 22.15]
m3 = [df for df in DF if df.m == 72.36]
dfM = [m1,m2,m3]

#### [1] (NO USADO) ω vs θ inicial con propagación de errores
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


##### [2] Variación de ω con largos (+ ajuste polinómico)
##### m, θ constantes: L1, L2, L3 

def angularVsL():
    wValues = [df.w for df in DF2]
    wErr = [df.angularError for df in DF2]
    lValues = [30, 20, 10]  # l es el promedio de r en cada experimento
    lErr = [0.1, 0.1, 0.1]
    initialMean = np.mean([max(DF2[0]['θ']), max(DF2[1]['θ']), max(DF2[2]['θ'])])

    plt.figure(figsize=(12,8))

    # Puntos + error de ω
    plt.errorbar(lValues, wValues, xerr=lErr, yerr=wErr, fmt='o', markersize=8, color='#ff484d', markeredgecolor='black', ecolor="#000000", elinewidth=1, capsize=3, barsabove=False, label=rf'$\theta_0 \approx {initialMean}°, m \approx 72.36g$')

    # Ajuste polimonial
    grado = 2
    coef = np.polyfit(lValues, wValues, grado)
    p = np.poly1d(coef)
    l_fit = np.linspace(min(lValues), max(lValues), 100)
    w_fit = p(l_fit)
    plt.plot(l_fit, w_fit, linestyle='--', color="#bc1f25", label=f'Ajuste polinómico (grado {grado})')

    plt.xlabel(r'Longitud del péndulo $L$ [cm]')
    plt.ylabel(r'Frecuencia angular $\omega$ [rad/s]')
    plt.title(r'$\omega$ vs $L$ con ajuste polinómico')
    plt.grid(True)
    plt.legend(fontsize=15, loc='upper right')
    plt.savefig("../Graphs/wVsL.png", bbox_inches='tight')

##### [3] Variación de ω con masas 
##### L, θ constantes: m1, m2, m3 
def angularVsM():
    colors = {
        0 : "#ff1616",
        1 : "#ff4124",
        2: "#ff8010"
    }

    dfs1 = [DF[0], DF[3], DF[6]]
    dfs2 = [DF[9], DF[12], DF[15]]
    dfs3 = [DF[11], DF[14], DF[17]]
    DFS = [dfs1, dfs2, dfs3]

    plt.figure(figsize=(12,8))

    for i, dfs in enumerate(DFS):
        wValues = [df.w for df in dfs]
        wErr = [df.angularError for df in dfs]
        mValues = [df.m for df in dfs] 
        mErr = [0.1, 0.1, 0.1]
        angleMean = round(np.mean([max(dfs[0]['θ']), max(dfs[1]['θ']), max(dfs[2]['θ'])]), 2)
        lMean = round(np.mean([max(dfs[0]['r']), max(dfs[1]['θ']), max(dfs[2]['θ'])]), 2)

        # Puntos + error de ω y m
        plt.errorbar(mValues, wValues, xerr=mErr, yerr=wErr, fmt='o', markersize=8, color=colors[i], markeredgecolor='black', ecolor="#000000", elinewidth=1, capsize=3, barsabove=False, label=rf'$\theta_0 \approx {angleMean}°, L \approx {lMean}cm$')

        # Interpolación
        l_fit = np.linspace(min(mValues), max(mValues), 200)
        spline = make_interp_spline(mValues, wValues, k=2)  # k=3 para cúbico
        w_fit = spline(l_fit)
        plt.plot(l_fit, w_fit, linestyle='--', color=colors[i], label='Interpolación (spline grado 2)')

    plt.xlabel(r'Masa $m$ [g]')
    plt.ylabel(r'Frecuencia angular $\omega$ [rad/s]')
    plt.title(r'$\omega$ vs $m$')
    plt.grid(True)
    plt.legend(fontsize=13, loc='center left', bbox_to_anchor=(0.6, 0.95))
    plt.savefig("../Graphs/wVsM.png", bbox_inches='tight')


##### [4] Todas las trayectorias
##### θ(t)

def trajectoriesGrid():
    colors = [ "#b0009e", "#1772bc", "#0ccf43"]
    
    DF_sub = DF[:18]
    
    masses = [6.07, 22.15, 72.36]
    styles = ['-', '--', ':']
    markers = ['o', 's', '^']
    
    fig1, axs1 = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    fig1.suptitle(r"Trayectorias angulares ($\theta$ vs $t$) de $L \approx 32.8$cm", fontsize=16)

    for i, m in enumerate(masses):
        df_mass = [df for df in DF_sub if df.m == m]
        long_l = [df for df in df_mass if df['r'][0] >= 30]
        ax = axs1[i]
        for k, df in enumerate(long_l):
            style = styles[k % len(styles)]
            marker = markers[k % len(markers)]
            col = colors[k % len(colors)]
            label = fr"$\theta_0 \approx {max(df['θ']):.1f}$°"
            ax.plot(df['t'], df['θ'], linestyle=style, marker=marker,
                    label=label, markevery=15, color=col)
        ax.set_title(rf"$m \approx$ {m} g")
        ax.grid(True)
        ax.set_ylabel(r"$\theta(t)$ [$°$]")
        if i == 2:
            ax.set_xlabel(r"Tiempo [$s$]")
        ax.legend(fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig("../Graphs/TrajectoriesGridLargos.png", bbox_inches='tight')

    fig2, axs2 = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    fig2.suptitle(r"Trayectorias angulares ($\theta$ vs $t$) de $L \approx 21.2$cm", fontsize=16)

    for i, m in enumerate(masses):
        df_mass = [df for df in DF_sub if df.m == m]
        long_c = [df for df in df_mass if df['r'][0] < 30]
        ax = axs2[i]
        for k, df in enumerate(long_c):
            style = styles[k % len(styles)]
            marker = markers[k % len(markers)]
            col = colors[k % len(colors)]
            label = fr"$\theta_0 \approx {max(df['θ']):.1f}$°"
            ax.plot(df['t'], df['θ'], linestyle=style, marker=marker,
                    label=label, markevery=15, color=col)
        ax.set_title(rf"$m \approx$ {m} g")
        ax.grid(True)
        ax.set_ylabel(r"$\theta(t)$ [$°$]")
        if i == 2:
            ax.set_xlabel(r"Tiempo [$s$]")
        ax.legend(fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig("../Graphs/TrajectoriesGridCortos.png", bbox_inches='tight')


###### (NO USADO) [5] ω en función de L y M 
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

###### [6] Máximos con ajuste lineal
###### θ(t)

def gridMaxLinealRegression():
    masas = [6.07, 22.15, 72.36]
    colores = ["#0C1CAB", "#750567", "#2c94a0"] 

    slopes = []

    # --- Largos mayores ---
    fig1, axs1 = plt.subplots(3, 1, figsize=(8, 14), sharex=True)
    fig1.suptitle(r"Amplitud de oscilación de $\theta(t)$° con ajuste lineal - $L \approx 32.8$cm", fontsize=17)

    for i, m in enumerate(masas):
        dfs_masa = [df for df in DF if df.m == m]
        dfs_larga = [df for df in dfs_masa if df['r'].iloc[0] >= 30]
        ax = axs1[i]
        for k, df in enumerate(dfs_larga):
            t = df['t'].values
            theta = df['θ'].values

            # Encontrar máximos
            peaks, _ = find_peaks(theta)
            t_peaks = t[peaks]
            theta_peaks = theta[peaks]

            color = colores[k % len(colores)]
            # Curva original
            ax.plot(t, theta, '-', color=color, label=rf'$\theta_0 \approx$ {max(theta):.1f}°')
            # Máximos
            ax.plot(t_peaks, theta_peaks, color=color, marker="o", markersize=8, linestyle='None',
                    label=rf'Máximos $\theta_0 \approx$ {max(theta):.1f}°')
            # Ajuste lineal
            if len(t_peaks) > 1:
                slope, intercept, _, _, _ = linregress(t_peaks, theta_peaks)
                
                if max(theta) <= 20: slopes.append(slope)

                t_fit = np.linspace(t_peaks.min(), t_peaks.max(), 100)
                theta_fit = slope * t_fit + intercept
                ax.plot(t_fit, theta_fit, linestyle='--', color=color,
                        label=rf'Ajuste max $\theta_0 \approx${max(theta):.1f}°')
        ax.set_title(rf'$m$ = {m} g, Long. Mayor')
        ax.grid(True)
        ax.set_ylabel(r'Ángulo $\theta$ [°]')
        if i == 2:
            ax.set_xlabel(r'Tiempo [$s$]')
        ax.legend(fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig("../Graphs/MaximosConFitLargos.png", bbox_inches='tight')

    # --- Largos menores ---
    fig2, axs2 = plt.subplots(3, 1, figsize=(8, 14), sharex=True)
    fig2.suptitle(r"Amplitud de oscilación de $\theta(t)$° con ajuste lineal - $L \approx 21.2$cm", fontsize=17)

    for i, m in enumerate(masas):
        dfs_masa = [df for df in DF if df.m == m]
        dfs_corta = [df for df in dfs_masa if df['r'].iloc[0] < 30]
        ax = axs2[i]
        for k, df in enumerate(dfs_corta):
            t = df['t'].values
            theta = df['θ'].values

            # Encontrar máximos
            peaks, _ = find_peaks(theta)
            t_peaks = t[peaks]
            theta_peaks = theta[peaks]

            color = colores[k % len(colores)]
            # Curva original
            ax.plot(t, theta, '-', color=color, label=rf'$\theta_0 \approx$ {max(theta):.1f}°')
            # Máximos
            ax.plot(t_peaks, theta_peaks, color=color, marker="o", markersize=8, linestyle='None',
                    label=rf'Máximos $\theta_0 \approx$ {max(theta):.1f}°')
            # Ajuste lineal
            if len(t_peaks) > 1:
                slope, intercept, _, _, _ = linregress(t_peaks, theta_peaks)
                
                if max(theta) <= 20: slopes.append(slope)

                t_fit = np.linspace(t_peaks.min(), t_peaks.max(), 100)
                theta_fit = slope * t_fit + intercept
                ax.plot(t_fit, theta_fit, linestyle='--', color=color,
                        label=rf'Ajuste max $\theta_0 \approx${max(theta):.1f}°')
        ax.set_title(rf'$m$ = {m} g, Long. Menor')
        ax.grid(True)
        ax.set_ylabel(r'Ángulo $\theta$ [°]')
        if i == 2:
            ax.set_xlabel(r'Tiempo [$s$]')
        ax.legend(fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig("../Graphs/MaximosConFitCortos.png", bbox_inches='tight')
    print(fr"La pendiente (promedio) del ajuste lineal para pequeñas oscilaciones ($\theta_0 \approx 10º$) es $a = {round(np.mean(slopes), 2)}$")

###### [7] (TERMINAR, DA CUALQ COSA) T vs L

def perdiodVsL():
    periods = [df.period for df in DF]
    periodErr = [df.periodError for df in DF]
    lengths = [radius(df)/100 for df in DF] 
    lengthsErr = [0.1/100 for _ in range(len(lengths))]

    plt.figure(figsize=(12,8))
    plt.errorbar(x=periods, y=lengths, xerr=periodErr, yerr=lengthsErr, fmt='o', markersize=10, color='grey', markeredgecolor='black', ecolor="#000000", elinewidth=1, capsize=4, barsabove=False)

    coef = np.polyfit(periods, lengths, 1)
    a, b = coef
    T_fit = np.linspace(min(periods), max(periods), 100)
    L_fit = a * T_fit + b
    plt.plot(T_fit, L_fit, 'r--', label=fr'Ajuste lineal: $L = {a:.2f}T + {b:.2f}$')

    plt.title(r'$T$ vs $L$')
    plt.grid()
    plt.xlabel(r'Período de oscilación $T$ [s]')
    plt.ylabel(r'Longitud del péndulo $L$ [m]')
    plt.tight_layout()

    print(a)

    plt.show()

###### [8] θ(t) θ chico vs θ(t) teórico

def oscilationComparative(dfIndex: int):
    ''' Teórico -> θ(t) = θ_0 cos ( sqrt(g/l) * t)'''
    df = DF[dfIndex] # m3, L1
    x = df['t']
    l = round(np.mean(df['r']), 2)

    y1 = [np.deg2rad(theta)+0.025 for theta in df['θ']] # + correción del ángulo
    theta0rad = max(y1)
    theta0deg = np.rad2deg(theta0rad)
    
    y2 = theta0rad* np.cos( np.sqrt( constants.g / ((l-5.8)/100) ) * x) # l a metros y ajuste del error 
    # y2 = theta0rad* np.cos( df.w * x) # l a metros y ajuste del error 

    plt.figure(figsize=(12,8))

    plt.plot(x, y2, '--', label=rf'Trayectoria Péndulo Teórico', color="#ff066a")
    plt.plot(x, y1, '-', label=rf'Trayectoria Péndulo Físico', color="#0380fc", linewidth=3, alpha=0.6)
    
    plt.grid()
    plt.xlabel('Tiempo [s]')
    plt.ylabel(r'$\theta(t)$ [rad]')

    plt.title(rf"$\theta(t)$ Teórico (Armónico) vs $\theta(t)$ Real ($\theta_0 \approx {round(theta0rad, 2)}$rad, $m \approx {df.m}g, L \approx {round(l-5.8, 2)}cm$)", fontsize=16)
    
    plt.legend(fontsize=14, loc='upper right')
    plt.tight_layout()
    plt.savefig(f'../Graphs/comparativeOsc(df[{dfIndex}]).png')