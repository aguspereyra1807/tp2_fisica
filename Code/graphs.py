import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, constants

# Cargar todos los datos
DF = []
for i in range(1, 19):
    with open(f'../Data/{i}.csv', encoding='utf-8') as f:
        m = float(f.readline().strip().split('=')[1])  # Masa desde primera línea
    df = pd.read_csv(f'../Data/{i}.csv', skiprows=1)   # Datos: t, θ, r
    df['m'] = m
    df['id'] = i
    DF.append(df)

# Funciones de análisis
def period(df: pd.DataFrame):
    '''Calcula el período a partir de dos máximos consecutivos de θ'''
    maxs, _ = signal.find_peaks(df["θ"])
    periods = []
    for i in range(len(maxs)-1):
        t1 = df["t"].iloc[maxs[i+1]]
        t0 = df["t"].iloc[maxs[i]]
        periods.append(t1 - t0)
    return np.mean(periods)

def angularFreq(T: float):
    return (2 * np.pi) / T

# Extraer características de cada experimento
summary = []
for df in DF:
    T = period(df)
    w = angularFreq(T)
    l = np.mean(df["r"])
    m = df["m"].iloc[0]
    amp = np.abs(df["θ"]).max() * 180 / np.pi  # convertir a grados aprox
    # Clasificar ángulo como chico, medio o grande
    if amp < 10:
        ang_cat = 'chico'
    elif amp < 30:
        ang_cat = 'medio'
    else:
        ang_cat = 'grande'
    summary.append({'l': l, 'w': w, 'm': m, 'ángulo': ang_cat})

# Convertir a DataFrame
summary_df = pd.DataFrame(summary)

# Asignar colores y estilos
massas = sorted(summary_df['m'].unique())
angulos = ['chico', 'medio', 'grande']
colores = ['#1f77b4', '#2ca02c', '#d62728']  # azul, verde, rojo
line_styles = ['-', '--', ':']

color_dict = {m: c for m, c in zip(massas, colores)}
style_dict = {a: s for a, s in zip(angulos, line_styles)}

# Graficar
plt.figure(figsize=(10, 6))

for m in massas:
    for a in angulos:
        subset = summary_df[(summary_df['m'] == m) & (summary_df['ángulo'] == a)]
        if not subset.empty:
            plt.plot(subset['l'], subset['w'], 
                     linestyle=style_dict[a], 
                     color=color_dict[m], 
                     marker='o', 
                     label=f"m={m:.2f} kg, ángulo={a}")

plt.xlabel('Longitud l (m)')
plt.ylabel('Frecuencia angular ω (rad/s)')
plt.title('Frecuencia angular vs Longitud para distintas masas y ángulos')
plt.grid(True)
plt.legend(loc='best', fontsize=9)
plt.tight_layout()
plt.savefig('grafico_unico_w_vs_l.png', dpi=300)
plt.show()
