# HACER

- Ponerle a cada DataFrame una m asociada.

```python
# Leer la masa de la primera línea
with open('../Data/1.csv', encoding='utf-8') as f:
    first_line = f.readline()
    m = float(first_line.strip().split('=')[1])

# Leer el DataFrame saltando la primera línea
df = pd.read_csv('../Data/1.csv', skiprows=1)

# Guardar la masa como atributo del DataFrame
df.m = m
```
