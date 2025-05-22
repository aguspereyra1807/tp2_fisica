import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

# Read CSVs
DF = [pd.read_csv(f'../Data/{i}.csv') for i in range(1,22)]