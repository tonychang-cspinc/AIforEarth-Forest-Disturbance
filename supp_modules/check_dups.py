import pandas as pd
import numpy as np

final=pd.read_csv('example_2states.csv')
values=final['CN'].values
print(len(np.unique(values)),len(values))
