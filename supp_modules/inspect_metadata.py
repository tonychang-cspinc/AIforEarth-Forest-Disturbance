import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


meta = pd.read_csv('csv_data/ALL_TREES.csv', low_memory=False).fillna(0)
nonwdlmeta = np.unique(meta[meta['DRYBIO_BOLE'] == 0]['PLT_CN'].tolist())
nonwdlmeta = pd.DataFrame(pd.Series(nonwdlmeta, index=range(len(nonwdlmeta))))
nonwdlmeta.to_csv('Non-woodland_CNs.csv', index=False)
print(len(meta), len(np.unique(nonwdlmeta)))
