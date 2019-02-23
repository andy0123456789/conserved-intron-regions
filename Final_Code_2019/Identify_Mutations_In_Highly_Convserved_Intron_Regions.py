'''dataframe'''
import pandas as pd
from pandas import DataFrame

sns.set_style("white")
Mutations = pd.read_table("data/mc3.v0.2.8.PUBLIC.maf",sep="\s+", header=0)
Newmutations = Mutations.iloc[0:10000000000000,[4,5,6]]
