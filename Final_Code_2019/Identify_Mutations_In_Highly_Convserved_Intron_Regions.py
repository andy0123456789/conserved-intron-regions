'''dataframe'''
import pandas as pd
from pandas import DataFrame

sns.set_style("white")
Mutations = pd.read_table("data/mc3.v0.2.8.PUBLIC.maf",sep="\s+", header=0)
Newmutations = Mutations.iloc[0:10000000000000,[4,5,6]]

data = pd.read_table("data/peak_df2.csv", header = None)

new = new

new = data[0].str.split("_", n = 1, expand = True) 
  
new["Numbers3"]= new[0] 
  

new["Huandy3"]= new[1] 

new = new.drop(new.columns[[0,2,3]], axis=1)

new2 = new

new2 = new[1].str.split("_", n = 1, expand = True) 
  
new2["Numbers3"]= new2[0] 
  

new2["Huandy3"]= new2[1] 

new2 = new2.drop(new2.columns[[0,2,3]], axis=1)
new3 = new2

new3 = new2[1].str.split("_", n = 1, expand = True) 
  
new3["Numbers3"]= new3[0] 
  

new3["Huandy3"]= new3[1] 
new3 = new3.drop(new3.columns[[2,3]], axis=1)

new3

new3.dropna(inplace = True) 

new4 = new3

new4 = new3[1].str.split("_", n = 1, expand = True) 
  
new4["Numbers3"]= new4[1] 
  

new5 = new3[[0]]
new6 = new4[[0,1]]
new5
new7 = pd.concat([new5, new6], axis=1)
new7.columns = ['chromosome', 'start', 'end']


new8 = new7['chromosome'].str.split("r", n = 3, expand = True) 
  
new8 = new8.drop(new8.columns[[0]], axis = 1)
new9 = pd.concat([new8, new6], axis = 1)
new9.columns = ['Chromosome', 'start', 'end']
new9

y = 0
filename = 'conserved_intron_regions/data/test0.csv'
Match_Peaks_Mutations7 = new9.iloc[0]
while y < 266:
    for index, row in Newmutations.iterrows():
        if str(row['Chromosome']) == 'X' and str(row['Start_Position']) == '129498779' and str(row['End_Position']) == '129498779':
            y += 1
            print(y, '  ',end='')
        if str(row['Chromosome']) == str(new9.iloc[y]['Chromosome']) and int(row["Start_Position"]) >= int(new9.iloc[y]['start']) and int(row["End_Position"]) <= int(new9.iloc[y]['end']):
            print("-------------","\n",row,new9.iloc[y, :], "\n","-------------")
            Match_Peaks_Mutations7 = Match_Peaks_Mutations7.append(row)
            Match_Peaks_Mutations7 = Match_Peaks_Mutations7.append(new9.iloc[y]
Match_Peaks_Mutations7.to_csv('Test11.csv', sep='\t')
Match_Peaks_Mutations2 = pd.read_table("Test11.csv",sep="\s+", header=0)
Match_Peaks_Mutations2
i = 1
while i < 10000:
    if i % 8 == 0:
        Match_Peaks_Mutations2[i - 6:i + 2].to_csv('data/ConservationMutations' + str(i) + '.csv', sep='\t')
        i += 1
    else:
        i += 1
                                                                   
