# -*- coding: utf-8 -*-
"""Untitled12.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1sdZBfYNLZnZo8mNbSzqZHhm-HG0emYSe
"""

import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
al_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/Al_Annealed_1150_T0.csv"
al_hard_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/Al_Work_hardened_5005_H34.csv"
cu_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/Copper_Annealed.csv"
ldpe_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/LDPE.csv"
ldpe_cold_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/LDPE_cold_down.csv"
steel_dat = "https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/all_material_instron_data/low_carbon_steel.csv"

recorded_data_raw = 'https://raw.githubusercontent.com/Tclack88/MatlabSucks/main/mechanics_and_materials/materials_practical_data/instrondata_machine3_15sep2022_2pm/machine_settings.csv'
old_cols = ['Specimen', ' Cross-Head Speed (mm/min)', ' Initial Guage Length (mm)',
       ' Final Gauge Length (mm)', ' thickness before (mm)',
       ' width before (mm)', ' thickness after (mm)', ' width after(mm)']
new_cols = ['specimen', 'speed', 'length_o', 'length_f', 'thickness_o',
       'width_o', 'thickness_f', 'width_f']
col_rename = dict(zip(old_cols,new_cols))
data_df = pd.read_csv(recorded_data_raw).rename(columns=col_rename)
data_df['spec'] = 'cu al al_hard steel ldpe ldpe_cold'.split() # shortened name
data_df['area_o'] = data_df.thickness_o*data_df.width_o # square mm
data_df['area_f'] = data_df.thickness_f*data_df.width_f # square mm
data_df

old_cols = ['Time measurement', 'Displacement', 'Primary force measurement',
       'Strain 1']
new_cols = ['time','displacement','force','strain']
col_rename = dict(zip(old_cols,new_cols))
# Units are s, mm, N, %

# make df
al_df = pd.read_csv(al_dat).rename(columns=col_rename)[1:].astype(float)
al_hard_df = pd.read_csv(al_hard_dat).rename(columns=col_rename)[1:].astype(float)
cu_df = pd.read_csv(cu_dat).rename(columns=col_rename)[1:].astype(float)
ldpe_df = pd.read_csv(ldpe_dat).rename(columns=col_rename)[1:].astype(float) # only has 3 columns, no extensiometer was used
ldpe_cold_df = pd.read_csv(ldpe_cold_dat).rename(columns=col_rename)[1:].astype(float)
steel_df = pd.read_csv(steel_dat).rename(columns=col_rename)[1:].astype(float)

material_dfs = [al_df, al_hard_df, cu_df, ldpe_df, ldpe_cold_df, steel_df]
material_dfs_labels = 'al_df al_hard_df cu_df ldpe_df ldpe_cold_df steel_df'.split()
df_labels = dict(zip(material_dfs_labels,material_dfs))

# Attach final and initial areas (mm^2) to dfs
for label in df_labels:
  spec = ('_'.join(str(label).split('_')[:-1]))
  area_o = data_df.area_o[data_df.spec == spec].iloc[0]
  area_f = data_df.area_f[data_df.spec == spec].iloc[0]
  length_o = data_df.length_o[data_df.spec == spec].iloc[0]
  df_labels[label]['area_o'] = area_o
  df_labels[label]['area_f'] = area_f
  df_labels[label]['stress'] = df_labels[label]['force']/df_labels[label]['area_o']
  df_labels[label]['strain_calc'] = df_labels[label]['displacement']/length_o
  try:
    df_labels[label]['strain'] = df_labels[label]['strain']/100 # convert percent to fraction for equal comparison
  except:
    pass # no strain exists

al_df

ldpe_df.plot('strain_calc','stress')

sb.set_style("whitegrid")
for label in df_labels:
  plt.figure()
  # plt.figure(figsize=(10,7))
  df = df_labels[label]
  spec = ('_'.join(str(label).split('_')[:-1]))
  material_name =  data_df.specimen[data_df.spec == spec].iloc[0]
  # print(material_name)
  scatterplot = sb.scatterplot(data=df,x='strain_calc',y='stress')
  try:
    sb.scatterplot(data=df,x='strain',y='stress')
  except:
    pass # no strain available (LDPE)
  scatterplot.set(title=f'stress-strain of {material_name}',xlabel='Strain (mm/mm)', ylabel='Engineering Stress (MPa)')
  plt.legend(labels=['calculated strain','extensometer strain'], loc='lower right')
  plt.show()

for label in df_labels:
  plt.figure()
  df = df_labels[label]
  spec = ('_'.join(str(label).split('_')[:-1]))
  material_name =  data_df.specimen[data_df.spec == spec].iloc[0]
  print(material_name)


