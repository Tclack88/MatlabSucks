import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

steel_data = {"energy":[45.1,26.5,25.5,32.3,8.8,12.7],'temp':[50,35,25,0,-15,-30]}
pvc_data = {'energy':[.728,.867,.762,.786,.346,.355,.37,.357],'temp':[25,25,25,25,-30,-30,-30,-30],'isavg':[0,0,0,1,0,0,0,1]}
hdpe_data = {'energy':[1.278,2.156,1.289,1.57,.909,.802,1.047,.919],'temp':[25,25,25,25,-30,-30,-30,-30],'isavg':[0,0,0,1,0,0,0,1]}

steel_df = pd.DataFrame(steel_data)
pvc_df = pd.DataFrame(pvc_data)
hdpe_df = pd.DataFrame(hdpe_data)

sb.set_style("whitegrid")
plt.figure(figsize=(10,8))
scatter = sb.scatterplot(data=pvc_df,x='temp',y='energy',hue='isavg', s=500, marker='x', legend=False)
scatter.set_title('PVC energy absorbed by temperature',size=20)
scatter.set_xlabel('Temperature (C)',size=15)
scatter.set_ylabel('Energy Absorbed (J)',size=15)
plt.show()

plt.figure(figsize=(10,8))
scatter = sb.scatterplot(data=hdpe_df,x='temp',y='energy',hue='isavg', s=500, marker='x', legend=False)
scatter.set_title('HDPE energy absorbed by temperature',size=20)
scatter.set_xlabel('Temperature (C)',size=15)
scatter.set_ylabel('Energy Absorbed (J)',size=15)
plt.show()

plt.figure(figsize=(10,8))
scatter = sb.scatterplot(data=steel_df,x='temp',y='energy', s=500, marker='x', legend=False)
scatter.set_title('Low Carbon Steel energy absorbed by temperature',size=20)
scatter.set_xlabel('Temperature (C)',size=15)
scatter.set_ylabel('Energy Absorbed (J)',size=15)
plt.show()

