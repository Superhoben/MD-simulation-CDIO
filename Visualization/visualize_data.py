import os
import json
import matplotlib.pyplot as plt

all_data = {}
path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/Cu_Ag_mixed_with_Au'
for file in os.listdir(path):
    if file[-4:] == ".txt":
        
        opened_file = open(path + "/" + file, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)
        
        file = file[:-4]
        data = file.split("_")
        
        all_data[data[-1]] = material_data_dict
        
Cu_data = all_data['mp-30']
Ag_data = all_data['mp-124']

percentages = []
for percentage in Cu_data['origin_files']:
    percentages.append(int(percentage[0:2]))

plt.plot(percentages, Cu_data['bulk_modulus'], "r^", label="Cu")
plt.plot(percentages, Ag_data['bulk_modulus'], "b^", label="Ag")
plt.legend(loc="upper left")
plt.xlabel("Percentage")
plt.ylabel("Bulk modulus")
plt.show()