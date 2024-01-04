import os
import json
import matplotlib.pyplot as plt


def prep_visualization(data_path, x_attribute, y_attribute):
    # If user chooses summarization file
    if data_path[-4:] == ".txt":
        opened_file = open(data_path, 'r')
        data = opened_file.readline()
        opened_file.close()
        material_data_dict = json.loads(data)

        file = data_path[:-4]
        label = file.split("/")[-1]

        percentages = []
        for percentage in material_data_dict['origin_files']:
            if percentage[1] == ".":
                percentages.append(int(percentage[0]))
            else:
                percentages.append(int(percentage[0:2]))
   
        mod_x_attribute = x_attribute.lower().replace(" ", "_")
        mod_y_attribute = y_attribute.lower().replace(" ", "_")

        x_vals = 0
        y_vals = 0
        if mod_x_attribute == "mix_percentage":
            x_vals = percentages
            y_vals = material_data_dict[mod_y_attribute]
        elif mod_y_attribute == "mix_percentage":
            x_vals = material_data_dict[mod_x_attribute]
            y_vals = percentages
        else:
            x_vals = material_data_dict[mod_x_attribute]
            y_vals = material_data_dict[mod_y_attribute]
        
        plt.plot(x_vals, y_vals, "r^", label=label.split("_")[-1])
        plt.legend(loc="best")
        plt.title(label)
        plt.xlabel(x_attribute)
        plt.ylabel(y_attribute)
        plt.show()

    # If user chooses folder to summarize
    else:
        label = data_path.split("_")
        plot_file_label = f"{label[-1]} mixed with {label[-2][-2:]}"
        all_data = {}
        for file in os.listdir(data_path):
            if file[-4:] == ".txt":
                opened_file = open(data_path + "/" + file, 'r')
                data = opened_file.readline()
                opened_file.close()
                material_data_dict = json.loads(data)
    
                material = file.split("_")
                material = material[-1]
                material = material[:-4]
                all_data[material] = material_data_dict
 
        mod_x_attribute = x_attribute.lower().replace(" ", "_")
        mod_y_attribute = y_attribute.lower().replace(" ", "_")
        plt.clf()

        for element in all_data:
            percentages = []
            for percentage in all_data[element]["origin_files"]:
                if percentage[1] == ".":
                    percentages.append(int(percentage[0]))
                else:
                    percentages.append(int(percentage[0:2]))
            
            x_vals = 0
            y_vals = 0
            if mod_x_attribute == "mix_percentage":
                x_vals = percentages
                y_vals = all_data[element][mod_y_attribute]
            elif mod_y_attribute == "mix_percentage":
                x_vals = all_data[element][mod_x_attribute]
                y_vals = percentages
            else:
                x_vals = all_data[element][mod_x_attribute]
                y_vals = all_data[element][mod_y_attribute]

            plt.scatter(x_vals, y_vals, label=element)

        plt.legend(loc="best")
        plt.title(plot_file_label)
        plt.xlabel(x_attribute)
        plt.ylabel(y_attribute)
        plt.show()
        
        # In case you want to save plots
        #element = data_path.split("_")[-1]
        #path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/Plots/' + element + "_base/"
        #file_name = f"{plot_file_label}_{mod_y_attribute}_{mod_x_attribute}.png" 
        #plt.savefig(path + file_name)

    
if __name__ == "__main__":
    prep_visualization()