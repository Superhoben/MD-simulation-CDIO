"""This module defines the entire GUI.

At the time of writing, all functionality is not implemented
"""
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
from ase.build import molecule
from ase.io.trajectory import Trajectory
from asap3 import EMT
from ase.visualize import view
from pathlib import Path
from User_interface.plot_in_gui import *
import Gather_data.configuration_file_script as cfs
import Gather_data.download_data 
from Simulation.run_md_simulation import run_single_md_simulation
#from User_API_key.start_program import *
from configparser import ConfigParser
from os import listdir
from os.path import isfile
from API_key.start_program import *
import json



def initiate_gui():
    """Create a GUI and defines funtcionality.

    The gui consists of two main windows. The left one is data inputs aswell
    as starting the simulation while the right one is graphic visualisations.
    """
    # Define Tkinter window
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("1200x900")
    gui.resizable(width=False, height=False)

    # 
    # Define the 2 different main sections of the gui.
    data_frame = Frame(gui, background="medium aquamarine")
    data_frame.pack_propagate(False)
    data_frame.grid_propagate(False)
    data_frame.pack(padx=10, pady=10, side=LEFT, expand=True, fill=BOTH)
    data_frame.columnconfigure(0, weight=2)
    data_frame.columnconfigure(1, weight=1)

    plot_frame = Frame(gui, background="light blue")
    plot_frame.pack_propagate(False)
    plot_frame.grid_propagate(False)
    plot_frame.pack(padx=10, pady=10, side=RIGHT, expand=True, fill=BOTH)

    tabs = ttk.Notebook(plot_frame)
    tabs.pack()
    tabframe1 = Frame(tabs, height=800, width=600, bg="SkyBlue1")
    tabframe1.pack_propagate(False)
    tabframe1.grid_propagate(False)
    tabframe1.pack(expand=True, fill=BOTH)

    tabframe2 = Frame(tabs, height=800, width=600, bg="SkyBlue1")
    tabframe2.pack_propagate(False)
    tabframe2.grid_propagate(False)
    tabframe2.pack(expand=True, fill=BOTH)

    tabs.add(tabframe1, text="Visualisation")
    tabs.add(tabframe2, text="Data output")

    text_box_frame = Frame(tabframe2, bg="lemon chiffon")
    text_box_frame.pack_propagate(False)
    text_box_frame.grid_propagate(False)
    text_box_frame.pack(fill="both", expand=True, padx=5, pady=225) 


    text_box = Text(text_box_frame, bg="lemon chiffon")
    text_box.pack(expand=True)
    text_box.config(state="disabled")

    ax_canvas, plot_title = plot_backbone(tabframe1)

    # Frame 1
    # Config settings
    # inits different data fields and buttons.
    rownumber = 0
    potential_label = Label(data_frame, text="Choose potential", width=20)
    potential_label.grid(row=rownumber, column=0)
    potential_list = ["EMT", "LennardJones"]
    value_inside_potential_list = StringVar(gui)
    value_inside_potential_list.set("Select an Option")
    potential_menu = OptionMenu(data_frame, value_inside_potential_list,
                                *potential_list)
    potential_menu.grid(row=rownumber, column=2, padx=15)

    rownumber += 1
    ensamble_label = Label(data_frame, text="Ensemble", width=20)
    ensamble_label.grid(row=rownumber, column=0)
    ensamble_list = ["NVE", "NVT"]
    value_inside_ensemble_list = StringVar(gui)
    value_inside_ensemble_list.set("Select an Option")
    ensemble_menu = OptionMenu(data_frame, value_inside_ensemble_list,
                               *ensamble_list)
    ensemble_menu.grid(row=rownumber, column=2)

    rownumber += 1
    temperature_label = Label(data_frame, text="Temperature (K)", width=20)
    temperature_label.grid(row=rownumber, column=0)
    temperature_entry = Entry(data_frame)
    temperature_entry.grid(row=rownumber, column=2)

    rownumber += 1
    steps_label = Label(data_frame, text="Number of steps", width=20)
    steps_label.grid(row=rownumber, column=0)
    steps_entry = Entry(data_frame)
    steps_entry.grid(row=rownumber, column=2)

    rownumber += 1
    time_steps_label = Label(data_frame, text="Time step (fs)", width=20)
    time_steps_label.grid(row=rownumber, column=0)
    time_steps_entry = Entry(data_frame)
    time_steps_entry.grid(row=rownumber, column=2)

    rownumber += 1
    friction_label = Label(data_frame, text="Friction", width=20)
    friction_label.grid(row=rownumber, column=0)
    friction_entry = Entry(data_frame)
    friction_entry.grid(row=rownumber, column=2)

    rownumber += 1
    rec_cohesive_label = Label(data_frame, text="Intervals for recording attributes", width=35)
    rec_cohesive_label.grid(row=rownumber, column=1)

    rownumber += 1    
    rec_energy_label = Label(data_frame, text="Energy", width=20)
    rec_energy_label.grid(row=rownumber, column=0)
    rec_energy_entry = Entry(data_frame)
    rec_energy_entry.grid(row=rownumber, column=2)

    rownumber += 1    
    rec_temp_label = Label(data_frame, text="Temperature", width=20)
    rec_temp_label.grid(row=rownumber, column=0)
    rec_temp_entry = Entry(data_frame)
    rec_temp_entry.grid(row=rownumber, column=2)

    rownumber += 1    
    rec_pressure_label = Label(data_frame, text="Pressure", width=20)
    rec_pressure_label.grid(row=rownumber, column=0)
    rec_pressure_entry = Entry(data_frame)
    rec_pressure_entry.grid(row=rownumber, column=2)

    rownumber += 1
    rec_config_label = Label(data_frame, text="Configuration", width=20)
    rec_config_label.grid(row=rownumber, column=0)
    rec_config_entry = Entry(data_frame)
    rec_config_entry.grid(row=rownumber, column=2)

    rownumber += 1
    rec_bulk_label = Label(data_frame, text="Bulk modulus", width=20)
    rec_bulk_label.grid(row=rownumber, column=0)
    rec_bulk_entry = Entry(data_frame)
    rec_bulk_entry.grid(row=rownumber, column=2)

    rownumber += 1
    rec_scaling_label = Label(data_frame, text="Optimal scaling", width=20)
    rec_scaling_label.grid(row=rownumber, column=0)
    rec_scaling_entry = Entry(data_frame)
    rec_scaling_entry.grid(row=rownumber, column=2)
    
    rownumber += 1
    rec_elastic_label = Label(data_frame, text="Elastic tensor", width=20)
    rec_elastic_label.grid(row=rownumber, column=0)
    rec_elastic_entry = Entry(data_frame)
    rec_elastic_entry.grid(row=rownumber, column=2)

    rownumber += 1
    config_name_label = Label(data_frame, text="Config file name", width=20)
    config_name_label.grid(row=rownumber, column=0)
    config_name_entry = Entry(data_frame)
    config_name_entry.grid(row=rownumber, column=2)

    rownumber += 1
    config_button = Button(data_frame, text='Write to config file',
                           command=lambda: write_to_config(
                               config_name_entry.get(),
                               value_inside_ensemble_list.get(),
                               temperature_entry.get() or "500",
                               value_inside_potential_list.get(),
                               steps_entry.get() or "5000",
                               time_steps_entry.get() or "5",
                               friction_entry.get() or "0.005",
                               rec_energy_entry.get() or "0",
                               "0",
                               rec_temp_entry.get() or "0",
                               rec_pressure_entry.get() or "0",
                               rec_config_entry.get() or "0",
                               rec_bulk_entry.get() or "0",
                               rec_scaling_entry.get() or "0",
                               rec_elastic_entry.get() or "0"
                               )
                           )
    config_button.grid(row=rownumber, column=1, pady=10)

    rownumber += 1
    sep_label1 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label1.grid(row=rownumber, column = 0, columnspan = 3)


    # Gather data
    rownumber += 1
    materialID_label = Label(data_frame, text="Material ID", width=20)
    materialID_label.grid(row=rownumber, column=0)
    materialID_entry = Entry(data_frame)
    materialID_entry.grid(row=rownumber, column=2)
    
    rownumber += 1
    cell_size_label = Label(data_frame, text="Supercell size", width=20)
    cell_size_label.grid(row=rownumber, column=0)
    cell_size_entry = Entry(data_frame)
    cell_size_entry.grid(row=rownumber, column=2)

    rownumber += 1
    gather_data_button = Button(data_frame, text='Gather material data',
                                command=lambda: send_mat_id_to_gather_data(
                                materialID_entry.get(), cell_size_entry.get()))
    gather_data_button.grid(row=rownumber, column=2, pady=10, padx=10)

    update_api_key = Button(data_frame, text='Update API key ',
                                command=lambda: prompt_for_api_key())
    update_api_key.grid(row=rownumber, column=0, pady=10)
    
    rownumber += 1
    sep_label2 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label2.grid(row=rownumber, column = 0, columnspan = 3)

    # Simulation
    rownumber += 1
    config_files_label = Label(data_frame, text="Config files", width=20)
    config_files_label.grid(row=rownumber, column=0)

    # Creates a list which contains the file names of all the files in a directory
    # In this case the directory is the one we're standing in (since listdir has no input)
    config_files = ["Initializing list"]
    value_inside_config_files_list = StringVar(gui)
    value_inside_config_files_list.set("Select a config file")

    config_files_menu = OptionMenu(data_frame, value_inside_config_files_list,
                                *config_files)

    # Argument trick to allow us to send in more than one argument to the event handler, config_files_menu.bind('<Button-1>', config_handler)
    # For the curious, see: https://anzeljg.github.io/rin2/book2/2405/docs/tkinter/extra-args.html 
    def config_handler(event, config_files_menu=config_files_menu, value_inside_config_files_list=value_inside_config_files_list):   
        return update_input_config_list(event, config_files_menu, value_inside_config_files_list)
    config_files_menu.bind('<Button-1>', config_handler)    

    config_files_menu.grid(row=rownumber, column=2)

    rownumber += 1
    ensamble_label = Label(data_frame, text=".traj files", width=20)
    ensamble_label.grid(row=rownumber, column=0)

    traj_list = [file for file in listdir() if isfile(file)]
    value_inside_traj_list = StringVar(gui)
    value_inside_traj_list.set("Select a .traj file")

    input_traj_menu = OptionMenu(data_frame, value_inside_traj_list,
                               *traj_list)

    # Same as for config handler
    def input_traj_handler(event, traj_menu=input_traj_menu, value_inside_traj_list=value_inside_traj_list):   
        return update_input_traj_list(event, traj_menu, value_inside_traj_list)
    input_traj_menu.bind('<Button-1>', input_traj_handler)  

    input_traj_menu.grid(row=rownumber, column=2)

    rownumber += 1
    # Load traj file data
    output_name_label = Label(data_frame, text="Output name", width=20)
    output_name_label.grid(row=rownumber, column=0)

    output_name_entry = Entry(data_frame)
    output_name_entry.grid(row=rownumber, column=2)

    # Start sim button 
    md_sim_button = Button(data_frame, text='Start Simulation',
                           command=lambda: run_single_md_simulation(value_inside_config_files_list.get(), 
                                                                    value_inside_traj_list.get(),
                                                                    output_name_entry.get()))
    rownumber += 1
    md_sim_button.grid(row=rownumber, column=1)

    rownumber += 1
    sep_label3 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label3.grid(row=rownumber, column = 0, columnspan = 3)

    rownumber += 1
    plottable_attributes = ["Total Energy", "Kinetic Energy", "Potential Energy", 
                            "Temperature", "Pressure", "Bulk Modulus",
                            "Optimal Scaling", "Elastic Tensor"]

    value_inside_plottable_list = StringVar(gui)
    value_inside_plottable_list.set("Attribute to plot")

    attributes_menu = OptionMenu(data_frame, value_inside_plottable_list,
                                *plottable_attributes)
    
    #rownumber += 1
    attributes_menu.grid(row=rownumber, column=0)

    output_data = [file for file in listdir() if isfile(file)]
    value_inside_output_data = StringVar(gui)
    value_inside_output_data.set("Select an output file")

    output_data_menu = OptionMenu(data_frame, value_inside_output_data,
                               *output_data)

    def output_data_handler(event, output_data_menu=output_data_menu, value_inside_output_data=value_inside_output_data):   
        return update_output_txt_list(event, output_data_menu, value_inside_output_data)
    output_data_menu.bind('<Button-1>', output_data_handler)  

    output_data_menu.grid(row=rownumber, column=2)
    


    # Visualise results button
    plot_button = Button(data_frame, text='Plot data',
                            command=lambda: visualise_2D(value_inside_plottable_list.get(), value_inside_output_data.get(), ax_canvas, text_box, tabframe1 ,False,plot_title))

    rownumber += 1
    plot_button.grid(row=rownumber, column=1)
    
    rownumber += 1
    Button(data_frame, text="Open plot in new window", command=lambda: visualise_2D(value_inside_plottable_list.get(), value_inside_output_data.get(), ax_canvas, text_box, tabframe1, True, plot_title="Attribute")).grid(row=rownumber,column=1,pady=10)

    rownumber += 1

    output_traj_files = [file for file in listdir() if isfile(file)]
    value_inside_trajoutput_data = StringVar(gui)
    value_inside_trajoutput_data.set("Select a trajectory file")

    output_traj_menu = OptionMenu(data_frame, value_inside_trajoutput_data,
                               *output_traj_files)

    def output_traj_handler(event, output_traj_menu=output_traj_menu, value_inside_trajoutput_data=value_inside_trajoutput_data):   
        return update_output_traj_list(event, output_traj_menu, value_inside_trajoutput_data)
    output_traj_menu.bind('<Button-1>', output_traj_handler)  

    output_traj_menu.grid(row=rownumber, column=0)
    


    trajbutton = Button(data_frame, text="Look at trajectory", command=lambda: animate_traj(value_inside_trajoutput_data.get()))
    trajbutton.grid(row=rownumber, column=2)
    
    
    # Quit
    quit_button = Button(data_frame, text="Exit Program", command=gui.quit)
    quit_button.grid(pady=50)


    Button(tabframe1, text="Clear Graph", command=lambda: clear_canvas(ax_canvas[0], ax_canvas[1], plot_title)).pack(pady=10)

    return gui


def update_input_config_list(event, config_files_menu, value_inside_config_files_list):
    """Updates the config selection dropdown menu.

    Args:
        event(tkinter.Event): Handles input events, in our case mouse left click
        config_files_menu(tkinter.OptionMenu): Dropdown menu which will be updated
        value_inside_config_files_list(tkinter.StringVar): A variable which is set
            when the user selects an item in the dropdown menu. The variable is a 
            string which specifies the filename of the selected item

    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    all_files = [file for file in listdir(path) if isfile(path+file)]
    config_files = [file for file in all_files if file[-4:] == ".ini"]
    menu = config_files_menu["menu"]
    menu.delete(0, "end")

    for config_file in config_files:
        menu.add_command(label=config_file, command=lambda value=config_file: value_inside_config_files_list.set(value))


def update_input_traj_list(event, traj_menu, value_inside_traj_list):
    """Updates the trajectory selection dropdown menu.

    Args:
        event(tkinter.Event): Handles input events, in our case mouse left click
        traj_menu(tkinter.OptionMenu): Dropdown menu which will be updated
        value_inside_traj_list(tkinter.StringVar): A variable which is set
            when the user selects an item in the dropdown menu. The variable is a 
            string which specifies the filename of the selected item

    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    all_files = [file for file in listdir(path) if isfile(path+file)]
    traj_files = [file for file in all_files if file[-5:] == ".traj"]
    menu = traj_menu["menu"]
    menu.delete(0, "end")

    for traj_file in traj_files:
        menu.add_command(label=traj_file, command=lambda value=traj_file: value_inside_traj_list.set(value))    


def update_output_txt_list(event, output_data_menu, value_inside_output_data):
    """Updates the output selection dropdown menu.

    Args:
        event(tkinter.Event): Handles input events, in our case mouse left click
        output_data_menu(tkinter.OptionMenu): Dropdown menu which will be updated
        value_inside_output_data(tkinter.StringVar): A variable which is set
            when the user selects an item in the dropdown menu. The variable is a 
            string which specifies the filename of the selected item
ax.set_xlabel("Time [femto seconds]")
    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/'
    print(path)
    all_files = [file for file in listdir(path) if isfile(path+file)]
    output_files = [file for file in all_files if file[-4:] == ".txt"]
    menu = output_data_menu["menu"]
    menu.delete(0, "end")

    for output_file in output_files:
        menu.add_command(label=output_file, command=lambda value=output_file: value_inside_output_data.set(value))    


def update_output_traj_list(event, traj_menu, value_inside_traj_list):
    """Updates the trajectory output selection dropdown menu.

    Args:
        event(tkinter.Event): Handles input events, in our case mouse left click
        traj_menu(tkinter.OptionMenu): Dropdown menu which will be updated
        value_inside_traj_list(tkinter.StringVar): A variable which is set
            when the user selects an item in the dropdown menu. The variable is a 
            string which specifies the filename of the selected item

    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/'
    all_files = [file for file in listdir(path) if isfile(path+file)]
    traj_files = [file for file in all_files if file[-5:] == ".traj"]
    menu = traj_menu["menu"]
    menu.delete(0, "end")

    for traj_file in traj_files:
        menu.add_command(label=traj_file, command=lambda value=traj_file: value_inside_traj_list.set(value))


def write_to_config(file_name='default_config', value_inside_ensemble_list='NVE', temperature=500, value_inside_potential_list='EMT',
                steps=5000, time_steps=5, friction=0.005, rec_energy = 0,
                rec_coh_e = 0, rec_temp = 0, rec_pressure = 0, 
                rec_config = 0, rec_bulk = 0, rec_scaling = 0,
                record_elastic = 0):
    """Create the configuration file

    Args:
    	file_name(string): Name of the file to create
        ensemble(string): Ensemble to use in simulation
        temperature(int): Initial temperature in simulation
        potential(string): Potential to use in simulation
        step_number(int): Number of steps to use in simulation
        time_step(int): Time step in fs to use in simulation
        friction(float): Friction for NVT simulation
        record_energy(int): Interval to record energy in simulation
        record_cohesive_energy(int): Interval to record cohesive energy in simulation
        record_temperature(int): Interval to record temperature in simulation
        record_pressure(int): Interval to record pressure in simulation
        record_configuration(int): Interval to record configuration in simulation
        record_bulk_modulus(int): Interval to record bulk modulus in simulation
        record_optimal_scaling(int): Interval to record optimal scaling in simulation
        record_elastic(int): Interval to record elastic properties in simulation

    Returns:
        None
    """
    
    # Check if temperature is valid
    if temperature.isdigit():
        if int(temperature) <= 0:
            messagebox.showerror("Value error", "Invalid temperature") 
            return None
        elif (int(temperature) > 10000):
            if not messagebox.askokcancel("Value doubt", "This value might result in an unstable simulation, you sure you'd like to continue?"):
                return None
    else:
        messagebox.showerror("Value error", "Invalid temperature")
        return None
   
    # Check if number of steps is valid
    if steps.isdigit():
        if (int(steps) <= 0):
            messagebox.showerror("Value error", "Invalid number of steps") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid number of steps")
        return None
    
    # Check if time step is valid
    if time_steps.isdigit():
        if (int(time_steps) <= 0):
            messagebox.showerror("Value error", "Invalid time step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid time step") 
        return None
    
    # Check if cohesive energy step is valid
    if rec_coh_e.isdigit():
        if (0 < int(rec_coh_e) < int(steps)):
            messagebox.showerror("Value error", "Invalid cohesive energy step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid cohesive energy step") 
        return None
    
    # Check if energy step is valid
    if rec_energy.isdigit():
        if int(rec_energy) < 0 or int(rec_energy) > int(steps):
            messagebox.showerror("Value error", "Invalid energy step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid energy step") 
        return None
    
    # Check if temperature step is valid
    if rec_temp.isdigit():
        if int(rec_temp) < 0 or int(rec_temp) > int(steps):
            messagebox.showerror("Value error", "Invalid temperature step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid temperature step") 
        return None
    
    # Check if pressure step is valid
    if rec_pressure.isdigit():
        if int(rec_pressure) < 0 or int(rec_pressure) > int(steps):
            messagebox.showerror("Value error", "Invalid pressure step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid pressure step") 
        return None
    
    # Check if configuration step is valid
    if rec_config.isdigit():
        if int(rec_config) < 0 or int(rec_config) > int(steps):
            messagebox.showerror("Value error", "Invalid configuration step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid configuration step") 
        return None
    
    # Check if bulk step is valid
    if rec_bulk.isdigit():
        if int(rec_bulk) < 0 or int(rec_bulk) > int(steps):
            messagebox.showerror("Value error", "Invalid bulk step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid bulk step") 
        return None
    
    # Check if scaling step is valid
    if rec_scaling.isdigit():
        if int(rec_scaling) < 0 or int(rec_scaling) > int(steps):
            messagebox.showerror("Value error", "Invalid scaling step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid scaling step") 
        return None
    
    cfs.config_file(file_name, value_inside_ensemble_list, temperature, value_inside_potential_list,
                steps, time_steps, friction, rec_energy,
                rec_coh_e, rec_temp, rec_pressure, 
                rec_config, rec_bulk, rec_scaling,
                record_elastic)



def send_mat_id_to_gather_data(materialID, cell_size):
    """Write user input data to config file.

    Args:
        materialID(string): specifies which material is to be downloaded from database


    Returns:
        None
    """
    
    #messagebox.showinfo("Information", "Please Check the terminal")

    api_key = get_api_key()
    if api_key == "API key doesnt exist":
        api_key = prompt_for_api_key()
    
    try:
        Gather_data.download_data.make_traj_from_material_id(materialID, api_key, int(cell_size))
    except:
        messagebox.showerror("Could not download data", "Please check if the material ID is correct or if the API-key is correct.")


def visualise_2D(attribute_to_plot, file_to_plot, ax_canvas, text_box, frame, boolean, plot_title):
    """Visualize data in 2D.

    Args:
        At the time of writing not fully clear.

    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/' + file_to_plot
    if attribute_to_plot == "Attribute to plot":
        messagebox.showerror("Missing attribute", "Choose attribute")

    try:
        opened_file = open(path, 'r')
    except FileNotFoundError:
        messagebox.showerror("Missing text file", "Choose text file")
        return None
    data = opened_file.readline()
    opened_file.close()
    material_data_dict = json.loads(data)

    config_file = material_data_dict['config_file']
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/' + config_file
    config_data = ConfigParser()
    config_data.read(path)

    config_ensemble = config_data['SimulationSettings']['ensemble']
    config_temperature = config_data['SimulationSettings']['temperature']
    config_potential = config_data['SimulationSettings']['potential']
    config_step_number = config_data['SimulationSettings']['step_number']
    config_time_step = config_data['SimulationSettings']['time_step']
    config_friction = config_data['SimulationSettings']['friction']

    attribute = attribute_to_plot.lower()

    if attribute == "bulk modulus":
        attribute = "bulk_modulus"
    elif attribute == "cohesive energy":
        attribute = "cohesive_energy"
    elif attribute == "optimal scaling":
        attribute = "optimal_scaling"
    elif attribute == "elastic tensor":
        attribute = "elastic_tensor"
    elif attribute == "total energy":
        attribute = "total_energy"
    elif attribute == "kinetic energy":
        attribute = "kinetic_energy"
    elif attribute == "potential energy":
        attribute = "potential_energy"

    try:
        material_data_dict[attribute]
    except KeyError:
        messagebox.showerror("Missing data", "Attribute not recorded")
        return None
    
    average_attribute = round(sum(material_data_dict[attribute])/len(material_data_dict[attribute]))
    max_attribute = round(max(material_data_dict[attribute]),4)
    min_attribute = round(min(material_data_dict[attribute]),4)
    data_points = len(material_data_dict[attribute])
    
    message = f"""
    Simulation inputs:
        Config file: {config_file}
        Ensemble: {config_ensemble}
        Temperature: {config_temperature}
        Potential: {config_potential}
        Step number: {config_step_number}
        Time step: {config_time_step}
        Friction: {config_friction}

    Simulation results:
        Average {attribute}: {average_attribute}
        Max {attribute}: {max_attribute}
        Min {attribute}: {min_attribute}
        Data points: {data_points}"""
    
    text_box.config(state="normal")
    text_box.delete('1.0', END)
    text_box.insert("end", message)
    text_box.config(font=("bitstream charter", 15))
    text_box.update()
    text_box.config(state="disabled")
    
    record_attribute = attribute
    if attribute in ["total_energy", "kinetic_energy", "potential_energy"]:
        record_attribute = "energy"
    elif attribute == "elastic_tensor":
        record_attribute = "elastic"
 
    x_values = []
    i = 0
    while i < len(material_data_dict[attribute]):
        x_values.append(i * int(config_data['SimulationSettings']['time_step']) * 
                      int(config_data['RecordingIntervals']['record_' + record_attribute]))
        i += 1
    
    x_lim = int(config_data['SimulationSettings']['time_step']) * int(config_data['SimulationSettings']['step_number'])

    if boolean:
        open_window(x_values, material_data_dict[attribute], x_lim, attribute_to_plot)
    else:
        plot_title.config(text = "Plotted Attribute: \n" + attribute_to_plot)
        plot(ax_canvas[0],ax_canvas[1], x_values, material_data_dict[attribute], x_lim, attribute_to_plot)


def animate_traj(traj_file):
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/'
    os.system("ase gui" + " " + path + traj_file)

if __name__ == "__main__":
    main_program = initiate_gui()
    main_program.mainloop()
