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
from configparser import ConfigParser
from os import listdir
from os.path import isfile
import json


def initiate_gui():
    """Create a GUI and defines funtcionality.

    The gui consists of two main windows. The left one is data inputs aswell
    as starting the simulation while the right one is graphic visualisations.
    """
    # Define Tkinter window
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("1200x800")
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


    # Frame 2
    # inits different visualization buttons aswell as calling the 2D-fig.
    #three_dim_vis_button = Button(tabframe1, text="Visualize 3D Atom "
    #
    #                              "(Basic water at the moment)",
    #                              command=visualise_3D)
    #plot_label = Label(tabframe1, text="Plot of data")
    #plot_label.pack(padx=20, pady=20)

    #import_data_button = Button(tabframe2, text="Import Data",
     #                           command=lambda: load_data(gui))
    #import_data_button.pack(padx=20, pady=5)

    ax_canvas = plot_backbone(tabframe1)

    # Frame 1
    # Config settings
    # inits different data fields and buttons.
    potential_label = Label(data_frame, text="Choose potential", width=20)
    potential_label.grid(row=0, column=0)
    potential_list = ["EMT", "LennardJones"]

    value_inside_potential_list = StringVar(gui)
    value_inside_potential_list.set("Select an Option")

    potential_menu = OptionMenu(data_frame, value_inside_potential_list,
                                *potential_list)
    potential_menu.grid(row=0, column=2, padx=15)

    ensamble_label = Label(data_frame, text="Ensemble", width=20)
    ensamble_label.grid(row=1, column=0)

    #ensamble_list = ["NVE (Microcanonical)", "NVT (Canonical)",
    #                 "muVT (Grand canonical)", "NpT (Isothermal-Isobaric)"]
    ensamble_list = ["NVE", "NVT"]

    value_inside_ensemble_list = StringVar(gui)
    value_inside_ensemble_list.set("Select an Option")

    ensemble_menu = OptionMenu(data_frame, value_inside_ensemble_list,
                               *ensamble_list)
    ensemble_menu.grid(row=1, column=2)

    temperature_label = Label(data_frame, text="Temperature (K)", width=20)
    temperature_label.grid(row=2, column=0)

    temperature_entry = Entry(data_frame)
    temperature_entry.grid(row=2, column=2)

    steps_label = Label(data_frame, text="Number of steps", width=20)
    steps_label.grid(row=3, column=0)

    steps_entry = Entry(data_frame)
    steps_entry.grid(row=3, column=2)

    time_steps_label = Label(data_frame, text="Time step", width=20)
    time_steps_label.grid(row=4, column=0)

    time_steps_entry = Entry(data_frame)
    time_steps_entry.grid(row=4, column=2)

    friction_label = Label(data_frame, text="Friction", width=20)
    friction_label.grid(row=5, column=0)

    friction_entry = Entry(data_frame)
    friction_entry.grid(row=5, column=2)

    rec_cohesive_label = Label(data_frame, text="Intervals for recording attributes", width=35)
    rec_cohesive_label.grid(row=6, column=1)

    rec_temp_label = Label(data_frame, text="Temperature", width=20)
    rec_temp_label.grid(row=7, column=0)

    rec_temp_entry = Entry(data_frame)
    rec_temp_entry.grid(row=7, column=2)

    rec_pressure_label = Label(data_frame, text="Pressure", width=20)
    rec_pressure_label.grid(row=8, column=0)

    rec_pressure_entry = Entry(data_frame)
    rec_pressure_entry.grid(row=8, column=2)

    rec_config_label = Label(data_frame, text="Configuration", width=20)
    rec_config_label.grid(row=9, column=0)

    rec_config_entry = Entry(data_frame)
    rec_config_entry.grid(row=9, column=2)

    rec_bulk_label = Label(data_frame, text="Bulk modulus", width=20)
    rec_bulk_label.grid(row=10, column=0)

    rec_bulk_entry = Entry(data_frame)
    rec_bulk_entry.grid(row=10, column=2)

    rec_scaling_label = Label(data_frame, text="Optimal scaling", width=20)
    rec_scaling_label.grid(row=11, column=0)

    rec_scaling_entry = Entry(data_frame)
    rec_scaling_entry.grid(row=11, column=2)

    config_name_label = Label(data_frame, text="Config file name", width=20)
    config_name_label.grid(row=12, column=0)

    config_name_entry = Entry(data_frame)
    config_name_entry.grid(row=12, column=2)

    config_button = Button(data_frame, text='Write to config file',
                           command=lambda: write_to_config(
                               config_name_entry.get(),
                               value_inside_ensemble_list.get(),
                               temperature_entry.get() or "500",
                               value_inside_potential_list.get(),
                               steps_entry.get() or "5000",
                               time_steps_entry.get() or "5",
                               friction_entry.get() or "0.005",
                               "0",
                               rec_temp_entry.get() or "0",
                               rec_pressure_entry.get() or "0",
                               rec_config_entry.get() or "0",
                               rec_bulk_entry.get() or "0",
                               rec_scaling_entry.get() or "0"
                               )
                           )

    config_button.grid(row=13, column=1, pady=10)

    sep_label1 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label1.grid(row=14, column = 0, columnspan = 3)


    # Gather data
    materialID_label = Label(data_frame, text="Material ID", width=20)
    materialID_label.grid(row=15, column=0)

    materialID_entry = Entry(data_frame)
    materialID_entry.grid(row=15, column=2)

    gather_data_button = Button(data_frame, text='Gather material data from Material ID',
                                command=lambda: send_mat_id_to_gather_data(
                                materialID_entry.get()))
    gather_data_button.grid(row=16, column=1, pady=10)

    sep_label2 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label2.grid(row=17, column = 0, columnspan = 3)


    # Simulation
    config_files_label = Label(data_frame, text="Config files", width=20)
    config_files_label.grid(row=18, column=0)

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

    config_files_menu.grid(row=18, column=2)

    traj_file_label = Label(data_frame, text=".traj files", width=20)
    traj_file_label.grid(row=19, column=0)

    traj_list = [file for file in listdir() if isfile(file)]
    value_inside_traj_list = StringVar(gui)
    value_inside_traj_list.set("Select a .traj file")

    input_traj_menu = OptionMenu(data_frame, value_inside_traj_list,
                               *traj_list)

    # Same as for config handler
    def input_traj_handler(event, traj_menu=input_traj_menu, value_inside_traj_list=value_inside_traj_list):   
        return update_input_traj_list(event, traj_menu, value_inside_traj_list)
    input_traj_menu.bind('<Button-1>', input_traj_handler)  

    input_traj_menu.grid(row=19, column=2)

    # Load traj file data
    output_name_label = Label(data_frame, text="Output name", width=20)
    output_name_label.grid(row=20, column=0)

    output_name_entry = Entry(data_frame)
    output_name_entry.grid(row=20, column=2)

    # Start sim button 
    md_sim_button = Button(data_frame, text='Start Simulation',
                           command=lambda: run_single_md_simulation(value_inside_config_files_list.get(), 
                                                                    value_inside_traj_list.get(),
                                                                    output_name_entry.get()))

    md_sim_button.grid(row=21, column=1)

    sep_label3 = Label(data_frame, text="-"*100, bg = "medium aquamarine")
    sep_label3.grid(row=22, column = 0, columnspan = 3)


    # For a future state
    
    # Load output file data
    #output_label = Label(data_frame, text="output files", width=20)
    #output_label.grid(row=13, column=0)

    output_data = [file for file in listdir() if isfile(file)]
    value_inside_output_data = StringVar(gui)
    value_inside_output_data.set("Select an output file")

    output_data_menu = OptionMenu(data_frame, value_inside_output_data,
                               *output_data)

    def output_data_handler(event, output_data_menu=output_data_menu, value_inside_output_data=value_inside_output_data):   
        return update_output_txt_list(event, output_data_menu, value_inside_output_data)
    output_data_menu.bind('<Button-1>', output_data_handler)  

    output_data_menu.grid(row=23, column=2)
    
    plottable_attributes = ["Temperature", "Pressure", "Bulk Modulus",
                            "Cohesive energy", "Optimal Scaling"]
    value_inside_plottable_list = StringVar(gui)
    value_inside_plottable_list.set("Attribute to plot")

    attributes_menu = OptionMenu(data_frame, value_inside_plottable_list,
                                *plottable_attributes)
    
    attributes_menu.grid(row=23, column=0)
    

    # Visualise results button
    plot_button = Button(data_frame, text='Plot data',
                            command=lambda: visualise_2D(value_inside_plottable_list.get(), value_inside_output_data.get(), ax_canvas, text_box, tabframe1 ,False))

    plot_button.grid(row=24, column=1)
    
    Button(data_frame, text="Open plot in new window", command=lambda: visualise_2D(value_inside_plottable_list.get(), value_inside_output_data.get(), ax_canvas, text_box, tabframe1, True)).grid(row=25,column=1,pady=30)

    # Quit
    quit_button = Button(data_frame, text="Exit Program", command=gui.quit)
    quit_button.grid(pady=50)

    Button(tabframe1, text="Clear Graph", command=lambda: clear_canvas(ax_canvas[0], ax_canvas[1])).pack()

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


def write_to_config(config_name, value_inside_ensemble_list, 
                    temperature, value_inside_potential_list, steps, 
                    time_steps, friction, rec_coh_e, rec_temp, 
                    rec_pressure, rec_config, rec_bulk, rec_scaling):
    """Create the configuration file

    Args:
        ensemble(string): Ensemble to use in simulation
        temperature(int): Initial temperature in simulation
        potential(string): Potential to use in simulation
        step_number(int): Number of steps to use in simulation
        time_step(int): Time step in fs to use in simulation
        friction(float): Friction for NVT simulation
        interval(int): Interval for which to calculate properties
        show_properties(bool): Show properties or not

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
    print(time_steps)
    if time_steps.isdigit():
        if (int(time_steps) <= 0):
            messagebox.showerror("Value error", "Invalid time step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid time step here?") 
        return None
    
    # Check if cohesive energy step is valid
    if rec_coh_e.isdigit():
        if (0 < int(rec_coh_e) < int(steps)):
            messagebox.showerror("Value error", "Invalid cohesive energy step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid cohesive energy step") 
        return None
    
    # Check if temperature step is valid
    if rec_temp.isdigit():
        if not (0 < int(rec_temp) < int(steps)):
            messagebox.showerror("Value error", "Invalid temperature step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid temperature step") 
        return None
    
    # Check if pressure step is valid
    if rec_pressure.isdigit():
        if not (0 < int(rec_pressure) < int(steps)):
            messagebox.showerror("Value error", "Invalid pressure step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid pressure step") 
        return None
    
    # Check if configuration step is valid
    if rec_config.isdigit():
        if not (0 < int(rec_config) < int(steps)):
            messagebox.showerror("Value error", "Invalid configuration step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid configuration step") 
        return None
    
    # Check if bulk step is valid
    if rec_bulk.isdigit():
        if not (0 < int(rec_bulk) < int(steps)):
            messagebox.showerror("Value error", "Invalid bulk step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid bulk step") 
        return None
    
    # Check if scaling step is valid
    if rec_scaling.isdigit():
        if not (0 < int(rec_scaling) < int(steps)):
            messagebox.showerror("Value error", "Invalid scaling step") 
            return None
    else:
        messagebox.showerror("Value error", "Invalid scaling step") 
        return None
    
    cfs.config_file(config_name, value_inside_ensemble_list, temperature,
                        value_inside_potential_list, steps, time_steps,
                        friction, 0, rec_temp, rec_pressure, rec_config, rec_bulk,
                        rec_scaling)


def send_mat_id_to_gather_data(materialID):
    """Write user input data to config file.

    Args:
        materialID(string): specifies which material is to be downloaded from database

    Returns:
        None
    """
    try:
        Gather_data.download_data.make_traj_from_material_id(materialID)
    except:
        messagebox.showerror("Invalid id", "Please enter a valid id")


def visualise_2D(attribute_to_plot, file_to_plot, ax_canvas, text_box, frame, boolean):
    """Visualize data in 2D.

    Args:
        At the time of writing not fully clear.

    Returns:
        None
    """
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/' + file_to_plot
    opened_file = open(path, 'r')
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
    try:
        material_data_dict[attribute]
    except KeyError:
        messagebox.showerror("Missing data", "Attribute not recorded")
        return None
    
    average_attribute = round(sum(material_data_dict[attribute])/len(material_data_dict[attribute]))
    max_attribute = round(max(material_data_dict[attribute]),4)
    min_attribute = round(min(material_data_dict[attribute]),4)
    
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
        Min {attribute}: {min_attribute}"""
    
    text_box.config(state="normal")
    text_box.delete('1.0', END)
    text_box.insert("end", message)
    text_box.config(font=("bitstream charter", 15))
    text_box.update()
    text_box.config(state="disabled")
 
    x_values = []
    i = 0
    while i < len(material_data_dict[attribute_to_plot.lower()]):
        x_values.append(i * int(config_data['SimulationSettings']['time_step']) * 
                      int(config_data['RecordingIntervals']['record_' + attribute_to_plot.lower()]))
        i += 1
    
    x_lim = int(config_data['SimulationSettings']['time_step']) * int(config_data['SimulationSettings']['step_number'])

    graph_in_plot = plot(ax_canvas[0], ax_canvas[1], x_values, material_data_dict[attribute_to_plot.lower()], x_lim, attribute_to_plot)
    if boolean:
        open_window(graph_in_plot[0], material_data_dict[attribute_to_plot.lower()], x_lim, attribute_to_plot)
    else:
        plot(ax_canvas[0],ax_canvas[1], x_values, material_data_dict[attribute_to_plot.lower()], x_lim, attribute_to_plot)


def visualise_3D(value_inside_output_data, value_inside_traj_data):
    """Visualize data in 3D.

    Args:
        At the time of writing not fully clear.

    Returns:
        None
    """
    h2o = molecule("H2O")
    view(h2o)


def load_data(gui):
    """Load atom data to plot.

    Args:
        At the time of writing not fully clear.

    Returns:
        None
    """
    gui.filename = filedialog.askopenfilename()
    print(gui.filename)


if __name__ == "__main__":
    main_program = initiate_gui()
    main_program.mainloop()
