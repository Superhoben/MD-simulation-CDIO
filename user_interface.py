"""This module defines the entire GUI.

At the time of writing, all functionality is not implemented
"""

from plot_in_gui import *
from tkinter import *
from tkinter import filedialog
from ase.build import molecule
from ase.visualize import view
from run_md_simulation import run_md_simulation
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')



def initiate_gui():
    """Create a GUI and defines funtcionality.

    The gui consists of two main windows. The left one is data inputs aswell
    as starting the simulation while the right one is graphic visualisations.
    """
    # Define Tkinter window
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("1200x800")

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

    # inits different visualization buttons aswell as calling the 2D-fig.
    three_dim_vis_button = Button(plot_frame, text="Visualize 3D Atom "
                                  "(Basic water at the moment)", command=visualise_3D)
    three_dim_vis_button.pack(padx=20, pady=20)

    import_data_button = Button(plot_frame, text="Import Data",
                                command=lambda: load_data(gui))
    import_data_button.pack(padx=20, pady=5)

    plot_backbone(plot_frame)

    # inits different data fields and buttons.
    potential_label = Label(data_frame, text="Choose potential", width=20)
    potential_label.grid(row=0, column=0)

    potential_list = ["Lennard Jones", "EMT", "List goes on"]

    value_inside_potential_list = StringVar(gui)
    value_inside_potential_list.set("Select an Option")

    potential_menu = OptionMenu(data_frame, value_inside_potential_list, *potential_list)
    potential_menu.grid(row=0, column=1)

    ensamble_label = Label(data_frame, text="Ensamble", width=20)
    ensamble_label.grid(row=1, column=0)

    ensamble_list = ["NVE (Microcanonical)", "NVT (Canonical)",
                     "muVT (Grand canonical)", "NpT (Isothermal-Isobaric)"]

    value_inside_ensemble_list = StringVar(gui)
    value_inside_ensemble_list.set("Select an Option")

    ensemble_menu = OptionMenu(data_frame, value_inside_ensemble_list, *ensamble_list)
    ensemble_menu.grid(row=1, column=1)

    materialID_label = Label(data_frame, text="Material ID", width=20)
    materialID_label.grid(row=2, column=0)

    materialID_entry = Entry(data_frame)
    materialID_entry.grid(row=2, column=1)

    temperature_label = Label(data_frame, text="Temperature (K)", width=20)
    temperature_label.grid(row=3, column=0)

    temperature_entry = Entry(data_frame)
    temperature_entry.grid(row=3, column=1)

    steps_label = Label(data_frame, text="Number of steps", width=20)
    steps_label.grid(row=4, column=0)

    steps_entry = Entry(data_frame)
    steps_entry.grid(row=4, column=1)

    config_button = Button(data_frame, text='Write to config file',
                           command=lambda: write_to_config(value_inside_potential_list.get(),
                           value_inside_ensemble_list.get(), materialID_entry,
                           temperature_entry, steps_entry))
    config_button.grid(row=5, column=0, pady=10)


    md_sim_button = Button(data_frame, text='Start Simulation',
                            command=lambda: run_md_simulation("filepath to config ini"))
    md_sim_button.grid(row=5, column=1)

    quit_button = Button(data_frame, text="Exit Program", command=gui.quit)
    quit_button.grid(pady=500)

    return gui

    

def write_to_config(potential, ensamble, materialID, temperature, steps):
    """Write user input data to config file.

    Args:
        At the time of writing not fully clear.

    Returns:
        None
    """
    print("Selected Option: {} and {} and {} and {} and {}".format(potential, ensamble, materialID.get(),
                                                     temperature.get(), steps.get()))


def visualise_3D():
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
    
