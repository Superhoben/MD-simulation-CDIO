"""This module defines the entire GUI.

At the time of writing, all functionality is not implemented
"""
#from plot_in_gui import * 
from tkinter import *
from tkinter import ttk
#from PIL import ImageTk, Image
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from functools import partial


def initiate_gui():
    """Create a GUI and defines funtcionality."""
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("1200x600")
    width = 30
    data_frame = Frame(gui, background="red", width = 400, height = 800)
    data_frame.pack_propagate(0)
    data_frame.pack(side = LEFT, fill=BOTH, expand = True)
    data_frame.columnconfigure(0,weight=2)
    data_frame.columnconfigure(1,weight=1)
    
    plot_frame = Frame(gui,background="blue", width = 800, height = 800, name="plot frame")
    plot_frame.pack_propagate(0)
    plot_frame.pack(side = RIGHT, fill=BOTH, expand = True)
    
    
    ensamble_label = Label(data_frame, text = "Ensamble", width=width)
    ensamble_label.grid(row=0, column=0)
    
    global variable
    variable = StringVar(data_frame, "Choose ensamble")
    variable.trace('w', chosen_ensamble)
    
    lst = ["NVE (Microcanonical)", "NVT (Canonical)", "mjuVT (Grand canonical)"]
    listbox = ttk.Combobox(data_frame, values=lst, width=width, textvariable=variable)
    listbox.grid(row=0,column=1)
    
    
    lattice_constant_label = Label(data_frame, text = "Lattice constant", width=width)
    lattice_constant_label.grid(row=1,column=0)
    global lattice_constant_enter
    lattice_constant_enter = Entry(data_frame, width=width, state="normal")
    lattice_constant_enter.grid(row=1,column=1)
    
    temperature_label = Label(data_frame, text = "Temperature", width=width)
    temperature_label.grid(row=2,column=0)
    
    temperature_enter = Entry(data_frame, width=width)
    temperature_enter.grid(row=2,column=1)
    
    start_sim_button = Button(data_frame, text="Start simulation", command=lambda: get_sim_values(listbox, lattice_constant_enter, temperature_enter))
    start_sim_button.grid(row=4,column=0)
    
    
    
    plot_button = Button(plot_frame, text="Plot", command=partial(display_plot, plot_frame))
    plot_button.pack()
    return gui


def start_simulation():
    """Start the simulation by calling another function."""
    print("Test Simulation_button")
    
    
def chosen_ensamble(*args):
    if (variable.get() == "NVE (Microcanonical)"):
        lattice_constant_enter.config(state="disabled")
    elif (variable.get() == "NVT (Canonical)"):
        lattice_constant_enter.config(state="normal")
    
def get_sim_values(ensamble, lattice, temperature):
    print("Ensamble: " + ensamble.get() + ", Lattice: " + lattice.get() + ", Temperature: " + temperature.get())
    
def display_plot(frame):
    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)
    fig = Figure(figsize=(5,4), dpi=100)
    plot = fig.add_subplot(1,1,1)
    plt.plot(x,y)
    
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    
    print(frame.winfo_children())
    
def on_closing():
    frame2

if __name__ == "__main__":
    gui = initiate_gui()
    plot_frame = gui.nametowidget("plot frame")
    plot_frame.configure(bg="green")
    gui.mainloop()
    

    