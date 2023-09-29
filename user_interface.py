"""This module defines the entire GUI.

At the time of writing, all functionality is not implemented
"""

from tkinter import *
from PIL import ImageTk, Image


def initiate_gui():
    """Create a GUI and defines funtcionality."""
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("800x500")

    start_sim_button = Button(gui, text="Start Simulation",
                              command=start_simulation)
    start_sim_button.grid(row=0)

    gui.mainloop()


def start_simulation():
    """Start the simulation by calling another function."""
    print("Test Simulation_button")


if __name__ == "__main__":
    initiate_gui()
