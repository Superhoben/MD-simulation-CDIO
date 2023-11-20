"""This file makes a 2D plot with functionality in the main gui."""
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import *
import numpy as np
import matplotlib.figure
import matplotlib.pyplot as plt

def plot_backbone(frame):
    """Define the figure and plot button in the right frame of the gui.

    Args:
        frame: Which side of the GUI it's supposed to be fixated on.

    Returns:
        None
    """
    fig = plt.figure()
    ax = fig.add_subplot()

    label = Label(frame, text="Data Analysis")
    label.config(font=("Courier", 12))
    label.pack(pady=40)
    
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.get_tk_widget().pack()

    toolbar = NavigationToolbar2Tk(canvas, frame)
    toolbar.update()
    toolbar.pack()
    ax.set_xlabel("Time [femto seconds]")
    
    canvas.draw()
    

    return [ax, canvas]


"""This file makes a 2D plot with functionality in the main gui."""
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import *
import numpy as np
import matplotlib.figure
import matplotlib.pyplot as plt

def plot_backbone(frame):
    """Define the figure and plot button in the right frame of the gui.

    Args:
        frame: Which side of the GUI it's supposed to be fixated on.

    Returns:
        None
    """
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot()

    label = Label(frame, text="Attribute to be plotted")
    label.config(font=("Courier", 12))
    label.pack(pady=40)

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.get_tk_widget().pack()

    toolbar = NavigationToolbar2Tk(canvas, frame)
    toolbar.update()
    toolbar.pack()
    ax.set_xlabel("Time [femto seconds]")
    
    canvas.draw()
    
    return [ax, canvas], label


def plot(ax, canvas, x, y, x_lim, y_title):
    """Plot values in figure.
    
    Args:
        ax and canvas: .

    Returns:
        None
    """
    ax.clear()
    ax.plot(x, y, marker="o")
    ax.set_xlim(0, x_lim)
    if y_title == "Temperature":
        ax.set_ylabel(y_title + " [" + "K" "]")
    elif y_title == "Pressure":
        ax.set_ylabel(y_title + " [" + "Pa" "]")
    ax.set_xlabel("Time [femto seconds]")
    canvas.draw()


def clear_canvas(ax, canvas, title):
    title.config(text="Attribute to be plotted")
    ax.clear()
    canvas.draw()
    ax.set_xlabel("Time [femto seconds]")

def open_window(x, y, x_lim ,y_title):

    plt.figure(y_title)
    plt.plot(x,y, marker = "o")
    plt.xlim(0, x_lim)
    plt.xlabel("Time [femto seconds]")
    if y_title == "Temperature":
        plt.ylabel(y_title + " [" + "K" "]")
    elif y_title == "Pressure":
        plt.ylabel(y_title + " [" + "Pa" "]")
    plt.show()
