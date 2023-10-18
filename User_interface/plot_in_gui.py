"""This file makes a 2D plot with functionality in the main gui."""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import *
import numpy as np
import matplotlib.figure


def plot(ax, canvas):
    """Plot values in figure.
    
    Args:
        ax and canvas: .

    Returns:
        None
    """
    ax.clear()
    x = np.random.randint(0, 10, 10)
    y = np.random.randint(0, 10, 10)
    ax.plot(x, y)
    canvas.draw()


def plot_backbone(frame):
    """Define the figure and plot button in the right frame of the gui.

    Args:
        frame: Which side of the GUI it's supposed to be fixated on.

    Returns:
        None
    """
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot()

    label = Label(frame, text="TITLE OF PLOT (Random Values at the moment is plotted)")
    label.config(font=("Courier", 12))
    label.pack()

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.get_tk_widget().pack()

    toolbar = NavigationToolbar2Tk(canvas, frame)
    toolbar.update()
    toolbar.pack()
    Button(frame, text="Plot Graph", command=lambda: plot(ax, canvas)).pack(pady=10)
