from tkinter import *
from PIL import ImageTk, Image
from ase.build import molecule
from ase.visualize import view
import test

    

def initiate_gui():
#Creates gui
    gui = Tk()
    gui.title("Molecular Dynamics simulations")
    gui.geometry("800x500")
    start_sim_button = Button(gui, text="Start Simulation", command = test.testa)
    start_sim_button.grid(row=0)
    img = ImageTk.PhotoImage(Image.open("<filnamn>"))
    panel = Label(gui, image=img)
    #panel.pack(side = "bottom", fill="both", expand="yes")
    panel.grid(row=0, column=3)
    gui.mainloop()

def start_simulation():
    print("BANGER")


if __name__ == "__main__":
    initiate_gui()