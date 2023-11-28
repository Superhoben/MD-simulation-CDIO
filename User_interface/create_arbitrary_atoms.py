from tkinter import *
from tkinter import ttk
from ase import Atoms
from ase.visualize import view
from ase.lattice.triclinic import Triclinic
from ase.lattice.monoclinic import SimpleMonoclinic, BaseCenteredMonoclinic
from ase.lattice.orthorhombic import SimpleOrthorhombic, BaseCenteredOrthorhombic, BodyCenteredOrthorhombic, FaceCenteredOrthorhombic
from ase.lattice.tetragonal import SimpleTetragonal, CenteredTetragonal
from ase.lattice.cubic import SimpleCubic, BodyCenteredCubic, FaceCenteredCubic


def inbetweener(system, window, **kwargs):
    window.destroy()
    print(system)
    create_view_and_save_crystal_guided(system, **kwargs)


def create_view_and_save_crystal_guided(system,a=0,b=0,c=0,alpha=0,beta=0,gamma=0,element=0,elements=0,lattice_constant=0, lattice_constants=0):
    primitive_cell = 0
    if system == "triclinic":
        primitive_cell = Triclinic(symbol=element, latticeconstant={'a':a,'b':b,'c':c,'alpha':alpha,'beta':beta,'gamma':gamma})
        
    elif system == "monoclinic":
        primitive_cell = SimpleMonoclinic(symbol=element, latticeconstant={'a':a,'b':b,'c':c,'alpha':90,'beta':beta,'gamma':gamma})
        
    elif system == "base_centered_monoclinic":
        elements = elements.split()
        primitive_cell = BaseCenteredMonoclinic(symbol=elements[0], latticeconstant={'a':a,'b':b,'c':c,'alpha':90,'beta':beta,'gamma':gamma})
        primitive_cell.set_chemical_symbols(elements)
    
    elif system == "simple_orthonormic":
        primitive_cell = SimpleOrthorhombic(symbol=element, latticeconstant={'a':a,'b':b,'c':c})
    
    elif system == "base_centered_orthonormic":
        elements = elements.split()
        primitive_cell = BaseCenteredOrthorhombic(symbol=elements[0], latticeconstant={'a':a,'b':b,'c':c})
        primitive_cell.set_chemical_symbols(elements)
        
    elif system == "body_centered_orthonormic":
        elements = elements.split()
        primitive_cell = BodyCenteredOrthorhombic(symbol=elements[0], latticeconstant={'a':a,'b':b,'c':c})
        primitive_cell.set_chemical_symbols(elements)

    elif system == "face_centered_orthonormic":
        primitive_cell = FaceCenteredOrthorhombic(symbol=element, latticeconstant={'a':a,'b':b,'c':c}) 
    
    elif system == "simple_tetragonal":
        primitive_cell = SimpleTetragonal(symbol=element, latticeconstant=lattice_constants)
        
    elif system == "centered_tetragonal":
        elements = elements.split()
        primitive_cell = CenteredTetragonal(symbol=elements[0], latticeconstant=lattice_constants)
        primitive_cell.set_chemical_symbols(elements)
        
    elif system == "simple_cubic":
        primitive_cell = SimpleCubic(symbol=element, latticeconstant=lattice_constant)
    
    elif system == "body_centered_cubic":
        elements = elements.split()
        primitive_cell = BodyCenteredCubic(symbol=elements[0], latticeconstant=lattice_constant)
        primitive_cell.set_chemical_symbols(elements)

    elif system == "face_centered_cubic":
        elements = elements.split()
        primitive_cell = FaceCenteredCubic(symbol=elements[0], latticeconstant=lattice_constant)
        primitive_cell.set_chemical_symbols(elements)

    number_atoms_primitive = len(primitive_cell)
    view(primitive_cell, block = False)


def create_atom():
    new_window = Tk()
    new_window.title("Arbitrary atom configuration")
    new_window.geometry("700x450")
    
    tabs = ttk.Notebook(new_window)
    tabs.pack()
    
    tabframe0 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe0.pack_propagate(False)
    tabframe0.grid_propagate(False)
    tabframe0.pack(expand=True, fill=BOTH)
    
    tabframe1 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe1.pack_propagate(False)
    tabframe1.grid_propagate(False)
    tabframe1.pack(expand=True, fill=BOTH)

    tabframe2 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe2.pack_propagate(False)
    tabframe2.grid_propagate(False)
    tabframe2.pack(expand=True, fill=BOTH)
    
    tabframe3 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe3.pack_propagate(False)
    tabframe3.grid_propagate(False)
    tabframe3.pack(expand=True, fill=BOTH)
    
    tabframe4 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe4.pack_propagate(False)
    tabframe4.grid_propagate(False)
    tabframe4.pack(expand=True, fill=BOTH)
    
    tabframe5 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe5.pack_propagate(False)
    tabframe5.grid_propagate(False)
    tabframe5.pack(expand=True, fill=BOTH)
    
    tabframe6 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe6.pack_propagate(False)
    tabframe6.grid_propagate(False)
    tabframe6.pack(expand=True, fill=BOTH)
    
    tabframe7 = Frame(tabs, height=400, width=600, bg="SkyBlue1")
    tabframe7.pack_propagate(False)
    tabframe7.grid_propagate(False)
    tabframe7.pack(expand=True, fill=BOTH)

    tabs.add(tabframe0, text="Homepage??")
    tabs.add(tabframe1, text="Triclinic")
    tabs.add(tabframe2, text="Monoclinic")
    tabs.add(tabframe3, text="Orthorhombic")
    tabs.add(tabframe4, text="Tetragonal")
    tabs.add(tabframe5, text="Trigonal")
    tabs.add(tabframe6, text="Hexagonal")
    tabs.add(tabframe7, text="Cubic")
    
    text_label = Label(tabframe0, text="Choose crystal system. System names and description of the lattice vectors follows below. \n" +
          "1 - Triclinic, arbitrary lengths and directions \n"+
          "2 - Monoclinic, arbitrary lengths, two vectors are orthagonal \n"+
          "3 - Orthorhombic, arbitrary lengths, all vectors are orthagonal \n"+
          "4 - Tetragonal, two equal lenghts, all vectors are orthagonol \n"+
          "5 - Trigonal, not yet implemented \n"+
          "6 - Hexagonal, not yet implemented \n"+
          "7 - Cubic, equal lenghts and all vectors are orthagonal")
    text_label.grid(row=0, column=1)

    
    label1 = Label(tabframe1, text="Input length of the first lattice translation vector in Å:")
    label1.grid(row=1,column=0)
    label2 = Label(tabframe1, text="Input length of the second lattice translation vector in Å:")
    label2.grid(row=2,column=0)
    label3 = Label(tabframe1, text="Input length of the third lattice translation vector in Å:")
    label3.grid(row=3,column=0)
    label4 = Label(tabframe1, text="Input angle between the first and second translation vector in degrees:")
    label4.grid(row=4,column=0)
    label5 = Label(tabframe1, text="Input angle between the first and third translation vector in degrees:")
    label5.grid(row=5,column=0)
    label6 = Label(tabframe1, text="Input angle between the second and third translation vector in degrees:")
    label6.grid(row=6,column=0)
    label7 = Label(tabframe1, text="Input the chemical symbol for the base atom:")
    label7.grid(row=7,column=0)
    entry1 = Entry(tabframe1)
    entry1.grid(row=1,column=1)
    entry2 = Entry(tabframe1)
    entry2.grid(row=2,column=1)
    entry3 = Entry(tabframe1)
    entry3.grid(row=3,column=1)
    entry4 = Entry(tabframe1)
    entry4.grid(row=4,column=1)
    entry5 = Entry(tabframe1)
    entry5.grid(row=5,column=1)
    entry6 = Entry(tabframe1)
    entry6.grid(row=6,column=1)
    entry7 = Entry(tabframe1)
    entry7.grid(row=7,column=1)
    
    triclinic_button = Button(tabframe1, text="Continue", command = lambda: inbetweener("triclinic", new_window,
                                                                                        a=float(entry1.get()),
                                                                                        b=float(entry2.get()),
                                                                                        c=float(entry3.get()),
                                                                                        alpha=float(entry4.get()),
                                                                                        beta=float(entry5.get()),
                                                                                        gamma=float(entry6.get()),
                                                                                        element=entry7.get()))
    triclinic_button.grid(row=8,column=0)
    
    
    monoclinic_label1 = Label(tabframe2, text="1 - Simple monoclinic (also known as primitive monoclinic)")
    monoclinic_label1.grid(row=0,column=0)
    monoclinic_label2 = Label(tabframe2, text="Input length of the first lattice translation vector in Å:")
    monoclinic_label2.grid(row=1,column=0)
    monoclinic_label3 = Label(tabframe2, text="Input length of the second lattice translation vector in Å:")
    monoclinic_label3.grid(row=2,column=0)
    monoclinic_label4 = Label(tabframe2, text="Input length of the third lattice translation vector in Å:")
    monoclinic_label4.grid(row=3,column=0)
    monoclinic_label5 = Label(tabframe2, text="First and second translation vectors are orthogonal")
    monoclinic_label5.grid(row=4,column=0)
    monoclinic_label6 = Label(tabframe2, text="Input angle between the first and third translation vector:")
    monoclinic_label6.grid(row=5,column=0)
    monoclinic_label7 = Label(tabframe2, text="Input angle between the second and third translation vector:")
    monoclinic_label7.grid(row=6,column=0)
    monoclinic_label8 = Label(tabframe2, text="Input one element in chemical notation:")
    monoclinic_label8.grid(row=7,column=0)
    monoclinic_entry1 = Entry(tabframe2)
    monoclinic_entry1.grid(row=1,column=1)
    monoclinic_entry2 = Entry(tabframe2)
    monoclinic_entry2.grid(row=2,column=1)
    monoclinic_entry3 = Entry(tabframe2)
    monoclinic_entry3.grid(row=3,column=1)
    monoclinic_entry4 = Entry(tabframe2)
    monoclinic_entry4.grid(row=5,column=1)
    monoclinic_entry5 = Entry(tabframe2)
    monoclinic_entry5.grid(row=6,column=1)
    monoclinic_entry6 = Entry(tabframe2)
    monoclinic_entry6.grid(row=7,column=1)
    monoclinic_button1 = Button(tabframe2, text="Continue", command=lambda:inbetweener("monoclinic", new_window,
                                                                                       a = float(monoclinic_entry1.get()),
                                                                                       b = float(monoclinic_entry2.get()),
                                                                                       c = float(monoclinic_entry3.get()),
                                                                                       beta = float(monoclinic_entry4.get()),
                                                                                       gamma = float(monoclinic_entry5.get()),
                                                                                       element = monoclinic_entry6.get()))
    monoclinic_button1.grid(row=8,column=1)
    
    centered_monoclinic_label1 = Label(tabframe2, text="2 - Base centered monoclinic")
    centered_monoclinic_label1.grid(row=9,column=0)
    centered_monoclinic_label2 = Label(tabframe2, text="Input length of the first lattice translation vector in Å:")
    centered_monoclinic_label2.grid(row=10,column=0)
    centered_monoclinic_label3 = Label(tabframe2, text="Input length of the second lattice translation vector in Å:")
    centered_monoclinic_label3.grid(row=11,column=0)
    centered_monoclinic_label4 = Label(tabframe2, text="Input length of the third lattice translation vector in Å:")
    centered_monoclinic_label4.grid(row=12,column=0)
    centered_monoclinic_label5 = Label(tabframe2, text="First and second translation vectors are orthogonal")
    centered_monoclinic_label5.grid(row=13,column=0)
    centered_monoclinic_label6 = Label(tabframe2, text="Input angle between the first and third translation vector:")
    centered_monoclinic_label6.grid(row=14,column=0)
    centered_monoclinic_label7 = Label(tabframe2, text="Input angle between the second and third translation vector:")
    centered_monoclinic_label7.grid(row=15,column=0)
    centered_monoclinic_label8 = Label(tabframe2, text="Input the chemical symbols for the two base atoms separated by space:")
    centered_monoclinic_label8.grid(row=16,column=0)
    centered_monoclinic_entry1 = Entry(tabframe2)
    centered_monoclinic_entry1.grid(row=10,column=1)
    centered_monoclinic_entry2 = Entry(tabframe2)
    centered_monoclinic_entry2.grid(row=11,column=1)
    centered_monoclinic_entry3 = Entry(tabframe2)
    centered_monoclinic_entry3.grid(row=12,column=1)
    centered_monoclinic_entry4 = Entry(tabframe2)
    centered_monoclinic_entry4.grid(row=14,column=1)
    centered_monoclinic_entry5 = Entry(tabframe2)
    centered_monoclinic_entry5.grid(row=15,column=1)
    centered_monoclinic_entry6 = Entry(tabframe2)
    centered_monoclinic_entry6.grid(row=16,column=1)
    centered_monoclinic_button1 = Button(tabframe2, text="Continue", command=lambda:inbetweener("base_centered_monoclinic", new_window,
                                                                                       a = float(centered_monoclinic_entry1.get()),
                                                                                       b = float(centered_monoclinic_entry2.get()),
                                                                                       c = float(centered_monoclinic_entry3.get()),
                                                                                       beta = float(centered_monoclinic_entry4.get()),
                                                                                       gamma = float(centered_monoclinic_entry5.get()),
                                                                                       elements = centered_monoclinic_entry6.get()))
    centered_monoclinic_button1.grid(row=17,column=1)
    
    
    simple_orthorhomic_label1 = Label(tabframe3, text="Input x-direction lattice constant in Ångström: ")
    simple_orthorhomic_label1.grid(row=1,column=0)
    simple_orthorhomic_label2 = Label(tabframe3, text="Input y-direction lattice constant in Ångström: ")
    simple_orthorhomic_label2.grid(row=2,column=0)
    simple_orthorhomic_label3 = Label(tabframe3, text="Input z-direction lattice constant in Ångström: ")
    simple_orthorhomic_label3.grid(row=3,column=0)
    simple_orthorhomic_label4 = Label(tabframe3, text="1 - Simple orthorhomic (also known as primitive orthorhombic)")
    simple_orthorhomic_label4.grid(row=4,column=0)
    simple_orthorhomic_label5 = Label(tabframe3, text="Input the chemical symbol for the base atom:")
    simple_orthorhomic_label5.grid(row=5,column=0)
    simple_orthorhomic_button1 = Button(tabframe3, text="Continue", command=lambda: inbetweener("simple_orthonormic", new_window, 
                                                                                                a = float(simple_orthorhomic_entry1.get()),
                                                                                                b = float(simple_orthorhomic_entry2.get()),
                                                                                                c = float(simple_orthorhomic_entry3.get()),
                                                                                                element = simple_orthorhomic_entry4.get()))
    simple_orthorhomic_button1.grid(row=6,column=1)
    simple_orthorhomic_label6 = Label(tabframe3, text="2 - Base centered orthorhomic")
    simple_orthorhomic_label6.grid(row=7,column=0)
    simple_orthorhomic_label7 = Label(tabframe3, text="Input the chemical symbol for the two base atoms separated by space:")
    simple_orthorhomic_label7.grid(row=8,column=0)
    simple_orthorhomic_button2 = Button(tabframe3, text="Continue", command=lambda: inbetweener("base_centered_orthonormic", new_window, 
                                                                                                a = float(simple_orthorhomic_entry1.get()),
                                                                                                b = float(simple_orthorhomic_entry2.get()),
                                                                                                c = float(simple_orthorhomic_entry3.get()),
                                                                                                elements = simple_orthorhomic_entry5.get()))
    simple_orthorhomic_button2.grid(row=9,column=1)
    simple_orthorhomic_label8 = Label(tabframe3, text="3 - Body centered orthorhomic")
    simple_orthorhomic_label8.grid(row=10,column=0)
    simple_orthorhomic_label9 = Label(tabframe3, text="Input the chemical symbol for the two base atoms separated by space:")
    simple_orthorhomic_label9.grid(row=11,column=0)
    simple_orthorhomic_button3 = Button(tabframe3, text="Continue", command=lambda: inbetweener("body_centered_orthonormic", new_window, 
                                                                                                a = float(simple_orthorhomic_entry1.get()),
                                                                                                b = float(simple_orthorhomic_entry2.get()),
                                                                                                c = float(simple_orthorhomic_entry3.get()),
                                                                                                elements = simple_orthorhomic_entry6.get()))
    simple_orthorhomic_button3.grid(row=12,column=1)
    simple_orthorhomic_label10 = Label(tabframe3, text="4 - Face centered orthorhomic")
    simple_orthorhomic_label10.grid(row=13,column=0)
    simple_orthorhomic_label11 = Label(tabframe3, text="Input the chemical symbol for the base atom:")
    simple_orthorhomic_label11.grid(row=14,column=0)
    simple_orthorhomic_button4 = Button(tabframe3, text="Continue", command=lambda: inbetweener("face_centered_orthonormic", new_window, 
                                                                                                a = float(simple_orthorhomic_entry1.get()),
                                                                                                b = float(simple_orthorhomic_entry2.get()),
                                                                                                c = float(simple_orthorhomic_entry3.get()),
                                                                                                element = simple_orthorhomic_entry7.get()))
    simple_orthorhomic_button4.grid(row=15,column=1)
    simple_orthorhomic_entry1 = Entry(tabframe3)
    simple_orthorhomic_entry1.grid(row=1, column=1)
    simple_orthorhomic_entry2 = Entry(tabframe3)
    simple_orthorhomic_entry2.grid(row=2, column=1)
    simple_orthorhomic_entry3 = Entry(tabframe3)
    simple_orthorhomic_entry3.grid(row=3, column=1)
    simple_orthorhomic_entry4 = Entry(tabframe3)
    simple_orthorhomic_entry4.grid(row=5, column=1)
    simple_orthorhomic_entry5 = Entry(tabframe3)
    simple_orthorhomic_entry5.grid(row=8, column=1)
    simple_orthorhomic_entry6 = Entry(tabframe3)
    simple_orthorhomic_entry6.grid(row=11, column=1)
    simple_orthorhomic_entry7 = Entry(tabframe3)
    simple_orthorhomic_entry7.grid(row=14, column=1)
    
    
    simple_tetragonal_label1 = Label(tabframe4, text="1 - Simple tetragonal")
    simple_tetragonal_label1.grid(row=0,column=0)
    simple_tetragonal_label2 = Label(tabframe4, text="Input x-, y-direction lattice constant in Ånström:")
    simple_tetragonal_label2.grid(row=1,column=0)
    simple_tetragonal_label3 = Label(tabframe4, text="Input z-direction lattice constant in Ånström:")
    simple_tetragonal_label3.grid(row=2,column=0)
    simple_tetragonal_label4 = Label(tabframe4, text="Input the chemical symbol for the base atom:")
    simple_tetragonal_label4.grid(row=3,column=0)
    simple_tetragonal_button1 = Button(tabframe4, text="Continue", command=lambda: inbetweener("simple_tetragonal", new_window, 
                                                                                                lattice_constants = [float(simple_tetragonal_entry1.get()),
                                                                                                float(simple_tetragonal_entry2.get())],
                                                                                                element = simple_tetragonal_entry3.get()))
    simple_tetragonal_button1.grid(row=4,column=1)
    simple_tetragonal_label5 = Label(tabframe4, text="2 - Centered tetragonal (also known as body centered tetragonal)")
    simple_tetragonal_label5.grid(row=5,column=0)
    simple_tetragonal_label6 = Label(tabframe4, text="Input x-, y-direction lattice constant in Ånström:")
    simple_tetragonal_label6.grid(row=6,column=0)
    simple_tetragonal_label7 = Label(tabframe4, text="Input z-direction lattice constant in Ånström:")
    simple_tetragonal_label7.grid(row=7,column=0)
    simple_tetragonal_label8 = Label(tabframe4, text="Input the chemical symbols for the two base atoms separated by space:")
    simple_tetragonal_label8.grid(row=8,column=0)
    simple_tetragonal_button2 = Button(tabframe4, text="Continue", command=lambda: inbetweener("centered_tetragonal", new_window, 
                                                                                                lattice_constants = [float(simple_tetragonal_entry4.get()),
                                                                                                float(simple_tetragonal_entry5.get())],
                                                                                                elements = simple_tetragonal_entry6.get()))
    simple_tetragonal_button2.grid(row=9,column=1)
    simple_tetragonal_entry1 = Entry(tabframe4)
    simple_tetragonal_entry1.grid(row=1,column=1)
    simple_tetragonal_entry2 = Entry(tabframe4)
    simple_tetragonal_entry2.grid(row=2,column=1)
    simple_tetragonal_entry3 = Entry(tabframe4)
    simple_tetragonal_entry3.grid(row=3,column=1)
    simple_tetragonal_entry4 = Entry(tabframe4)
    simple_tetragonal_entry4.grid(row=6,column=1)
    simple_tetragonal_entry5 = Entry(tabframe4)
    simple_tetragonal_entry5.grid(row=7,column=1)
    simple_tetragonal_entry6 = Entry(tabframe4)
    simple_tetragonal_entry6.grid(row=8,column=1)
    
    
    cubic_label1 = Label(tabframe7, text = "Input lattice constant in Ångström:")
    cubic_label1.grid(row=0,column=0)
    cubic_label2 = Label(tabframe7, text = "1 - Simple cubic")
    cubic_label2.grid(row=1,column=0)
    cubic_label3 = Label(tabframe7, text = "Input one element in element notation:")
    cubic_label3.grid(row=2,column=0)
    cubic_button1 = Button(tabframe7, text = "Continue", command=lambda: inbetweener("simple_cubic", new_window,
                                                                                     lattice_constant = float(cubic_entry1.get()),
                                                                                     element = cubic_entry2.get()))
    cubic_button1.grid(row=3,column=1)
    cubic_label2 = Label(tabframe7, text = "2 - Body centered cubic")
    cubic_label2.grid(row=4,column=0)
    cubic_label3 = Label(tabframe7, text = "Input the chemical symbols for the two base atoms separated by space:")
    cubic_label3.grid(row=5,column=0)
    cubic_button2 = Button(tabframe7, text = "Continue", command=lambda: inbetweener("body_centered_cubic", new_window,
                                                                                     lattice_constant = float(cubic_entry1.get()),
                                                                                     elements = cubic_entry3.get()))
    cubic_button2.grid(row=6,column=1)
    cubic_label2 = Label(tabframe7, text = "3 - Face centered cubic")
    cubic_label2.grid(row=7,column=0)
    cubic_label3 = Label(tabframe7, text = "Input the chemical symbols for the four base atoms separated by space:")
    cubic_label3.grid(row=8,column=0)
    cubic_button3 = Button(tabframe7, text = "Continue", command=lambda: inbetweener("face_centered_cubic", new_window,
                                                                                     lattice_constant = float(cubic_entry1.get()),
                                                                                     elements = cubic_entry4.get()))
    cubic_button3.grid(row=9,column=1)
    cubic_entry1 = Entry(tabframe7)
    cubic_entry1.grid(row=0, column=1)
    cubic_entry2 = Entry(tabframe7)
    cubic_entry2.grid(row=2, column=1)
    cubic_entry3 = Entry(tabframe7)
    cubic_entry3.grid(row=5, column=1)
    cubic_entry4 = Entry(tabframe7)
    cubic_entry4.grid(row=8, column=1)
    
    new_window.mainloop()