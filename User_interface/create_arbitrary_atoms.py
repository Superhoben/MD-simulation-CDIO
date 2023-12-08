import os
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.lattice.triclinic import Triclinic
from ase.lattice.monoclinic import SimpleMonoclinic, BaseCenteredMonoclinic
from ase.lattice.orthorhombic import SimpleOrthorhombic, BaseCenteredOrthorhombic, BodyCenteredOrthorhombic, FaceCenteredOrthorhombic
from ase.lattice.tetragonal import SimpleTetragonal, CenteredTetragonal
from ase.lattice.cubic import SimpleCubic, BodyCenteredCubic, FaceCenteredCubic
from math import floor
from numpy import cbrt
from ase.build import make_supercell


def inbetweener(system, window, **kwargs):
    """Transition from to create_view_and_save_crystal_guide and destroy window"""
    
    window.destroy()
    create_view_and_save_crystal_guided(system, **kwargs)


def create_view_and_save_crystal_guided(system, a=0, b=0, c=0, alpha=0, beta=0, gamma=0, element=0,
                                        elements=0, lattice_constant=0, lattice_constants=0):
    """Is called by create_atom() to as the next atep in creating the crystal

    Args:
        system(str): Name of the crystal system the user is creating
        a, b, c (int): Length of the translations vectors, if two or more vectors always share lengths
                       in a crystal system the only the nessecary vectors should be added. For example
                       for the cubic cell only the parameter 'a' shoud be entered
        alpha, beta, gamma: The angles between the translation vectors, if two or more vectors always have
                            perpendicular angles in crystal system the only the nessecary vectors should be 
                            added. For example for the cubic cell no parameter should be entered
        elements(str or List(str)): The chemical symbol for the base atom as a string or if there are many
                                    base atoms as a list of strings.
        lattice_constant: To be adjusted
        lattice_constant(dict): To be adjusted

    Returns:
        None
    """
    primitive_cell = 0
    if system == "triclinic":
        primitive_cell = Triclinic(symbol=element, latticeconstant={'a': a, 'b': b, 'c': c, 'alpha': alpha,
                                                                    'beta': beta, 'gamma': gamma})

    elif system == "monoclinic":
        primitive_cell = SimpleMonoclinic(symbol=element, latticeconstant={'a': a, 'b': b, 'c': c,
                                                                           'beta': beta, 'gamma': gamma})

    elif system == "base_centered_monoclinic":
        elements = elements.split()
        primitive_cell = BaseCenteredMonoclinic(symbol=elements[0], latticeconstant={'a': a, 'b': b, 'c': c,
                                                                                     'beta': beta, 'gamma': gamma})
        primitive_cell.set_chemical_symbols(elements)

    elif system == "simple_orthonormic":
        primitive_cell = SimpleOrthorhombic(symbol=element, latticeconstant={'a': a, 'b': b, 'c': c})

    elif system == "base_centered_orthonormic":
        elements = elements.split()
        primitive_cell = BaseCenteredOrthorhombic(symbol=elements[0], latticeconstant={'a': a, 'b': b, 'c': c})
        primitive_cell.set_chemical_symbols(elements)

    elif system == "body_centered_orthonormic":
        elements = elements.split()
        primitive_cell = BodyCenteredOrthorhombic(symbol=elements[0], latticeconstant={'a': a, 'b': b, 'c': c})
        primitive_cell.set_chemical_symbols(elements)

    elif system == "face_centered_orthonormic":
        primitive_cell = FaceCenteredOrthorhombic(symbol=element, latticeconstant={'a': a, 'b': b, 'c': c})

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

    view(primitive_cell, block=False)

    new_window = Tk()
    new_window.title("Create trajectory file")

    info_label = Label(new_window, text=f"Number of atoms in primitive cell: {len(primitive_cell)}").grid(row=0, column=0)
    input_label = Label(new_window, text="Input the target number of atoms in the super cell:").grid(row=1, column=0)
    input_entry = Entry(new_window)
    input_entry.grid(row=1, column=1, padx=(40, 0))
    traj_label = Label(new_window, text="Input name of .traj file:").grid(row=2, column=0)
    traj_entry = Entry(new_window)
    traj_entry.grid(row=2, column=1, padx=(40, 0))
    final_button = Button(new_window, text="Create .traj file", 
                          command=lambda: create_traj(len(primitive_cell), input_entry.get(),
                                                      traj_entry.get(), primitive_cell, new_window))
    final_button.grid(row=3, column=1, padx=(40, 0))

    new_window.mainloop()
    new_window.quit()


def create_atom():
    """Makes a window from which the user can create their own crystal"""

    new_window = Tk()
    new_window.title("Create and save input atoms objects")
    new_window.geometry("770x500")

    tabs = ttk.Notebook(new_window)
    tabs.pack()

    frame_heights = 500
    frame_widths = 770

    # Tabs
    tabframe0 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe0.pack_propagate(False)
    tabframe0.grid_propagate(False)
    tabframe0.pack(expand=True, fill=BOTH)

    tabframe1 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe1.pack_propagate(False)
    tabframe1.grid_propagate(False)
    tabframe1.pack(expand=True, fill=BOTH)

    tabframe2 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe2.pack_propagate(False)
    tabframe2.grid_propagate(False)
    tabframe2.pack(expand=True, fill=BOTH)

    tabframe3 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe3.pack_propagate(False)
    tabframe3.grid_propagate(False)
    tabframe3.pack(expand=True, fill=BOTH)

    tabframe4 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe4.pack_propagate(False)
    tabframe4.grid_propagate(False)
    tabframe4.pack(expand=True, fill=BOTH)

    tabframe5 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe5.pack_propagate(False)
    tabframe5.grid_propagate(False)
    tabframe5.pack(expand=True, fill=BOTH)

    tabframe6 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe6.pack_propagate(False)
    tabframe6.grid_propagate(False)
    tabframe6.pack(expand=True, fill=BOTH)

    tabframe7 = Frame(tabs, height=frame_heights, width=frame_widths, bg="SkyBlue1")
    tabframe7.pack_propagate(False)
    tabframe7.grid_propagate(False)
    tabframe7.pack(expand=True, fill=BOTH)

    tabs.add(tabframe0, text="Guide")
    tabs.add(tabframe1, text="Triclinic")
    tabs.add(tabframe2, text="Monoclinic")
    tabs.add(tabframe3, text="Orthorhombic")
    tabs.add(tabframe4, text="Tetragonal")
    tabs.add(tabframe5, text="Trigonal")
    tabs.add(tabframe6, text="Hexagonal")
    tabs.add(tabframe7, text="Cubic")

    # General description of crystal systems
    text_label = Label(tabframe0, text="Select a tab corresponding to the crystal system you want to create.\n"+
                       "System names and description of the lattice vectors follows below. \n" +
                       "1 - Triclinic, arbitrary lengths and directions \n"+
                       "2 - Monoclinic, arbitrary lengths, two vectors are orthagonal \n"+
                       "3 - Orthorhombic, arbitrary lengths, all vectors are orthagonal \n"+
                       "4 - Tetragonal, two equal lenghts, all vectors are orthagonol \n"+
                       "5 - Trigonal, not yet implemented \n"+
                       "6 - Hexagonal, not yet implemented \n"+
                       "7 - Cubic, equal lenghts and all vectors are orthagonal")
    text_label.place(relx=0.5, rely=0.5, anchor=S)

    # Triclinc cell inputs
    label0 = Label(tabframe1, text="Triclinic cell, arbitrary lengths and angles of translation vectors")
    label0.grid(row=0, column=0, columnspan=2)
    label1 = Label(tabframe1, text="Input length of the first lattice translation vector in Å:")
    label1.grid(row=1, column=0)
    entry1 = Entry(tabframe1)
    entry1.grid(row=1, column=1, padx=(40, 0))
    label2 = Label(tabframe1, text="Input length of the second lattice translation vector in Å:")
    label2.grid(row=2, column=0)
    entry2 = Entry(tabframe1)
    entry2.grid(row=2, column=1, padx=(40, 0))
    label3 = Label(tabframe1, text="Input length of the third lattice translation vector in Å:")
    label3.grid(row=3, column=0)
    entry3 = Entry(tabframe1)
    entry3.grid(row=3, column=1, padx=(40, 0))
    label4 = Label(tabframe1, text="Input angle between the first and second translation vector in degrees:")
    label4.grid(row=4, column=0)
    entry4 = Entry(tabframe1)
    entry4.grid(row=4, column=1, padx=(40, 0))
    label5 = Label(tabframe1, text="Input angle between the first and third translation vector in degrees:")
    label5.grid(row=5, column=0)
    entry5 = Entry(tabframe1)
    entry5.grid(row=5, column=1, padx=(40, 0))
    label6 = Label(tabframe1, text="Input angle between the second and third translation vector in degrees:")
    label6.grid(row=6, column=0)
    entry6 = Entry(tabframe1)
    entry6.grid(row=6, column=1, padx=(40, 0))
    label7 = Label(tabframe1, text="Input the chemical symbol for the base atom:")
    label7.grid(row=7, column=0)
    entry7 = Entry(tabframe1)
    entry7.grid(row=7, column=1, padx=(40, 0))

    triclinic_button = Button(tabframe1, text="Continue",
                              command=lambda: inbetweener("triclinic", new_window,
                                                          a=float(entry1.get()),
                                                          b=float(entry2.get()),
                                                          c=float(entry3.get()),
                                                          alpha=float(entry4.get()),
                                                          beta=float(entry5.get()),
                                                          gamma=float(entry6.get()),
                                                          element=entry7.get()))
    triclinic_button.grid(row=8, column=1, padx=(40, 0))

    # Monoclinic cell inputs
    monoclinic_label0 = Label(tabframe2, text="Monoclinic cell, arbitrary lengths, first and second translation vectors are orthagonal")
    monoclinic_label0.grid(row=0, column=0, columnspan=2, pady=(10, 5))
    monoclinic_label1 = Label(tabframe2, text="Input length of a1, the first translation vector, in Å:")
    monoclinic_label1.grid(row=1, column=0)
    monoclinic_entry1 = Entry(tabframe2)
    monoclinic_entry1.grid(row=1, column=1, padx=(40, 0))
    monoclinic_label2 = Label(tabframe2, text="Input length of a2, the second translation vector, in Å:")
    monoclinic_label2.grid(row=2, column=0)
    monoclinic_entry2 = Entry(tabframe2)
    monoclinic_entry2.grid(row=2, column=1, padx=(40, 0))
    monoclinic_label3 = Label(tabframe2, text="Input length of a3, the third translation vector, in Å:")
    monoclinic_label3.grid(row=3,column=0)
    monoclinic_entry3 = Entry(tabframe2)
    monoclinic_entry3.grid(row=3, column=1, padx=(40, 0))
    monoclinic_label4 = Label(tabframe2, text="Input angle between the first and third translation vector:")
    monoclinic_label4.grid(row=4, column=0)
    monoclinic_entry4 = Entry(tabframe2)
    monoclinic_entry4.grid(row=6, column=1, padx=(40, 0))
    monoclinic_label5 = Label(tabframe2, text="Input angle between the second and third lattice translation vector:")
    monoclinic_label5.grid(row=6, column=0)
    monoclinic_entry5 = Entry(tabframe2)
    monoclinic_entry5.grid(row=6, column=1, padx=(40, 0))

    # Basis inputs simple monoclinic
    simple_monoclinic_label0 = Label(tabframe2, text="Simple monoclinic, base atom in 0")
    simple_monoclinic_label0.grid(row=7, column=0, columnspan=2, pady=(20, 5))
    simple_monoclinic_label1 = Label(tabframe2, text="Input the chemical symbol for the base atom: ")
    simple_monoclinic_label1.grid(row=8, column=0)
    simple_monoclinic_entry1 = Entry(tabframe2)
    simple_monoclinic_entry1.grid(row=8, column=1, padx=(40, 0))
    simple_monoclinic_button1 = Button(tabframe2, text="Continue",
                                       command=lambda:inbetweener("monoclinic", new_window,
                                                                  a=float(monoclinic_entry1.get()),
                                                                  b=float(monoclinic_entry2.get()),
                                                                  c=float(monoclinic_entry3.get()),
                                                                  beta=float(monoclinic_entry4.get()),
                                                                  gamma=float(monoclinic_entry5.get()),
                                                                  element=simple_monoclinic_entry1.get()))
    simple_monoclinic_button1.grid(row=9, column=1, padx=(40, 0))

    # Basis inputs base centered monoclinc
    centered_monoclinic_label1 = Label(tabframe2, text="Base centered monoclinic, base atoms in 0 and a1+a2")
    centered_monoclinic_label1.grid(row=10, column=0, columnspan=2, pady=(20, 5))
    centered_monoclinic_label2 = Label(tabframe2, text="Input the chemical symbols for the two base atoms separated by space:")
    centered_monoclinic_label2.grid(row=11, column=0)
    centered_monoclinic_entry2 = Entry(tabframe2)
    centered_monoclinic_entry2.grid(row=11, column=1, padx=(40, 0))
    centered_monoclinic_button = Button(tabframe2, text="Continue", 
                                        command=lambda:inbetweener("base_centered_monoclinic", new_window,
                                                                   a=float(monoclinic_entry1.get()),
                                                                   b=float(monoclinic_entry2.get()),
                                                                   c=float(monoclinic_entry3.get()),
                                                                   beta=float(monoclinic_entry4.get()),
                                                                   gamma=float(monoclinic_entry5.get()),
                                                                   elements=centered_monoclinic_entry2.get()))
    centered_monoclinic_button.grid(row=12, column=1, padx=(40, 0))

    # Orthorhombic cell inputs
    orthorhombic_label0 = Label(tabframe3, text="Orthorhombic cell, arbitrary lengths but all translation vectors are orthagonal")
    orthorhombic_label0.grid(row=0, column=0, columnspan=2, pady=(10, 5))
    orthorhombic_label1 = Label(tabframe3, text="Input x-direction lattice constant in Ångström: ")
    orthorhombic_label1.grid(row=1, column=0)
    orthorhombic_entry1 = Entry(tabframe3)
    orthorhombic_entry1.grid(row=1, column=1, padx=(40, 0))
    orthorhombic_label2 = Label(tabframe3, text="Input y-direction lattice constant in Ångström: ")
    orthorhombic_label2.grid(row=2, column=0)
    orthorhombic_entry2 = Entry(tabframe3)
    orthorhombic_entry2.grid(row=2, column=1, padx=(40, 0))
    orthorhombic_label3 = Label(tabframe3, text="Input z-direction lattice constant in Ångström: ")
    orthorhombic_label3.grid(row=3, column=0)
    orthorhombic_entry3 = Entry(tabframe3)
    orthorhombic_entry3.grid(row=3, column=1, padx=(40, 0))

    # Basis inputs simple orthorhombic
    simple_orthorhombic_label0 = Label(tabframe3, text="Simple orthorhombic (also known as primitive orthorhombic, base atom in 0")
    simple_orthorhombic_label0.grid(row=4, column=0, columnspan=2, pady=(20, 5))
    simple_orthorhombic_label1 = Label(tabframe3, text="Input the chemical symbol for the base atom:")
    simple_orthorhombic_label1.grid(row=5, column=0)
    simple_orthorhombic_entry1 = Entry(tabframe3)
    simple_orthorhombic_entry1.grid(row=5, column=1, padx=(40, 0))
    simple_orthorhombic_button1 = Button(tabframe3, text="Continue", 
                                         command=lambda: inbetweener("simple_orthonormic", new_window, 
                                                                     a = float(orthorhombic_entry1.get()),
                                                                     b = float(orthorhombic_entry2.get()),
                                                                     c = float(orthorhombic_entry3.get()),
                                                                     element = simple_orthorhombic_entry1.get()))
    simple_orthorhombic_button1.grid(row=6, column=1, padx=(40, 0))

    # Basis inputs base centered orthorhombic 
    base_centered_orthorhombic_label0 = Label(tabframe3, text="Base centered orthorhombic, base atoms in 0 and a1+a2")
    base_centered_orthorhombic_label0.grid(row=7, column=0, columnspan=2, pady=(20, 5))
    base_centered_orthorhombic_label1 = Label(tabframe3, text="Input the chemical symbol for the two base atoms separated by space:")
    base_centered_orthorhombic_label1.grid(row=8, column=0)
    base_centered_orthorhombic_entry1 = Entry(tabframe3)
    base_centered_orthorhombic_entry1.grid(row=8, column=1, padx=(40, 0))
    base_centered_orthorhombic_button1 = Button(tabframe3, text="Continue",
                                                command=lambda: inbetweener("base_centered_orthonormic", new_window, 
                                                                            a=float(orthorhombic_entry1.get()),
                                                                            b=float(orthorhombic_entry2.get()),
                                                                            c=float(orthorhombic_entry3.get()),
                                                                            elements=base_centered_orthorhombic_entry1.get()))
    base_centered_orthorhombic_button1.grid(row=9, column=1, padx=(40, 0))

    # Basis inputs body centered orthorhombic
    body_centered_orthorhombic_label8 = Label(tabframe3, text="Body centered orthorhombic, base atoms in 0 and a1+a2+a3")
    body_centered_orthorhombic_label8.grid(row=10, column=0, columnspan=2, pady=(20, 5))
    body_centered_orthorhombic_label9 = Label(tabframe3, text="Input the chemical symbol for the two base atoms separated by space:")
    body_centered_orthorhombic_label9.grid(row=11, column=0)
    body_centered_orthorhombic_entry1 = Entry(tabframe3)
    body_centered_orthorhombic_entry1.grid(row=11, column=1, padx=(40, 0))
    body_centered_orthorhombic_button3 = Button(tabframe3, text="Continue", 
                                                command=lambda: inbetweener("body_centered_orthonormic", new_window, 
                                                                            a=float(orthorhombic_entry1.get()),
                                                                            b=float(orthorhombic_entry2.get()),
                                                                            c=float(orthorhombic_entry3.get()),
                                                                            elements=body_centered_orthorhombic_entry1.get()))
    body_centered_orthorhombic_button3.grid(row=12, column=1, padx=(40, 0))

    # Basis inputs face centered orthorhombic
    face_centered_orthorhombic_label10 = Label(tabframe3, text="Face centered orthorhombic")
    face_centered_orthorhombic_label10.grid(row=13, column=0, columnspan=2, pady=(20, 5))
    face_centered_orthorhombic_label11 = Label(tabframe3, text="Input the chemical symbol for the base atom:")
    face_centered_orthorhombic_label11.grid(row=14, column=0)
    face_centered_orthorhombic_entry1 = Entry(tabframe3)
    face_centered_orthorhombic_entry1.grid(row=14, column=1, padx=(40, 0))
    face_centered_orthorhombic_button1 = Button(tabframe3, text="Continue", 
                                         command=lambda: inbetweener("face_centered_orthonormic", new_window, 
                                                                     a=float(orthorhombic_entry1.get()),
                                                                     b=float(orthorhombic_entry2.get()),
                                                                     c=float(orthorhombic_entry3.get()),
                                                                     element=face_centered_orthorhombic_entry1.get()))
    face_centered_orthorhombic_button1.grid(row=15, column=1, padx=(40, 0))

    # Tetragonal cell inputs
    tetragonal_label0 = Label(tabframe4, text="Tetragonal, two equal lenghts, all vectors are orthagonol")
    tetragonal_label0.grid(row=0, column=0, columnspan=2, pady=(10, 5))
    tetragonal_label1 = Label(tabframe4, text="Input x-, y-direction lattice constant in Ångström:")
    tetragonal_label1.grid(row=1, column=0)
    tetragonal_entry1 = Entry(tabframe4)
    tetragonal_entry1.grid(row=1, column=1, padx=(40, 0))
    tetragonal_label2 = Label(tabframe4, text="Input z-direction lattice constant in Ångström:")
    tetragonal_label2.grid(row=2, column=0)
    tetragonal_entry2 = Entry(tabframe4)
    tetragonal_entry2.grid(row=2, column=1, padx=(40, 0))

    # Basis input simple tetragonal
    simple_tetragonal_label1 = Label(tabframe4, text="Simple tetragonal, base atom in 0")
    simple_tetragonal_label1.grid(row=3, column=0, columnspan=2, pady=(20, 5))
    simple_tetragonal_label2 = Label(tabframe4, text="Input the chemical symbol for the base atom:")
    simple_tetragonal_label2.grid(row=4, column=0)
    simple_tetragonal_entry2 = Entry(tabframe4)
    simple_tetragonal_entry2.grid(row=4, column=1, padx=(40, 0))
    simple_tetragonal_button1 = Button(tabframe4, text="Continue",
                                       command=lambda: inbetweener("simple_tetragonal", new_window,
                                                                   lattice_constants=[float(tetragonal_entry1.get()),
                                                                   float(tetragonal_entry2.get())],
                                                                   element = simple_tetragonal_entry2.get()))
    simple_tetragonal_button1.grid(row=5, column=1, padx=(40, 0))

    # Basis input centered tetragonal (also known as body centered tetragonal)
    centered_tetragonal_label1 = Label(tabframe4, text="Centered tetragonal (also known as body centered tetragonal)")
    centered_tetragonal_label1.grid(row=4, column=0)
    centered_tetragonal_label2 = Label(tabframe4, text="Input the chemical symbols for the two base atoms separated by space:")
    centered_tetragonal_label2.grid(row=5, column=0)
    centered_tetragonal_entry2 = Entry(tabframe4)
    centered_tetragonal_entry2.grid(row=5, column=1, padx=(40, 0))
    centered_tetragonal_button2 = Button(tabframe4, text="Continue",
                                         command=lambda: inbetweener("centered_tetragonal", new_window, 
                                                                     lattice_constants=[float(tetragonal_entry1.get()),
                                                                     float(tetragonal_entry2.get())],
                                                                     elements=centered_tetragonal_entry2.get()))
    centered_tetragonal_button2.grid(row=6, column=1, padx=(40, 0))

    # Cubic cell inputs
    cubic_label0 = Label(tabframe7, text = "Input lattice constant in Ångström:")
    cubic_label0.grid(row=0, column=0, columnspan=2, pady=(10, 5))
    cubic_label1 = Label(tabframe7, text = "Input lattice constant in Ångström:")
    cubic_label1.grid(row=1, column=0)
    cubic_entry1 = Entry(tabframe7)
    cubic_entry1.grid(row=1, column=1, padx=(40, 0))

    # Basis inputs simple cubic
    simple_cubic_label0 = Label(tabframe7, text="Simple cubic, base atom in 0")
    simple_cubic_label0.grid(row=2, column=0, columnspan=2, pady=(20, 5))
    simple_cubic_label1 = Label(tabframe7, text="Input the chemical symbol for the basis atom:")
    simple_cubic_label1.grid(row=3, column=0)
    simple_cubic_entry1 = Entry(tabframe7)
    simple_cubic_entry1.grid(row=3, column=1, padx=(40, 0))
    simple_cubic_button1 = Button(tabframe7, text="Continue", 
                           command=lambda: inbetweener("simple_cubic", new_window,
                                                       lattice_constant=float(cubic_entry1.get()),
                                                       element=simple_cubic_entry1.get()))
    simple_cubic_button1.grid(row=4, column=1, padx=(40, 0))

    # Basis inputs body centered cubic
    body_centered_cubic_label0 = Label(tabframe7, text="Body centered cubic, base atoms in 0 and a1+a2+a3")
    body_centered_cubic_label0.grid(row=5, column=0, columnspan=2, pady=(20, 5))
    body_centered_cubic_label1 = Label(tabframe7, text="Input the chemical symbols for the two base atoms separated by space:")
    body_centered_cubic_label1.grid(row=6, column=0)
    body_centered_cubic_entry1 = Entry(tabframe7)
    body_centered_cubic_entry1.grid(row=6, column=1, padx=(40, 0))
    body_centered_cubic_button1 = Button(tabframe7, text="Continue",
                                         command=lambda: inbetweener("body_centered_cubic", new_window,
                                                                     lattice_constant=float(cubic_entry1.get()),
                                                                     elements=body_centered_cubic_entry1.get()))
    body_centered_cubic_button1.grid(row=7, column=1, padx=(40, 0))


    face_centered_cubic_label0 = Label(tabframe7, text="Face centered cubic, base atoms 0, a1+a2, a1+a3 and a2+a3")
    face_centered_cubic_label0.grid(row=8, column=0, columnspan=2, pady=(20, 5))
    face_centered_cubic_label1 = Label(tabframe7, text="Input the chemical symbols for the four base atoms separated by space:")
    face_centered_cubic_label1.grid(row=9, column=0)
    face_centered_cubic_entry1 = Entry(tabframe7)
    face_centered_cubic_entry1.grid(row=9, column=1, padx=(40, 0))
    face_centered_cubic_button1 = Button(tabframe7, text="Continue",
                                         command=lambda: inbetweener("face_centered_cubic", new_window,
                                                                     lattice_constant=float(cubic_entry1.get()),
                                                                     elements=face_centered_cubic_entry1.get()))
    face_centered_cubic_button1.grid(row=10, column=1, padx=(40,0))

    new_window.mainloop()
    new_window.quit()
    return new_window


def create_traj(n_atoms_prim, target_n_atoms, file_name, primitive_cell, new_window):
    if target_n_atoms.isdigit():
        target_n_atoms = int(target_n_atoms)
        if target_n_atoms < n_atoms_prim:
            messagebox.showerror("Invalid data", "Number of atoms in supercell has to be greater than number in primitive cell")
        else:
            n = floor(cbrt(target_n_atoms/n_atoms_prim))
            M = [[n, 0, 0], [0, n, 0], [0, 0, n]]
            atoms = make_supercell(primitive_cell, M)
            path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
            location_and_name = path_to_traj_folder + file_name + ".traj"
            traj = Trajectory(location_and_name, "w")
            traj.write(atoms)
            new_window.destroy()
            messagebox.showinfo("Success!", "A trajectory file has been created!")
    else:
        messagebox.showerror("Invalid data", "Please put a proper value")


if __name__ == '__main__':
    create_atom_window = create_atom()
#    create_atom_window.mainloop()
