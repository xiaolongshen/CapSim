#! /usr/bin/env python
#
#This file is used to read data from a batch script to build a "System" object
#for CapSim to execute.

import tkFileDialog as tkfd, tkMessageBox as tkmb, cPickle as pickle, sys, os, csv
import _winreg as wreg

if sys.path[0][-3:] == 'zip': 
    os.chdir(sys.path[0][:-12])
    path = sys.path[0][:-12]
else: path = sys.path[0]

CapSimReg = wreg.ConnectRegistry(None, wreg.HKEY_CURRENT_USER)
CapSimKey = wreg.OpenKey(CapSimReg, r'Software\CapSim')
Filepath =wreg.QueryValueEx(CapSimKey, 'FilePath')[0]
CapSimKey.Close()

sys.path.append(path + r'/solvers')

from Tkinter             import Frame, Label, Button
from capsim_object_types import CapSimWindow, System, Chemical, Matrix, MatrixComponent, Sorption, Reaction, Layer, Coefficient, BC, IC, SolidIC
from reactioneditor      import Reactant, Product
from capsim_functions    import formula_converter
from random              import random
from numpy               import exp, cos, sin, log, pi

class BatchFile:
    """A class that creates a window to execute batch CapSim files."""

    def __init__(self, master, system):

        self.master   = master
        self.system   = system
        self.fonttype = system.fonttype
        self.version  = system.version
        self.tframe   = Frame(master.tframe)
        self.frame    = Frame(master.frame)
        self.bframe   = Frame(master.bframe)
        self.top      = None

    def make_widgets(self):

        self.label      = Label(self.frame, text = 'Please select the batch ' +
                                'file for CapSim to use to build the system:')
        self.loadbutton = Button(self.frame, text = 'Load File', width = 20,
                                 command = self.get_filename)
        self.runbutton  = Button(self.frame, text = 'Run Batch', width = 20,
                                 command = self.run)
        self.blank      = Label(self.frame, text = '')

        self.label.grid(row = 0, padx = 10)
        self.blank.grid(row = 1)
        self.loadbutton.grid(row = 2)

        self.loadbutton.bind('<Return>', self.get_filename)
        self.runbutton.bind('<Return>', self.run)
        self.focusbutton = self.loadbutton
        
    def get_filename(self, event = None):

        self.filename = tkfd.askopenfilename(initialdir = Filepath +r'\batch_files',  filetypes = [('CapSim batch files','*.txt')])
        if self.filename != '':

            systems       = []

            commands      = open(self.filename, 'r')
            run = 0
            n = 1

            exec(commands)

            isotherms      = ['Linear--Kd specified', 'Linear--Kocfoc', 'Freundlich', 'Langmuir']
            if len(systems) > 0:
                for i in range(len(systems)):

                    if systems[i].layers[0].number == 0:
                        systems[i].dep  = 'Deposition'
                        systems[i].Vdep = systems[i].layers[0].h
                        systems[i].ptotal = sum(systems[i].players[1:])
                    else:
                        systems[i].dep  = 0
                        systems[i].Vdep = 0
                        systems[i].ptotal = sum(systems[i].players)

                    matrix_list = []
                    for matrix in systems[i].matrices:
                        matrix_list.append(matrix.name)
                        volume = 0
                        matrix.e = 0
                        matrix.rho = 0
                        moc = 0
                        for component in matrix.components:
                            volume = volume + component.mfraction / component.rho
                        for component in matrix.components:
                            component.fraction = component.mfraction / component.rho / volume
                            matrix.e = matrix.e + component.fraction * component.e
                            matrix.rho = matrix.rho + component.fraction * component.rho
                            moc = moc + component.fraction * component.foc
                        matrix.foc = moc/matrix.rho

                        matrix.e   = round(matrix.e,6)
                        matrix.rho = round(matrix.rho,6)
                        matrix.foc = round(matrix.foc,6)

                    for component in systems[i].components:
                        for chemical in systems[i].chemicals:
                            systems[i].sorptions[component.name][chemical.name].thalf   = 1
                            if systems[i].sorptions[component.name][chemical.name].isotherm == isotherms[0]:
                                systems[i].sorptions[component.name][chemical.name].kdesorp = component.e*systems[i].sorptions[component.name][chemical.name].ksorp/component.rho/systems[i].sorptions[component.name][chemical.name].K
                                if component.e*systems[i].sorptions[component.name][chemical.name].ksorp > 0:
                                    systems[i].sorptions[component.name][chemical.name].thalf   = round(0.693/component.e*systems[i].sorptions[component.name][chemical.name].ksorp/(1+component.e/component.rho/systems[i].sorptions[component.name][chemical.name].K), 4)
                            if systems[i].sorptions[component.name][chemical.name].isotherm == isotherms[1]:
                                systems[i].sorptions[component.name][chemical.name].kdesorp = component.e*systems[i].sorptions[component.name][chemical.name].ksorp/component.rho/(10**systems[i].sorptions[component.name][chemical.name].Koc)/component.foc
                                if component.e*systems[i].sorptions[component.name][chemical.name].ksorp > 0:
                                    systems[i].sorptions[component.name][chemical.name].thalf   = round(0.693/component.e*systems[i].sorptions[component.name][chemical.name].ksorp/(1+component.e/component.rho/10**systems[i].sorptions[component.name][chemical.name].Koc/component.foc), 4)
                            if systems[i].sorptions[component.name][chemical.name].isotherm == isotherms[2]:
                                systems[i].sorptions[component.name][chemical.name].kdesorp = component.e*systems[i].sorptions[component.name][chemical.name].ksorp/component.rho/(systems[i].sorptions[component.name][chemical.name].Kf**(1/systems[i].sorptions[component.name][chemical.name].N))
                                if component.e*systems[i].sorptions[component.name][chemical.name].ksorp > 0:
                                    systems[i].sorptions[component.name][chemical.name].thalf   = round(0.693/component.e*systems[i].sorptions[component.name][chemical.name].ksorp/(1+component.e/component.rho/systems[i].sorptions[component.name][chemical.name].Kf), 4)
                            if systems[i].sorptions[component.name][chemical.name].isotherm == isotherms[3]:
                                systems[i].sorptions[component.name][chemical.name].kdesorp = component.e*systems[i].sorptions[component.name][chemical.name].ksorp/component.rho/systems[i].sorptions[component.name][chemical.name].b

                    systems[i].delz = []
                    for j in range(systems[i].nlayers):
                        systems[i].layers[j].type_index = matrix_list.index(systems[i].layers[j].type)
                        delz_temp = systems[i].layers[j].h/ systems[i].players[j]
                        if delz_temp > 1:  systems[i].delz.append(round(delz_temp,3))
                        else:
                            j = 2
                            while delz_temp/100 < 0.1**j:
                                j = j + 1
                            systems[i].delz.append(round(delz_temp,j))

                    if systems[i].adv == 'None':
                        systems[-1].Vdar   = 0
                        systems[-1].Vtidal = 0
                        systems[-1].ptidal = 10
                    if systems[-1].bio == 'None':
                        systems[-1].hbio   = 0
                        systems[-1].sigma  = 0
                        systems[-1].Dbiop  = 0
                        systems[-1].Dbiopw = 0
                    if systems[-1].con == 'None':
                        systems[-1].hcon   = 0
                        systems[-1].t90    = 0

                    systems[i].bltype             = 'River'
                    systems[i].blcoefs            = {}
                    systems[i].blcoefs['vx']      = 1.
                    systems[i].blcoefs['n']       = 0.02
                    systems[i].blcoefs['hriver']  = 5.
                    systems[i].blcoefs['rh']      = 5.
                    systems[i].blcoefs['nu']      = 1e-6

                    systems[i].blcoefs['rhoair']  = 1.
                    systems[i].blcoefs['rhowater']= 1000.
                    systems[i].blcoefs['vwind']   = 5.
                    systems[i].blcoefs['hlake']   = 10.
                    systems[i].blcoefs['llake']   = 1000.

            self.systems = systems

            self.loadbutton.grid_forget()
            if run == 0:    self.runbutton.grid(row = 2)
            self.runbutton.focus_set()

    def run(self, event = None):

        self.master.tk.quit()

    def show_error(self, error, n, event = None):

        tkmb.showerror(title = 'Run Error', message = 'Unable to process' + 'Line %d, error: %r' % (n, error))

        self.runbutton.grid_forget()
        self.loadbutton.grid(row = 2)
        self.loadbutton.focus_set()

        return 2

def make_batch(system):

    root = CapSimWindow(buttons = 2)
    root.make_window(BatchFile(root, system))
    root.mainloop()

    if root.main.get() == 1:
        root.window.systems     = None
        root.window.type        = None
        root.window.filenames   = None

    root.destroy()

    return root.window.systems, root.window.type, root.window.filenames, root.main.get()

def get_outputfile(filename, version, fonttype, formulatype):

    content       = []
    file          = open(filename, 'r')
    file_content  = csv.reader(file)
    for row in file_content:    content.append(row)
    rows     = len(content)

    row_file            = []
    row_units           = []
    row_chemicals       = []
    row_matrices        = []
    row_components      = []
    row_sorptions       = []
    row_layers          = []
    row_reactions       = []
    row_coefficients    = []
    row_system          = []
    row_conditions      = []
    row_solver          = []
    row_ow              = []

    # Determine the rows of each section
    row = 0

    while row < rows:
        if len(content[row]) > 0:
            if content[row][0] == 'System Units':                       row_units        = row
            if content[row][0] == 'Chemicals':                          row_chemicals    = row
            if content[row][0] == 'Matrices':                           row_matrices     = row
            if content[row][0] == 'Matrix Components':                  row_components   = row
            if content[row][0] == 'Sorptions':                          row_sorptions    = row
            if content[row][0] == 'Layers':                             row_layers       = row
            if content[row][0] == 'Reactions':                          row_reactions    = row
            if content[row][0] == 'Reaction Coefficients':              row_coefficients = row
            if content[row][0] == 'System Properties':                  row_system       = row
            if content[row][0] == 'Auxiliary Conditions':               row_conditions   = row
            if content[row][0] == 'Solver Options':                     row_solver       = row
            if content[row][0] == 'Overlying Water Column Properties':  row_ow           = row
        row = row + 1

    system = System(version, fonttype, formulatype)
    # Read information in each section

    system.cpsmfilename    = 'input'
    system.batchfileoption = 'None'
    system.batchfilename   = 'Batch_file'

    lengthunits = ['um', 'cm', 'm']
    concunits   = ['ug/L', 'mg/L', 'g/L', 'umol/L', 'mmol/L', 'mol/L']
    timeunits   = ['s', 'min', 'hr', 'day', 'yr']
    diffunits   = ['cm^2/s', 'cm^2/yr']

    system.lengthunits = [u'\u03BCm', 'cm', 'm']
    system.concunits   = [u'\u03BCg/L', 'mg/L', 'g/L', u'\u03BCmol/L', 'mmol/L', 'mol/L']
    system.timeunits   = ['s', 'min', 'hr', 'day', 'yr']
    system.diffunits   = [u'cm\u00B2/s', u'cm\u00B2/yr']

    system.lengthunit     = system.lengthunits[lengthunits.index(str(content[row_units + 1][2]))]
    system.concunit       = system.concunits[concunits.index(str(content[row_units + 2][2]))]
    system.timeunit       = system.timeunits[timeunits.index(str(content[row_units + 3][2]))]
    system.diffunit       = system.diffunits[diffunits.index(str(content[row_units + 4][2]))]

    system.nchemicals = int(content[row_chemicals][1])
    system.chemicals  = []
    for i in range(system.nchemicals):
        system.chemicals.append(Chemical(i+1, 1))
        system.chemicals[-1].name    = content[row_chemicals + 2 + i][1]
        system.chemicals[-1].formula = formula_converter(content[row_chemicals + 2 + i][2])
        system.chemicals[-1].MW      = float(content[row_chemicals + 2 + i][3])
        system.chemicals[-1].temp    = float(content[row_chemicals + 2 + i][4])
        system.chemicals[-1].Dw      = float(content[row_chemicals + 2 + i][5])
        system.chemicals[-1].Koc     = float(content[row_chemicals + 2 + i][6])
        system.chemicals[-1].Kdoc    = float(content[row_chemicals + 2 + i][7])
        system.chemicals[-1].Ref     = (content[row_chemicals + 2 + i][8])
        system.chemicals[-1].Kf      = 0
        system.chemicals[-1].N       = 0

    system.nmatrices  = int(content[row_matrices][1])
    system.matrices   = []
    current_row_components         = row_components + 2
    for i in range(system.nmatrices):
        system.matrices.append(Matrix(i+1))
        system.matrices[-1].name    = content[row_matrices + 2 + i][1]
        system.matrices[-1].e       = float(content[row_matrices + 2 + i][3])
        system.matrices[-1].rho     = float(content[row_matrices + 2 + i][4])
        system.matrices[-1].foc     = float(content[row_matrices + 2 + i][5])
        system.matrices[-1].model   = content[row_matrices + 2 + i][6]
        system.matrices[-1].components = []
        for j in range (int(content[row_matrices + 2 + i][2])):
            system.matrices[-1].components.append(MatrixComponent(j))
            system.matrices[-1].components[-1].name     =       content[current_row_components+ j][2]
            system.matrices[-1].components[-1].e        = float(content[current_row_components+ j][4])
            system.matrices[-1].components[-1].rho      = float(content[current_row_components+ j][5])
            system.matrices[-1].components[-1].foc      = float(content[current_row_components+ j][6])
            system.matrices[-1].components[-1].sorp     = 'Linear--Kd specified'
            system.matrices[-1].components[-1].tort     = 'None'
            if content[row_components+1][3] == 'Weight fraction':
                system.matrices[-1].components[-1].mfraction = float(content[current_row_components+ j][3])
                system.matrices[-1].components[-1].fraction = system.matrices[-1].components[-1].mfraction/system.matrices[-1].components[-1].rho * system.matrices[-1].rho
            else:
                system.matrices[-1].components[-1].fraction  = float(content[current_row_components+ j][3])
                system.matrices[-1].components[-1].mfraction = system.matrices[-1].components[-1].fraction * system.matrices[-1].components[-1].rho / system.matrices[-1].rho

        current_row_components = current_row_components + int(content[row_matrices + 2 + i][2])

    system.component_list      = []
    system.components          = []
    for matrix in system.matrices:
        for component in matrix.components:
            try: system.component_list.index(component.name)
            except:
                 system.components.append(component)
                 system.component_list.append(component.name)

    system.sorptions  = {}
    for i in range(len(system.components)):
        system.sorptions[system.components[i].name]={}
        for j in range(system.nchemicals):
            system.sorptions[system.components[i].name][system.chemicals[j].name] = Sorption(system.components[i], system.chemicals[j])
            system.sorptions[system.components[i].name][system.chemicals[j].name].isotherm = content[row_sorptions + 2 + i * system.nchemicals+ j][3]
            system.sorptions[system.components[i].name][system.chemicals[j].name].kinetic  = content[row_sorptions + 2 + i * system.nchemicals+ j][4]
            if system.sorptions[system.components[i].name][system.chemicals[j].name].isotherm == 'Linear--Kd specified':
                system.sorptions[system.components[i].name][system.chemicals[j].name].K  = float(content[row_sorptions + 2 + i * system.nchemicals+ j][5])
            elif system.sorptions[system.components[i].name][system.chemicals[j].name].isotherm == 'Linear--Kocfoc':
                system.sorptions[system.components[i].name][system.chemicals[j].name].Koc  = float(content[row_sorptions + 2 + i * system.nchemicals+ j][5])
            elif system.sorptions[system.components[i].name][system.chemicals[j].name].isotherm == 'Freundlich':
                system.sorptions[system.components[i].name][system.chemicals[j].name].Kf  = float(content[row_sorptions + 2 + i * system.nchemicals+ j][5])
                system.sorptions[system.components[i].name][system.chemicals[j].name].N   = float(content[row_sorptions + 2 + i * system.nchemicals+ j][6])
            else:
                system.sorptions[system.components[i].name][system.chemicals[j].name].qmax  = float(content[row_sorptions + 2 + i * system.nchemicals+ j][5])
                system.sorptions[system.components[i].name][system.chemicals[j].name].b     = float(content[row_sorptions + 2 + i * system.nchemicals+ j][6])
            if system.sorptions[system.components[i].name][system.chemicals[j].name].kinetic == 'Transient':
                system.sorptions[system.components[i].name][system.chemicals[j].name].ksorp   = float(content[row_sorptions + 2 + i * system.nchemicals+ j][7])
                system.sorptions[system.components[i].name][system.chemicals[j].name].kdesorp = float(content[row_sorptions + 2 + i * system.nchemicals+ j][8])
                system.sorptions[system.components[i].name][system.chemicals[j].name].thalf   = float(content[row_sorptions + 2 + i * system.nchemicals+ j][9])

    matrix_list = [matrix.name for matrix in system.matrices]
    system.nlayers  = int(content[row_layers][1])
    system.layers   = []
    for i in range(system.nlayers):
        if content[row_layers + 2][1] == 'Deposition':  system.layers.append(Layer(i))
        else:                                           system.layers.append(Layer(i + 1))
        system.layers[i].name   = content[row_layers + 2 + i][1]
        system.layers[i].type   = content[row_layers + 2 + i][2]
        system.layers[i].h      = float(content[row_layers + 2 + i][3])
        system.layers[i].tort   = content[row_layers + 2 + i][4]
        system.layers[i].alpha  = float(content[row_layers + 2 + i][5])
        system.layers[i].doc    = float(content[row_layers + 2 + i][6])
        system.layers[i].type_index = matrix_list.index(system.layers[i].type)

    system.nreactions  = int(content[row_reactions][1])
    system.reactions   = []
    for i in range(system.nreactions):
        system.reactions.append(Reaction(i+1))
        system.reactions[i].name     = content[row_reactions + 1 + i * 7][2]
        system.reactions[i].model    = content[row_reactions + 1 + i * 7][3]
        system.reactions[i].equation = formula_converter(content[row_reactions + 1 + i * 7][4])
        system.reactions[i].reactants = []
        system.reactions[i].products  = []
        column = 2
        for j in range(int(content[row_reactions + 2 + i * 7][3])):
            system.reactions[i].reactants.append(Reactant(j+1))
            system.reactions[i].reactants[j].name   = content[row_reactions + 3 + i * 7][column]
            system.reactions[i].reactants[j].formula = formula_converter(content[row_reactions + 4 + i * 7][column])
            system.reactions[i].reactants[j].MW = float(content[row_reactions + 5 + i * 7][column])
            system.reactions[i].reactants[j].coef = float(content[row_reactions + 6 + i * 7][column])
            system.reactions[i].reactants[j].index = float(content[row_reactions + 7 + i * 7][column])
            column = column + 1
        column = column + 1

        for j in range(int(content[row_reactions + 2 + i * 7][column+1])):
            system.reactions[i].products.append(Product(j+1))
            system.reactions[i].products[j].name   = content[row_reactions + 3 + i * 7][column]
            system.reactions[i].products[j].formula = formula_converter(content[row_reactions + 4 + i * 7][column])
            system.reactions[i].products[j].MW = float(content[row_reactions + 5 + i * 7][column])
            system.reactions[i].products[j].coef = float(content[row_reactions + 6 + i * 7][column])
            column = column + 1

    system.coefficients   = {}
    for i in range(system.nlayers):
        system.coefficients[system.layers[i].name]={}
        for j in range(system.nreactions):
            system.coefficients[system.layers[i].name][system.reactions[j].name]=Coefficient(system.layers[i], system.reactions[j])
            system.coefficients[system.layers[i].name][system.reactions[j].name].lam = float(content[row_coefficients + 2 + i * system.nreactions + j][3])

    system.adv    = content[row_system + 1][2]
    system.Vdar   = float(content[row_system + 2][2])
    system.Vtidal = float(content[row_system + 3][2])
    system.ptidal = float(content[row_system + 4][2])
    system.bio    = content[row_system + 5][2]
    system.hbio   = float(content[row_system + 6][2])
    system.sigma  = float(content[row_system + 7][2])
    system.Dbiop  = float(content[row_system + 8][2])
    system.Dbiopw = float(content[row_system + 9][2])
    system.con    = content[row_system + 10][2]
    system.hcon   = float(content[row_system + 11][2])
    system.t90    = float(content[row_system + 12][2])

    system.biomix = 0
    if system.bio <> 'None' and system.Dbiop <> 0:
        if system.layers[0].name == 'Deposition':
            system.biomix = 1
        elif system.hbio >  system.layers[0].h:
            system.biomix = 1

    system.BCs = {}
    row = row_conditions + 2
    system.topBCtype   = content[row][2]
    column = 4
    for chemical in system.chemicals:
        system.BCs[chemical.name] = BC(chemical.name)
        if system.topBCtype == 'Fixed Concentration':
            system.BCs[chemical.name].Co = float(content[row][column])
        elif system.topBCtype == 'Mass transfer':
            system.BCs[chemical.name].k   = float(content[row  ][column])
            system.BCs[chemical.name].Cw  = float(content[row+1][column])
        else:
            system.BCs[chemical.name].k   = float(content[row  ][column])
            system.BCs[chemical.name].tau = float(content[row+1][column])
        column = column + 1

    system.ICs      = {}
    system.SolidICs = {}

    if system.topBCtype == 'Fixed Concentration': row = row_conditions + 4
    else:                                         row = row_conditions + 5
    for layer in system.layers:
        layer.ICtype = content[row][2]
        system.ICs[layer.name] = {}
        system.SolidICs[layer.name] = {}
        column = 4
        if layer.ICtype == 'Uniform':
            for chemical in system.chemicals:
                system.ICs[layer.name][chemical.name] = IC(layer.name, chemical.name)
                system.ICs[layer.name][chemical.name].uniform = float(content[row][column])
                column = column + 1

            for component in system.matrices[layer.type_index].components:
                system.SolidICs[layer.name][component.name] = {}
                for chemical in system.chemicals:
                    system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(component.name,layer.name,chemical.name)
                    if system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                        system.SolidICs[layer.name][component.name][chemical.name].uniform = float(content[row][column])
                        column = column + 1
            row = row + 1
        else:
            for chemical in system.chemicals:
                system.ICs[layer.name][chemical.name] = IC(layer.name, chemical.name)
                system.ICs[layer.name][chemical.name].top = float(content[row][column])
                system.ICs[layer.name][chemical.name].bot = float(content[row+1][column])
                column = column + 1
                for component in system.matrices[layer.type_index].components:
                    system.SolidICs[layer.name][component.name] = {}
                    for chemical in system.chemicals:
                        system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(component.name,layer.name,chemical.name)
                        if system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            system.SolidICs[layer.name][component.name][chemical.name].top = float(content[row][column])
                            system.SolidICs[layer.name][component.name][chemical.name].bot = float(content[row + 1][column])
                            column = column + 1

            row = row + 2

    row     = row + 1
    column  = 4
    system.botBCtype   = content[row][2]
    for chemical in system.chemicals:
        if system.botBCtype == 'Fixed Concentration' or system.botBCtype == 'Flux-matching':
            system.BCs[chemical.name].Cb = float(content[row][column])
        column = column + 1

    system.tfinal       = float(content[row_solver + 1][2])
    system.outputsteps  = int(content[row_solver + 2][2])
    system.timeoption   = content[row_solver + 3][2]
    system.discrete     = content[row_solver + 4][2]
    system.ptotal       = content[row_solver + 5][2]
    system.delt         = float(content[row_solver + 6][2])
    system.ptype        = content[row_solver + 7][2]
    system.tvariable    = content[row_solver + 8][2]
    system.delz         = []
    system.players      = []
    column = 2
    for layer in system.layers:
        system.delz.append(float(content[row_solver + 10][column]))
        system.players.append(float(content[row_solver + 11][column]))
        column = column + 1
    system.tidalsteps    = float(content[row_solver + 12][2])
    system.nonlinear     = content[row_solver + 13][2]
    system.nlerror       = float(content[row_solver + 14][2])
    system.depgrid      = float(content[row_solver + 15][2])
    system.averageoption= content[row_solver + 16][2]
    system.massoption   = content[row_solver + 17][2]

    system.taucoefs           = {}
    system.taucoefs['Q']      = float(content[row_ow + 1][2])
    system.taucoefs['V']      = float(content[row_ow + 2][2])
    system.taucoefs['h']      = float(content[row_ow + 3][2])
    system.taucoefs['DOC']    = float(content[row_ow + 4][2])
    system.taucoefs['Qevap']  = float(content[row_ow + 5][2])
    for n in range(system.nchemicals):
        system.BCs[system.chemicals[n].name].kdecay = float(content[row_ow + 7][2 + n])
        system.BCs[system.chemicals[n].name].kevap  = float(content[row_ow + 8][2 + n])

    # Fill the rest of the system properties using default properties
    system.bltype             = 'River'
    system.blcoefs            = {}
    system.blcoefs['vx']      = 1.
    system.blcoefs['n']       = 0.02
    system.blcoefs['hriver']  = 5.
    system.blcoefs['rh']      = 5.
    system.blcoefs['nu']      = 1e-6

    system.blcoefs['rhoair']  = 1.
    system.blcoefs['rhowater']= 1000.
    system.blcoefs['vwind']   = 5.
    system.blcoefs['hlake']   = 10.
    system.blcoefs['llake']   = 1000.

    if system.nlayers > 0:
        if system.layers[0].name == 'Deposition':
            system.dep  = 'Deposition'
            system.Vdep = system.layers[0].h
        else:
            system.dep  = 0
            system.Vdep = 0

    return system