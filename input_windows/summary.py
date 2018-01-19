#! /usr/bin/env python
#
#This is the Tkinter GUI for the summary window for the simulation.

import tkMessageBox as tkmb, cPickle as pickle, sys, os
import _winreg as wreg

if sys.path[0][-3:] == 'zip': 
    os.chdir(sys.path[0][:-12])
    path = sys.path[0][:-12]
else: path = sys.path[0]

CapSimReg = wreg.ConnectRegistry(None, wreg.HKEY_CURRENT_USER)
CapSimKey = wreg.OpenKey(CapSimReg, r'Software\CapSim')
Filepath =wreg.QueryValueEx(CapSimKey, 'FilePath')[0]
CapSimKey.Close()

from Tkinter              import Tk, Toplevel, Canvas, Frame, Label, Entry, Text, Button, Scrollbar, OptionMenu, StringVar, DoubleVar, IntVar, FLAT, RAISED, Checkbutton
from capsim_object_types  import CapSimWindow, Chemical, Matrix, MatrixComponent, Layer, Sorption, Coefficient, IC, BC, SolidIC
from chemicalproperties   import ChemicalProperties
from matrixproperties     import MatrixProperties
from layerproperties      import LayerProperties
from sorptionproperties   import SorptionProperties
from reactionproperties   import ReactionProperties
from reactioncoefficients import ReactionCoefficients
from systemproperties     import SystemProperties
from layerconditions      import LayerConditions
from solidlayerconditions import SolidLayerConditions
from solveroptions        import SolverOptions
from inputoptions         import InputOptions

from capsim_functions     import text_converter


class Summary:
    """Summary window of all the input file simulation."""

    def __init__(self, master, system, database, materials):
        """The parameters that are to be displayed by the GUI are defined here
        first, then set by the main program."""

        self.system    = system           #stores the system data
        self.version   = system.version   #CapSim version
        self.fonttype  = system.fonttype  #font for summary window
        self.database  = database
        self.materials = materials
        
        self.systemrow = max(14, len(self.system.layers) + 6)
        self.chemrow   = max(10, len(self.system.chemicals) + 2)
        
        self.master    = master
        self.bframe    = Frame(master.tframe)
        self.frame     = Frame(master.frame)
        self.tframe    = Frame(master.bframe)

        self.bgcolor   = self.frame.cget('bg')
        self.top       = None             #existence of toplevel widget flag

    def make_widgets(self):
        """Make and grid all the widgets for the summary GUI window."""
        self.frame_width = 80

        self.intro        = Label(self.frame, text = 'The following summarizes the information you have provided for ' +
                                  'this simulation. Please verify that it is correct and update as necessary.\n', justify = 'left')

        self.rowcolumn     = Label(self.frame, text = '', font = 'courier 10', width = 2)
        self.titelcolumn   = Label(self.frame, text = '', font = 'courier 10', width = 10)
        self.framecolumn   = Label(self.frame, text = '', font = 'courier 10', width = self.frame_width)
        self.blank2column  = Label(self.frame, text = '', font = 'courier 10', width = 2)
        self.buttoncolumn  = Label(self.frame, text = '', font = 'courier 10', width = 18)
        self.blank3column  = Label(self.frame, text = '', font = 'courier 10', width = 2)

        self.chemicallabel = Label(self.frame, text = 'Chemicals:')
        self.reactionlabel = Label(self.frame, text = 'Reactions:')
        self.layerlabel    = Label(self.frame, text = 'Layers:')
        self.systemlabel   = Label(self.frame, text = 'System:')
        self.solverlabel   = Label(self.frame, text = 'Solver:')

        self.chemicalbutton    = Button(self.frame, text = 'Edit Chemical Properties',   command = self.edit_chemicalproperties)
        self.reactionbutton    = Button(self.frame, text = 'Edit Reaction Properties',   command = self.edit_reactionproperties)
        self.matrixbutton      = Button(self.frame, text = 'Edit Material Properties',   command = self.edit_matrixproperties)
        self.sorptionbutton    = Button(self.frame, text = 'Edit Sorption Properties',   command = self.edit_sorptionproperties)
        self.layerbutton       = Button(self.frame, text = 'Edit Layer Properties',      command = self.edit_layerproperties)
        self.coefficientbutton = Button(self.frame, text = 'Edit Reaction Coefficients', command = self.edit_reactioncoefficients)
        self.systembutton      = Button(self.frame, text = 'Edit System Parameters',     command = self.edit_systemproperties)
        self.conditionbutton   = Button(self.frame, text = 'Edit Auxiliary Conditions',  command = self.edit_layerconditions)
        self.solverbutton      = Button(self.frame, text  = 'Edit Solver Options',       command = self.edit_solveroptions)
        self.inputbutton       = Button(self.frame, text  = 'Edit File Options',         command = self.edit_inputoptions)

        row_number = 21
        self.rowlabels          = []
        for i in range(row_number):
            self.rowlabels.append(Label(self.frame, text = '', width = 1))

        row = 2
        for rowlabel in self.rowlabels:
            rowlabel.grid(row = row,  column = 0, sticky = 'WE', padx = 1, pady = 1)
            row = row + 1

        self.intro.grid(row = 0, columnspan = 6, pady = 2, sticky = 'W', padx = 10)
        
        self.rowcolumn.grid(     row = 1,  column = 0,  pady = 1, sticky = 'WE', padx = 1)
        self.titelcolumn.grid(   row = 1,  column = 1,  pady = 1, sticky = 'WE', padx = 1)
        self.framecolumn.grid(   row = 1,  column = 2,  pady = 1, sticky = 'WE', padx = 1)
        self.blank2column.grid(  row = 1,  column = 3,  pady = 1, sticky = 'WE', padx = 1)
        self.buttoncolumn.grid(  row = 1,  column = 4,  pady = 1, sticky = 'WE', padx = 1)
        self.blank3column.grid(  row = 1,  column = 5,  pady = 1, sticky = 'WE', padx = 1)

        self.chemical_row = 2
        self.reaction_row = 3
        self.layer_row    = 8
        self.system_row   = 4
        self.solver_row   = 4

        row = 2
        self.chemicallabel.grid(        row = row,      column = 1,  sticky = 'WE', padx = 1,  pady = 1 )
        self.chemicalbutton.grid(       row = row,      column = 4,  sticky = 'WE', padx = 1)
        row = row + self.chemical_row
        self.reactionlabel.grid(        row = row,      column = 1,  sticky = 'WE', padx = 1,  pady = 1)
        self.reactionbutton.grid(       row = row,      column = 4,  sticky = 'WE', padx = 1)
        row = row + self.reaction_row
        self.layerlabel.grid(           row = row,      column = 1,  sticky = 'WE', padx = 1,  pady = 1)
        self.matrixbutton.grid(         row = row    ,  column = 4,  sticky = 'WE', padx = 1)
        self.sorptionbutton.grid(       row = row + 2,  column = 4,  sticky = 'WE', padx = 1)
        self.layerbutton.grid(          row = row + 4,  column = 4,  sticky = 'WE', padx = 1)
        self.coefficientbutton.grid(    row = row + 6,  column = 4,  sticky = 'WE', padx = 1)
        row = row + self.layer_row
        self.systemlabel.grid(          row = row,      column = 1,  sticky = 'WE', padx = 1,  pady = 1)
        self.systembutton.grid(         row = row    ,  column = 4,  sticky = 'WE', padx = 1)
        self.conditionbutton.grid(      row = row + 2,  column = 4,  sticky = 'WE', padx = 1)

        row = row + self.system_row
        self.solverlabel.grid(          row = row,      column = 1,  sticky = 'WE', padx = 1,  pady = 1)
        self.solverbutton.grid(         row = row,      column = 4,  sticky = 'WE', padx = 1)
        self.inputbutton.grid(          row = row + 2,  column = 4,  sticky = 'WE', padx = 1)

        #bind the "Return" key to the appropriate methods (listed next)

        self.chemicalbutton.bind('<Return>',    self.edit_chemicalproperties)
        self.reactionbutton.bind('<Return>',    self.edit_reactionproperties)
        self.matrixbutton.bind('<Return>',      self.edit_matrixproperties)
        self.sorptionbutton.bind('<Return>',    self.edit_sorptionproperties)
        self.layerbutton.bind('<Return>',       self.edit_layerproperties)
        self.coefficientbutton.bind('<Return>', self.edit_reactioncoefficients)
        self.systembutton.bind('<Return>',      self.edit_systemproperties)
        self.conditionbutton.bind('<Return>',   self.edit_layerconditions)
        self.solverbutton.bind('<Return>',      self.edit_solveroptions)

        self.updatesummary()

    def updatesummary(self):

        self.update_chemicals_widgets(row = 2)
        self.update_reactions_widgets(row = 2 + self.chemical_row)
        self.update_layers_widgets(row = 2 + self.chemical_row+ self.reaction_row)
        self.update_system_widgets(row = 2 + self.chemical_row+ self.reaction_row+ self.layer_row)
        self.update_solver_widgets(row = 2 + self.chemical_row+ self.reaction_row+ self.layer_row+ self.system_row)

        self.focusbutton = None
        self.master.geometry()
        self.master.center()

    def edit_chemicalproperties(self, event = None):
        """Makes a window to edit simulation parameters."""

        if self.top is None:

            chemicals = [chemical.copy() for chemical in self.system.chemicals]
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(ChemicalProperties(self.top, self.system, self.database))
            self.top.mainloop()

            if self.top is not None:

                self.system.get_chemicalproperties(self.top.window)
                for chemical in self.system.chemicals: chemical.remove_chemicalwidgets()
                self.top.destroy()
                self.top = None
                
            #update the summary screen
            #need to update reactions, sorptions and auxiliary conditions if the chemicals of layers changes
            property_check = 0
            
            add_list    = []
            remove_list = []
            
            for chemical in chemicals: remove_list.append(chemical.number)
            for new_chemical in self.system.chemicals: add_list.append(new_chemical.number)
            
            for chemical in chemicals:
                for new_chemical in self.system.chemicals:
                    if new_chemical.name == chemical.name:
                        remove_list.remove(chemical.number)
                        add_list.remove(new_chemical.number)

                        if chemical.number != new_chemical.number:  property_check = 1
                        if chemical.temp   != new_chemical.temp:    property_check = 1
                        if chemical.Dw     != new_chemical.Dw:      property_check = 1
                        if chemical.Koc    != new_chemical.Koc:     property_check = 1
                        if chemical.Kdoc   != new_chemical.Kdoc:    property_check = 1
                        
            for add_chemical in add_list:
                chemical = self.system.chemicals[add_chemical - 1]
                for component in self.system.components:
                    self.system.sorptions[component.name][chemical.name] = Sorption(component, chemical)

                self.system.BCs[chemical.name] = BC(chemical.name)
                for layer in self.system.layers:
                    self.system.ICs[layer.name][chemical.name] = IC(layer.name, chemical.name)
                    for component in self.system.matrices[layer.type_index].components:
                        self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)

            for remove_chemical in remove_list:
                chemical = chemicals[remove_chemical - 1]
                for component in self.system.components:
                    del self.system.sorptions[component.name][chemical.name]
                reactions_remove = []
                for reaction in self.system.reactions:
                    reactant_list = [reactant.name for reactant in reaction.reactants]
                    product_list  = [product.name for product in reaction.products]
                    if (reactant_list.count(chemical.name)+product_list.count(chemical.name)) > 0:
                        reactions_remove.append(self.system.reactions.index(reaction))
                        for layer in self.system.layers:
                            del self.system.coefficients[layer.name][reaction.name]
                if len(reactions_remove) > 0:
                    reactions_remove.reverse()
                    for reaction_index in reactions_remove:
                        self.system.reactions.remove(self.system.reactions[reaction_index])

                del self.system.BCs[chemical.name]

                for layer in self.system.layers:
                    del self.system.ICs[layer.name][chemical.name]
                    for component in self.system.matrices[layer.type_index].components:
                        del self.system.SolidICs[layer.name][component.name][chemical.name]

            if len(add_list) != 0:
                self.edit_sorptionproperties(flag = 1)
                self.edit_layerconditions()

            elif property_check != 0:
                self.edit_sorptionproperties(flag = 1)

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

    def edit_reactionproperties(self, flag = None, event = None):
        """Makes a window to edit flow and system parameters."""

        if self.top is None:
            
            reactions = [reaction.copy() for reaction in self.system.reactions] 
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(ReactionProperties(self.top, self.system))
            self.top.mainloop()

            if self.top is not None: 
                
                self.system.get_reactionproperties(self.top.window)
                for reaction in self.system.reactions:
                    reaction.remove_propertieswidgets()
                self.top.destroy()
                self.top = None
                              
            add_list    = []
            remove_list = []
            
            for reaction in reactions: remove_list.append(reaction.number)
            for new_reaction in self.system.reactions: add_list.append(new_reaction.number)

            for reaction in reactions:
                for new_reaction in self.system.reactions:
                    reaction_check = 0
                    if new_reaction.name            != reaction.name:           reaction_check = 1
                    if new_reaction.equation        != reaction.equation:       reaction_check = 1
                    if new_reaction.model           != reaction.model:          reaction_check = 1
                    if len(new_reaction.reactants)  != len(reaction.reactants): reaction_check = 1
                    else:
                        for r_index in range(len(new_reaction.reactants)):
                            if new_reaction.reactants[r_index].coef    != reaction.reactants[r_index].coef:     reaction_check = 1
                            if new_reaction.reactants[r_index].name    != reaction.reactants[r_index].name:     reaction_check = 1
                            if new_reaction.reactants[r_index].formula != reaction.reactants[r_index].formula:  reaction_check = 1
                            if new_reaction.reactants[r_index].index   != reaction.reactants[r_index].index:    reaction_check = 1
                    if len(new_reaction.products)   != len(reaction.products):  reaction_check = 1
                    else:
                        for r_index in range(len(new_reaction.products)):
                            if new_reaction.products[r_index].coef     != reaction.products[r_index].coef:      reaction_check = 1
                            if new_reaction.products[r_index].name     != reaction.products[r_index].name:      reaction_check = 1
                            if new_reaction.products[r_index].formula  != reaction.products[r_index].formula:   reaction_check = 1
                                        
                    if reaction_check == 0:
                        remove_list.remove(reaction.number)
                        add_list.remove(new_reaction.number)
                        
            for remove_reaction in remove_list:
                reaction = reactions[remove_reaction - 1]
                for layer in self.system.layers:
                    del self.system.coefficients[layer.name][reaction.name]
                    
            for add_reaction in add_list:
                reaction = self.system.reactions[add_reaction - 1]
                for layer in self.system.layers:
                    self.system.coefficients[layer.name][reaction.name] = Coefficient(layer, reaction)

            if len(add_list) != 0 or len(remove_list) != 0 or flag == 1:

                self.edit_reactioncoefficients()

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

    def edit_matrixproperties(self, event = None):

        
        if self.top is None:

            matrices   = [matrix.copy()     for matrix      in self.system.matrices]
            components = [component.copy()  for component   in self.system.components]
            component_list = self.system.component_list

            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(MatrixProperties(self.top, self.system, self.materials))
            self.top.mainloop()

            if self.top is not None: #update the screen
                
                self.system.get_matricesproperties(self.top.window)
                for matrix in self.system.matrices:
                    matrix.remove_propertieswidgets()
                self.top.destroy()
                self.top = None
                
            property_check = 0
            layer_check    = 0

            layer_list = []
            matrix_list = [matrix.name for matrix in self.system.matrices]
            for layer in self.system.layers:
                if matrix_list.count(layer.type) == 0:
                    layer.type = matrix_list[0]
                    layer_list.append(self.system.layers.index(layer))
                    layer_check = 1
                layer.type_index = matrix_list.index(layer.type)

            add_list    = []
            remove_list = []

            for component_name in component_list:                remove_list.append(component_list.index(component_name))
            for component_name in self.system.component_list:    add_list.append(self.system.component_list.index(component_name))

            indeces     = []
            new_indeces = []
            
            for index in remove_list:
                for new_index in add_list:
                    if component_list[index] == self.system.component_list[new_index]:
                        indeces.append(index)
                        new_indeces.append(new_index)
            for index in indeces:
                remove_list.remove(index)
            for new_index in new_indeces:
                add_list.remove(new_index)

            for remove_index in remove_list:
                component = components[remove_index]
                del self.system.sorptions[component.name]
                for layer in self.system.layers:
                    for layer_component in matrices[layer.type_index].components:
                        if layer_component.name == component.name:
                            del self.system.SolidICs[layer.name][component.name]


            for add_index in add_list:
                component = self.system.components[add_index]
                self.system.sorptions[component.name] = {}
                for chemical in self.system.chemicals:
                    self.system.sorptions[component.name][chemical.name] = Sorption(component, chemical)

                for layer in self.system.layers:
                    for layer_component in self.system.matrices[layer.type_index].components:
                        if layer_component.name == component.name:
                            self.system.SolidICs[layer.name][component.name] = {}
                            for chemical in self.system.chemicals:
                                self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)

            for layer_index in layer_list:
                for layer_component in matrices[self.system.layers[layer_index].type_index].components:
                    self.system.SolidICs[layer.name][layer_component.name] = {}
                    for chemical in self.system.chemicals:
                        self.system.SolidICs[layer.name][layer_component.name][chemical.name] = SolidIC(layer.name, layer_component.name, chemical.name)

            for component in self.system.components:
                for chemical in self.system.chemicals:
                    self.system.sorptions[component.name][chemical.name].matrix = component
                    self.system.sorptions[component.name][chemical.name].foc    = component.foc
                    self.system.sorptions[component.name][chemical.name].update_material()

            if len(add_list) != 0 or len(remove_list) != 0:

                self.edit_sorptionproperties(flag = 1)
                
            if layer_check == 1:

                self.edit_layerproperties()

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

    def edit_layerproperties(self, flag = None, event = None):
        """Makes a window to get some basic parameters for each layer."""

        
        if self.top is None:
            
            layers          = [layer.copy() for layer in self.system.layers]
            players         = [player for player in self.system.players]
            coefficients    = self.system.coefficients
            ICs             = self.system.ICs
            SolidICs        = self.system.SolidICs
            dep             = self.system.dep

            num_record = [j for j in range(self.system.nlayers)]
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(LayerProperties(self.top, self.system, editflag = 1))
            self.top.mainloop()

            if self.top is not None:
                self.system.get_layerproperties(self.top.window)
                num_record = self.top.window.num_record
                for layer in self.system.layers:
                    layer.remove_propertywidgets()
                self.top.destroy()
                self.top = None
                
            #update the summary screen
            self.updatesummary()
            self.master.geometry()

            #Name_check = 0
            type_check  = 0
            num_check   = 0
            h_check     = 0
            dep_check   = 0

            if len(num_record) != len(layers): num_check = 1
            else:
                for j in range(len(layers)):
                    if num_record[j] != j:
                        num_check = 1
                    else:
                         if layers[j].h <> self.system.layers[j].h: h_check = 1

            if num_check == 1:

                self.system.coefficients = {}
                self.system.ICs          = {}
                self.system.SolidICs     = {}
                self.system.players      = []

                for layer in self.system.layers:
                    if num_record[self.system.layers.index(layer)] >= 0:
                        self.system.coefficients[layer.name] = coefficients[layers[num_record[self.system.layers.index(layer)]].name]
                        self.system.ICs[layer.name]          = ICs[layers[num_record[self.system.layers.index(layer)]].name]
                        self.system.SolidICs[layer.name]     = {}
                        for component in self.system.matrices[layer.type_index].components:
                            if [M_component.name for M_component in self.system.matrices[layers[num_record[self.system.layers.index(layer)]].type_index].components].count(component.name) > 0:
                                self.system.SolidICs[layer.name][component.name]     = SolidICs[layers[num_record[self.system.layers.index(layer)]].name][component.name].copy()
                            else:
                                self.system.SolidICs[layer.name][component.name]     = {}
                                for chemical in self.system.chemicals:
                                    self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)
                        self.system.players.append(players[num_record[self.system.layers.index(layer)]])
                    else:
                        self.system.coefficients[layer.name]    = {}
                        self.system.ICs[layer.name]             = {}
                        self.system.SolidICs[layer.name]        = {}
                        for reaction in self.system.reactions:
                            self.system.coefficients[layer.name][reaction.name] = Coefficient(layer, reaction)
                        for chemical in self.system.chemicals:
                            self.system.ICs[layer.name][chemical.name]          = IC(layer.name, chemical.name)
                        for component in self.system.matrices[layer.type_index].components:
                            self.system.SolidICs[layer.name][component.name] = {}
                            for chemical in self.system.chemicals:
                                self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)
                        self.system.players.append(10)
            else:
                self.system.SolidICs     = {}
                for layer in self.system.layers:
                    #component_list = [M_component.name for M_component in self.system.matrices[layers[self.system.layers.index(layer)].type_index].components]
                    self.system.SolidICs[layer.name] = {}
                    for component in self.system.matrices[layer.type_index].components:
                        try:
                            self.system.SolidICs[layer.name][component.name]     = SolidICs[layers[num_record[self.system.layers.index(layer)]].name][component.name].copy()
                        except:
                            self.system.SolidICs[layer.name][component.name]     = {}
                            for chemical in self.system.chemicals:
                                self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)
                            type_check = 1

            if self.system.bio <> 'None' and (self.system.layers[0].name == 'Deposition' or self.system.layers[0].h < self.system.hbio or self.system.bio == 'Depth-dependent') and self.system.Dbiop> 0 and self.system.nlayers > 1:     self.system.biomix = 1
            else:                                                                                                                                                                                                                         self.system.biomix = 0


            if dep != self.system.dep: dep_check = 1

            if num_check == 1:
                self.edit_reactioncoefficients()
                self.edit_layerconditions()
                self.edit_solveroptions(editflag = 1)
            else:
                if type_check == 1:                 self.edit_layerconditions()
                if h_check == 1 or dep_check == 1:  self.edit_solveroptions(editflag = 1)

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()
        
    def edit_sorptionproperties(self, flag = None, event = None):
        
        #coefficients = self.system.coefficients
        old_dynamic_flags = []
        dynamic_flags     = []


        for component in self.system.components:
            for chemical in self.system.chemicals:
                if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                    old_dynamic_flags.append(1)
                else:
                    old_dynamic_flags.append(0)

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 1)
            self.top.make_window(SorptionProperties(self.top, self.system))
            self.top.mainloop()

            if self.top is not None: #update the screen
                
                self.system.get_sorptionproperties(self.top.window)
                for component in self.system.components:
                    component.remove_sorptionwidgets()
                for chemical in self.system.chemicals:
                    if chemical.soluable == 1:
                        chemical.remove_sorptionwidgets()
                for component in self.system.components:
                    for chemical in self.system.chemicals:
                        if chemical.soluable == 1:
                            self.system.sorptions[component.name][chemical.name].remove_propertieswidgets()
                        if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            dynamic_flags.append(1)
                        else:
                            dynamic_flags.append(0)
                self.top.destroy()
                self.top = None

            dynamic_check = 0
            for flag_index in range (len(dynamic_flags)):
                if dynamic_flags[flag_index] != old_dynamic_flags[flag_index]:
                    dynamic_check = 1

            if dynamic_check == 1:
                self.system.SolidICs = {}
                for layer in self.system.layers:
                    self.system.SolidICs[layer.name] = {}
                    for component in self.system.matrices[layer.type_index].components:
                        self.system.SolidICs[layer.name][component.name] = {}
                        for chemical in self.system.chemicals:
                            self.system.SolidICs[layer.name][component.name][chemical.name] = SolidIC(layer.name, component.name, chemical.name)

                self.edit_layerconditions()
                self.edit_solveroptions()

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

        
    def edit_reactioncoefficients(self, event = None):
        
        #coefficients = self.system.coefficients

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(ReactionCoefficients(self.top, self.system, editflag = 1))
            self.top.mainloop()

            if self.top is not None: #update the screen
                
                self.system.get_reactioncoefficients(self.top.window)
                for layer in self.system.layers:
                    layer.remove_reactionwidgets()
                    for reaction in self.system.reactions:
                        self.system.coefficients[layer.name][reaction.name].remove_propertieswidgets()
                self.top.destroy()
                self.top = None

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()
        

    def edit_systemproperties(self, flag = None, event = None):
        """Makes a window to get models for the layers."""

        adv     = self.system.adv
        bio     = self.system.bio
        con     = self.system.con
        hbio    = self.system.hbio
        ptidal  = self.system.ptidal


        if self.top is None:

            self.adv    = self.system.adv
            self.bio    = self.system.bio
            self.con    = self.system.con
            
            self.Vdar   = self.system.Vdar
            self.Vtidal = self.system.Vtidal
            self.ptidal = self.system.ptidal        
            self.hbio   = self.system.hbio
            self.Dbiop  = self.system.Dbiop
            self.Dbiopw = self.system.Dbiopw
            self.hcon   = self.system.hcon
            self.t90    = self.system.t90
            
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(SystemProperties(self.top, self.system))
            self.top.mainloop()

            if self.top is not None:
                self.system.get_systemproperties(self.top.window)
                self.top.destroy()
                self.top = None

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

        if self.system.adv != 'Period oscillation': self.system.averageoption = 'Instaneous'

        if adv != self.system.adv or self.Vdar != self.system.adv:  self.edit_layerconditions()

        if adv != self.system.adv or bio != self.system.bio or con != self.system.con or hbio != self.system.hbio or ptidal != self.system.ptidal:

            self.edit_solveroptions()


    def edit_layerconditions(self, event = None):
        
        
        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(LayerConditions(self.top, self.system))
            self.top.mainloop()

            if self.top is not None: #update the screen
                
                self.system.get_layerconditions(self.top.window)
                
                for layer in self.system.layers:
                    layer.get_layerconditions()
                    layer.remove_ICwidgets()
                    for chemical in self.system.chemicals:
                        if chemical.soluable == 1:
                            self.system.ICs[layer.name][chemical.name].remove_widgets()

                for chemical in self.system.chemicals:
                    if chemical.soluable == 1:
                        chemical.remove_ICwidgets()
                        self.system.BCs[chemical.name].remove_widgets()

                self.top.destroy()
                self.top = None

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

        solid_check = 0

        for chemical in self.system.chemicals:
            for component in self.system.components:
                if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                    solid_check = 1

        if solid_check > 0:

            if self.top is None:
                self.top = CapSimWindow(master = self.master, buttons = 3)
                self.top.make_window(SolidLayerConditions(self.top, self.system))
                self.top.mainloop()

                if self.top is not None: #update the screen

                    self.system.get_solidlayerconditions(self.top.window)

                    for layer in self.system.layers:
                        layer.remove_SolidICwidgets()
                        for component in self.system.matrices[layer.type_index].components:
                            component.remove_SolidICswidgets()
                            for chemical in self.system.chemicals:
                                self.system.SolidICs[layer.name][component.name][chemical.name].remove_widgets()

                    for chemical in self.system.chemicals:
                        try:    chemical.remove_ICwidgets()
                        except: pass

                    self.top.destroy()
                    self.top = None

                    self.updatesummary()
                    self.master.geometry()
                    self.master.center()


            else: self.master.open_toplevel()

    def edit_solveroptions(self, event = None, editflag = 0):
        """Makes a window to edit output options."""

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(SolverOptions(self.top, self.system, editflag))
            self.top.mainloop()

            if self.top is not None:
                self.system.get_solveroptions(self.top.window)
                self.top.destroy()
                self.top = None
            
            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

    def edit_inputoptions(self, event = None):
        """Makes a window to edit output options."""

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 3)
            self.top.make_window(InputOptions(self.top, self.system))
            self.top.mainloop()

            if self.top is not None:
                self.system.get_inputoptions(self.top.window)
                self.top.destroy()
                self.top = None

            self.updatesummary()
            self.master.geometry()
            self.master.center()

        else: self.master.open_toplevel()

    def update_chemicals_widgets(self, row):

        try:
            for chemicalnamelabel in self.chemicalnamelabels: chemicalnamelabel.grid_forget()
            self.chemical_frame.grid_forget()
            self.chemical_canvas.grid_forget()
            self.chemical_vscrollbar.grid_forget()
        except: pass

        self.chemical_frame      = Frame(self.frame)
        self.chemical_canvas     = Canvas(self.chemical_frame, height = 22 * self.chemical_row)

        self.chemical_frame.grid_rowconfigure(0, weight=1)
        self.chemical_frame.grid_columnconfigure(0, weight=1)

        self.chemical_vscrollbar = Scrollbar(self.chemical_frame)
        self.chemical_canvas.config(yscrollcommand = self.chemical_vscrollbar.set)
        self.chemical_vscrollbar.config(command = self.chemical_canvas.yview)

        self.chemical_canvas.grid(row     = 0, column = 0, sticky = 'NSWE')

        self.chemical_content   = Frame(self.chemical_canvas)

        self.chemicalnamelabels = [Label(self.chemical_content, text = chemical.name, width = 18) for chemical in self.system.chemicals]
        rows = int(self.system.nchemicals/4)
        lastrows = self.system.nchemicals-rows*4
        for i in range(rows):
            for column in range(4):
                self.chemicalnamelabels[i*4+column].grid(row = i, column = column, pady = 1)
        for column in range(lastrows):
            self.chemicalnamelabels[rows*4+column].grid(row = rows, column = column, pady = 1)

        if (self.system.nchemicals) > self.chemical_row*4:
            self.chemical_vscrollbar.grid(row = 0, column = 1, sticky = 'NS')
        self.chemical_canvas.create_window(0, 0, anchor = 'nw', window = self.chemical_content)
        self.chemical_canvas.config(scrollregion = (0,0,59,(rows+(lastrows > 0)*1)*23))

        self.chemical_frame.grid( row     = row, column = 2, sticky = 'NWE', rowspan = self.chemical_row)

        self.chemical_frame.update()


    def update_reactions_widgets(self, row):

        try:
            for i in range(len(self.system.reactions)):
                self.reactionnamelabels[i].grid_forget()
                self.reactionequalabels[i].grid_forget()
            self.reaction_frame.grid_forget()
            self.reaction_canvas.grid_forget()
            self.reaction_vscrollbar.grid_forget()
            self.reaction_hscrollbar.grid_forget()
        except: pass

        self.reaction_frame      = Frame(self.frame)
        self.reaction_canvas     = Canvas(self.reaction_frame, height = 22 * self.reaction_row)

        self.reaction_frame.grid_rowconfigure(0, weight=1)
        self.reaction_frame.grid_columnconfigure(0, weight=1)

        self.reaction_vscrollbar = Scrollbar(self.reaction_frame)
        self.reaction_canvas.config(yscrollcommand = self.reaction_vscrollbar.set)
        self.reaction_vscrollbar.config(command = self.reaction_canvas.yview)

        self.reaction_hscrollbar = Scrollbar(self.reaction_frame, orient = 'horizontal')
        self.reaction_canvas.config(xscrollcommand = self.reaction_hscrollbar.set)
        self.reaction_hscrollbar.config(command = self.reaction_canvas.xview)

        self.reaction_canvas.grid(row     = 0, column = 0, sticky = 'NSWE')

        self.reaction_content  = Frame(self.reaction_canvas)

        #self.reactionnumberlabels = [Label(self.reaction_content, text = reaction.number,   width = 3)                  for reaction in self.system.reactions]
        self.reactionnamelabels   = [Label(self.reaction_content, text = reaction.name,  width = 18, justify = 'left')              for reaction in self.system.reactions]
        self.reactionequalabels   = [Label(self.reaction_content, text = reaction.equation, font = 'Calibri 11')  for reaction in self.system.reactions]

        for i in range(len(self.system.reactions)):
            #self.reactionnumberlabels[i].grid(row = i, column = 0, pady = 1)
            self.reactionnamelabels[i].grid(row = i, column = 0, pady = 1, sticky = 'W')
            self.reactionequalabels[i].grid(row = i, column = 1, pady = 1)

        vscrollbar_check = 0
        vscrollbar_width = 0
        for i in range(len(self.reactionequalabels)):
            if self.reactionequalabels[i].winfo_reqwidth() + self.reactionnamelabels[i].winfo_reqwidth() + 16  > self.frame_width * 8:
                self.reaction_hscrollbar.grid(row = 1, column = 0, sticky = 'WE')
                vscrollbar_check = 1
            if (self.reactionequalabels[i].winfo_reqwidth() + self.reactionnamelabels[i].winfo_reqwidth()) > vscrollbar_width:
                vscrollbar_width = self.reactionequalabels[i].winfo_reqwidth() + self.reactionnamelabels[i].winfo_reqwidth()

        if (len(self.system.reactions) + vscrollbar_check) >= self.reaction_row:
            self.reaction_vscrollbar.grid(row = 0, column = 1, sticky = 'NS')

        self.reaction_canvas.create_window(0, 0, anchor = 'nw', window = self.reaction_content)
        self.reaction_canvas.config(scrollregion = (0,0,vscrollbar_width+16,sum([reactionequalabel.winfo_reqheight() for reactionequalabel in self.reactionequalabels])+16))

        self.reaction_frame.grid( row     = row, column = 2, sticky = 'NWE', rowspan = self.reaction_row)

        self.reaction_frame.update()

    def update_layers_widgets(self, row):

        try:
            for i in range(len(self.system.layers)):
                self.layernamelabels[i].grid_forget()
                self.layertypelabels[i].grid_forget()
                self.layerthicklabels[i].grid_forget()
            self.layer_frame.grid_forget()
            self.layer_canvas.grid_forget()
            self.layer_vscrollbar.grid_forget()
            self.layer_hscrollbar.grid_forget()

        except: pass

        self.layer_frame     = Frame(self.frame)
        self.layer_canvas    = Canvas(self.layer_frame, height = 22 * self.layer_row)

        self.layer_frame.grid_rowconfigure(0, weight=1)
        self.layer_frame.grid_columnconfigure(0, weight=1)

        self.layer_vscrollbar = Scrollbar(self.layer_frame)
        self.layer_canvas.config(yscrollcommand = self.layer_vscrollbar.set)
        self.layer_vscrollbar.config(command = self.layer_canvas.yview)

        self.layer_hscrollbar = Scrollbar(self.layer_frame, orient = 'horizontal')
        self.layer_canvas.config(xscrollcommand = self.layer_hscrollbar.set)
        self.layer_hscrollbar.config(command = self.layer_canvas.xview)

        self.layer_canvas.grid(row     = 0, column = 0, sticky = 'NSWE')

        self.layer_content  = Frame(self.layer_canvas)

        self.layernamelabels        = [Label(self.layer_content, text = layer.name,     width = 14)            for layer in self.system.layers]
        self.layertypelabels        = [Label(self.layer_content, text = layer.type)                            for layer in self.system.layers]
        self.layerthicklabels       = [Label(self.layer_content, text = str(layer.h) + ' ' + self.system.lengthunit, width = 15, justify = 'center') for layer in self.system.layers]

        if self.system.layers[0].name == 'Deposition':
            self.layerthicklabels[0] = Label(self.layer_content, text = str(self.system.layers[0].h) + ' ' + self.system.lengthunit + '/' + self.system.timeunit, width = 15, justify = 'center')

        for i in range(len(self.system.layers)):
            self.layernamelabels[i].grid(       row = i, column = 0)
            self.layertypelabels[i].grid(       row = i, column = 1)
            self.layerthicklabels[i].grid(      row = i, column = 2)

        vscrollbar_check = 0
        vscrollbar_width = 0
        for i in range(len(self.layernamelabels)):
            if self.layernamelabels[i].winfo_reqwidth() + self.layertypelabels[i].winfo_reqwidth()+ self.layerthicklabels[i].winfo_reqwidth()  > self.frame_width * 8:
                self.layer_hscrollbar.grid(row = 1, column = 0, sticky = 'WE')
                vscrollbar_check = 1
            if (self.layernamelabels[i].winfo_reqwidth() + self.layertypelabels[i].winfo_reqwidth()+self.layerthicklabels[i].winfo_reqwidth()) > vscrollbar_width:
                vscrollbar_width = self.layernamelabels[i].winfo_reqwidth() + self.layertypelabels[i].winfo_reqwidth() + self.layerthicklabels[i].winfo_reqwidth()

        if (len(self.system.layers) + vscrollbar_check) > self.layer_row:  self.layer_vscrollbar.grid(row = 0, column = 1, sticky = 'NS')

        self.layer_canvas.create_window(0, 0, anchor = 'nw', window = self.layer_content)
        self.layer_canvas.config(scrollregion = (0,0,vscrollbar_width,sum([layernamelabel.winfo_reqheight() for layernamelabel in self.layernamelabels])+16))

        self.layer_frame.grid( row     = row, column = 2, sticky = 'NWE', rowspan = self.layer_row)

        self.layer_frame.update()


    def update_system_widgets(self, row):

        try:
            self.system_frame.grid_forget()
            self.system_canvas.grid_forget()
            self.system_content.grid_forget()

            self.NoVdarwidget.grid_forget()
            self.UniVdarlabel.grid_forget()
            self.UniVdarwidget.grid_forget()
            self.TidalVdarlabel.grid_forget()
            self.TidalVdarwidget.grid_forget()
            self.ptidallabel.grid_forget()
            self.ptidalwidget.grid_forget()

            self.Nobiolabel.grid_forget()
            self.biolabel.grid_forget()
            self.tBClabel.grid_forget()
            self.tBCwidget.grid_forget()
            self.bBClabel.grid_forget()
            self.bBCwidget.grid_forget()

        except: pass

        self.system_frame     = Frame(self.frame)
        self.system_canvas    = Canvas(self.system_frame, height = 23 * self.system_row)

        self.system_frame.grid_rowconfigure(0, weight=1)
        self.system_frame.grid_columnconfigure(0, weight=1)

        self.system_canvas.grid(row     = 0, column = 0, sticky = 'NSWE')

        self.system_content  = Frame(self.system_canvas)

        self.blanklabel         = Label(self.system_content,  text = '       ')

        self.NoVdarwidget       = Label(self.system_content,  text = 'No advection flow')
        self.UniVdarlabel       = Label(self.system_content,  text = 'Darcy velocity:')
        self.UniVdarwidget      = Label(self.system_content,  text = '%.1f %s/%s' %(self.system.Vdar, self.system.lengthunit,self.system.timeunit))
        self.TidalVdarlabel     = Label(self.system_content,  text = 'Darcy velocity:')
        self.TidalVdarwidget    = Label(self.system_content,  text = '%.1f' %(self.system.Vdar) + u'\u00B1' + '%.1f %s/%s' %(self.system.Vtidal, self.system.lengthunit,self.system.timeunit))
        self.ptidallabel        = Label(self.system_content,  text = 'Oscillation period:', width = 25)
        self.ptidalwidget       = Label(self.system_content,  text = '%.4f %s' % (self.system.ptidal, self.system.timeunit))

        self.Nobiolabel         = Label(self.system_content,  text = 'Ignore the bioturbation')
        self.biolabel           = Label(self.system_content,  text = 'Model the bioturbation')

        self.tBClabel           = Label(self.system_content, text = 'Benthic surface boundary:  ')
        self.tBCwidget          = Label(self.system_content, text = str(self.system.topBCtype))
        
        self.bBClabel           = Label(self.system_content, text = 'Underlying sediment boundary:  ')
        self.bBCwidget          = Label(self.system_content, text = str(self.system.botBCtype))

        self.blanklabel.grid(row      = 0,      column = 0, sticky= 'W')

        if self.system.adv == 'None':
            
            self.NoVdarwidget.grid(row      = 0,      column = 1, sticky= 'W',    columnspan = 2)

        elif self.system.adv == 'Steady flow':
            
            self.UniVdarlabel.grid(row      = 0,      column = 1, sticky= 'W')
            self.UniVdarwidget.grid(row     = 0,      column = 2, sticky= 'W')

        elif self.system.adv == 'Period oscillation':
            
            self.TidalVdarlabel.grid(row    = 0,      column = 1, sticky= 'W')
            self.TidalVdarwidget.grid(row   = 0,      column = 2, sticky= 'W')

            self.ptidallabel.grid(row       = 0 ,     column = 3, sticky= 'W')
            self.ptidalwidget.grid(row      = 0 ,     column = 4, sticky= 'W')

        if self.system.bio == 'None': self.Nobiolabel.grid( row       = 1,  column = 1, sticky = 'W', columnspan = 2, pady = 1)
        else:                         self.biolabel.grid( row         = 1,  column = 1, sticky = 'W', columnspan = 2, pady = 1)

        self.tBClabel.grid( row             = 2  ,  column = 1, sticky = 'W', pady = 1)
        self.tBCwidget.grid( row            = 2  ,  column = 2, sticky = 'W', pady = 1)

        self.bBClabel.grid( row             = 3  ,  column = 1, sticky = 'W', pady = 1)
        self.bBCwidget.grid( row            = 3  ,  column = 2, sticky = 'W', pady = 1)
        
        self.system_canvas.create_window(0, 0, anchor = 'nw', window = self.system_content)

        self.system_frame.grid( row     = row, column = 2, sticky = 'WE', rowspan = self.system_row)

        self.system_frame.update()


    def update_solver_widgets(self, row):

        try:
            self.solver_frame.grid_forget()
            self.solver_canvas.grid_forget()
            self.solver_content.grid_forget()

        except: pass

        self.solver_frame     = Frame(self.frame, height = 23 *self.solver_row)
        self.solver_canvas    = Canvas(self.solver_frame, height = 23 *self.solver_row)

        self.solver_frame.grid_rowconfigure(0, weight=1)
        self.solver_frame.grid_columnconfigure(0, weight=1)

        self.solver_canvas.grid(row     = 0, column = 0, sticky = 'NSWE')

        self.solver_content  = Frame(self.solver_canvas)

        self.blanksolverlabel = Label(self.solver_content,  text = '       ')

        self.durationwidget  = Label(self.solver_content,  text  = 'Run a simulation of %.1f %s' %(self.system.tfinal, self.system.timeunit))
        self.nonlinearlabel  = Label(self.solver_content,  text  = 'Non-linear solver')
        self.nonlineariwdget = Label(self.solver_content,  text  = self.system.nonlinear)

        self.blanksolverlabel.grid(row    = 0,        column = 0,     sticky=  'W',  pady = 1)
        self.durationwidget.grid(  row    = 0,        column = 1,     sticky = 'W', pady = 1)
        self.nonlinearlabel.grid(  row    = 1,        column = 1,     sticky = 'W', pady = 1)
        self.nonlineariwdget.grid( row    = 1,        column = 2,     sticky = 'W', pady = 1)

        self.solver_canvas.create_window(0, 0, anchor = 'nw', window = self.solver_content)

        self.solver_frame.grid( row     = row, column = 2, sticky = 'WE', rowspan = self.solver_row)

        self.solver_frame.update()

    def make_batch_file(self, Filepath):

        lengthunits = ["'um'", "'cm'", "'m'"]
        concunits   = ["'ug/L'", "'mg/L'", "'g/L'", "'umol/L'", "'mmol/L'", "'mol/L'"]
        timeunits   = ["'s'", "'min'", "'hr'", "'day'", "'yr'"]
        diffunits   = ["'cm^2/s'", "'cm^2/yr'"]

        lengthunit  = lengthunits[self.system.lengthunits.index(self.system.lengthunit)]
        concunit    = concunits[self.system.concunits.index(self.system.concunit)]
        timeunit    = timeunits[self.system.timeunits.index(self.system.timeunit)]
        diffunit    = diffunits[self.system.diffunits.index(self.system.diffunit)]

        file = open(Filepath + "/batch_files/"+ self.system.batchfilename + ".txt", "w")            # Create a csv file

        file.write("# Batch function \n")
        file.write("\n")

        file.write("# CapSim," + self.system.version + "\n\n")
        file.write("\n")

        file.write("# Batch simulation available functions\n")
        file.write("# random() - returns a random value between 0 to 1\n")
        file.write("# exp()    - exponential function\n")
        file.write("# log()    - natural logarithm \n")
        file.write("# sin()    - sine function \n")
        file.write("# cos()    - cosine function \n")
        file.write("# pi       - the constant pi \n")
        file.write("# abs()    - return absolute value \n")

        file.write("\n")
        file.write("# Batch simulation steps\n")
        file.write("#    1. Define the type of batch function (self.type)\n")
        file.write("#        1) 'Multiple Senarios' generates multiple output files for each senario\n")
        file.write("#        2) 'Monte Carlo Simulation' generates one output file with the mean values for all senarios.\n")
        file.write("#    2. Determine the number of senarios (i_num) in the batch function\n")
        file.write("#    3. Define the output files in the batch function\n")
        file.write("#        1) 'Multiple Senarios': the number of names in the vector(self.filenames) should be the same as i_num\n")
        file.write("#        2) 'Monte Carlo Simulation': the number of names in (self.filenames) should be 1\n")
        file.write("#    4. Assign the given values to the target parameter in each scenario\n")
        file.write("#       1) Develop a vector with the same number of desired values\n")
        file.write("#        2) Assign those values to the parameter\n")
        file.write("#    5. Assign random values to the target parameter in each scenario\n")
        file.write("#        1) use the function random() to generate a random number between 0 and 1\n")
        file.write("\n\n")

        file.write("# Batch simulation parameters\n")
        file.write("self.type = '" + self.system.batchfileoption + "'\n")
        if self.system.batchfileoption == "Monte Carlo Simulation": file.write("self.filenames = ['" + self.system.batchfilename + "']\n")
        else:                                                       file.write("self.filenames = ['"+self.system.batchfilename+"']\n")
        file.write("\n")
        file.write("i_num = 1\n\n")
        file.write("for i in range(i_num):\n\n")
        file.write("    systems.append(self.system.copy())\n")
        file.write("    systems[i].executablelines = []\n")
        file.write("    systems[i].iterativelines  = []\n")
        file.write("\n\n")

        file.write("    # System Properties\n")
        file.write("    #    Parameter: \n")
        file.write("    #       adv:    advection option\n")
        if self.system.adv == 'Steady flow' or self.system.adv == 'Period oscillation':
            file.write("    #       Vdar:   average Darcy velocity\n")
        if self.system.adv == 'Period oscillation':
            file.write("    #       Vtidal: maximum velocity of the periodical flow\n")
            file.write("    #       ptidal: period of the periodical flow\n")
        file.write("    #       bio:    bioturbation option\n")
        if self.system.bio == 'Uniform':
            file.write("    #       hbio:   bioturbation depth in uniform bioturbation option\n")
        elif self.system.bio == 'Depth-dependent':
            file.write("    #       sigma:  Gaussian function coefficient for depth-dependent bioturbation\n")
        if self.system.bio == 'Uniform' or self.system.bio == 'Depth-dependent':
            file.write("    #       Dbiop:  biodiffusion coefficient of solid particles at the benthic surface(z[0])\n")
            file.write("    #       Dbiopw: biodiffusion coefficient of solutes at the benthic surface(z[0])\n")
        file.write("    #       con:    consolidation option\n")
        if self.system.con == 'Consolidation':
            file.write("    #       hcon:   consolidation thickness\n")
            file.write("    #       t90:    time to reach a 90% consolidation option\n")
        file.write("    #   Batch function guide:\n")
        if self.system.adv == 'Steady flow':
            file.write("    #       U = Vdar \n")
        elif self.system.adv == 'Period oscillation':
            file.write("    #       U = Vdar + Vtidal*sin(2*pi*t/ptidal) \n")
        if self.system.bio == 'Uniform':
            file.write("    #       Dbiop = Dbiop (z < hbio) \n")
        elif self.system.bio == 'Depth-dependent':
            file.write("    #       Dbio(z)/Dbio(0) = exp(-((z-z[0])/sigma)^2) \n")
        if self.system.con == 'Consolidation':
            file.write("    #       Vcon = hcon * (2.3026 / system.t90) * exp(-(2.3026 / system.t90) * t) \n")
        file.write("    systems[i].adv      = '"+self.system.adv        +"'\n")
        if self.system.adv == 'Steady flow' or self.system.adv == 'Period oscillation':
            file.write("    systems[i].Vdar     = "+str(self.system.Vdar)   +"\n")
        if self.system.adv =='Period oscillation':
            file.write("    systems[i].Vtidal   = "+str(self.system.Vtidal) +"\n")
            file.write("    systems[i].ptidal   = "+str(self.system.ptidal) +"\n")
        file.write("    systems[i].bio      = '"+self.system.bio        +"'\n")
        if self.system.bio == 'Uniform':
            file.write("    systems[i].hbio     = "+str(self.system.hbio)   +"\n")
        elif self.system.bio == 'Depth-dependent':
            file.write("    systems[i].sigma    = "+str(self.system.sigma)  +"\n")
        if self.system.bio == 'Uniform' or self.system.bio == 'Depth-dependent':
            file.write("    systems[i].Dbiop    = "+str(self.system.Dbiop)  +"\n")
            file.write("    systems[i].Dbiopw   = "+str(self.system.Dbiopw) +"\n")
        file.write("    systems[i].con      = '"+self.system.con        +"'\n")
        if self.system.con == 'Consolidation':
            file.write("    systems[i].hcon     = "+str(self.system.hcon)   +"\n")
            file.write("    systems[i].t90      = "+str(self.system.t90)    +"\n")
        file.write("\n\n")

        file.write("    # Chemical Properties\n")
        file.write("    #   chemicals: solute chemicals\n")
        file.write("    #   nchemicals: the number of the chemicals\n")
        file.write("    #   Parameters: \n")
        file.write("    #       Dw:     water diffusivity\n")
        file.write("    #       Kdoc:   log value of DOC-water partioning coefficient\n")
        file.write("    systems[i].nchemicals = " + str(self.system.nchemicals) + "\n")
        file.write("    systems[i].chemicals = []\n")
        for n in range(self.system.nchemicals):
            file.write("    systems[i].chemicals.append(Chemical("+ str(n+1) +", 1))\n")
            file.write("    systems[i].chemicals["+str(n)+"].name    ='" + self.system.chemicals[n].name     + "'\n")
            file.write("    systems[i].chemicals["+str(n)+"].Dw      =" + str(self.system.chemicals[n].Dw)   + "\n")
            file.write("    systems[i].chemicals["+str(n)+"].Kdoc    =" + str(self.system.chemicals[n].Kdoc) + "\n")
        file.write("\n\n")

        file.write("    # Material Properties\n")
        file.write("    #   matrices:   solid matrices\n")
        file.write("    #   components: solid materials\n")
        file.write("    #   nmatrices:  the number of the matrices\n")
        file.write("    #   Parameters: \n")
        file.write("    #       e:         porosity\n")
        file.write("    #       rho:       bulk density\n")
        file.write("    #       foc:       organic carbon fraction\n")
        file.write("    #       mfraction: volumetric fraction of a component in the matrix\n")
        file.write("    #   Batch function guide:\n")
        file.write("    #       1. the same component that show up in various matrices should be consistent\n")
        file.write("    #       2. keep the sum of the mfractions to be 1\n")
        file.write("    systems[i].nmatrices = " + str(self.system.nmatrices) + "\n")
        file.write("    systems[i].matrices = []\n")
        for m in range(self.system.nmatrices):
            file.write("    systems[i].matrices.append(Matrix("+str(m+1)+"))\n")
            file.write("    systems[i].matrices["+str(m)+"].name     ='" + self.system.matrices[m].name      + "'\n")
            file.write("    systems[i].matrices["+str(m)+"].components = []\n")
            for mm in range(len(self.system.matrices[m].components)):
                file.write("    systems[i].matrices["+str(m)+"].components.append(MatrixComponent("+str(mm)+"))\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].name      ='" + self.system.matrices[m].components[mm].name          + "'\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].e         =" + str(self.system.matrices[m].components[mm].e)         + "\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].rho       =" + str(self.system.matrices[m].components[mm].rho)       + "\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].foc       =" + str(self.system.matrices[m].components[mm].foc)       + "\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].mfraction =" + str(self.system.matrices[m].components[mm].mfraction) + "\n")
        file.write("\n\n")

        file.write("    # Matrix Component Properties\n")
        file.write("    systems[i].component_list = []\n")
        file.write("    systems[i].components = []\n")
        file.write("    for matrix in systems[i].matrices:\n")
        file.write("        for component in matrix.components:\n")
        file.write("            try:    systems[i].component_list.index(component.name)\n")
        file.write("            except: \n")
        file.write("                systems[i].component_list.append(component.name)\n")
        file.write("                systems[i].components.append(component)\n")
        file.write("\n\n")

        file.write("    # Sorption Properties\n")
        file.write("    #   sorptions:  sorption\n")
        file.write("    #   Parameters: \n")
        file.write("    #       isotherm: isotherm option\n")
        file.write("    #       kinetic:  kinetic option\n")
        file.write("    #       K:        solid-water partitioning coefficient\n")
        file.write("    #       Koc:      organic carbon-water partitioning coefficient \n")
        file.write("    #       Kf, N:    Freundlich isotherm coefficients\n")
        file.write("    #       qmax,b:   Langmuir isotherm coefficients\n")
        file.write("    #       ksorp:    first-order sorption rate coefficient\n")
        file.write("    #   Batch function guide:\n")
        file.write("    #       1. Available isotherms:['Linear--Kd specified','Linear--Kocfoc','Freundlich','Langmuir']\n")
        file.write("    #       2. Available kinetic options:['Equilibrium', 'Transient']\n")
        file.write("    #       3. Assign the sorption parameters under the selected isotherm, the rest will be ignored\n")
        file.write("    systems[i].sorptions = {}\n")
        file.write("    for component in systems[i].components: \n")
        file.write("        systems[i].sorptions[component.name]={}\n")
        file.write("        for chemical in systems[i].chemicals: \n")
        file.write("            systems[i].sorptions[component.name][chemical.name]= Sorption(component, chemical)\n")
        for component in self.system.components:
            for chemical in self.system.chemicals:
                file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].isotherm = '"+   self.system.sorptions[component.name][chemical.name].isotherm +"'\n")
                file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].kinetic  = '"+   self.system.sorptions[component.name][chemical.name].kinetic +"'\n")
                if self.system.sorptions[component.name][chemical.name].isotherm == 'Linear--Kd specified':
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].K  = "+         str(self.system.sorptions[component.name][chemical.name].K) +"\n")
                elif self.system.sorptions[component.name][chemical.name].isotherm == 'Linear--Kocfoc':
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].Koc  = "+       str(self.system.sorptions[component.name][chemical.name].Koc) +"\n")
                elif self.system.sorptions[component.name][chemical.name].isotherm == 'Freundlich':
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].Kf  = "+        str(self.system.sorptions[component.name][chemical.name].Kf) +"\n")
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].N  = "+         str(self.system.sorptions[component.name][chemical.name].N) +"\n")
                elif self.system.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].qmax  = "+      str(self.system.sorptions[component.name][chemical.name].qmax) +"\n")
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].b  = "+         str(self.system.sorptions[component.name][chemical.name].b) +"\n")
                if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                    file.write("    systems[i].sorptions['"+component.name+"']['"+chemical.name+"'].ksorp  = "+     str(self.system.sorptions[component.name][chemical.name].ksorp) +"\n")
        file.write("\n\n")

        file.write("    # Layers Properties\n")
        file.write("    #   nlayers: number of layer\n")
        file.write("    #   Parameters: \n")
        file.write("    #       type:  solid matrix of the layer\n")
        file.write("    #       h:     layer thickness\n")
        file.write("    #       tort:  layer tortuosity model\n")
        file.write("    #       alpha: hydrodynamic dispersion coefficient nomalized by Darcy velocity\n")
        file.write("    #       doc:   concentration of dissovled organic carbon [mg/L]\n")
        file.write("    #   Batch function guide:\n")
        file.write("    #       hydrodynamic dispersion coefficient: Ddisp =  alpha * U\n")
        file.write("    systems[i].nlayers = " + str(self.system.nlayers) + "\n")
        file.write("    systems[i].layers = []\n")
        for j in range(self.system.nlayers):
            file.write("    systems[i].layers.append(Layer("+str(self.system.layers[j].number)+"))\n")
            file.write("    systems[i].layers["+str(j)+"].name       ='" + self.system.layers[j].name        + "'\n")
            file.write("    systems[i].layers["+str(j)+"].type       ='" + self.system.layers[j].type        + "'\n")
            file.write("    systems[i].layers["+str(j)+"].h          =" + str(self.system.layers[j].h)       + "\n")
            file.write("    systems[i].layers["+str(j)+"].alpha      =" + str(self.system.layers[j].alpha)   + "\n")
            file.write("    systems[i].layers["+str(j)+"].doc        =" + str(self.system.layers[j].doc)     + "\n")
        file.write("\n\n")

        file.write("    # Reaction Properties\n")
        file.write("    #   Parameter: \n")
        file.write("    #       coef:      stoichiochemical coefficients of reactants/products\n")
        file.write("    #       index:     rate index coefficients of reactants\n")
        file.write("    #   Batch function guide:\n")
        file.write("    #       1. the rate of the equation is expressed as r = coef * reactant^index \n")
        file.write("    systems[i].reactions = []\n")
        for l in range(len(self.system.reactions)):
            file.write("    systems[i].reactions.append(Reaction("+str(l+1)+"))\n")
            file.write("    systems[i].reactions["+str(l)+"].name      ='" + self.system.reactions[l].name     + "'\n")
            file.write("    systems[i].reactions["+str(l)+"].reactants = []\n")
            file.write("    systems[i].reactions["+str(l)+"].products  = []\n")
            for j in range(len(self.system.reactions[l].reactants)):
                file.write("    systems[i].reactions["+str(l)+"].reactants.append(Reactant("+str(j+1)+"))\n")
                file.write("    systems[i].reactions["+str(l)+"].reactants["+str(j)+"].name    ='" + self.system.reactions[l].reactants[j].name + "'\n")
                file.write("    systems[i].reactions["+str(l)+"].reactants["+str(j)+"].coef    =" + str(self.system.reactions[l].reactants[j].coef) + "\n")
                file.write("    systems[i].reactions["+str(l)+"].reactants["+str(j)+"].index   =" + str(self.system.reactions[l].reactants[j].index) + "\n")
            for j in range(len(self.system.reactions[l].products)):
                file.write("    systems[i].reactions["+str(l)+"].products.append(Product("+str(j+1)+"))\n")
                file.write("    systems[i].reactions["+str(l)+"].products["+str(j)+"].name    ='" + self.system.reactions[l].products[j].name + "'\n")
                file.write("    systems[i].reactions["+str(l)+"].products["+str(j)+"].coef    =" + str(self.system.reactions[l].products[j].coef) + "\n")
        file.write("\n\n")

        file.write("    #Reaction Rate Coefficients\n")
        file.write("    #   Parameter: \n")
        file.write("    #       lam:      rate coefficient of the reaction\n")
        file.write("    systems[i].coefficients = {}\n")
        file.write("    for layer in systems[i].layers: \n")
        file.write("        systems[i].coefficients[layer.name]={}\n")
        file.write("        for reaction in systems[i].reactions: \n")
        file.write("            systems[i].coefficients[layer.name][reaction.name]= Coefficient(layer, reaction)\n")
        for layer in self.system.layers:
            for reaction in self.system.reactions:
                file.write("    systems[i].coefficients['"+layer.name+"']['"+reaction.name+"'].lam = "+ str(self.system.coefficients[layer.name][reaction.name].lam) +"\n")
        file.write("\n\n")

        file.write("    # Boundary Conditions\n")
        file.write("    #   Parameter: \n")
        if self.system.topBCtype == 'Fixed Concentration':
            file.write("    #       Co:        fixed concentration at the top boundary\n")
        elif self.system.topBCtype == 'Mass transfer':
            file.write("    #       k:         benthic surface mass transfer coefficient\n")
            file.write("    #       Cw:        concentration in the overlying water body\n")
        elif self.system.topBCtype == 'Finite mixed water column':
            file.write("    #       k:         benthic surface mass transfer coefficient\n")
            file.write("    #       Cw:        concentration in the overlying water body\n")
            file.write("    #       kdecay:    first-order decay coefficient of the contaminant in the overlying water body\n")
            file.write("    #       kevap:     first-order evaporation coefficient of the contaminant in the overlying water body\n")
            file.write("    #       taucoefs:  parameters to estimate the mean retention time tau\n")
        file.write("    #       Cb:        fixed concentration or flux matched concentration at the bottom boundary\n")
        file.write("    #   Batch function guide:\n")
        if self.system.topBCtype == 'Fixed Concentration':
            file.write("    #       'Fixed concentration' top boundary:  C = Co\n")
        elif self.system.topBCtype == 'Mass transfer':
            file.write("    #       'Mass transfer' top boundary:  dC/dz = k*(C - Cw)\n")
        elif self.system.topBCtype == 'Finite mixed water column':
            file.write("    #       'Finite mixed water column' top boundary:  dC/dz = k*C*(1 - (U+k)/h/(1/tau+k/h))\n")
        if self.system.botBCtype == 'Fixed Concentration':
            file.write("    #       'Fixed concentration' bottom boundary:  C = Cb\n")
        else:
            file.write("    #       'Fixed concentration' bottom boundary:  F = U * Cb\n")
        file.write("    systems[i].BCs = {}\n")
        file.write("    systems[i].topBCtype = '"+self.system.topBCtype+"'\n")
        file.write("    systems[i].botBCtype = '"+self.system.botBCtype+"'\n")
        for chemical in self.system.chemicals:
            file.write("    systems[i].BCs['"+chemical.name+"'] = BC('"+chemical.name+"')\n")
            if self.system.topBCtype == 'Fixed Concentration':
                file.write("    systems[i].BCs['"+chemical.name+"'].Co     = "+str(self.system.BCs[chemical.name].Co)+"\n")
            else:
                file.write("    systems[i].BCs['"+chemical.name+"'].k      = "+str(self.system.BCs[chemical.name].k)+"\n")
                file.write("    systems[i].BCs['"+chemical.name+"'].Cw     = "+str(self.system.BCs[chemical.name].Cw)+"\n")
            if self.system.topBCtype == 'Finite mixed water column':
                file.write("    systems[i].BCs['"+chemical.name+"'].kdecay = "+str(self.system.BCs[chemical.name].kdecay)+"\n")
                file.write("    systems[i].BCs['"+chemical.name+"'].kevap  = "+str(self.system.BCs[chemical.name].kevap)+"\n")
            file.write("    systems[i].BCs['"+chemical.name+"'].Cb     = "+str(self.system.BCs[chemical.name].Cb)+"\n")
        file.write("\n")
        if self.system.topBCtype == 'Finite mixed water column':
            file.write("    systems[i].taucoefs = {}\n")
            file.write("    systems[i].taucoefs['Q'] = "    +str(self.system.taucoefs["Q"])+"\n")
            file.write("    systems[i].taucoefs['V'] = "    +str(self.system.taucoefs["V"])+"\n")
            file.write("    systems[i].taucoefs['h'] = "    +str(self.system.taucoefs["h"])+"\n")
            file.write("    systems[i].taucoefs['DOC'] = "  +str(self.system.taucoefs["DOC"])+"\n")
            file.write("    systems[i].taucoefs['Qevap'] = "+str(self.system.taucoefs["Qevap"])+"\n")
        file.write("\n\n")

        file.write("    # Initial Conditions\n")
        file.write("    #   Parameter: \n")
        file.write("    #       ICtype:  intial condition type\n")
        file.write("    #       uniform: the layer intital concentration when ICtype is 'Uniform'\n")
        file.write("    #       top:     the initial concentration at top of a layer when ICtype is 'Linear'\n")
        file.write("    #       bot:     the initial concentration at bottom of a layer when ICtype is 'Linear'\n")
        file.write("    #   Batch function guide:\n")
        file.write("    #       initial uniform concentration C = uniform \n")
        file.write("    #       initial linear concentration C = top + (bot-top)*(z-ztop)/(zbot-ztop) \n")
        file.write("    systems[i].ICs      = {}\n")
        file.write("    for layer in systems[i].layers: \n")
        file.write("        systems[i].ICs[layer.name]={}\n")
        file.write("        for chemical in systems[i].chemicals: \n")
        file.write("            systems[i].ICs[layer.name][chemical.name]= IC(layer.name, chemical.name)\n")
        for layer in self.system.layers:
            for chemical in self.system.chemicals:
                file.write("    systems[i].ICs['"+layer.name+"']['"+chemical.name+"'].ICtype  ='" + self.system.ICs[layer.name][chemical.name].ICtype    +"'\n")
                if self.system.ICs[layer.name][chemical.name].ICtype == 'Uniform':
                    file.write("    systems[i].ICs['"+layer.name+"']['"+chemical.name+"'].uniform =" + str(self.system.ICs[layer.name][chemical.name].uniform) +"\n")
                else:
                    file.write("    systems[i].ICs['"+layer.name+"']['"+chemical.name+"'].top     =" + str(self.system.ICs[layer.name][chemical.name].top)     +"\n")
                    file.write("    systems[i].ICs['"+layer.name+"']['"+chemical.name+"'].bot     =" + str(self.system.ICs[layer.name][chemical.name].bot)     +"\n")
        file.write("    systems[i].SolidICs = {}\n")
        file.write("    for layer in systems[i].layers: \n")
        file.write("        systems[i].SolidICs[layer.name]={}\n")
        file.write("        for component in systems[i].components: \n")
        file.write("            systems[i].SolidICs[layer.name][component.name]={}\n")
        file.write("            for chemical in systems[i].chemicals: \n")
        file.write("                systems[i].SolidICs[layer.name][component.name][chemical.name]= SolidIC(layer.name, component.name, chemical.name)\n")
        for layer in self.system.layers:
            for component in self.system.matrices[layer.type_index].components:
                for chemical in self.system.chemicals:
                    file.write("    systems[i].SolidICs['"+layer.name+"']['"+component.name+"']['"+chemical.name+"'].uniform  =" + str(self.system.SolidICs[layer.name][component.name][chemical.name].uniform) +"\n")
        file.write("\n\n")

        file.write("    # Solver Options\n")
        file.write("    #   Parameter: \n")
        file.write("    #       tfinal:        total simulation time\n")
        file.write("    #       outputsteps:   time steps in output files\n")
        file.write("    #       delt:          time step size\n")
        file.write("    #       players:       grid number in each layer\n")
        file.write("    #       averageoption: average options for periodical advection or deposition\n")
        file.write("    #       nlerror:       error tolerance in each time step\n")
        file.write("    systems[i].tfinal        = "+str(self.system.tfinal)      +"\n")
        file.write("    systems[i].outputsteps   = "+str(self.system.outputsteps) +"\n")
        file.write("    systems[i].delt          = "+str(self.system.delt)        +"\n")
        file.write("    systems[i].players       = [" + str(self.system.players[0]))
        for j in range(self.system.nlayers-1):
            file.write( ", " + str(self.system.players[j+1]))
        file.write("]\n")
        file.write("    systems[i].nlerror       = "+str(self.system.nlerror)     +"\n")
        file.write("    systems[i].averageoption = '"+self.system.averageoption   +"'\n")
        file.write("\n\n\n")

        file.write("#  The following lines define additional non-assignable parameters for the system\n")
        file.write("    # System Units\n")
        file.write("    systems[i].lengthunit = " + lengthunit + "\n")
        file.write("    systems[i].concunit = "   + concunit + "\n")
        file.write("    systems[i].timeunit = "   + timeunit + "\n")
        file.write("    systems[i].diffunit = "   + diffunit + "\n")
        file.write("\n")
        file.write("    # Chemical Properties\n")
        for n in range(self.system.nchemicals):
            file.write("    systems[i].chemicals["+str(n)+"].MW      =" + str(self.system.chemicals[n].MW)   + "\n")
            file.write("    systems[i].chemicals["+str(n)+"].formula ='" + text_converter(self.system.chemicals[n].formula)   + "'\n")
            file.write("    systems[i].chemicals["+str(n)+"].temp    =" + str(self.system.chemicals[n].temp) + "\n")
            file.write("    systems[i].chemicals["+str(n)+"].Ref     ='" + self.system.chemicals[n].Ref       + "'\n")
            file.write("    systems[i].chemicals["+str(n)+"].Koc     =" + str(self.system.chemicals[n].Koc)  + "\n")
            file.write("    systems[i].chemicals["+str(n)+"].Kf      =" + str(self.system.chemicals[n].Kf)   + "\n")
            file.write("    systems[i].chemicals["+str(n)+"].N       =" + str(self.system.chemicals[n].N)    + "\n")
        file.write("\n")
        file.write("    # Solid Material Properties\n")
        for m in range(self.system.nmatrices):
            file.write("    systems[i].matrices["+str(m)+"].model    ='" + self.system.matrices[m].model     + "'\n")
            for mm in range(len(self.system.matrices[m].components)):
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].sorp      ='" + self.system.matrices[m].components[mm].sorp          + "'\n")
                file.write("    systems[i].matrices["+str(m)+"].components["+str(mm)+"].tort      ='" + self.system.matrices[m].components[mm].tort          + "'\n")
        file.write("    # Solid Layer Properties\n")
        for j in range(self.system.nlayers):
            file.write("    systems[i].layers["+str(j)+"].tort       ='" + self.system.layers[j].tort        + "'\n")
        for l in range(len(self.system.reactions)):
            file.write("    systems[i].reactions["+str(l)+"].model     ='" + self.system.reactions[l].model    + "'\n")
            file.write("    systems[i].reactions["+str(l)+"].equation  ='" + text_converter(self.system.reactions[l].equation) + "'\n")
            for j in range(len(self.system.reactions[l].reactants)):
                file.write("    systems[i].reactions["+str(l)+"].reactants["+str(j)+"].formula ='" + text_converter(self.system.reactions[l].reactants[j].formula) + "'\n")
                file.write("    systems[i].reactions["+str(l)+"].reactants["+str(j)+"].MW      =" + str(self.system.reactions[l].reactants[j].MW) + "\n")
            for j in range(len(self.system.reactions[l].products)):
                file.write("    systems[i].reactions["+str(l)+"].products["+str(j)+"].formula ='" + text_converter(self.system.reactions[l].products[j].formula) + "'\n")
                file.write("    systems[i].reactions["+str(l)+"].products["+str(j)+"].MW      =" + str(self.system.reactions[l].products[j].MW) + "\n")
        for chemical in self.system.chemicals:
                file.write("    systems[i].BCs['"+chemical.name+"'].tau    = "+str(self.system.BCs[chemical.name].tau)+"\n")

        file.write("    systems[i].discrete      = '"+self.system.discrete        +"'\n")
        file.write("    systems[i].ptype         = '"+self.system.ptype           +"'\n")
        file.write("    systems[i].tvariable     = '"+self.system.tvariable       +"'\n")
        file.write("    systems[i].nonlinear     = '"+self.system.nonlinear       +"'\n")
        file.write("    systems[i].massoption    = '"+self.system.massoption      +"'\n")
        file.write("    systems[i].depgrid       = "+str(self.system.depgrid)     +"\n")
        file.write("    systems[i].tidalsteps    = "+str(self.system.tidalsteps)  +"\n")
        file.write("    systems[i].timeoption    = '"+self.system.timeoption      +"'\n")

def get_summary(system, database, materials):
    """Takes the system, shows the summary and allows user to modify."""

    root = CapSimWindow(buttons = 3)
    root.make_window(Summary(root, system, database, materials))
    root.mainloop()

    if root.main.get() == 0:
        check = 1
        try:
            pickle.dump(root.window.system, open(Filepath + '/temp_test_file.cpsm', 'w'))
        except:
            check = 0

        if check == 1:
            pickle.dump(root.window.system, open(Filepath + '/input_cpsm_files/%s.cpsm' %system.cpsmfilename, 'w'))
            if root.window.system.batchfileoption <> 'None':
                root.window.make_batch_file(Filepath)

        else:
            report_check = 0

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.chemicals, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to chemical properties in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.matrices, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to material properties in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.sorptions, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to sorption properties in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.layers, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to layer properties in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.reactions, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to reaction properties in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.coefficients, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to reaction coefficients in saving your input file. Please re-try it.')

            if report_check == 0:
                try:
                    pickle.dump(root.window.system.ICs, open(Filepath + '/temp_test_file.cpsm', 'w'))
                    pickle.dump(root.window.system.BCs, open(Filepath + '/temp_test_file.cpsm', 'w'))
                    pickle.dump(root.window.system.SolidICs, open(Filepath + '/temp_test_file.cpsm', 'w'))
                except:
                    report_check = 1
                    tkmb.showwarning('Administrator Error', 'It appears you have a problem related to auxiliary conditions in saving your input file. Please re-try it.')

            if report_check == 0:
                tkmb.showwarning('Administrator Error', 'It appears you have a problem related to system parameters or solver options in saving your input file. Please re-try it.')

        if check == 0:   system = None

    else: system = None

    root.destroy()

    return system
