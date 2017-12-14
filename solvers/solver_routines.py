#! /usr/bin/env python
#
#This file contains a variety of functions that are used collectively to
#simulate transport through a sediment cap system.  In essence the "System"
#object type that is built using the GUI or batch file is used to create an 
#instance of the "Parameters" class that lumps together the data for 
#discretization.  Based on the appropriate solver for the system, the 
#simulation is performed to ensure an accurate and stable solution.  The 
#details of each function are presented below.

import tkMessageBox as tkmb, cPickle as pickle, sys, os
import time as timer
import _winreg as wreg

if sys.path[0][-3:] == 'zip': 
    os.chdir(sys.path[0][:-12])
    path = sys.path[0][:-12]
else: path = sys.path[0]

CapSimReg = wreg.ConnectRegistry(None, wreg.HKEY_CURRENT_USER)
CapSimKey = wreg.OpenKey(CapSimReg, r'Software\CapSim')
Filepath =wreg.QueryValueEx(CapSimKey, 'FilePath')[0]
CapSimKey.Close()

import math, cPickle as pickle

from numpy               import matrix, array, linalg, zeros, transpose, interp, ceil, exp
from capsim_object_types import System, Layer, Chemical, BC, Reaction, Coefficient, SolidIC
from reactioneditor      import Reactant, Product
from capsim_functions    import consolidation, tidal, tauwater

class Parameters:
    """This object type creates the variables used for the finite 
    differencing to solve the transport equations for simulating a sediment
    cap.  The discretized matrix form of the differential equations used to
    solve for the unknown concentrations at the succeeding time step "Cn+1" 
    using the known concentration "Cn" is: 
        
    A * Cn+1 + a = B * Cn + b.

    Parameter instances can obtain the following attributes:

    z    -- the array of the grid points
    A    -- the matrix defined above representing the coefficients for the
            unknown values at the next time step
    B    -- the matrix defined above that is multiplied by the 
            concentrations at the previous time step
    a    -- the vector defined above containing boundary values and forcing
            functions, it also contains the non-linear terms at the next time step
    b    -- the vector defined above containing boundary values and forcing 
            functions, it also contains the non-linear terms at the next time step

    p    -- a list of the boundary points of each of the layers
    U    -- the Darcy velocity at the time step (scalar)
    D    -- an array of the diffusion coefficients at each grid point
    R    -- an array of the retardation factors at each grid point
    elam -- an array of the product of porosity and decay rate at each 
            grid point
    delz -- an array of the grid spacing in each layer
    delt -- the time step size
    """
    
    def __init__(self, system):
        """Constructor method.  Gets parameters from a "System" instance that
        is either built by the GUI or the user (for batch processing) and uses
        them to build a finite difference approximation for cap simulation.
        Variable definitions can be found in the "capsim_object_types" file."""
        
        self.bio            = 0
        self.con            = 0
        self.sorp           = 0
        self.tidal          = 0
        self.reac           = 0
        self.nsolidchemicals= 0
        self.varit          = 0
        self.mass           = 0
        self.topwater       = 0

        self.Vdar           = system.Vdar
        
        self.nlayers        = system.nlayers
        self.nchemicals     = system.nchemicals
        self.component_list = system.component_list
        self.topBCtype      = system.topBCtype
        self.botBCtype      = system.botBCtype

        if self.topBCtype == 'Finite mixed water column':
            self.topwater = 1
            self.taucoefs       = system.taucoefs

        try:    self.tstart     = system.tstart
        except: self.tstart     = 0

        try:    self.timeoption = system.timeoption
        except: self.timeoption = 'Implicit'

        self.outputsteps    = system.outputsteps
        self.discrete       = system.discrete
        self.delt           = system.delt
        self.ptype          = system.ptype
        self.tvariable      = system.tvariable
        self.delz           = system.delz
        self.players        = system.players
        self.tidalsteps     = system.tidalsteps
        self.nloption       = system.nonlinear
        self.nlerror        = system.nlerror
        self.averageoption  = system.averageoption
        self.depsteps       = system.depgrid
        self.depgrid        = 1

        if system.layers[0].name == 'Deposition':  self.steps      = self.depsteps
        else:                                      self.steps      = 0

        if self.averageoption == 'Average' and self.depsteps <= self.tidalsteps:
            self.steps      = self.tidalsteps

        if system.massoption == 'Track':           self.mass      = 1

        self.tfinal         = system.tfinal + self.delt * (1 + self.steps)
        self.tfinal_ori     = system.tfinal


        self.lengthunits = [u'\u03BCm', 'cm', 'm']
        self.concunits   = [u'\u03BCg/L', 'mg/L', 'g/L', u'\u03BCmol/L', 'mmol/L', 'mol/L']
        self.timeunits   = ['s', 'min', 'hr', 'day', 'yr']
        self.diffunits   = [u'cm\u00B2/s', u'cm\u00B2/yr']

        lengthunits = ['um', 'cm', 'm']
        concunits   = ['ug/L', 'mg/L', 'g/L', 'umol/L', 'mmol/L', 'mol/L']
        diffunits   = ['cm^2/s', 'cm^2/yr']

        if self.lengthunits.count(system.lengthunit) > 0:   self.lengthunit     = system.lengthunit
        else:                                               self.lengthunit     = self.lengthunits[lengthunits.index(system.lengthunit)]
        if self.concunits.count(system.concunit) > 0:       self.concunit       = system.concunit
        else:                                               self.concunit       = self.concunits[concunits.index(system.concunit)]
        if self.diffunits.count(system.diffunit) > 0:       self.diffunit       = system.diffunit
        else:                                               self.diffunit       = self.diffunits[diffunits.index(system.diffunit)]

        self.timeunit       = system.timeunit

        self.diff_factor = 1.
        self.k_factor    = 1.
        self.flux_factor = 0.001
        self.h_factor    = 100

        if self.lengthunit == self.lengthunits[0]:
            self.diff_factor = self.diff_factor * (10000**2)
            self.k_factor    = self.k_factor * 10000
            self.flux_factor = self.flux_factor / (10000**2)
            self.h_factor    = self.h_factor * 10000
        elif self.lengthunit == self.lengthunits[2]:
            self.diff_factor = self.diff_factor / (100**2)
            self.k_factor    = self.k_factor / 100
            self.flux_factor = self.flux_factor * (100**2)
            self.h_factor    = self.h_factor/100

        if self.diffunit == self.diffunits[0]:
            if self.timeunit == self.timeunits[1]:
                self.diff_factor = self.diff_factor * 60
            elif self.timeunit == self.timeunits[2]:
                self.diff_factor = self.diff_factor * 60 * 60
            elif self.timeunit == self.timeunits[3]:
                self.diff_factor = self.diff_factor * 60 * 60 * 24
            elif self.timeunit == self.timeunits[4]:
                self.diff_factor = self.diff_factor * 60 * 60 * 24 * 365.25
        else:
            if self.timeunit == self.timeunits[0]:
                self.diff_factor = self.diff_factor /365.25/24/60/60
            elif self.timeunit == self.timeunits[1]:
                self.diff_factor = self.diff_factor /365.25/24/60
            elif self.timeunit == self.timeunits[2]:
                self.diff_factor = self.diff_factor /365.25/24
            elif self.timeunit == self.timeunits[3]:
                self.diff_factor = self.diff_factor /365.25

        if self.timeunit == self.timeunits[0]:
            self.k_factor = self.k_factor / 60 / 60
        elif self.timeunit == self.timeunits[1]:
            self.k_factor = self.k_factor / 60
        elif self.timeunit == self.timeunits[3]:
            self.k_factor = self.k_factor * 24
        elif self.timeunit == self.timeunits[4]:
            self.k_factor = self.k_factor * 24 * 365.25

        if system.layers[0].name == 'Deposition':
            self.dep            = 1
            self.layers         = [layer.copy() for layer   in system.layers[1:]]
            self.delzdep        = self.delz[0]
            self.delz           = self.delz[1:]
            self.depplayer      = self.players[0]
            self.players        = self.players[1:]
            self.deplayer       = [system.layers[0].copy()]
            self.Vdep           = system.Vdep
            self.deplayer[0].h  = self.Vdep * self.tfinal
        else:
            self.dep            = 0
            self.layers         = [layer.copy() for layer   in system.layers]
            self.deplayer       = []
            self.Vdep           = 0
            self.depplayer      = 0

        self.depcheck     = 'None'
        self.toppoints    = 1
        if self.dep == 1:
            self.depsize = self.delzdep

        self.nlayers = len(self.layers)

        for i in range(len(self.players)):
            self.players[i] = int(self.players[i])

        self.layertot           = self.layers

        self.chemicals          = [chemical.copy() for chemical  in system.chemicals]
        self.solidchemicals     = []
        self.matrices       = [matrix.copy()   for matrix    in system.matrices]
        self.components     = [component.copy()for component in system.components]
        self.reactions      = [reaction.copy() for reaction  in system.reactions]

        self.ncomponents    = len(self.components)

        if self.concunits[:3].count(self.concunit) == 0:
            for chemical in self.chemicals: chemical.MW = 1

        for chemical in self.chemicals:
            chemical.Dw = chemical.Dw * self.diff_factor

        self.layerchemicals     = []
        self.deplayerchemicals  = []

        self.sorptions      = {}
        self.coefficients   = {}
        self.BCs            = {}
        self.ICs            = {}

        # Convert the transient sorption to kinetic reactions
        for component in self.components:
            self.sorptions[component.name] = {}
            for chemical in self.chemicals:
                self.sorptions[component.name][chemical.name] = system.sorptions[component.name][chemical.name].copy()
                if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                    for layer in self.layers:
                        for layer_component in self.matrices[layer.type_index].components:
                            if layer_component.name == component.name:   self.sorp = 1

                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                    if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich':
                        self.sorptions[component.name][chemical.name].ksorp   = self.sorptions[component.name][chemical.name].ksorp
                        self.sorptions[component.name][chemical.name].kdesorp = self.sorptions[component.name][chemical.name].kdesorp/(chemical.MW**(1.-1./self.sorptions[component.name][chemical.name].N))

                    elif self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                        self.sorptions[component.name][chemical.name].ksorp = self.sorptions[component.name][chemical.name].ksorp/chemical.MW

                    self.solidchemicals.append(Chemical(self.nchemicals+self.nsolidchemicals, soluable = 0))
                    self.solidchemicals[-1].get_solid_chemical(component.name, chemical.name, chemical.MW, self.chemicals.index(chemical))
                    self.nsolidchemicals = self.nsolidchemicals + 1

                    self.reactions.append(Reaction(len(self.reactions)))
                    reactants = Reactant(1)
                    products =  Product(1)

                    if self.sorptions[component.name][chemical.name].isotherm == 'Linear--Kd specified' or self.sorptions[component.name][chemical.name].isotherm == 'Linear--Kocfoc':
                        reactants.get_dynamic_sorption(chemical, 1)
                        products.get_dynamic_sorption(self.solidchemicals[-1])
                        self.reactions[-1].get_dynamic_sorption(component.name +  '_' + chemical.name + '_sorption', component.name, reactants, products, 'Fundamental')

                    elif self.sorptions[component.name][chemical.name].isotherm == 'Freundlich':
                        reactants.get_dynamic_sorption(chemical, 1)
                        products.get_dynamic_sorption(self.solidchemicals[-1])
                        self.reactions[-1].get_dynamic_sorption(component.name +  '_' + chemical.name + '_sorption', component.name, reactants, products, 'Fundamental')

                    elif self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                        reactants.get_dynamic_sorption(chemical, self.sorptions[component.name][chemical.name].qmax/chemical.MW)
                        products.get_dynamic_sorption(self.solidchemicals[-1])
                        self.reactions[-1].get_dynamic_sorption(component.name +  '_' + chemical.name + '_sorption', component.name, reactants, products, 'Langmuir')

                    self.reactions.append(Reaction(len(self.reactions)))
                    reactants = Reactant(1)
                    products =  Product(1)
                    if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich':
                        reactants.get_dynamic_desorption(self.solidchemicals[-1], 1./self.sorptions[component.name][chemical.name].N)
                        products.get_dynamic_desorption(chemical)
                        self.reactions[-1].get_dynamic_sorption(component.name +  '_' + chemical.name + '_desorption', component.name, reactants, products, 'User-defined')
                    else:
                        reactants.get_dynamic_desorption(self.solidchemicals[-1], 1)
                        products.get_dynamic_desorption(chemical)
                        self.reactions[-1].get_dynamic_sorption(component.name +  '_' + chemical.name + '_desorption', component.name, reactants, products, 'Fundamental')

        for reaction in self.reactions:
            for reactant in reaction.reactants:
                reactant.coef = - reactant.coef

        if system.adv == 'Period oscillation':
            self.tidal      = 1
            self.U          = self.Vdar
            self.Vtidal     = system.Vtidal
            self.ptidal     = system.ptidal

        elif system.adv == 'Steady flow':
            self.U = self.Vdar
            self.tidalsteps = 0
        else:
            self.U = 0
            self.tidalsteps = 0

        if system.bio <> 'None' and (system.layers[0].name == 'Deposition' or system.bio == 'Depth-dependent' or system.layers[0].h < system.hbio) and system.Dbiop> 0 and system.nlayers > 1:     self.biomix = 1
        else:                                                                                                                                                                                      self.biomix = 0

        self.biolayerinit           = 0
        self.mixcomponents          = []
        self.mixcomponentslist      = []
        self.nmixcomponents         = 0
        self.mixsolidchemicalslist  = []
        self.mixsolidchemicals      = []
        self.solidcomponents        = []

        if system.bio == 'Uniform' or system.bio == 'Depth-dependent':
            self.bio     = 1
            self.biotype = system.bio
            self.sigma   = system.sigma
            if system.bio == 'Depth-dependent':
                if self.dep == 1:
                    self.hbio    = sum([layer.h for layer in self.layers]) + self.tfinal*self.deplayer[0].h
                else:
                    self.hbio    = sum([layer.h for layer in self.layers])
            else:
                self.hbio    = system.hbio

            if self.diffunit == self.diffunits[0]:
                self.Dbiop   = system.Dbiop  /86400/365 * self.diff_factor
                self.Dbiopw  = system.Dbiopw /86400/365 * self.diff_factor
            else:
                self.Dbiop   = system.Dbiop  * self.diff_factor
                self.Dbiopw  = system.Dbiopw * self.diff_factor


            if self.biomix == 1:
                H = 0
                for j in range(len(self.layers)):
                    layer = self.layers[j]
                    if self.hbio > H:
                        for component in self.matrices[layer.type_index].components:
                            if self.mixcomponentslist.count(component.name) == 0:
                                self.mixcomponentslist.append(component.name)
                                self.mixcomponents.append(component.copy())
                                for solidchemical in self.solidchemicals:
                                    if component.name == solidchemical.component_name:
                                        if self.mixsolidchemicalslist.count(solidchemical.name) == 0:
                                            self.mixsolidchemicals.append(solidchemical.copy())
                                            self.mixsolidchemicalslist.append(solidchemical.name)
                    if self.hbio > H + layer.h:
                        if j < len(self.layers) - 1:
                            self.biolayerinit = self.biolayerinit + 1
                    H = H + layer.h

                if self.dep == 1:
                    for component in self.matrices[self.deplayer[0].type_index].components:
                        if self.mixcomponentslist.count(component.name) == 0:
                            self.mixcomponentslist.append(component.name)
                            self.mixcomponents.append(component.copy())
                            for solidchemical in self.solidchemicals:
                                if component.name == solidchemical.component_name:
                                    if self.mixsolidchemicalslist.count(solidchemical.name) == 0:
                                        self.mixsolidchemicals.append(solidchemical.copy())
                                        self.mixsolidchemicalslist.append(solidchemical.name)
        else:
            self.Dbiop   = 0

        if system.con == 'Consolidation':
            self.con    = 1
            self.kcon   = 2.3026 / system.t90
            self.Vcon0  = system.hcon * self.kcon
            self.U      = self.Vdar + self.Vcon0

        self.U_plus_1 = self.U

        self.chemical_list = [chemical.name for chemical in (self.chemicals + self.solidchemicals)]
        self.solidchemical_list = [chemical.name for chemical in (self.solidchemicals)]

        # Load the parameters of chemical informations in layers

        layer_num = 0
        if self.biolayerinit > 0:
            for i in range(self.biolayerinit+1):
                layer = self.layers[i]
                layer.nchemicals            = self.nchemicals
                layer.nsolidchemicals       = 0
                layer.nsolidcomponents      = 0
                layer.chemicals             = [chemical.copy() for chemical in self.chemicals]
                layer.solidchemicals        = []
                layer.solidcomponents       = []
                layer.components            = [component.copy() for component in self.mixcomponents]
                layer.component_list        = [component.name   for component in self.mixcomponents]
                layer.init_components       = [component.copy() for component in self.matrices[layer.type_index].components]
                layer.init_component_list   = [component.name   for component in self.matrices[layer.type_index].components]
                layer.lreactions            = []
                layer.nlreactions           = []

                for solidchemical in self.mixsolidchemicals:
                    layer.solidchemicals.append(solidchemical.copy())
                    layer.solidchemicals[-1].component_name = solidchemical.component_name
                    layer.solidchemicals[-1].chemical_name  = solidchemical.chemical_name
                    layer.solidchemicals[-1].component_index = self.mixcomponentslist.index(layer.solidchemicals[-1].component_name)
                    layer.chemicals.append(solidchemical.copy())
                    layer.nsolidchemicals = layer.nsolidchemicals + 1
                    layer.nchemicals      = layer.nchemicals + 1

                layer_num = layer_num + 1

        for i in range(layer_num, self.nlayers):
            layer = self.layers[i]
            layer.nchemicals            = self.nchemicals
            layer.nsolidchemicals       = 0
            layer.chemicals             = [chemical.copy() for chemical in self.chemicals]
            layer.solidchemicals        = []
            layer.components            = [component.copy() for component in self.matrices[layer.type_index].components]
            layer.component_list        = [component.name for component in layer.components]
            layer.init_components       = [component.copy() for component in self.matrices[layer.type_index].components]
            layer.init_component_list   = [component.name for component in layer.components]
            layer.lreactions            = []
            layer.nlreactions           = []

            for solidchemical in self.solidchemicals:
                for component in layer.components:
                     if component.name == solidchemical.component_name:
                        layer.solidchemicals.append(solidchemical.copy())
                        layer.solidchemicals[-1].component_name = solidchemical.component_name
                        layer.solidchemicals[-1].chemical_name  = solidchemical.chemical_name
                        layer.solidchemicals[-1].component_index = layer.components.index(component)
                        layer.chemicals.append(solidchemical.copy())
                        layer.nsolidchemicals = layer.nsolidchemicals + 1
                        layer.nchemicals      = layer.nchemicals + 1

        if self.dep == 1:
            if self.biomix == 1:
                self.deplayer[0].nchemicals             = self.nchemicals
                self.deplayer[0].nsolidchemicals        = 0
                self.deplayer[0].chemicals              = [chemical.copy() for chemical in self.chemicals]
                self.deplayer[0].solidchemicals         = []
                self.deplayer[0].components             = [component.copy() for component in self.mixcomponents]
                self.deplayer[0].component_list         = [component.name   for component in self.mixcomponents]
                self.deplayer[0].init_components        = [component.copy() for component in self.matrices[self.deplayer[0].type_index].components]
                self.deplayer[0].init_component_list    = [component.name   for component in self.matrices[self.deplayer[0].type_index].components]
                self.deplayer[0].lreactions             = []
                self.deplayer[0].nlreactions            = []

                for solidchemical in self.mixsolidchemicals:
                    self.deplayer[0].solidchemicals.append(solidchemical.copy())
                    self.deplayer[0].solidchemicals[-1].component_name = solidchemical.component_name
                    self.deplayer[0].solidchemicals[-1].chemical_name  = solidchemical.chemical_name
                    self.deplayer[0].solidchemicals[-1].component_index = self.mixcomponentslist.index(self.deplayer[0].solidchemicals[-1].component_name)
                    self.deplayer[0].chemicals.append(solidchemical.copy())
                    self.deplayer[0].nsolidchemicals = self.deplayer[0].nsolidchemicals + 1
                    self.deplayer[0].nchemicals      = self.deplayer[0].nchemicals + 1

                self.deplayer[0].chemical_list      = [chemical.name for chemical in self.deplayer[0].chemicals]
                self.deplayer[0].solidchemical_list = [chemical.name for chemical in self.deplayer[0].solidchemicals]

            else:
                self.deplayer[0].nchemicals             = self.nchemicals
                self.deplayer[0].nsolidchemicals        = 0
                self.deplayer[0].chemicals              = [chemical.copy() for chemical in self.chemicals]
                self.deplayer[0].solidchemicals         = []
                self.deplayer[0].components             = [component.copy() for component in self.matrices[self.deplayer[0].type_index].components]
                self.deplayer[0].component_list         = [component.name   for component in self.deplayer[0].components]
                self.deplayer[0].init_components        = [component.copy() for component in self.matrices[self.deplayer[0].type_index].components]
                self.deplayer[0].init_component_list    = [component.name   for component in self.deplayer[0].components]
                self.deplayer[0].lreactions             = []
                self.deplayer[0].nlreactions            = []

                for solidchemical in self.solidchemicals:
                    for component in self.deplayer[0].components:
                        if component.name == solidchemical.component_name:
                            self.deplayer[0].solidchemicals.append(solidchemical.copy())
                            self.deplayer[0].solidchemicals[-1].component_name = solidchemical.component_name
                            self.deplayer[0].solidchemicals[-1].chemical_name  = solidchemical.chemical_name
                            self.deplayer[0].solidchemicals[-1].component_index = self.deplayer[0].components.index(component)
                            self.deplayer[0].nsolidchemicals = self.deplayer[0].nsolidchemicals + 1
                            self.deplayer[0].nchemicals      = self.deplayer[0].nchemicals + 1
                            self.deplayer[0].chemicals.append(solidchemical.copy())

                self.deplayer[0].chemical_list      = [chemical.name for chemical in self.deplayer[0].chemicals]
                self.deplayer[0].solidchemical_list = [chemical.name for chemical in self.deplayer[0].solidchemicals]

        for layer in self.deplayer + self.layers:
            self.coefficients[layer.name] = {}
            for reaction in self.reactions:
                if reaction.name.count('_sorption') == 1:
                    self.coefficients[layer.name][reaction.name] = Coefficient(layer, reaction)
                    for component in layer.components:
                        if component.name == reaction.component:
                            self.coefficients[layer.name][reaction.name].get_dynamic_sorption(self.sorptions[component.name][reaction.reactants[0].name].ksorp)
                elif reaction.name.count('_desorption') == 1:
                    self.coefficients[layer.name][reaction.name] = Coefficient(layer, reaction)
                    for component in layer.components:
                        if component.name == reaction.component:
                            self.coefficients[layer.name][reaction.name].get_dynamic_sorption(self.sorptions[component.name][reaction.products[0].name].kdesorp)
                else:
                    self.coefficients[layer.name][reaction.name] = system.coefficients[layer.name][reaction.name].copy()


        for layer in self.layers:
            layer.chemical_list         = [chemical.name for chemical in layer.chemicals]
            layer.solidchemical_list    = [chemical.name for chemical in layer.solidchemicals]

            for reaction in self.reactions:
                if self.coefficients[layer.name][reaction.name].lam > 0:
                    indeces    = [reactant.index for reactant in reaction.reactants]
                    if sum(indeces) == 1 and indeces.count(1.) == 1 and reaction.model != 'Langmuir':
                        layer.lreactions.append(reaction.copy())
                        for chemical in (layer.lreactions[-1].reactants + layer.lreactions[-1].products):
                            chemical.count = layer.chemical_list.index(chemical.name)
                        for chemical in (layer.lreactions[-1].reactants):
                            if chemical.index == 1: layer.lreactions[-1].key = chemical
                    else:
                        layer.nlreactions.append(reaction.copy())
                        for chemical in (layer.nlreactions[-1].reactants + layer.nlreactions[-1].products):
                            chemical.count = layer.chemical_list.index(chemical.name)
                        self.reac = 1

        if self.dep == 1:
            for reaction in self.reactions:
                if self.coefficients[self.deplayer[0].name][reaction.name].lam > 0:
                    indeces    = [reactant.index for reactant in reaction.reactants]
                    if sum(indeces) == 1 and indeces.count(1.) == 1 and reaction.model != 'Langmuir':
                        self.deplayer[0].lreactions.append(reaction.copy())
                        for chemical in (self.deplayer[0].lreactions[-1].reactants + self.deplayer[0].lreactions[-1].products):
                            chemical.count = self.deplayer[0].chemical_list.index(chemical.name)
                        for chemical in (self.deplayer[0].lreactions[-1].reactants):
                            if chemical.index == 1: self.deplayer[0].lreactions[-1].key = chemical
                    else:
                        self.deplayer[0].nlreactions.append(reaction.copy())
                        for chemical in (self.deplayer[0].nlreactions[-1].reactants + self.deplayer[0].nlreactions[-1].products):
                            chemical.count = self.deplayer[0].chemical_list.index(chemical.name)
                        self.reac = 1

        for j in range(self.nlayers):
            layer = self.layers[j]
            layer.nlowersolidchemicals      = 0
            layer.lowersolidchemicalnums    = []
            layer.lowersolidchemicals       = []

            if j == self.nlayers - 1:
                layer.lowersolidchemicalnums   = [n + self.nchemicals for n in range(layer.nsolidchemicals)]

            else:
                for solidchemical in layer.solidchemicals:
                    if self.layers[j+1].solidchemical_list.count(solidchemical.name) == 1:
                        layer.lowersolidchemicalnums.append(self.nchemicals + self.layers[j+1].solidchemical_list.index(solidchemical.name))
                    else:
                        layer.lowersolidchemicalnums.append(self.layers[j+1].nchemicals + layer.nlowersolidchemicals)
                        layer.lowersolidchemicals.append(solidchemical)
                        layer.nlowersolidchemicals = layer.nlowersolidchemicals + 1

        if self.dep == 1:
            self.deplayer[0].nlowersolidchemicals   = 0
            self.deplayer[0].lowersolidchemicalnums = []
            self.deplayer[0].lowersolidchemicals    = []

            for solidchemical in self.deplayer[0].solidchemicals:
                if self.layers[0].solidchemical_list.count(solidchemical.name) == 1:
                    self.deplayer[0].lowersolidchemicalnums.append(self.nchemicals + self.layers[0].solidchemical_list.index(solidchemical.name))
                else:
                    self.deplayer[0].lowersolidchemicalnums.append(self.layers[0].nchemicals + self.deplayer[0].nlowersolidchemicals)
                    self.deplayer[0].lowersolidchemicals.append(solidchemical)
                    self.deplayer[0].nlowersolidchemicals = self.deplayer[0].nlowersolidchemicals + 1

        for chemical in self.chemicals:
            self.BCs[chemical.name]           = system.BCs[chemical.name].copy()
            self.BCs[chemical.name].k         = self.BCs[chemical.name].k * self.k_factor
            self.BCs[chemical.name].Co        = self.BCs[chemical.name].Co/chemical.MW
            self.BCs[chemical.name].Cb        = self.BCs[chemical.name].Cb/chemical.MW
            self.BCs[chemical.name].Cw        = self.BCs[chemical.name].Cw/chemical.MW
            self.BCs[chemical.name].Cw_plus_1 = self.BCs[chemical.name].Cw/chemical.MW

        for chemical in self.solidchemicals:
            self.BCs[chemical.name] = BC(chemical.name, 0)

        for layer in self.deplayer + self.layers:
            self.ICs[layer.name] = {}
            for chemical in self.chemicals:
                self.ICs[layer.name][chemical.name] = system.ICs[layer.name][chemical.name].copy()
                self.ICs[layer.name][chemical.name].uniform = self.ICs[layer.name][chemical.name].uniform / chemical.MW
                self.ICs[layer.name][chemical.name].top = self.ICs[layer.name][chemical.name].top / chemical.MW
                self.ICs[layer.name][chemical.name].bot = self.ICs[layer.name][chemical.name].bot / chemical.MW

            for chemical in self.solidchemicals:
                if [component.name for component in system.matrices[layer.type_index].components].count(chemical.component_name) >= 1:
                    self.ICs[layer.name][chemical.name] = system.SolidICs[layer.name][chemical.component_name][chemical.chemical_name].copy()
                    self.ICs[layer.name][chemical.name].uniform = self.ICs[layer.name][chemical.name].uniform / chemical.MW
                else:
                    self.ICs[layer.name][chemical.name] = SolidIC(layer.name, chemical.component_name, chemical.chemical_name)
                    self.ICs[layer.name][chemical.name].uniform = 0 / chemical.MW


        self.Cmax = {}
        for chemical in (self.chemicals + self.solidchemicals):
            self.Cmax[chemical.name] = max(self.BCs[chemical.name].Co, self.BCs[chemical.name].Cw, self.BCs[chemical.name].Cb)

        for layer in (self.deplayer + self.layers):
            for chemical in layer.chemicals:
                try:
                    if self.Cmax[chemical.name] < max(self.ICs[layer.name][chemical.name].uniform, self.ICs[layer.name][chemical.name].top, self.ICs[layer.name][chemical.name].bot):
                        self.Cmax[chemical.name] = max(self.ICs[layer.name][chemical.name].uniform, self.ICs[layer.name][chemical.name].top, self.ICs[layer.name][chemical.name].bot)
                except: pass

        for layer in (self.deplayer + self.layers):
            for chemical in self.chemicals:
                Solid_C = 0
                for component in self.matrices[layer.type_index].components:
                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                        Solid_C = Solid_C + self.sorptions[component.name][chemical.name].get_C(system.SolidICs[layer.name][component.name][chemical.name].uniform) * component.mfraction
                Liquid_C = max(self.ICs[layer.name][chemical.name].uniform, self.ICs[layer.name][chemical.name].top, self.ICs[layer.name][chemical.name].bot)
                if  Solid_C * self.matrices[layer.type_index].rho+ self.matrices[layer.type_index].e * Liquid_C > self.Cmax[chemical.name]:
                    self.Cmax[chemical.name] = Solid_C * self.matrices[layer.type_index].rho+ self.matrices[layer.type_index].e * Liquid_C


    def make_uniform_grid(self):

        """Creates a uniform grid for assessing contaminant transport in a cap
        using the depth in each layer with a spacing of approximately delzmax.
        The actual spacing is set up such that the boundary of each layer falls
        exactly on a grid point.  Creates a list with the depth of each grid
        point "z," the actual grid spacing in each layer "delz," the list number
        of the boundary points "p."
        """

        self.z      = []
        self.p      = []
        self.ptot   = []
        self.cp     = []
        self.cptot  = []
        self.cpL    = []
        self.cptotL = []

        self.z.append(0)
        self.p.append(0)
        self.ptot.append(0)
        self.cp.append(0)
        self.cpL.append(0)
        self.cptot.append(0)
        self.cptotL.append(0)

        for i in range(len(self.layers)):
            layer = self.layers[i]
            self.p.append(self.p[-1] + self.players[i])
            self.ptot.append(self.ptot[-1] + self.players[i])
            for j in range (self.p[-2], self.p[-1]):
                self.z.append(self.z[-1]+self.delz[i])
            if i == 0:
                self.cp.append(self.cp[-1] + self.players[i]* layer.nchemicals)
                self.cptot.append(self.cptot[-1] + self.players[i]* layer.nchemicals)
            else:
                self.cp.append(self.cp[-1] + self.players[i]* layer.nchemicals+self.layers[i-1].nlowersolidchemicals)
                self.cptot.append(self.cptot[-1] + self.players[i]* layer.nchemicals+self.layers[i-1].nlowersolidchemicals)

            self.cpL.append(self.cpL[-1] + self.players[i]* layer.nchemicals+self.layers[i].nlowersolidchemicals)
            self.cptotL.append(self.cptotL[-1] + self.players[i]* layer.nchemicals+self.layers[i].nlowersolidchemicals)

        if self.dep == 1:
            points       = int(round(self.depplayer*self.tfinal,10))-1
            self.zdep    = [-self.delzdep]
            for i in range(0, points):
                self.zdep.append(self.zdep[-1] - self.delzdep)
            self.zdep.reverse()
            self.pdep = [len(self.zdep)]
            self.bdep = [len(self.zdep) * self.deplayer[0].nchemicals]

        else:
            self.zdep = []
            self.pdep = [0]
            self.bdep = [0]

        self.ptotbio     = [ptot    for ptot    in self.ptot]
        self.cptotbio    = [cptot   for cptot   in self.cptot]
        self.cptotLbio   = [cptotL  for cptotL  in self.cptotL]
        self.layertotbio = [layer   for layer   in self.layertot]


        if self.bio <> 1:   self.pbio = -1

    def get_initial_concentrations(self):
        """Uses the list of "Layer" objects and their concentrations to fill in
        the values of concentration at each grid point "C0" """

        C0L = []

        if self.layertot[0].ICtype == 'Uniform':
            for chemical in (self.chemicals):
                C0L.append(self.ICs[self.layertot[0].name][chemical.name].uniform)
            for solidchemical in self.layertot[0].solidchemicals:
                try:    C0L.append(self.ICs[self.layertot[0].name][solidchemical.name].uniform)
                except: C0L.append(0)

        elif self.layertot[0].ICtype == 'Linear':
            for chemical in (self.chemicals):
                C0L.append(self.ICs[self.layertot[0].name][chemical.name].top)
            for solidchemical in self.layertot[0].solidchemicals:
                try:    C0L.append(self.ICs[self.layertot[0].name][solidchemical.name].uniform)
                except: C0L.append(0)

        for j in range(len(self.ptot) - 1):
            layer = self.layertot[j]
            if layer.ICtype == 'Uniform':
                for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                    for chemical in (self.chemicals):
                        C0L.append(self.ICs[layer.name][chemical.name].uniform)
                    for solidchemical in layer.solidchemicals:
                        try: C0L.append(self.ICs[layer.name][solidchemical.name].uniform)
                        except: C0L.append(0)

            elif layer.ICtype == 'Linear':
                for i in range (self.ptot[j] + 1,self.ptot[j + 1]):
                    for chemical in (self.chemicals):
                        top = self.ICs[layer.name][chemical.name].top
                        bot = self.ICs[layer.name][chemical.name].bot
                        C0L.append(top + (bot-top)/(self.z[self.ptot[j + 1]] - self.z[self.ptot[j]])*(self.z[i]-self.z[self.ptot[j]]))
                    for solidchemical in layer.solidchemicals:
                        try: C0L.append(self.ICs[layer.name][solidchemical.name].uniform)
                        except: C0L.append(0)

            if j < len(self.ptot) - 2:
                if self.layertot[j+1].ICtype == 'Uniform':
                    for chemical in (self.chemicals):
                        C0L.append(self.ICs[self.layertot[j+1].name][chemical.name].uniform)
                    for solidchemical in self.layertot[j+1].solidchemicals:
                        try: C0L.append(self.ICs[self.layertot[j+1].name][solidchemical.name].uniform)
                        except: C0L.append(0)
                    for solidchemical in layer.lowersolidchemicals:
                        C0L.append(0)

                elif self.layertot[j+1].ICtype == 'Linear':
                    for chemical in (self.chemicals):
                        C0L.append(self.ICs[self.layertot[j+1].name][chemical.name].top)
                    for solidchemical in self.layertot[j+1].solidchemicals:
                        try: C0L.append(self.ICs[self.layertot[j+1].name][solidchemical.name].uniform)
                        except: C0L.append(0)
                    for solidchemical in layer.lowersolidchemicals:
                        C0L.append(0)
            else:
                if self.layertot[-1].ICtype == 'Uniform':
                    for chemical in (self.chemicals):
                        C0L.append(self.ICs[self.layertot[-1].name][chemical.name].uniform)
                    for solidchemical in self.layertot[-1].solidchemicals:
                        try:    C0L.append(self.ICs[self.layertot[-1].name][solidchemical.name].uniform)
                        except: C0L.append(0)

                elif self.layertot[-1].ICtype == 'Linear':
                    for chemical in (self.chemicals):
                        C0L.append(self.ICs[self.layertot[-1].name][chemical.name].bot)
                    for solidchemical in self.layertot[-1].solidchemicals:
                        try:    C0L.append(self.ICs[self.layertot[-1].name][solidchemical.name].uniform)
                        except: C0L.append(0)

        C0 = [C for C in C0L]

        for jj in range(len(self.ptot) - 2):
            j = jj + 1
            if self.layertot[j-1].ICtype == 'Uniform':
                for n in range(self.nchemicals):
                    C0[self.cptot[j]+n] = self.ICs[self.layertot[j-1].name][self.chemicals[n].name].uniform

            elif self.layertot[j-1].ICtype == 'Linear':
                for n in range(self.nchemicals):
                    C0[self.cptot[j]+n] = self.ICs[self.layertot[j-1].name][self.chemicals[n].name].bot

            for n in range(self.layertot[j].nsolidchemicals):
                solidchemical = self.layertot[j].solidchemicals[n]
                if self.layertot[j-1].solidchemical_list.count(solidchemical.name) > 0:
                    C0[self.cptot[j]+n+self.nchemicals] = self.ICs[self.layertot[j-1].name][solidchemical.name].uniform
                else:
                    C0[self.cptot[j]+n+self.nchemicals] = 0

            for solidchemical in self.layertot[j-1].lowersolidchemicals:
                n = self.layertot[j-1].lowersolidchemicalnums[self.layertot[j-1].solidchemical_list.index(solidchemical.name)]
                C0[self.cptot[j]+n] = self.ICs[self.layertot[j-1].name][solidchemical.name].uniform

        return C0, C0L


    def get_initial_boundary_concentration(self, Cn, CnL):

        C0 = [C for C in Cn]


        for jj in range(len(self.ptotbio) - 2):
            j = jj + 1
            for n in range(self.nchemicals):
                C0[self.cptotbio[j]+n] = ((Cn[self.cptotbio[j]+n]*self.Rs_plus_1[self.cptotbio[j]+n]*(self.z[self.ptotbio[j]]-self.z[self.ptotbio[j]-1]) + CnL[self.cptotbio[j]+n]*self.RsL_plus_1[self.cptotbio[j]+n]*(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]]))/
                                    (self.Rs_plus_1[self.cptotbio[j]+n]*(self.z[self.ptotbio[j]]-self.z[self.ptotbio[j]-1]) + self.RsL_plus_1[self.cptotbio[j]+n]*(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])))

        return C0

    def get_initial_component_fractions(self):
        """Uses the list of "Layer" objects and their concentrations to fill in
        the values of concentration at each grid point "C0" """

        Fis0 = {}
        for component in (self.components):
            Fis0[component.name] = []
            if self.layers[0].init_component_list.count(component.name) == 0: Fis0[component.name].append(0)
            else:                                                             Fis0[component.name].append(self.layers[0].init_components[self.layers[0].init_component_list.index(component.name)].fraction)
            for j in range(len(self.p) - 1):
                layer = self.layers[j]
                for i in range(self.p[j]+1, self.p[j + 1]+1):
                    if layer.init_component_list.count(component.name) == 0:  Fis0[component.name].append(0)
                    else:                                                     Fis0[component.name].append(layer.init_components[layer.init_component_list.index(component.name)].fraction)

        Fis0L = {}
        for component in (self.components):
            Fis0L[component.name] = [Fi for Fi in Fis0[component.name]]
            for j in range(len(self.p) - 1):
                layer = self.layers[j]
                if layer.init_component_list.count(component.name) == 0:     Fis0L[component.name][self.p[j]] = 0
                else:                                                        Fis0L[component.name][self.p[j]] = layer.init_components[layer.init_component_list.index(component.name)].fraction

        Fis0M = {}
        for component in (self.components):
            Fis0M[component.name] = [Fi for Fi in Fis0[component.name]]
            for j in range(len(self.p) - 1):
                if self.bio == 1 and self.p[j] < self.pbio:
                    layer = self.layers[j]
                    if layer.init_component_list.count(component.name) == 0:     Fis0M[component.name][self.p[j]] = (Fis0M[component.name][self.p[j]]* (self.z[self.p[j]]-self.z[self.p[j]-1]))/(self.z[self.p[j]+1]-self.z[self.p[j]-1])
                    else:                                                        Fis0M[component.name][self.p[j]] = (Fis0M[component.name][self.p[j]]* (self.z[self.p[j]]-self.z[self.p[j]-1]) + layer.init_components[layer.init_component_list.index(component.name)].fraction*(self.z[self.p[j]+1]-self.z[self.p[j]]))/(self.z[self.p[j]+1]-self.z[self.p[j]-1])

        return Fis0, Fis0L, Fis0M

    def make_grid_es_rhos(self, Fis, FisL):

        """Makes the list of values of "R" for the finite difference equations
        for the grid "z" with boundary points "p" and concentrations "Cn." """

        self.es_plus_1      = []
        self.rhos_plus_1    = []

        self.es_plus_1.append(0)
        self.rhos_plus_1.append(0)
        for component in self.components:
            self.es_plus_1[-1]      = self.es_plus_1[-1]    + Fis[component.name][0] * component.e
            self.rhos_plus_1[-1]    = self.rhos_plus_1[-1]  + Fis[component.name][0] * component.rho

        for j in range(len(self.ptotbio) - 1):
            for i in range(self.ptotbio[j]+1, self.ptotbio[j + 1] + 1):
                self.es_plus_1.append(0)
                self.rhos_plus_1.append(0)
                for component in self.components:
                    self.es_plus_1[-1]   = self.es_plus_1[-1]   + Fis[component.name][i] * component.e
                    self.rhos_plus_1[-1] = self.rhos_plus_1[-1] + Fis[component.name][i] * component.rho

        self.esL_plus_1     = [e for e in self.es_plus_1]
        self.rhosL_plus_1   = [rho for rho in self.rhos_plus_1]

        for j in range(len(self.ptotbio) - 1):
            self.esL_plus_1[self.ptotbio[j]]   = 0
            self.rhosL_plus_1[self.ptotbio[j]] = 0
            for component in self.components:
                self.esL_plus_1[self.ptotbio[j]]   = self.esL_plus_1[self.ptotbio[j]]   + FisL[component.name][self.ptotbio[j]] * component.e
                self.rhosL_plus_1[self.ptotbio[j]] = self.rhosL_plus_1[self.ptotbio[j]] + FisL[component.name][self.ptotbio[j]] * component.rho

    def make_grid_Rs(self, Cn, CnL, Fis_plus_1, FisL_plus_1):

        """Makes the list of values of "R" for the finite difference equations
        for the grid "z" with boundary points "p" and concentrations "Cn." """


        num     = self.nchemicals

        self.RsL_plus_1  = []
        self.SsL_plus_1  = []
        self.KsL_plus_1  = []

        layer  = self.layertotbio[0]
        for n in range(num):
            chemical = self.chemicals[n]
            Kd_rho   = 0
            for component in layer.components:
                Kd_rho = Kd_rho + FisL_plus_1[component.name][0] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[self.cptotbio[0]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
            self.RsL_plus_1.append(self.esL_plus_1[0] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)) + Kd_rho )
            self.SsL_plus_1.append(self.esL_plus_1[0] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
            self.KsL_plus_1.append(Kd_rho)
        for solidchemical in layer.solidchemicals:
            self.RsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][0] * layer.components[solidchemical.component_index].rho)
            self.SsL_plus_1.append(0)
            self.KsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][0] * layer.components[solidchemical.component_index].rho)

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for ii in range(self.ptotbio[j + 1] - self.ptotbio[j] - 1):
                i = ii + 1
                for n in range(num):
                    chemical = self.chemicals[n]
                    Kd_rho   = 0
                    for component in layer.components:
                        Kd_rho = Kd_rho + FisL_plus_1[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[self.cptotLbio[j]+i*layer.nchemicals+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                    self.RsL_plus_1.append(self.esL_plus_1[self.ptotbio[j]+i] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)) + Kd_rho )
                    self.SsL_plus_1.append(self.esL_plus_1[self.ptotbio[j]+i] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
                    self.KsL_plus_1.append(Kd_rho)
                for solidchemical in layer.solidchemicals:
                    self.RsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j]+i] * layer.components[solidchemical.component_index].rho)
                    self.SsL_plus_1.append(0)
                    self.KsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j]+i] * layer.components[solidchemical.component_index].rho)

            if j < len(self.ptotbio) - 2:
                layer  = self.layertotbio[j+1]
                for n in range(num):
                    chemical = self.chemicals[n]
                    Kd_rho   = 0
                    for component in layer.components:
                        Kd_rho = Kd_rho + FisL_plus_1[component.name][self.ptotbio[j+1]] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[self.cptotbio[j+1]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                    self.RsL_plus_1.append(self.esL_plus_1[self.ptotbio[j+1]] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)) + Kd_rho )
                    self.SsL_plus_1.append(self.esL_plus_1[self.ptotbio[j+1]] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
                    self.KsL_plus_1.append(Kd_rho)
                for solidchemical in layer.solidchemicals:
                    self.RsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j+1]] * layer.components[solidchemical.component_index].rho)
                    self.SsL_plus_1.append(0)
                    self.KsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j+1]] * layer.components[solidchemical.component_index].rho)
                for solidchemical in self.layertotbio[j].lowersolidchemicals:
                    self.RsL_plus_1.append(Fis_plus_1[self.layertotbio[j].components[solidchemical.component_index].name][self.ptotbio[j+1]] * self.layertotbio[j].components[solidchemical.component_index].rho)
                    self.SsL_plus_1.append(0)
                    self.KsL_plus_1.append(Fis_plus_1[self.layertotbio[j].components[solidchemical.component_index].name][self.ptotbio[j+1]] * self.layertotbio[j].components[solidchemical.component_index].rho)
            else:
                layer  = self.layertotbio[-1]
                for n in range(num):
                    chemical = self.chemicals[n]
                    Kd_rho   = 0
                    for component in layer.components:
                        Kd_rho = Kd_rho + FisL_plus_1[component.name][-1] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[self.cptotbio[-1]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)

                    self.RsL_plus_1.append(self.esL_plus_1[-1] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)) + Kd_rho )
                    self.SsL_plus_1.append(self.esL_plus_1[-1] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
                    self.KsL_plus_1.append(Kd_rho)
                for solidchemical in layer.solidchemicals:
                    self.RsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][-1] * layer.components[solidchemical.component_index].rho)
                    self.SsL_plus_1.append(0)
                    self.KsL_plus_1.append(FisL_plus_1[layer.components[solidchemical.component_index].name][-1] * layer.components[solidchemical.component_index].rho)

        self.Rs_plus_1  = [R for R in self.RsL_plus_1]
        self.Ss_plus_1  = [S for S in self.SsL_plus_1]
        self.Ks_plus_1  = [K for K in self.KsL_plus_1]

        for j in range(1, len(self.ptotbio)-1):
            layer  = self.layertotbio[j-1]
            for n in range(num):
                chemical = self.chemicals[n]
                Kd_rho   = 0
                for component in layer.components:
                    Kd_rho = Kd_rho + Fis_plus_1[component.name][self.ptotbio[j]] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotbio[j]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                self.Rs_plus_1[self.cptotbio[j]+n]=   self.es_plus_1[self.ptotbio[j]] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc)) + Kd_rho
                self.Ss_plus_1[self.cptotbio[j]+n]=   self.es_plus_1[self.ptotbio[j]] *(1 + layer.doc/(10**6) * 10**(chemical.Kdoc))
                self.Ks_plus_1[self.cptotbio[j]+n]=   Kd_rho
            for n in range(self.layertotbio[j].nsolidchemicals):
                solidchemicals = self.layertotbio[j].solidchemicals[n]
                if layer.component_list.count(solidchemicals.component_name) > 0:
                    self.Rs_plus_1[self.cptotbio[j]+num+n] = Fis_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j]] * layer.components[solidchemical.component_index].rho
                    self.Ss_plus_1[self.cptotbio[j]+num+n] = 0
                    self.Ks_plus_1[self.cptotbio[j]+num+n] = Fis_plus_1[layer.components[solidchemical.component_index].name][self.ptotbio[j]] * layer.components[solidchemical.component_index].rho
                else:
                    self.Rs_plus_1[self.cptotbio[j]+num+n] = 0
                    self.Ss_plus_1[self.cptotbio[j]+num+n] = 0
                    self.Ks_plus_1[self.cptotbio[j]+num+n] = 0

    def make_grid_Us(self):
        """Makes the list of values of "D" for the finite difference equations
        for the grid "z" with boundary points "p." """

        self.Us_plus_1  = []

        layer  = self.layertotbio[0]
        for chemical in (self.chemicals):
            self.Us_plus_1.append(   self.U_plus_1 * (1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
        for n in range(self.layertotbio[0].nsolidchemicals):
            self.Us_plus_1.append(0)

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for i in range(self.ptotbio[j] + 1, self.ptotbio[j + 1]):
                for chemical in (self.chemicals):
                    self.Us_plus_1.append(self.U_plus_1 * (1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
                for n in range(self.layertotbio[j].nsolidchemicals):
                    self.Us_plus_1.append(0)
            if j < len(self.ptotbio) - 2:
                for chemical in (self.chemicals):
                    self.Us_plus_1.append(self.U_plus_1 * (1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
                for n in range(self.layertotbio[j+1].nsolidchemicals):
                    self.Us_plus_1.append(0)
                for n in range(self.layertotbio[j].nlowersolidchemicals):
                    self.Us_plus_1.append(0)

        layer  = self.layertotbio[-1]
        for chemical in (self.chemicals):
            self.Us_plus_1.append(self.U_plus_1 * (1 + layer.doc/(10**6) * 10**(chemical.Kdoc)))
        for n in range(self.layertotbio[-1].nsolidchemicals):
            self.Us_plus_1.append(0)

        self.UsL_plus_1  = [U for U in self.Us_plus_1]
        for j in range(len(self.ptotbio) - 1):
            if j < len(self.ptotbio) - 2:
                layer  = self.layertotbio[j+1]
                for chemical in (self.chemicals):
                    n = self.chemicals.index(chemical)
                    self.UsL_plus_1[self.cptotbio[j]+n]=self.U_plus_1 * (1 + layer.doc/(10**6) * 10**(chemical.Kdoc))


    def make_grid_Dbiops(self):
        """Makes the list of values of "D" for the finite difference equations
        for the grid "z" with boundary points "p." """

        self.Dbiops_plus_1  = []
        self.Dbiops_plus_1.append(0)

        for j in range(len(self.ptotbio) - 1):
            for i in range(self.ptotbio[j]+1, self.ptotbio[j + 1]+1):
                self.Dbiops_plus_1.append(0)

        if self.biotype == 'Depth-dependent':
            for n in range (self.pbio + 1):
                self.Dbiops_plus_1[n] = self.Dbiops_plus_1[n] + self.Dbiop * exp(-(self.z[n]-self.z[0])**2/2/self.sigma**2)
        else:
            for n in range (self.pbio + 1):
                self.Dbiops_plus_1[n] = self.Dbiops_plus_1[n] + self.Dbiop

        self.DbiopsL_plus_1  = [Dbiop for Dbiop in self.Dbiops_plus_1]
        self.DbiopsL_plus_1[self.pbio] = 0


    def make_grid_Ds(self):
        """Makes the list of values of "D" for the finite difference equations
        for the grid "z" with boundary points "p." """

        self.DsL_plus_1   = []

        Dn = 0

        layer  = self.layertotbio[0]
        for chemical in (self.chemicals):
            self.DsL_plus_1.append(layer.get_D(chemical.Dw, self.UsL_plus_1[Dn], self.esL_plus_1[0]))
            Dn = Dn + 1
        for n in range(self.layertotbio[0].nsolidchemicals):
            self.DsL_plus_1.append(layer.get_D(0, 0, self.esL_plus_1[0]))
            Dn = Dn + 1

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for i in range(self.ptotbio[j] + 1, self.ptotbio[j + 1]):
                for chemical in (self.chemicals):
                    self.DsL_plus_1.append(layer.get_D(chemical.Dw, self.UsL_plus_1[Dn], self.esL_plus_1[i]))
                    Dn = Dn + 1
                for n in range(self.layertotbio[j].nsolidchemicals):
                    self.DsL_plus_1.append(layer.get_D(0, 0, self.esL_plus_1[i]))
                    Dn = Dn + 1

            if j < len(self.ptotbio) - 2:
                layer  = self.layertotbio[j+1]
                for chemical in (self.chemicals):
                    self.DsL_plus_1.append(layer.get_D(chemical.Dw, self.UsL_plus_1[Dn], self.esL_plus_1[self.ptotbio[j+1]]))
                    Dn = Dn + 1
                for n in range(layer.nsolidchemicals):
                    self.DsL_plus_1.append(layer.get_D(0, 0, self.esL_plus_1[self.ptotbio[j+1]]))
                    Dn = Dn + 1
                for n in range(self.layertotbio[j].nlowersolidchemicals):
                    self.DsL_plus_1.append(layer.get_D(0, 0, self.esL_plus_1[self.ptotbio[j+1]]))
                    Dn = Dn + 1
            else:
                layer  = self.layertotbio[-1]
                for chemical in (self.chemicals):
                    self.DsL_plus_1.append(layer.get_D(chemical.Dw, self.UsL_plus_1[Dn], self.esL_plus_1[-1]))
                    Dn = Dn + 1
                for n in range(self.layertotbio[-1].nsolidchemicals):
                    self.DsL_plus_1.append(layer.get_D(0, 0, self.esL_plus_1[-1]))
                    Dn = Dn + 1

        if self.bio == 1:
            if self.biotype == 'Uniform':
                for n in range (self.cpbio):
                    self.DsL_plus_1[n] = self.DsL_plus_1[n] + self.Dbiopw * self.Ss_plus_1[n]
            else:
                for j in range(len(self.ptotbio) - 1):
                    for ii in range(self.ptotbio[j + 1] - self.ptotbio[j] - 1):
                        i = ii + 1
                        for n in range(len(self.layertotbio[j].chemicals)):
                            self.DsL_plus_1[self.cptotLbio[j]+i*self.layertotbio[j].nchemicals+n] = self.DsL_plus_1[self.cptotLbio[j]+i*self.layertotbio[j].nchemicals+n] + self.Dbiopw * self.SsL_plus_1[self.cptotLbio[j]+i*self.layertotbio[j].nchemicals+n]* exp(-(self.z[self.ptotbio[j]+i]-self.z[0])**2/2/self.sigma**2)
                    for n in range(len(self.layertotbio[j].chemicals)):
                        self.DsL_plus_1[self.cptotbio[j]+n] = self.DsL_plus_1[self.cptotbio[j]+n] + self.Dbiopw * self.SsL_plus_1[self.cptotbio[j]+n]* exp(-(self.z[self.ptotbio[j]]-self.z[0])**2/2/self.sigma**2)
                for n in range(len(self.layertotbio[-1].chemicals)):
                    self.DsL_plus_1[self.cptotbio[-1]+n] = self.DsL_plus_1[self.cptotbio[-1]+n] + self.Dbiopw * self.SsL_plus_1[self.cptotbio[-1]+n]* exp(-(self.z[self.ptotbio[-1]]-self.z[0])**2/2/self.sigma**2)

        self.Ds_plus_1  = [D for D in self.DsL_plus_1]

        for j in range(1, len(self.ptotbio) - 1):
            layer  = self.layertotbio[j-1]
            for n in range(self.nchemicals):
                chemical = self.chemicals[n]
                self.Ds_plus_1[self.cptotbio[j]+n] = layer.get_D(chemical.Dw, self.Us_plus_1[self.ptotbio[j]], self.es_plus_1[self.ptotbio[j]])
            for n in range(self.layertotbio[j].nsolidchemicals+self.layertotbio[j-1].nlowersolidchemicals):
                self.Ds_plus_1[self.cptotbio[j]+self.nchemicals+n] = layer.get_D(0, 0, 1)

        if self.bio == 1:
            if self.biotype == 'Uniform':
                for j in range(1, len(self.ptotbio) - 1):
                    if self.ptotbio[j] <= self.pbio:
                        for n in range(layer.nchemicals):
                            self.Ds_plus_1[self.cptotbio[j]+n] = self.Ds_plus_1[self.cptotbio[j]+n] + self.Dbiopw * self.Ss_plus_1[self.cptotbio[j]+n]
            else:
                for j in range(1, len(self.ptotbio) - 1):
                    for n in range(layer.nchemicals):
                        self.Ds_plus_1[self.cptotbio[j]+n] = self.Ds_plus_1[self.cptotbio[j]+n] + self.Dbiopw * self.Ss_plus_1[self.cptotbio[j]+n]* exp(-(self.z[self.ptotbio[j]]-self.z[0])**2/2/self.sigma**2)* exp(-(self.z[self.ptotbio[j]]-self.z[0])**2/2/self.sigma**2)


    def make_grid_elams(self, Fis, FisL):
        """Makes the list of values of "elam" for the finite difference 
        equations for the grid "z" with boundary points "p." """

        self.elams_plus_1  = zeros([self.cptotbio[-1] + self.layertotbio[-1].nchemicals, self.cptotbio[-1] + self.layertotbio[-1].nchemicals])

        j = 0
        i = 0
        layer  = self.layertotbio[0]
        for reaction in layer.lreactions:
            n = reaction.key.count
            if reaction.name.count('_sorption') == 1:
                component_name = reaction.component
                e = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j]+i]
                for chemical in (reaction.reactants + reaction.products):
                    nn = chemical.count
                    self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n] = (self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n]
                                                                                    +chemical.coef * self.coefficients[layer.name][reaction.name].lam * e* (self.z[1]-self.z[0])/2)
            elif reaction.name.count('_desorption') == 1:
                component_name = reaction.component
                rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j]+i]
                for chemical in (reaction.reactants + reaction.products):
                    nn = chemical.count
                    self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n] = (self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n]
                                                                                    + chemical.coef * self.coefficients[layer.name][reaction.name].lam * rho* (self.z[1]-self.z[0])/2)
            else:
                for chemical in (reaction.reactants + reaction.products):
                    nn = chemical.count
                    self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n] = self.elams_plus_1[self.cptotbio[j]+nn, self.cptotbio[j]+n] + chemical.coef * self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j]+i]* (self.z[1]-self.z[0])/2

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for ii in range(self.ptotbio[j + 1] - self.ptotbio[j] - 1):
                i = ii + 1
                for reaction in layer.lreactions:
                    n = reaction.key.count
                    if reaction.name.count('_sorption') == 1:
                        component_name = reaction.component
                        e = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j]+i]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n] = (self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n]
                                                                                                                                    +chemical.coef * self.coefficients[layer.name][reaction.name].lam * e)
                    elif reaction.name.count('_desorption') == 1:
                        component_name = reaction.component
                        rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j]+i]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n] = (self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n]
                                                                                                                                 + chemical.coef * self.coefficients[layer.name][reaction.name].lam * rho)
                    else:
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n] = self.elams_plus_1[self.cptotLbio[j]+i*layer.nchemicals+nn, self.cptotLbio[j]+i*layer.nchemicals+n] \
                                                                                                                                                    + chemical.coef * self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j]+i]
            if j < len(self.ptotbio) - 2:

                layer  = self.layertotbio[j]
                for reaction in layer.lreactions:
                    if reaction.name.count('_sorption') == 1:
                        n = reaction.key.count
                        component_name = reaction.component
                        e = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j+1]]
                        for chemical in (reaction.reactants + reaction.products):
                            if chemical.count < self.nchemicals:    nn = chemical.count
                            else:                                   nn = layer.lowersolidchemicalnums[chemical.count-self.nchemicals]
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n]
                                                                                              +chemical.coef * self.coefficients[layer.name][reaction.name].lam * e * (self.z[self.ptotbio[j+1]]-self.z[self.ptotbio[j+1]-1])/2)
                    elif reaction.name.count('_desorption') == 1:
                        n = layer.lowersolidchemicalnums[reaction.key.count-self.nchemicals]
                        component_name = reaction.component
                        rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j+1]]
                        for chemical in (reaction.reactants + reaction.products):
                            if chemical.count < self.nchemicals:    nn = chemical.count
                            else:                                   nn = layer.lowersolidchemicalnums[chemical.count-self.nchemicals]
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n]
                                                                                                + chemical.coef * self.coefficients[layer.name][reaction.name].lam * rho * (self.z[self.ptotbio[j+1]]-self.z[self.ptotbio[j+1]-1])/2)
                    else:
                        n = reaction.key.count
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] +
                                                                                              chemical.coef * self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j+1]]* (self.z[self.ptotbio[j+1]]-self.z[self.ptotbio[j+1]-1])/2)

                layer  = self.layertotbio[j+1]
                for reaction in layer.lreactions:
                    n = reaction.key.count
                    if reaction.name.count('_sorption') == 1:
                        component_name = reaction.component
                        e = self.components[self.component_list.index(component_name)].e * FisL[component_name][self.ptotbio[j+1]]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n]
                                                                                              +chemical.coef * self.coefficients[layer.name][reaction.name].lam * e * (self.z[self.ptotbio[j+1]+1]-self.z[self.ptotbio[j+1]])/2)
                    elif reaction.name.count('_desorption') == 1:
                        component_name = reaction.component
                        rho = self.components[self.component_list.index(component_name)].rho * FisL[component_name][self.ptotbio[j+1]]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n]
                                                                                                + chemical.coef * self.coefficients[layer.name][reaction.name].lam * rho * (self.z[self.ptotbio[j+1]+1]-self.z[self.ptotbio[j+1]])/2)
                    else:
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] = (self.elams_plus_1[self.cptotbio[j+1]+nn, self.cptotbio[j+1]+n] +
                                                                                              chemical.coef * self.coefficients[layer.name][reaction.name].lam * self.esL_plus_1[self.ptotbio[j+1]]* (self.z[self.ptotbio[j+1]+1]-self.z[self.ptotbio[j+1]])/2)

            else:
                layer  = self.layertotbio[-1]
                for reaction in layer.lreactions:
                    n = reaction.key.count
                    if reaction.name.count('_sorption') == 1:
                        component_name = reaction.component
                        e = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[-1]]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n] = (self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n]
                                                                                            +chemical.coef * self.coefficients[layer.name][reaction.name].lam * e* (self.z[self.ptotbio[-1]]-self.z[self.ptotbio[-1]-1])/2)
                    elif reaction.name.count('_desorption') == 1:
                        component_name = reaction.component
                        rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[-1]]
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n] = (self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n]
                                                                                        + chemical.coef * self.coefficients[layer.name][reaction.name].lam * rho* (self.z[self.ptotbio[-1]]-self.z[self.ptotbio[-1]-1])/2)
                    else:
                        for chemical in (reaction.reactants + reaction.products):
                            nn = chemical.count
                            self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n] = (self.elams_plus_1[self.cptotbio[-1]+nn, self.cptotbio[-1]+n]
                                                                                            + chemical.coef * self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[-1]]* (self.z[self.ptotbio[-1]]-self.z[self.ptotbio[-1]-1])/2)

    def make_grid_rates_plus_1(self, Cn, CnL, Fis, FisL):

        """Makes the list of values of "rates" for the finite difference
        equations for the grid "z" with boundary points "p." """

        self.rates_plus_1      = zeros(self.cptot[-1] + self.layertotbio[-1].nchemicals)

        layer  = self.layertotbio[0]
        for reaction in layer.nlreactions:
            if reaction.model == 'Langmuir':
                component_name = reaction.component
                e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][0]
                rate = (self.coefficients[layer.name][reaction.name].lam * Cn[reaction.reactants[0].count] * e
                        *(reaction.reactants[0].index - Cn[reaction.products[0].count]))
            else:
                if reaction.name.count('_desorption') == 1:
                    component_name = reaction.component
                    rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][0]
                    rate = self.coefficients[layer.name][reaction.name].lam * rho
                    for reactant in reaction.reactants:
                        if     Cn[reactant.count] <= 0:    rate = 0
                        else:                              rate = rate * Cn[reactant.count]**reactant.index
                else:
                    rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[0]
                    for reactant in reaction.reactants:
                        if     Cn[reactant.count] <= 0:    rate = 0
                        else:                              rate = rate * Cn[reactant.count]**reactant.index

            for chemical in (reaction.reactants + reaction.products):
                self.rates_plus_1[chemical.count] = self.rates_plus_1[chemical.count] + rate * chemical.coef * (self.z[1]-self.z[0])/2

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for ii in range(self.ptotbio[j + 1] - self.ptotbio[j]-1):
                i = ii + 1
                for reaction in layer.nlreactions:
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j]+i]
                        rate = (self.coefficients[layer.name][reaction.name].lam * Cn[self.cptotLbio[j]+i*layer.nchemicals+reaction.reactants[0].count] * e
                                *(reaction.reactants[0].index - Cn[self.cptotLbio[j]+i*layer.nchemicals+reaction.products[0].count]))
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j]+i]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho
                            for reactant in reaction.reactants:
                                if     Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count] <= 0:    rate = 0
                                else:                                                                   rate = rate * Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count]**reactant.index
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j]+i]
                            for reactant in reaction.reactants:
                                if     Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count] <= 0:    rate = 0
                                else:                                                                   rate = rate * Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count]**reactant.index

                    for chemical in (reaction.reactants + reaction.products):
                        self.rates_plus_1[self.cptotLbio[j]+i*layer.nchemicals+chemical.count]      = self.rates_plus_1[self.cptotLbio[j]+i*layer.nchemicals+chemical.count] + rate * chemical.coef

            if j < len(self.ptotbio) - 2:
                layer  = self.layertotbio[j+1]
                for reaction in layer.nlreactions:
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * FisL[component_name][self.ptotbio[j+1]]
                        rate = (self.coefficients[layer.name][reaction.name].lam * CnL[self.cptotbio[j+1]+reaction.reactants[0].count] * e
                                *(reaction.reactants[0].index - CnL[self.cptotbio[j+1]+reaction.products[0].count]))
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * FisL[component_name][self.ptotbio[j+1]]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho
                            for reactant in reaction.reactants:
                                if     CnL[self.cptotbio[j+1]+reactant.count] <= 0:    rate = 0
                                else:                                                  rate = rate * CnL[self.cptotbio[j+1]+reactant.count]**reactant.index
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j+1]]
                            for reactant in reaction.reactants:
                                if     CnL[self.cptotbio[j+1]+reactant.count] <= 0:    rate = 0
                                else:                                                  rate = rate * CnL[self.cptotbio[j+1]+reactant.count]**reactant.index
                    for chemical in (reaction.reactants + reaction.products):
                        self.rates_plus_1[self.cptotbio[j+1]+chemical.count]      = self.rates_plus_1[self.cptotbio[j+1]+chemical.count] + rate * chemical.coef * (self.z[self.ptotbio[j+1]+1]-self.z[self.ptotbio[j+1]])/2

                layer  = self.layertotbio[j]
                for reaction in layer.nlreactions:
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j+1]]
                        rate = (self.coefficients[layer.name][reaction.name].lam * Cn[self.cptotbio[j+1]+reaction.reactants[0].count] * e
                                *(reaction.reactants[0].index - Cn[self.cptotbio[j+1]+reaction.products[0].count]))
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j+1]]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho
                            for reactant in reaction.reactants:
                                nn = layer.lowersolidchemicalnums[reactant.count - self.nchemicals]
                                if     Cn[self.cptotbio[j+1]+reactant.count] <= 0:    rate = 0
                                else:                                                 rate = rate * Cn[self.cptotbio[j+1]+nn]**reactant.index
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j+1]]
                            for reactant in reaction.reactants:
                                if     Cn[self.cptotbio[j+1]+reactant.count] <= 0:    rate = 0
                                else:                                                 rate = rate * Cn[self.cptotbio[j+1]+reactant.count]**reactant.index

                    for chemical in (reaction.reactants + reaction.products):
                        self.rates_plus_1[self.cptotbio[j+1]+chemical.count]      = self.rates_plus_1[self.cptotbio[j+1]+chemical.count] + rate * chemical.coef * (self.z[self.ptotbio[j+1]]-self.z[self.ptotbio[j+1]-1])/2


        j = -1
        i = 0
        layer  = self.layertotbio[j]
        for reaction in layer.nlreactions:
            if reaction.model == 'Langmuir':
                component_name = reaction.component
                e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j]+i]
                rate = (self.coefficients[layer.name][reaction.name].lam * Cn[self.cptotbio[j]+i*layer.nchemicals+reaction.reactants[0].count] * e
                        *(reaction.reactants[0].index - Cn[self.cptotbio[j]+i*layer.nchemicals+reaction.products[0].count]))
            else:
                if reaction.name.count('_desorption') == 1:
                    component_name = reaction.component
                    rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j]+i]
                    rate = self.coefficients[layer.name][reaction.name].lam * rho
                    for reactant in reaction.reactants:
                        if     Cn[self.cptotbio[j]+i*layer.nchemicals+reactant.count] <= 0:    rate = 0
                        else:                                                               rate = rate * Cn[self.cptotbio[j]+i*layer.nchemicals+reactant.count]**reactant.index
                else:
                    rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j]+i]
                    for reactant in reaction.reactants:
                        if     Cn[self.cptotbio[j]+i*layer.nchemicals+reactant.count] <= 0:    rate = 0
                        else:                                                               rate = rate * Cn[self.cptotbio[j]+i*layer.nchemicals+reactant.count]**reactant.index

            for chemical in (reaction.reactants + reaction.products):
                self.rates_plus_1[self.cptotbio[j]+i*layer.nchemicals+chemical.count]      = self.rates_plus_1[self.cptotbio[j]+i*layer.nchemicals+chemical.count] + rate * chemical.coef * (self.z[-1]-self.z[-2])/2

    def make_grid_rates_diff(self, Cn, CnL, Fis, FisL):

        """Makes the list of values of "rates differentiation" for the newton's medthod finite difference
        equations for the grid "z" with boundary points "p." """

        self.rates_diff  = zeros([self.cptotbio[-1] + self.layertotbio[-1].nchemicals, self.cptotbio[-1] + self.layertotbio[-1].nchemicals])

        layer  = self.layertotbio[0]
        for reaction in layer.nlreactions:
            diff = zeros(layer.nchemicals)
            if reaction.model == 'Langmuir':
                component_name = reaction.component
                e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][0]
                diff[reaction.reactants[0].count] = diff[reaction.reactants[0].count] + self.coefficients[layer.name][reaction.name].lam * e * (reaction.reactants[0].index - Cn[reaction.products[0].count])

                diff[reaction.products[0].count]  = diff[reaction.products[0].count] - self.coefficients[layer.name][reaction.name].lam * e * Cn[reaction.reactants[0].count]
            else:
                if reaction.name.count('_desorption') == 1:
                    component_name = reaction.component
                    rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][0]
                    rate = self.coefficients[layer.name][reaction.name].lam * rho
                else:
                    rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[0]

                for reactant in reaction.reactants:
                    if  Cn[reactant.count] <= 0:    rate = 0
                    else:                           rate = rate * Cn[reactant.count]**reactant.index

                for reactant in reaction.reactants:
                    if  Cn[reactant.count] >0:
                        diff[reactant.count] = rate / Cn[reactant.count] * reactant.index

            for chemical in (reaction.reactants + reaction.products):
                for chemical_a in (reaction.reactants + reaction.products):
                    self.rates_diff[chemical.count, chemical_a.count] = self.rates_diff[chemical.count, chemical_a.count] + diff[chemical_a.count] * chemical.coef * (self.z[1]-self.z[0])/2

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            for ii in range(self.ptotbio[j + 1] - self.ptotbio[j] - 1):
                i = ii + 1
                for reaction in layer.nlreactions:
                    diff = zeros(layer.nchemicals)
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j]+i]
                        diff[reaction.reactants[0].count] = diff[reaction.reactants[0].count] + self.coefficients[layer.name][reaction.name].lam * e * (reaction.reactants[0].index - Cn[self.cptotLbio[j]+i*layer.nchemicals+reaction.products[0].count])

                        diff[reaction.products[0].count]  = diff[reaction.products[0].count] - self.coefficients[layer.name][reaction.name].lam * e * Cn[self.cptotLbio[j]+i*layer.nchemicals+reaction.reactants[0].count]
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j]+i]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j]+i]

                        for reactant in reaction.reactants:
                            if     Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count] <= 0:    rate = 0
                            else:                                                                   rate = rate * Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count]**reactant.index
                        for reactant in reaction.reactants:
                            if  Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count] > 0:
                                diff[reactant.count] = rate / Cn[self.cptotLbio[j]+i*layer.nchemicals+reactant.count] * reactant.index
                    for chemical in (reaction.reactants + reaction.products):
                        for chemical_a in (reaction.reactants + reaction.products):
                            self.rates_diff[self.cptotLbio[j]+i*layer.nchemicals+chemical.count, self.cptotLbio[j]+i*layer.nchemicals+chemical_a.count]      \
                                = self.rates_diff[self.cptotLbio[j]+i*layer.nchemicals+chemical.count, self.cptotLbio[j]+i*layer.nchemicals+chemical_a.count] + diff[chemical_a.count] * chemical.coef

            if j < len(self.ptotbio) - 2:
                layer  = self.layertotbio[j+1]
                for reaction in layer.nlreactions:
                    diff = zeros(layer.nchemicals)
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * FisL[component_name][self.ptotbio[j+1]]
                        diff[reaction.reactants[0].count] = diff[reaction.reactants[0].count] + self.coefficients[layer.name][reaction.name].lam * e * (reaction.reactants[0].index - CnL[self.cptotbio[j+1]+reaction.products[0].count])

                        diff[reaction.products[0].count]  = diff[reaction.products[0].count] - self.coefficients[layer.name][reaction.name].lam * e * CnL[self.cptotbio[j+1]+reaction.reactants[0].count]
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * FisL[component_name][self.ptotbio[j+1]]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j+1]]

                        for reactant in reaction.reactants:
                            if     CnL[self.cptotbio[j+1]+reactant.count] <= 0:                 rate = 0
                            else:                                                               rate = rate * CnL[self.cptotbio[j+1]+reactant.count]**reactant.index
                        for reactant in reaction.reactants:
                            if  CnL[self.cptotbio[j+1]+reactant.count] > 0:
                                diff[reactant.count] = rate / CnL[self.cptotbio[j+1]+reactant.count] * reactant.index
                    for chemical in (reaction.reactants + reaction.products):
                        for chemical_a in (reaction.reactants + reaction.products):
                            self.rates_diff[self.cptotbio[j+1]+chemical.count, self.cptotbio[j+1]+chemical_a.count]      \
                                = self.rates_diff[self.cptotbio[j+1]+chemical.count, self.cptotbio[j+1]+chemical_a.count] + diff[chemical_a.count] * chemical.coef * (self.z[self.ptotbio[j+1]+1]-self.z[self.ptotbio[j+1]])/2

                layer  = self.layertotbio[j]
                for reaction in layer.nlreactions:
                    diff = zeros(layer.nchemicals)
                    if reaction.model == 'Langmuir':
                        component_name = reaction.component
                        e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[j+1]]
                        diff[reaction.reactants[0].count] = diff[reaction.reactants[0].count] + self.coefficients[layer.name][reaction.name].lam * e * (reaction.reactants[0].index - Cn[self.cptotbio[j+1]+reaction.products[0].count])

                        diff[reaction.products[0].count]  = diff[reaction.products[0].count] - self.coefficients[layer.name][reaction.name].lam * e * Cn[self.cptotbio[j+1]+reaction.reactants[0].count]
                    else:
                        if reaction.name.count('_desorption') == 1:
                            component_name = reaction.component
                            rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[j+1]]
                            rate = self.coefficients[layer.name][reaction.name].lam * rho

                            for reactant in reaction.reactants:
                                nn = layer.lowersolidchemicalnums[reactant.count - self.nchemicals]
                                if     Cn[self.cptotbio[j+1]+nn] <= 0:                                 rate = 0
                                else:                                                                  rate = rate * Cn[self.cptotbio[j+1]+nn]**reactant.index
                            for reactant in reaction.reactants:
                                nn = layer.lowersolidchemicalnums[reactant.count - self.nchemicals]
                                if  Cn[self.cptotbio[j+1]+reactant.count] >0:
                                    diff[reactant.count] = rate / Cn[self.cptotbio[j+1]+nn] * reactant.index
                        else:
                            rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[j+1]]

                            for reactant in reaction.reactants:
                                nn = reactant.count
                                if     Cn[self.cptotbio[j+1]+nn] <= 0:                                 rate = 0
                                else:                                                                  rate = rate * Cn[self.cptotbio[j+1]+nn]**reactant.index
                            for reactant in reaction.reactants:
                                nn = reactant.count
                                if  Cn[self.cptotbio[j+1]+reactant.count] >0:
                                    diff[reactant.count] = rate / Cn[self.cptotbio[j+1]+nn] * reactant.index

                    for chemical in (reaction.reactants + reaction.products):
                        for chemical_a in (reaction.reactants + reaction.products):
                            self.rates_diff[self.cptotbio[j+1]+chemical.count, self.cptotbio[j+1]+chemical_a.count]      \
                                = self.rates_diff[self.cptotbio[j+1]+chemical.count, self.cptotbio[j+1]+chemical_a.count] + diff[chemical_a.count] * chemical.coef * (self.z[self.ptotbio[j+1]]-self.z[self.ptotbio[j+1]-1])/2

        layer  = self.layertotbio[-1]
        for reaction in layer.nlreactions:
            diff = zeros(layer.nchemicals)
            if reaction.model == 'Langmuir':
                component_name = reaction.component
                e   = self.components[self.component_list.index(component_name)].e * Fis[component_name][self.ptotbio[-1]]
                diff[reaction.reactants[0].count] = diff[reaction.reactants[0].count] + self.coefficients[layer.name][reaction.name].lam * e * (reaction.reactants[0].index - Cn[self.cptotbio[-1]+reaction.products[0].count])

                diff[reaction.products[0].count]  = diff[reaction.products[0].count] - self.coefficients[layer.name][reaction.name].lam * e * Cn[self.cptotbio[-1]+reaction.reactants[0].count]
            else:
                if reaction.name.count('_desorption') == 1:
                    component_name = reaction.component
                    rho = self.components[self.component_list.index(component_name)].rho * Fis[component_name][self.ptotbio[-1]]
                    rate = self.coefficients[layer.name][reaction.name].lam * rho
                else:
                    rate = self.coefficients[layer.name][reaction.name].lam * self.es_plus_1[self.ptotbio[-1]]

                for reactant in reaction.reactants:
                    if     Cn[self.cptotbio[-1]+reactant.count] <= 0:    rate = 0
                    else:                                                               rate = rate * Cn[self.cptotbio[-1]+reactant.count]**reactant.index
                for reactant in reaction.reactants:
                    if  Cn[self.cptotbio[-1]+reactant.count] >0:
                        diff[reactant.count] = rate / Cn[self.cptotbio[-1]+reactant.count] * reactant.index

            for chemical in (reaction.reactants + reaction.products):
                for chemical_a in (reaction.reactants + reaction.products):
                    self.rates_diff[self.cptotbio[-1]+chemical.count, self.cptotbio[-1]+chemical_a.count]      \
                        = self.rates_diff[self.cptotbio[-1]+chemical.count, self.cptotbio[-1]+chemical_a.count] + diff[chemical_a.count] * chemical_a.coef * (self.z[-1]-self.z[-2])/2

    def make_boundary_equations(self, Flag = None, Cn = None, CnL = None):
        """Makes the finite difference equations for the boundary conditions."""

        #top boundary
        num = self.layertotbio[0].nchemicals

        delz = self.z[1]-self.z[0]

        for n in range(self.nchemicals, self.layertotbio[0].nchemicals):
            solidchemical = self.layertotbio[0].chemicals[n]
            nn = solidchemical.kineticchemnum
            if self.Ks_plus_1[n+num] > 0 or self.Ks_plus_1[n] > 0:
                self.A[n,n]     = self.Ks_plus_1[n] * delz/2/self.delt
                self.B[n,n]     = self.Ks[n] * delz/2/self.delt
                if self.bio == 1 and self.Dbiop > 0:
                    self.A[n,n]     = self.A[n,n]     + self.Dbiops_plus_1[0]/delz*self.Ks_plus_1[n]
                    self.A[n,n+num] = self.A[n,n+num] - self.Dbiops_plus_1[0]/delz*self.Ks_plus_1[n+num]
                if self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Freundlich':
                    self.A[n,nn]    =   self.A[n,nn] - self.elams[n,nn]
                    self.a[n]       =   self.a[n]    - self.rates_plus_1[n]
                elif self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Langmuir':
                    self.A[n,n]     =   self.A[n,n]  - self.elams[n,n]
                    self.a[n]       =   self.a[n]    - self.rates_plus_1[n]
                else:
                    self.A[n,n]     =   self.A[n,n]  - self.elams[n,n]
                    self.A[n,nn]    =   self.A[n,nn] - self.elams[n,nn]
            else:
                self.A[n,n] = 1

        for n in range(self.nchemicals):
            chemical = self.layertotbio[0].chemicals[n]
            if self.topBCtype == 'Fixed Concentration':
                self.A[n,n], self.b[n] = 1, self.BCs[chemical.name].Co
            else:
                if self.bio == 1:
                    KDbiops_plus_1  = [self.DbiopsL_plus_1[0]  *self.Ks_plus_1[n], self.DbiopsL_plus_1[0]  *self.Ks_plus_1[num+n]]
                else:
                    KDbiops_plus_1 = 0

                [self.A[n,n], self.A[n,num+n]], self.b[n] = top_mass_transfer_boundary_kinetic(self.DsL_plus_1[n],  self.Us_plus_1[n], KDbiops_plus_1,  self.BCs[chemical.name].k, self.BCs[chemical.name].Cw_plus_1, self.z[0:2])

                self.A[n,n]     = self.A[n,n] + self.Rs_plus_1[n]   * delz/2/self.delt
                self.B[n,n]     = self.B[n,n] + self.Rs[n]          * delz/2/self.delt

                for np in range(num):
                    self.A[n,np]     =   self.A[n,np] - self.elams[n,np]

                if self.reac == 1:
                    self.a[n]       =   self.a[n] - self.rates_plus_1[n]

        #flux boundary conditions at interfaces
        for j in range(len(self.layertotbio) - 1):
            numa =self.layertotbio[j].nchemicals
            numb =self.layertotbio[j+1].nchemicals
            i  = self.cptotbio[j+1]
            iL = self.cptotLbio[j+1]
            p = self.ptotbio[j+1]
            delza = self.z[p]-self.z[p-1]
            #if self.bio == 1 and p == self.pbio:  delzb = 0
            delzb = self.z[p+1]-self.z[p]
            for n in range(self.layertotbio[j+1].nsolidchemicals):
                nb = n + self.nchemicals
                solidchemical = self.layertotbio[j+1].solidchemicals[n]
                nn = solidchemical.kineticchemnum
                if self.layertotbio[j].solidchemical_list.count(solidchemical.name) > 0 and self.bio == 1 and self.Dbiop > 0:
                    na = self.layertotbio[j].solidchemical_list.index(solidchemical.name) + self.nchemicals
                    KDbiops_plus_1  = [self.Dbiops[p]  *self.Ks_plus_1[i-numa+na],
                                       self.Dbiops[p]  *self.Ks_plus_1[i+nb]]
                    KDbiopsL_plus_1 = [self.DbiopsL[p] *self.KsL_plus_1[i+nb],
                                       self.DbiopsL[p] *self.KsL_plus_1[iL+numb+nb]]

                    [self.A[i+nb,i-numa+na], self.A[i+nb,i+nb], self.A[i+nb,iL+numb+nb]] = interface_boundary_kinetic(0,  0,  KDbiops_plus_1, 0,  0,  KDbiopsL_plus_1, self.z[p-1:p+2])

                self.A[i+nb,i+nb]      = self.A[i+nb,i+nb] + (self.Ks_plus_1[i+nb] * delza + self.KsL_plus_1[i+nb] * delzb)/2/self.delt
                if Flag == 1:
                    self.b[i+nb]           = self.b[i+nb] + (self.Ks[i+nb] * delza * Cn[i+nb] + self.KsL[i+nb] * delzb* CnL[i+nb])/2/self.delt
                else:
                    self.B[i+nb,i+nb]      = self.B[i+nb,i+nb] + (self.Ks[i+nb] * delza + self.KsL[i+nb] * delzb)/2/self.delt

                if self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Freundlich':
                    self.A[i+nb,i+nn]  =  self.A[i+nb,i+nn] - self.elams[i+nb,i+nn]
                    self.a[i+nb]       =  self.a[i+nb]      - self.rates_plus_1[i+nb]
                elif self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Langmuir':
                    self.A[i+nb,i+nb]  =  self.A[i+nb,i+nb] - self.elams[i+nb,i+nb]
                    self.a[i+nb]       =  self.a[i+nb]      - self.rates_plus_1[i+nb]
                else:
                    self.A[i+nb,i+nb]  =  self.A[i+nb,i+nb] - self.elams[i+nb,i+nb]
                    self.A[i+nb,i+nn]  =  self.A[i+nb,i+nn] - self.elams[i+nb,i+nn]

            for solidchemical in self.layertotbio[j].lowersolidchemicals:
                n  = self.layertotbio[j].solidchemicals.index(solidchemical)
                nb = self.layertotbio[j].lowersolidchemicalnums[n]
                if self.bio == 1 and self.Dbiop > 0:
                    KDbiops_plus_1  = [self.Dbiops[p]  *self.Ks_plus_1[i-numa+n],
                                       self.Dbiops[p]  *self.Ks_plus_1[i+nb]]
                    KDbiopsL_plus_1 = [0,0]

                    [self.A[i+nb,i-numa+n], self.A[i+nb,i+nb], self.A[i+nb,iL+numb+nb]] = interface_boundary_kinetic(self.Ds_plus_1[i+nb],  self.Us_plus_1[i+nb],  KDbiops_plus_1,
                                                                                                                     self.DsL_plus_1[i+nb], self.UsL_plus_1[i+nb], KDbiopsL_plus_1, self.z[p-1:p+2])

                self.A[i+nb,i+nb]      = self.A[i+nb,i+nb] + (self.Ks_plus_1[i+nb] * delza)/2/self.delt

                if Flag == 1:
                    self.b[i+nb]           = self.b[i+nb] + (self.Ks[i+nb] * delza * Cn[i+nb])/2/self.delt
                else:
                    self.B[i+nb,i+nb]      = self.B[i+nb,i+nb] + (self.Ks[i+nb] * delza)/2/self.delt

                if self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Freundlich':
                    self.A[i+nb,i+nn]  =  self.A[i+nb,i+nn] - self.elams[i+nb,i+nn]
                    self.a[i+nb]       =  self.a[i+nb]      - self.rates_plus_1[i+nb]
                elif self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Langmuir':
                    self.A[i+nb,i+nb]  =  self.A[i+nb,i+nb] - self.elams[i+nb,i+nb]
                    self.a[i+nb]       =  self.a[i+nb]      - self.rates_plus_1[i+nb]
                else:
                    self.A[i+nb,i+nb]  =  self.A[i+nb,i+nb] - self.elams[i+nb,i+nb]
                    self.A[i+nb,i+nn]  =  self.A[i+nb,i+nn] - self.elams[i+nb,i+nn]

            delzb = self.z[p+1]-self.z[p]
            for n in range(self.nchemicals):
                if self.bio == 1 and self.Dbiop > 0:
                    KDbiops_plus_1  = [self.Dbiops_plus_1[p]  *self.Ks_plus_1[ i-numa+n],
                                       self.Dbiops_plus_1[p]  *self.Ks_plus_1[ i+n]]
                    KDbiopsL_plus_1 = [self.DbiopsL_plus_1[p] *self.KsL_plus_1[i+n],
                                       self.DbiopsL_plus_1[p] *self.KsL_plus_1[iL+numb+n]]
                else:
                    KDbiops_plus_1  =  0
                    KDbiopsL_plus_1 =  0

                [self.A[i+n,i-numa+n], self.A[i+n,i+n], self.A[i+n,iL+numb+n]] = interface_boundary_kinetic( self.Ds_plus_1[i+n],  self.Us_plus_1[i+n],  KDbiops_plus_1,
                                                                                                             self.DsL_plus_1[i+n], self.UsL_plus_1[i+n], KDbiopsL_plus_1, self.z[p-1:p+2])

                self.A[i+n,i+n]      = self.A[i+n,i+n] + (self.Rs_plus_1[i+n] * delza + self.RsL_plus_1[i+n] *delzb)/2/self.delt
                if Flag == 1:
                    self.b[i+n]          = self.b[i+n] + (self.Rs[i+n] * delza * Cn[i+n] + self.RsL[i+n]*delzb* CnL[i+n])/2/self.delt
                else:
                    self.B[i+n,i+n]      = self.B[i+n,i+n] + (self.Rs[i+n] * delza        + self.RsL[i+n]        *delzb)/2/self.delt

                for np in range(numb):
                    self.A[i+n,i+np]  =  self.A[i+n,i+np] - self.elams[i+n,i+np]
                if self.reac == 1:
                    self.a[i+n] = self.a[i+n] - self.rates_plus_1[i+n]

        #bottom boundary
        num = self.layertotbio[-1].nchemicals
        i   = self.cptotbio[-1]
        p   = self.ptotbio[-1]
        for n in range(self.nchemicals, self.layertotbio[-1].nchemicals):
            solidchemical = self.layertotbio[-1].chemicals[n]
            nn = solidchemical.kineticchemnum
            self.A[i+n,i+n]     =  self.Ks_plus_1[i+n] * self.delz[-1]/2/self.delt
            self.B[i+n,i+n]     =  self.Ks[i+n] * self.delz[-1]/2/self.delt

            if self.bio == 1 and self.Dbiop > 0:
                self.A[i+n,i+n]     = self.A[i+n,i+n] +   self.Ks_plus_1[i+n] * self.Dbiops_plus_1[-1]/self.delz[-1]
                self.A[i+n,i-num+n] = self.A[i+n,i-num+n]-self.Ks_plus_1[i+n] * self.Dbiops_plus_1[-1]/self.delz[-1]

            if self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Freundlich':
                self.A[i+n,i+nn]    = self.A[i+n,i+nn]-self.elams[i+n,i+nn]
                self.a[i+n]         = self.a[i+n]     -self.rates_plus_1[i+n]
            elif self.sorptions[solidchemical.component_name][solidchemical.chemical_name].isotherm == 'Langmuir':
                self.A[i+n,i+n]     =  self.A[i+n,i+n] - self.elams[i+n,i+n]
                self.a[i+n]         = self.a[i+n] -self.rates_plus_1[i+n]
            else:
                self.A[i+n,i+n]     = self.A[i+n,i+n] - self.elams[i+n,i+n]
                self.A[i+n,i+nn]    = self.A[i+n,i+nn]- self.elams[i+n,i+nn]

        for n in range(self.nchemicals):
            chemical = self.chemicals[n]
            if self.botBCtype == 'Fixed Concentration':
                self.A[i+n,i+n], self.b[i+n] = 1, self.BCs[chemical.name].Cb
            else:
                if self.bio == 1 and self.ptotbio[-1] <= self.pbio:
                    KDbiops_plus_1  = [self.Dbiops_plus_1[p]  *self.Ks_plus_1[i-num+n],
                                       self.Dbiops_plus_1[p]  *self.Ks_plus_1[i+n]]
                else:
                    KDbiops_plus_1  =  0
                [self.A[i+n,i-num+n], self.A[i+n,i+n]], self.b[i+n] = flux_bottom_boundary_kinetic(self.Us_plus_1[i+n], self.Ds_plus_1[i+n], KDbiops_plus_1, self.BCs[chemical.name].Cb, self.z[-2:])

                self.A[i+n,i+n]     =  self.A[i+n,i+n] + self.Rs_plus_1[i+n] * self.delz[-1]/2/self.delt
                self.B[i+n,i+n]     =  self.B[i+n,i+n] + self.Rs[i+n] * self.delz[-1]/2/self.delt

                for np in range(num):     self.A[i+n,i+np]    =  self.A[i+n,i+np] - self.elams[i+n,i+np]

                if self.reac == 1:     self.a[i+n]         = self.a[i+n] - self.rates_plus_1[i+n]

    def make_components_equations(self):

        if self.timeoption == 'Implicit':
            for j in range(len(self.ptotbio) - 1):
                for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                    i = ii + 1
                    p = self.ptotbio[j] + i
                    if p < self.pbio:
                        Ds_plus_1      = [ self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4,
                                           self.DbiopsL_plus_1[p],
                                          -self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4]
                        [self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+(i-1)], self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+i], self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+(i+1)]] , \
                        [self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+(i-1)], self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+i], self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+(i+1)]] = \
                        get_3pt_adr_fde_imp(1, 1, 0, Ds_plus_1, [0, 0, 0], self.delt, self.z[p-1:p+2])
                    else:
                        self.Acomp[p, p], self.Bcomp[p, p] = 1, 1

                if j < len(self.ptotbio) - 2:
                    p = self.ptotbio[j+1]
                    if p < self.pbio:
                        self.Acomp[p,p-1], self.Acomp[p,p], self.Acomp[p,p+1] = interface_boundary_kinetic(self.Dbiops_plus_1[p],  0,  0, self.DbiopsL_plus_1[p], 0,  0, self.z[p-1:p+2])
                        self.Acomp[p,p]      = self.Acomp[p,p] + (self.z[p+1] - self.z[p-1])/2/self.delt
                        self.Bcomp[p,p]      = self.Bcomp[p,p] + (self.z[p+1] - self.z[p-1])/2/self.delt
                        #[self.Acomp[p,p-2], self.Acomp[p,p-1], self.Acomp[p,p], self.Acomp[p,p+1], self.Acomp[p,p+2]] = comp_interface_boundary(self.Dbiops_plus_1[p], 0, 0, self.DbiopsL_plus_1[p], 0, 0, self.z[p-2:p+3])

                    elif p == self.pbio:
                        [self.Acomp[p,p-1], self.Acomp[p,p]], self.bcomp[p] = flux_bottom_boundary_kinetic(0, self.Dbiops_plus_1[p], 0, 0, self.z[p-1:p+1])

                        self.Acomp[p,p] =  self.Acomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt
                        self.Bcomp[p,p] =  self.Bcomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt
                    else:
                        self.Acomp[p, p], self.Bcomp[p, p] = 1, 1
        else:
            for j in range(len(self.ptotbio) - 1):
                for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                    i = ii + 1
                    p = self.ptotbio[j] + i
                    if p < self.pbio:
                        Ds_plus_1      = [ self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4,
                                           self.DbiopsL_plus_1[p],
                                          -self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4]
                        Ds             = [ self.DbiopsL[p-1]/4+self.Dbiops[p] - self.Dbiops[p+1]/4,
                                           self.DbiopsL[p],
                                          -self.DbiopsL[p-1]/4+self.Dbiops[p] + self.Dbiops[p+1]/4]
                        [self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+(i-1)], self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+i], self.Acomp[self.ptotbio[j]+i,self.ptotbio[j]+(i+1)]] , \
                        [self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+(i-1)], self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+i], self.Bcomp[self.ptotbio[j]+i,self.ptotbio[j]+(i+1)]] = \
                        get_3pt_adr_fde_CN(1, 1, 0, 0, Ds_plus_1, Ds, [0, 0, 0], [0, 0, 0], self.delt, self.z[p-1:p+2])
                    else:
                        self.Acomp[p, p], self.Bcomp[p, p] = 1, 1

                if j < len(self.ptotbio) - 2:
                    p = self.ptotbio[j+1]
                    if p < self.pbio:
                        self.Acomp[p,p-1], self.Acomp[p,p], self.Acomp[p,p+1] = interface_boundary_kinetic(self.Dbiops_plus_1[p],  0,  0, self.DbiopsL_plus_1[p], 0,  0, self.z[p-1:p+2])
                        self.Acomp[p,p]      = self.Acomp[p,p] + (self.z[p+1] - self.z[p-1])/2/self.delt
                        self.Bcomp[p,p]      = self.Bcomp[p,p] + (self.z[p+1] - self.z[p-1])/2/self.delt
                    elif p == self.pbio:
                        [self.Acomp[p,p-1], self.Acomp[p,p]], self.bcomp[p] = flux_bottom_boundary_kinetic(0, self.Dbiops_plus_1[p], 0, 0, self.z[p-1:p+1])

                        self.Acomp[p,p] =  self.Acomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt
                        self.Bcomp[p,p] =  self.Bcomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt
                    else:
                        self.Acomp[p, p], self.Bcomp[p, p] = 1, 1

        [self.Acomp[0,0], self.Acomp[0,1]], self.bcomp[0] = top_mass_transfer_boundary_kinetic(self.DbiopsL_plus_1[0],  0, 0,  0, 0, self.z[0:2])
        self.Acomp[0,0] =  self.Acomp[0,0] + 1 * (self.z[1]-self.z[0])/2/self.delt
        self.Bcomp[0,0] =  self.Bcomp[0,0] + 1 * (self.z[1]-self.z[0])/2/self.delt

        p = self.ptot[-1]
        [self.Acomp[p,p-1], self.Acomp[p,p]], self.bcomp[p] = flux_bottom_boundary_kinetic(0, self.Dbiops_plus_1[p], 0, 0, self.z[p-1:p+1])
        self.Acomp[p,p] =  self.Acomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt
        self.Bcomp[p,p] =  self.Bcomp[p,p] + 1 * (self.z[p]-self.z[p-1])/2/self.delt

    def make_governing_equations(self):
        """Makes finite difference equations for the governing PDEs."""

        if self.timeoption == 'Implicit':
            for j in range(len(self.ptotbio) - 1):
                num = self.layertotbio[j].nchemicals
                for n in range(self.nchemicals):
                    chemical = self.layertotbio[j].chemicals[n]
                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 2)):
                        i = ii + 1
                        p = self.ptotbio[j] + i
                        if i == 1: cp = self.cptotbio[j]
                        else:      cp = self.cptotLbio[j]
                        if self.bio == 1:
                            KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4) * self.KsL_plus_1[cp+(i-1)*num+n],
                                              ( self.Dbiops_plus_1[p])                                                        * self.Ks_plus_1[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4) * self.Ks_plus_1[self.cptotLbio[j]+(i+1)*num+n]]
                            if self.biotype == 'Uniform':
                                Ds_plus_1      = [self.Ds_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n]]
                            else:
                                Ds_plus_1      = [ self.DsL_plus_1[cp+(i-1)*num+n]/4+self.Ds_plus_1[self.cptotLbio[j]+i*num+n]-self.Ds_plus_1[self.cptotLbio[j]+(i+1)*num+n]/4,
                                                                                     self.Ds_plus_1[self.cptotLbio[j]+i*num+n],
                                                  -self.DsL_plus_1[cp+(i-1)*num+n]/4+self.Ds_plus_1[self.cptotLbio[j]+i*num+n]+self.Ds_plus_1[self.cptotLbio[j]+(i+1)*num+n]/4]
                        else:
                            KDbiops_plus_1 = [0,0,0]
                            Ds_plus_1      = [self.DsL_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n]]

                        [self.A[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n]] , \
                        [self.B[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n]] = \
                        get_4pt_adr_fde_imp(self.RsL[self.cptotLbio[j]+i*num+n], self.RsL_plus_1[self.cptotLbio[j]+i*num+n], self.UsL_plus_1[self.cptotLbio[j]+i*num+n], Ds_plus_1, KDbiops_plus_1, self.delt, self.z[p-1:p+3])

                        for np in range(num):
                            self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] - self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]
                        self.a[self.cptotLbio[j]+i*num+n] = -self.rates_plus_1[self.cptotLbio[j]+i*num+n]

                        if self.RsL_plus_1[self.cptotLbio[j]+i*num+n] == 0: self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = 1

                    p        = self.ptotbio[j+1]
                    if self.bio == 1:
                        KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] - self.Dbiops_plus_1[p]/4) *self.KsL_plus_1[self.cptotbio[j+1]-2*num+n],
                                          ( self.Dbiops_plus_1[p-1])                                                      *self.Ks_plus_1[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] + self.Dbiops_plus_1[p]/4) *self.Ks_plus_1[self.cptotbio[j+1]+n]]

                        if self.biotype == 'Uniform':
                            Ds_plus_1      = [self.DsL_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n]]
                        else:
                            Ds_plus_1      = [ self.DsL_plus_1[self.cptotbio[j+1]-2*num+n]/4+self.Ds_plus_1[self.cptotbio[j+1]-num+n]-self.Ds_plus_1[self.cptotbio[j+1]+n]/4,
                                               self.Ds_plus_1[self.cptotbio[j+1]-num+n],
                                              -self.DsL_plus_1[self.cptotbio[j+1]-2*num+n]/4+self.Ds_plus_1[self.cptotbio[j+1]-num+n]+self.Ds_plus_1[self.cptotbio[j+1]+n]/4]
                    else:
                        KDbiops_plus_1 = [0,0,0]
                        Ds_plus_1      = [self.DsL_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n]]

                    [self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n]], \
                    [self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n]] = \
                    get_3pt_adr_fde_imp(self.Rs[self.cptotbio[j+1]-num+n], self.Rs_plus_1[self.cptotbio[j+1]-num+n],  self.Us_plus_1[self.cptotbio[j+1]-num+n], Ds_plus_1, KDbiops_plus_1, self.delt, self.z[self.ptotbio[j+1]-2:self.ptotbio[j+1]+1])
                    for np in range(num):
                        self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] - self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]
                    self.a[self.cptotbio[j+1]-num+n] = -self.rates_plus_1[self.cptotbio[j+1]-num+n]
                    if self.Rs_plus_1[self.cptotbio[j+1]-num+n] == 0: self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n] = 1

                for nn in range(self.layertotbio[j].nsolidchemicals):
                    solidchemical = self.layertotbio[j].solidchemicals[nn]
                    n = nn + self.nchemicals
                    n_lower = self.layertotbio[j].lowersolidchemicalnums[nn]
                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 2)):
                        i = ii + 1
                        p = self.ptotbio[j] + i
                        if i == 1: cp = self.cptotbio[j]
                        else:      cp = self.cptotLbio[j]
                        if i == self.ptotbio[j+1]-1: n_l = n_lower
                        else:                        n_l = n

                        if self.bio == 1:
                            KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4) *self.KsL_plus_1[cp+(i-1)*num+n],
                                              ( self.Dbiops_plus_1[p])                                                         *self.Ks_plus_1[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4) *self.Ks_plus_1[self.cptotLbio[j]+(i+1)*num+n]]

                        else:
                            KDbiops_plus_1 = [0,0,0]

                        if self.Rs_plus_1[self.cptotLbio[j]+i*num+n]/self.components[solidchemical.component_index].rho < 0.0000000001:
                            self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = 1
                        else:
                            [self.A[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n_l]] , \
                            [self.B[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n_l]] = \
                            get_4pt_adr_fde_imp(self.RsL[self.cptotLbio[j]+i*num+n], self.RsL_plus_1[self.cptotLbio[j]+i*num+n], self.UsL_plus_1[self.cptotLbio[j]+i*num+n],  [0,0,0], KDbiops_plus_1, self.delt, self.z[p-1:p+3])
                            for np in range(num):
                                self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] - self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]
                            self.a[self.cptotLbio[j]+i*num+n] = -self.rates_plus_1[self.cptotLbio[j]+i*num+n]


                    p        = self.ptotbio[j+1]
                    if self.bio == 1:
                        KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] - self.Dbiops_plus_1[p]/4) *self.KsL_plus_1[self.cptotbio[j+1]-2*num+n],
                                          ( self.Dbiops_plus_1[p-1])                                                      *self.Ks_plus_1[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] + self.Dbiops_plus_1[p]/4) *self.Ks_plus_1[self.cptotbio[j+1]+n_lower]]

                    else:
                        KDbiops_plus_1 = [0,0,0]

                    if self.Rs_plus_1[self.cptotbio[j+1]-num+n]/self.components[solidchemical.component_index].rho < 0.0000000001:
                        self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n] = 1
                    else:
                        [self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n_lower]], \
                        [self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n_lower]] = \
                        get_3pt_adr_fde_imp(self.Rs[self.cptotbio[j+1]-num+n], self.Rs_plus_1[self.cptotbio[j+1]-num+n], self.Us_plus_1[self.cptotbio[j+1]-num+n],[0,0,0], KDbiops_plus_1, self.delt, self.z[self.ptotbio[j+1]-2:self.ptotbio[j+1]+1])
                        for np in range(num):
                            self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] - self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]
                        self.a[self.cptotbio[j+1]-num+n] = -self.rates_plus_1[self.cptotbio[j+1]-num+n]


        if self.timeoption == 'Crank-Nicolson':

            for j in range(len(self.ptotbio) - 1):
                num = self.layertotbio[j].nchemicals
                for n in range(self.nchemicals):
                    chemical = self.layertotbio[j].chemicals[n]
                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 2)):
                        i = ii + 1
                        p = self.ptotbio[j] + i
                        if i == 1: cp = self.cptotbio[j]
                        else:      cp = self.cptotLbio[j]
                        if self.bio == 1:
                            KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4) *self.KsL_plus_1[cp+(i-1)*num+n],
                                              ( self.Dbiops_plus_1[p])                                                        *self.Ks_plus_1[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4) *self.Ks_plus_1[self.cptotLbio[j]+(i+1)*num+n]]
                            KDbiops        = [( self.DbiopsL[p-1]/4+self.DbiopsL[p] - self.Dbiops[p+1]/4) *self.KsL[cp+(i-1)*num+n],
                                              ( self.DbiopsL[p])                                           *self.Ks[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL[p-1]/4+self.DbiopsL[p] + self.Dbiops[p+1]/4) *self.Ks[self.cptotLbio[j]+(i+1)*num+n]]
                            if self.biotype == 'Uniform':
                                Ds_plus_1      = [self.DsL_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n]]
                                Ds             = [self.DsL[self.cptotLbio[j]+i*num+n], self.Ds[self.cptotLbio[j]+i*num+n], self.Ds[self.cptotLbio[j]+i*num+n]]
                            else:
                                Ds_plus_1      = [ self.DsL_plus_1[cp+(i-1)*num+n]/4+self.Ds_plus_1[self.cptotLbio[j]+i*num+n]-self.Ds_plus_1[self.cptotLbio[j]+(i+1)*num+n]/4,
                                                                                     self.Ds_plus_1[self.cptotLbio[j]+i*num+n],
                                                  -self.DsL_plus_1[cp+(i-1)*num+n]/4+self.Ds_plus_1[self.cptotLbio[j]+i*num+n]+self.Ds_plus_1[self.cptotLbio[j]+(i+1)*num+n]/4]
                                Ds             = [ self.DsL[cp+(i-1)*num+n]/4+self.Ds[self.cptotLbio[j]+i*num+n]-self.Ds[self.cptotLbio[j]+(i+1)*num+n]/4,
                                                                              self.Ds[self.cptotLbio[j]+i*num+n],
                                                  -self.DsL[cp+(i-1)*num+n]/4+self.Ds[self.cptotLbio[j]+i*num+n]+self.Ds[self.cptotLbio[j]+(i+1)*num+n]/4]

                        else:
                            KDbiops_plus_1 = [0,0,0]
                            KDbiops        = [0,0,0]
                            Ds_plus_1      = [self.DsL_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n], self.Ds_plus_1[self.cptotLbio[j]+i*num+n]]
                            Ds             = [self.DsL[self.cptotLbio[j]+i*num+n], self.Ds[self.cptotLbio[j]+i*num+n], self.Ds[self.cptotLbio[j]+i*num+n]]

                        [self.A[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotbio[j]+(i+2)*num+n]] , \
                        [self.B[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotbio[j]+(i+2)*num+n]] = \
                        get_4pt_adr_fde_CN(self.RsL[self.cptotLbio[j]+i*num+n], self.RsL_plus_1[self.cptotLbio[j]+i*num+n], self.UsL[self.cptotLbio[j]+i*num+n], self.UsL_plus_1[self.cptotLbio[j]+i*num+n], Ds, Ds_plus_1, KDbiops, KDbiops_plus_1, self.delt, self.z[p-1:p+3])

                        for np in range(num):
                            self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] - self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]/2
                            self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] + self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]/2

                        self.a[self.cptotLbio[j]+i*num+n] = -self.rates_plus_1[self.cptotLbio[j]+i*num+n]/2
                        self.b[self.cptotLbio[j]+i*num+n] =  self.rates_plus_1[self.cptotLbio[j]+i*num+n]/2

                        if self.RsL_plus_1[self.cptotLbio[j]+i*num+n] == 0: self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = 1

                    p        = self.ptotbio[j+1]
                    if self.bio == 1:
                        KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] - self.Dbiops_plus_1[p]/4) *self.KsL_plus_1[self.cptotbio[j+1]-2*num+n],
                                          ( self.Dbiops_plus_1[p-1])                                                      *self.Ks_plus_1[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] + self.Dbiops_plus_1[p]/4) *self.Ks_plus_1[self.cptotbio[j+1]+n]]
                        KDbiops         = [( self.DbiopsL[p-2]/4+self.Dbiops[p-1] - self.Dbiops[p]/4) *self.KsL[self.cptotbio[j+1]-2*num+n],
                                          ( self.Dbiops[p-1])                                        *self.Ks[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL[p-2]/4+self.Dbiops[p-1] + self.Dbiops[p]/4) *self.Ks[self.cptotbio[j+1]+n]]

                        if self.biotype == 'Uniform':
                            Ds_plus_1      = [self.DsL_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n]]
                            Ds             = [self.DsL[self.cptotbio[j+1]-num+n], self.Ds[self.cptotbio[j+1]-num+n], self.Ds[self.cptotbio[j+1]-num+n]]
                        else:
                            Ds_plus_1      = [ self.DsL_plus_1[self.cptotbio[j+1]-2*num+n]/4+self.Ds_plus_1[self.cptotbio[j+1]-num+n]-self.Ds_plus_1[self.cptotbio[j+1]+n]/4,
                                               self.Ds_plus_1[self.cptotbio[j+1]-num+n],
                                              -self.DsL_plus_1[self.cptotbio[j+1]-2*num+n]/4+self.Ds_plus_1[self.cptotbio[j+1]-num+n]+self.Ds_plus_1[self.cptotbio[j+1]+n]/4]
                            Ds             = [ self.DsL[self.cptotbio[j+1]-2*num+n]/4+self.Ds[self.cptotbio[j+1]-num+n]-self.Ds[self.cptotbio[j+1]+n]/4,
                                               self.Ds[self.cptotbio[j+1]-num+n],
                                              -self.DsL[self.cptotbio[j+1]-2*num+n]/4+self.Ds[self.cptotbio[j+1]-num+n]+self.Ds[self.cptotbio[j+1]+n]/4]
                    else:
                        KDbiops_plus_1 = [0,0,0]
                        KDbiops        = [0,0,0]
                        Ds_plus_1      = [self.DsL_plus_1[self.cptotbio[j+1]-num+n], self.Ds_plus_1[self.cptotbio[j+1]-num+n],  self.Ds_plus_1[self.cptotbio[j+1]-num+n]]
                        Ds             = [self.DsL[self.cptotbio[j+1]-num+n],        self.Ds[self.cptotbio[j+1]-num+n],         self.Ds[self.cptotbio[j+1]-num+n]]

                    [self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n]], \
                    [self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n]] = \
                    get_3pt_adr_fde_CN(self.Rs[self.cptotbio[j+1]-num+n], self.Rs_plus_1[self.cptotbio[j+1]-num+n],  self.Us[self.cptotbio[j+1]-num+n], self.Us_plus_1[self.cptotbio[j+1]-num+n], Ds, Ds_plus_1, KDbiops, KDbiops_plus_1, self.delt, self.z[self.ptotbio[j+1]-2:self.ptotbio[j+1]+1])
                    for np in range(num):
                        self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] - self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]/2
                        self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] + self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]/2
                    self.a[self.cptotbio[j+1]-num+n] = -self.rates_plus_1[self.cptotbio[j+1]-num+n]/2
                    self.b[self.cptotbio[j+1]-num+n] = -self.rates_plus_1[self.cptotbio[j+1]-num+n]/2
                    if self.Rs_plus_1[self.cptotbio[j+1]-num+n] == 0: self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n] = 1

                for nn in range(self.layertotbio[j].nsolidchemicals):
                    solidchemical = self.layertotbio[j].solidchemicals[nn]
                    n = nn + self.nchemicals
                    n_lower = self.layertotbio[j].lowersolidchemicalnums[nn]
                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 2)):
                        i = ii + 1
                        p = self.ptotbio[j] + i
                        if i == 1: cp = self.cptotbio[j]
                        else:      cp = self.cptotLbio[j]
                        if i == self.ptotbio[j+1]-1: n_l = n_lower
                        else:                n_l = n

                        if self.bio == 1:
                            KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] - self.Dbiops_plus_1[p+1]/4) *self.KsL_plus_1[self.cptotLbio[j]+(i-1)*num+n],
                                              ( self.Dbiops_plus_1[p])                                                         *self.Ks_plus_1[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL_plus_1[p-1]/4+self.Dbiops_plus_1[p] + self.Dbiops_plus_1[p+1]/4) *self.Ks_plus_1[self.cptotLbio[j]+(i+1)*num+n]]
                            KDbiops        = [( self.DbiopsL[p-1]/4+self.DbiopsL[p] - self.DbiopsL[p+1]/4) *self.KsL[self.cptotLbio[j]+(i-1)*num+n],
                                              ( self.DbiopsL[p])                                           *self.Ks[self.cptotLbio[j]+i*num+n],
                                              (-self.DbiopsL[p-1]/4+self.DbiopsL[p] + self.DbiopsL[p+1]/4) *self.Ks[self.cptotLbio[j]+(i+1)*num+n]]
                        else:
                            KDbiops_plus_1 = [0,0,0]
                            KDbiops        = [0,0,0]
                        if self.Rs_plus_1[self.cptotLbio[j]+i*num+n]/self.components[solidchemical.component_index].rho < 0.0000000001:
                            self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = 1
                        else:
                            [self.A[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n_l]] , \
                            [self.B[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n], self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+2)*num+n_l]] = \
                            get_4pt_adr_fde_CN(self.RsL[self.cptotLbio[j]+i*num+n], self.RsL_plus_1[self.cptotLbio[j]+i*num+n], self.UsL[self.cptotLbio[j]+i*num+n], self.UsL_plus_1[self.cptotLbio[j]+i*num+n], [0,0,0], [0,0,0], KDbiops, KDbiops_plus_1, self.delt, self.z[p-1:p+3])
                            for np in range(num):
                                self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.A[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] - self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]/2
                                self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] = self.B[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np] + self.elams[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+np]/2
                            self.a[self.cptotLbio[j]+i*num+n] = -self.rates_plus_1[self.cptotLbio[j]+i*num+n]/2
                            self.b[self.cptotLbio[j]+i*num+n] =  self.rates_plus_1[self.cptotLbio[j]+i*num+n]/2


                    p        = self.ptotbio[j+1]
                    if self.bio == 1:
                        KDbiops_plus_1 = [( self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] - self.Dbiops_plus_1[p]/4) *self.KsL_plus_1[self.cptotbio[j+1]-2*num+n],
                                          ( self.Dbiops_plus_1[p-1])                                                       *self.Ks_plus_1[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL_plus_1[p-2]/4+self.Dbiops_plus_1[p-1] + self.Dbiops_plus_1[p]/4) *self.Ks_plus_1[self.cptotbio[j+1]+n_lower]]
                        KDbiops        = [( self.DbiopsL[p-2]/4+self.DbiopsL[p-1] - self.Dbiops[p]/4) *self.KsL[self.cptotbio[j+1]-2*num+n],
                                          ( self.DbiopsL[p-1])                                        *self.Ks[self.cptotbio[j+1]-num+n],
                                          (-self.DbiopsL[p-2]/4+self.DbiopsL[p-1] + self.Dbiops[p]/4) *self.Ks[self.cptotbio[j+1]+n_lower]]
                    else:
                        KDbiops_plus_1 = [0,0,0]
                        KDbiops        = [0,0,0]

                    if self.Rs_plus_1[self.cptotbio[j+1]-num+n]/self.components[solidchemical.component_index].rho < 0.0000000001:
                        self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n] = 1
                    else:
                        [self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n_lower]], \
                        [self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-2*num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+n], self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]+n_lower]] = \
                        get_3pt_adr_fde_CN(self.Rs[self.cptotbio[j+1]-num+n], self.Rs_plus_1[self.cptotbio[j+1]-num+n], self.Us[self.cptotbio[j+1]-num+n], self.Us_plus_1[self.cptotbio[j+1]-num+n],[0,0,0], [0,0,0], KDbiops, KDbiops_plus_1, self.delt, self.z[self.ptotbio[j+1]-2:self.ptotbio[j+1]+1])
                        for np in range(num):
                            self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.A[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] - self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]/2
                            self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] = self.B[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np] + self.elams[self.cptotbio[j+1]-num+n,self.cptotbio[j+1]-num+np]/2
                        self.a[self.cptotbio[j+1]-num+n] = -self.rates_plus_1[self.cptotbio[j+1]-num+n]/2
                        self.b[self.cptotbio[j+1]-num+n] =  self.rates_plus_1[self.cptotbio[j+1]-num+n]/2


    def make_Newton_Raphson_equations(self, Cn, CnL, Fis, FisL):

        self.NR     = self.A.copy()

        if self.reac == 1:
            if self.topBCtype <> 'Fixed Concentration':
                for n in range(self.layertotbio[0].nchemicals):
                    for np in range(self.layertotbio[0].nchemicals):
                        self.NR[n,np]    = self.NR[n,np] -self.rates_diff[n,np]

            for j in range(len(self.layertotbio) - 1):
                i = self.cptotbio[j+1]
                p = self.ptotbio[j+1]
                for n in range (self.layertotbio[j+1].nchemicals):
                    for np in range(self.layertotbio[j+1].nchemicals+self.layertotbio[j].nlowersolidchemicals):
                        self.NR[i+n,i+np]    = self.NR[i+n,i+np]  - self.rates_diff[i+n,i+np]

            i   = self.cptotbio[-1]
            if self.botBCtype <> 'Fixed Concentration':
                for n in range(self.layertotbio[-1].nchemicals):
                    for np in range(self.layertotbio[-1].nchemicals):
                        self.NR[i+n,i+np]   = self.NR[i+n,i+np] -self.rates_diff[i+n, i+np]
        j = 0
        num = self.layertotbio[0].nchemicals
        if self.sorp == 1 and self.topBCtype <> 'Fixed Concentration':
            for n in range(self.nchemicals):
                chemical = self.chemicals[n]
                for component in self.layertotbio[0].components:
                    if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                        if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                            self.NR[n,n] = (self.NR[n,n]-Fis[component.name][0] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[1]-self.z[0])/2
                                                        +Fis[component.name][0] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[1]-self.z[0])/2)

                            if self.bio == 1 and self.Dbiop > 0:
                                self.NR[n,n]    = (self.NR[n,n]     - self.Dbiops_plus_1[0]*FisL[component.name][0]*component.rho*self.sorptions[component.name][chemical.name].get_K( component, Cn[n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                                                                    + self.Dbiops_plus_1[0]*FisL[component.name][0]*component.rho*self.sorptions[component.name][chemical.name].get_NR(Cn[n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW))
                                self.NR[n,n+num]= (self.NR[n,n+num] + self.Dbiops_plus_1[0]*FisL[component.name][1]*component.rho*self.sorptions[component.name][chemical.name].get_K( component, Cn[n+num]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                                                                    - self.Dbiops_plus_1[0]*FisL[component.name][1]*component.rho*self.sorptions[component.name][chemical.name].get_NR(Cn[n+num]*chemical.MW, self.Cmax[chemical.name]*chemical.MW))


        j = -1
        if self.sorp == 1 and self.botBCtype <> 'Fixed Concentration':
            for n in range(self.nchemicals):
                chemical = self.chemicals[n]
                for component in self.layertotbio[j].components:
                    if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                        if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                            self.NR[self.cptotbio[j]+n,self.cptotbio[j]+n] = (self.NR[self.cptotbio[j]+n,self.cptotbio[j]+n]-Fis[component.name][self.ptotbio[j]] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotbio[j]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* self.delz[-1]/2
                                                                                            +Fis[component.name][self.ptotbio[j]] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotbio[j]+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* self.delz[-1]/2)

                            if self.bio == 1 and self.Dbiops[-1] > 0:
                                self.NR[self.cptotbio[j]+n,self.cptotbio[j]+n]    = (self.NR[self.cptotbio[j]+n,self.cptotbio[j]+n]     + self.Dbiops_plus_1[-1]*FisL[component.name][-1]*component.rho*self.sorptions[component.name][chemical.name].get_K( component, Cn[self.cptotbio[j]+n]    *chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                                                                                                        - self.Dbiops_plus_1[-1]*FisL[component.name][-1]*component.rho*self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotbio[j]+n]    *chemical.MW, self.Cmax[chemical.name]*chemical.MW))
                                self.NR[self.cptotbio[j]+n,self.cptotbio[j]-num+n]= (self.NR[self.cptotbio[j]+n,self.cptotbio[j]-num+n] - self.Dbiops_plus_1[-1]*FisL[component.name][-2]*component.rho*self.sorptions[component.name][chemical.name].get_K( component, Cn[self.cptotbio[j]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)
                                                                                                        + self.Dbiops_plus_1[-1]*FisL[component.name][-2]*component.rho*self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotbio[j]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW))


        for j in range(1, len(self.ptotbio)-1):
            layer = self.layertotbio[j]
            num   = self.layertotbio[j].nchemicals
            numa =self.layertotbio[j-1].nchemicals
            numb =self.layertotbio[j].nchemicals
            i = self.cptotbio[j]
            iL = self.cptotLbio[j]
            p = self.ptotbio[j]
            for n in range(self.nchemicals):
                chemical = self.chemicals[n]
                for component in self.layertotbio[j].components:
                    if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                        if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                            self.NR[i+n,i+n] = (self.NR[i+n,i+n]-FisL[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[p+1]-self.z[p])/2
                                                                +FisL[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_NR(CnL[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[p+1]-self.z[p])/2)

                for component in self.layertotbio[j-1].components:
                    if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                        if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                            self.NR[i+n,i+n] = (self.NR[i+n,i+n]-Fis[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[p]-self.z[p-1])/2
                                                                +Fis[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt* (self.z[p]-self.z[p-1])/2)

            if self.bio == 1 and self.Dbiop > 0 and self.sorp == 1:
                for n in range(self.nchemicals):
                    chemical = self.chemicals[n]
                    for component in self.layertotbio[j-1].components:
                        if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                            if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                self.NR[i+n,i-numa+n] = (self.NR[i+n,i-numa+n]
                                                            +self.Dbiops_plus_1[p] * Fis[component.name][p-1] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[i-numa+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p]-self.z[p-1])
                                                            -self.Dbiops_plus_1[p] * Fis[component.name][p-1] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[i-numa+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p]-self.z[p-1]))
                                self.NR[i+n,i+n] = (self.NR[i+n,i+n]
                                                            -self.Dbiops_plus_1[p] * Fis[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p]-self.z[p-1])
                                                            +self.Dbiops_plus_1[p] * Fis[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p]-self.z[p-1]))
                    for component in self.layertotbio[j].components:
                        if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                            if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                self.NR[i+n,i+n] = (self.NR[i+n,i+n]
                                                            -self.DbiopsL_plus_1[p] * FisL[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p+1]-self.z[p])
                                                            +self.DbiopsL_plus_1[p] * FisL[component.name][p] * component.rho * self.sorptions[component.name][chemical.name].get_NR(CnL[i+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p+1]-self.z[p]))

                                self.NR[i+n,iL+numb+n] = (self.NR[i+n,iL+numb+n]
                                                            +self.DbiopsL_plus_1[p] * FisL[component.name][p+1] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[iL+numb+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p+1]-self.z[p])
                                                            -self.DbiopsL_plus_1[p] * FisL[component.name][p+1] * component.rho * self.sorptions[component.name][chemical.name].get_NR(CnL[iL+numb+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[p+1]-self.z[p]))


        if self.timeoption == 'Implicit':
            # Correct the non-linear sorption terms in Newton Raphson Matrix
            for j in range(len(self.ptotbio) - 1):
                layer = self.layertotbio[j]
                num = self.layertotbio[j].nchemicals
                if self.sorp == 1:
                    for n in range(self.nchemicals):
                        chemical = self.chemicals[n]
                        for component in layer.components:
                            if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                                if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                    for ii in range((self.ptotbio[j+1]-self.ptotbio[j]-2)):
                                        i = ii + 1
                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n]
                                                                    -Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt
                                                                    +Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt)

                                    self.NR[self.cptotLbio[j+1]-num+n,self.cptotLbio[j+1]-num+n] = (self.NR[self.cptotLbio[j+1]-num+n,self.cptotLbio[j+1]-num+n]
                                                                -Fis[component.name][self.ptotbio[j+1]-1] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j+1]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt
                                                                +Fis[component.name][self.ptotbio[j+1]-1] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j+1]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt)

                # Correct the non-linear reaction terms in Newton Raphson Matrix
                if self.reac == 1:
                    for n in range(num):
                        for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                            i = ii + 1
                            for nn in range(num):
                                self.NR[self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn] = self.NR[self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn] - self.rates_diff[self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn]

            if self.bio == 1 and self.Dbiop > 0 and self.sorp == 1:
                for j in range(self.layerbio+1):
                    layer = self.layertotbio[j]
                    num = self.layertotbio[j].nchemicals
                    for n in range(self.nchemicals):
                        chemical = self.chemicals[n]
                        for component in layer.components:
                            if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                                if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                                        i = ii + 1
                                        if i == 1:  cp = self.cptotbio[j]
                                        else:       cp = self.cptotLbio[j]
                                        self.NR[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n]
                                                                                            +(self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]-self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * FisL[component.name][self.ptotbio[j]+(i-1)] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[cp+(i-1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2
                                                                                            -(self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]-self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * FisL[component.name][self.ptotbio[j]+(i-1)] * component.rho * self.sorptions[component.name][chemical.name].get_NR(CnL[cp+(i-1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2)

                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n]
                                                                                        -(self.DbiopsL_plus_1[self.ptotbio[j]+i]) * Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2*2
                                                                                        +(self.DbiopsL_plus_1[self.ptotbio[j]+i]) * Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2*2)

                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n]
                                                                                        +(-self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]+self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * Fis[component.name][self.ptotbio[j]+(i+1)] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+(i+1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2
                                                                                        -(-self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]+self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * Fis[component.name][self.ptotbio[j]+(i+1)] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+(i+1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2)


        if self.timeoption == 'Crank-Nicolson':
            # Correct the non-linear sorption terms in Newton Raphson Matrix
            for j in range(len(self.ptotbio) - 1):
                layer = self.layertotbio[j]
                num = self.layertotbio[j].nchemicals
                if self.sorp == 1:
                    for n in range(self.nchemicals):
                        chemical = self.chemicals[n]
                        for component in layer.components:
                            if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                                if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                    for ii in range((self.ptotbio[j+1]-self.ptotbio[j]-2)):
                                        i = ii + 1
                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n]
                                                                    -Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt
                                                                    +Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt)

                                    self.NR[self.cptotLbio[j+1]-num+n,self.cptotLbio[j+1]-num+n] = (self.NR[self.cptotLbio[j+1]-num+n,self.cptotLbio[j+1]-num+n]
                                                                -Fis[component.name][self.ptotbio[j+1]-1] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j+1]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt
                                                                +Fis[component.name][self.ptotbio[j+1]-1] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j+1]-num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/self.delt)

                # Correct the non-linear reaction terms in Newton Raphson Matrix
                    if self.reac == 1:
                        for n in range(num):
                            for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                                i = ii + 1
                                for nn in range(num):
                                    self.NR[self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn] = self.NR[self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn] - self.rates_diff[self.cptotLbio[j]+i*num+n, self.cptotbio[j]+i*num+nn]/2

            if self.bio == 1 and self.Dbiop > 0 and self.sorp == 1:
                for j in range(self.layerbio+1):
                    layer = self.layertotbio[j]
                    num = self.layertotbio[j].nchemicals
                    for n in range(self.nchemicals):
                        chemical = self.chemicals[n]
                        for component in layer.components:
                            if self.sorptions[component.name][chemical.name].kinetic == 'Equilibrium':
                                if self.sorptions[component.name][chemical.name].isotherm == 'Freundlich' or self.sorptions[component.name][chemical.name].isotherm == 'Langmuir':
                                    for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                                        i = ii + 1
                                        if i == 1:  cp = self.cptotbio[j]
                                        else:       cp = self.cptotLbio[j]
                                        self.NR[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,cp+(i-1)*num+n]
                                                                                            +(self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]-self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * FisL[component.name][self.ptotbio[j]+(i-1)] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, CnL[cp+(i-1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2/2
                                                                                            -(self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]-self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * FisL[component.name][self.ptotbio[j]+(i-1)] * component.rho * self.sorptions[component.name][chemical.name].get_NR(CnL[cp+(i-1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2/2)

                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+i*num+n]
                                                                                        -(self.DbiopsL_plus_1[self.ptotbio[j]+i]) * Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2*2/2
                                                                                        +(self.DbiopsL_plus_1[self.ptotbio[j]+i]) * Fis[component.name][self.ptotbio[j]+i] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+i*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2*2/2)

                                        self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n] = (self.NR[self.cptotLbio[j]+i*num+n,self.cptotLbio[j]+(i+1)*num+n]
                                                                                        +(-self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]+self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * Fis[component.name][self.ptotbio[j]+(i+1)] * component.rho * self.sorptions[component.name][chemical.name].get_K(component, Cn[self.cptotLbio[j]+(i+1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2/2
                                                                                        -(-self.DbiopsL_plus_1[self.ptotbio[j]+(i-1)]/4+self.Dbiops_plus_1[self.ptotbio[j]+i]+self.Dbiops_plus_1[self.ptotbio[j]+(i+1)]/4) * Fis[component.name][self.ptotbio[j]+(i+1)] * component.rho * self.sorptions[component.name][chemical.name].get_NR(Cn[self.cptotLbio[j]+(i+1)*num+n]*chemical.MW, self.Cmax[chemical.name]*chemical.MW)/(self.z[self.ptotbio[j]+1]-self.z[self.ptotbio[j]])**2/2)



    def make_matrix_parameter_vectors(self, Cn, CnL, Fis, FisL):

        self.make_grid_es_rhos(Fis, FisL)
        self.make_grid_Rs(Cn, CnL, Fis, FisL)

    def make_transport_parameter_vectors(self):

        self.make_grid_Us()
        self.make_grid_Ds()

    def make_reaction_parameter_vectors(self, Cn, CnL, Fis, FisL):

        self.make_grid_elams(Fis, FisL)
        self.make_grid_rates_plus_1(Cn, CnL, Fis, FisL)

        self.rates = self.rates_plus_1.copy()

    def make_matrices(self, Flag = None, Cn = None, CnL= None):
        """Generates the discretized variables necessary to perform a finite
        difference analysis of the underlying governing differential equations
        and auxiliary conditions for a linear system.  The "A," "B," and "b"
        matrices are used to determine the concentrations at the succeeding
        time step from:

        A * Cn+1 = B * Cn + b

        Cn -- the concentrations in the grid at time n
        """

        self.A = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, self.cptot[-1] + self.layertot[-1].nchemicals]))
        self.B = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, self.cptot[-1] + self.layertot[-1].nchemicals]))
        self.a = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, 1]))
        self.b = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, 1]))

        self.make_boundary_equations(Flag, Cn, CnL)
        self.make_governing_equations()

    def make_components_matrices(self):

        self.Acomp = matrix(zeros([self.ptot[-1] + 1, self.ptot[-1] + 1]))
        self.Bcomp = matrix(zeros([self.ptot[-1] + 1, self.ptot[-1] + 1]))
        self.acomp = matrix(zeros([self.ptot[-1] + 1, 1]))
        self.bcomp = matrix(zeros([self.ptot[-1] + 1, 1]))

        self.make_components_equations()

    def update_time_dependents(self):

        self.U      = self.U_plus_1

        if self.bio == 1:
            self.Dbiops = [Dbiop for Dbiop in self.Dbiops_plus_1]
            self.DbiopsL= [Dbiop for Dbiop in self.DbiopsL_plus_1]

        self.es     = [e for e in self.es_plus_1]
        self.rhos   = [rho for rho in self.rhos_plus_1]
        self.Ks     = [R for R in self.Ks_plus_1]
        self.Rs     = [R for R in self.Rs_plus_1]
        self.Ss     = [S for S in self.Ss_plus_1]
        self.Us     = [U for U in self.Us_plus_1]
        self.Ds     = [D for D in self.Ds_plus_1]

        self.esL     = [e for e in self.esL_plus_1]
        self.rhosL   = [rho for rho in self.rhosL_plus_1]
        self.KsL     = [R for R in self.KsL_plus_1]
        self.SsL     = [S for S in self.SsL_plus_1]
        self.RsL     = [R for R in self.RsL_plus_1]
        self.UsL     = [U for U in self.UsL_plus_1]
        self.DsL     = [D for D in self.DsL_plus_1]

        self.elams  = self.elams_plus_1.copy()
        self.rates  = self.rates_plus_1.copy()

        for chemical in self.chemicals:
            self.BCs[chemical.name].Cw = self.BCs[chemical.name].Cw_plus_1

    def update_bioturbation(self):

        i                = 0
        self.pbio        = 0  # make sure self.pbio is not defined for the Gaussian model case
        flag             = 0
        while flag == 0:

            if i <= len(self.ptot) - 1:

                if (self.z[self.ptot[i]]-self.z[0]) < self.hbio:
                    i = i + 1
                else:
                    self.pbio   = self.ptot[i-1] + int(round((self.hbio - (self.z[self.ptot[i-1]]-self.z[0])) / (self.z[self.ptot[i-1]+1] -self.z[self.ptot[i-1]] ), 8))
                    self.cpbio  = self.cptotL[i-1] + (self.pbio - self.ptot[i-1]) * self.layertot[i-1].nchemicals
                    self.cpbioL = self.cptotL[i-1] + (self.pbio - self.ptot[i-1]) * self.layertot[i-1].nchemicals
                    self.layerbio = i - 1
                    flag = 1
            else:
                self.pbio   = self.ptot[i-1]
                self.cpbio  = self.cptotL[i-1]
                self.cpbioL = self.cptotL[i-1]
                self.layerbio = i - 2
                flag = 2

        if flag == 1:
            if self.ptot[self.layerbio + 1] - self.pbio   <= 1:
                self.pbio   = self.ptot[self.layerbio + 1]
                self.cpbio  = self.cptotL[self.layerbio + 1]
                self.cpbioL = self.cptotL[self.layerbio + 1]

            elif self.pbio - self.ptot[self.layerbio]     <= 1:
                if self.layerbio == 0:
                    self.pbio   = self.ptot[self.layerbio] + 2
                    self.cpbio  = self.cptotL[self.layerbio] + 2
                    self.cpbioL = self.cptotL[self.layerbio] + 2

                else:
                    self.pbio = self.ptot[self.layerbio]
                    self.cpbio = self.cptotL[self.layerbio]
                    self.cpbioL = self.cptotL[self.layerbio]
                    self.layerbio = self.layerbio - 1

        self.ptotbio     = [ptot               for ptot    in self.ptot]
        self.cptotbio    = [cptot              for cptot   in self.cptot]
        self.cptotLbio   = [cptotL             for cptotL  in self.cptotL]
        self.layertotbio = [layer.bio_copy()   for layer   in self.layertot]

        if self.bio == 1 and self.ptotbio.count(self.pbio) == 0:
            self.ptotbio.insert(self.layerbio + 1,   self.pbio)
            self.cptotbio.insert(self.layerbio + 1,  self.cpbio)
            self.cptotLbio.insert(self.layerbio + 1, self.cpbioL)
            self.layertotbio.insert(self.layerbio + 1, self.layertot[self.layerbio].bio_copy())

            self.layertotbio[self.layerbio].nlowersolidchemicals      = 0
            self.layertotbio[self.layerbio].lowersolidchemicals       = []
            self.layertotbio[self.layerbio].lowersolidchemicalnums    = [n + self.nchemicals for n in range(self.layertotbio[self.layerbio].nsolidchemicals)]


    def update_deposition(self, Cn, CnL, Fis_load, FisL_load, t):
        """Updates the grid for deposition."""

        updatecheck = 0

        Fis  = {}
        FisL = {}
        for component in self.components:
            Fis[component.name]  = [Fi for Fi in Fis_load[component.name]]
            FisL[component.name] = [Fi for Fi in FisL_load[component.name]]

        try: Cn  = list(Cn.copy())
        except: pass
        try: CnL = list(CnL.copy())
        except: pass

        if self.depsize >= self.delzdep:                                                                                                                self.depcheck = -1
        else:
            if self.z[0] == 0:
                if self.Vdep * t >= self.depsize * 2:                                                                                                   self.depcheck = 2
                else:                                                                                                                                   self.depcheck = 0
            elif round(self.Vdep*t-abs(self.z[0]), 8) >= round(self.depsize, 8):                                                                        self.depcheck = 1
            else:                                                                                                                                       self.depcheck = 0

        if self.depcheck < 0:
            if self.z[0] == 0:
                if self.Vdep*t >= 2 * self.delzdep:
                    points = int(round((self.Vdep * t)/self.delzdep, 8))
                    self.ptot   = [0] + [p +  points for p  in self.p]
                    self.cptot  = [0] + [cp + points * self.deplayer[0].nchemicals + self.deplayer[0].nlowersolidchemicals for cp in self.cp]
                    self.cptot[1] = self.cptot[1] - self.deplayer[0].nlowersolidchemicals
                    self.cptotL = [0] + [cp + points * self.deplayer[0].nchemicals + self.deplayer[0].nlowersolidchemicals for cp in self.cpL]
                    self.layertot = self.deplayer + self.layers
                    for component in self.components:
                        Fis[component.name][0] = 0
                    for component in self.deplayer[0].init_components:
                        Fis[component.name][0] = component.fraction
                else:
                    points = 0
            else:
                points = int(round((self.Vdep*t -abs(self.z[0]))/self.delzdep, 8))
                self.ptot     = [ptot + points for ptot in self.ptot]
                self.cptot    = [cptot + points * self.deplayer[0].nchemicals for cptot in self.cptot]
                self.cptotL   = [cptot + points * self.deplayer[0].nchemicals for cptot in self.cptotL]
                self.ptot[0]  = 0
                self.cptot[0] = 0
                self.cptotL[0] = 0

            if points > 0:
                try: Cn = list(Cn.copy())
                except: pass
                Cn_solid = []
                for solidchemical in self.deplayer[0].solidchemicals:
                    Cn_solid.append(self.ICs[self.deplayer[0].name][solidchemical.name].uniform)
                if self.z[0] == 0:
                    Lowersolidchemicals = []
                    for soldichemical in self.deplayer[0].lowersolidchemicals:
                        Lowersolidchemicals.append(self.ICs[self.deplayer[0].name][soldichemical.name].uniform)
                    CnL = CnL[:self.layers[0].nchemicals] + Lowersolidchemicals + CnL[self.layers[0].nchemicals:]
                    Cn = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn[self.nchemicals:self.layers[0].nchemicals] + Lowersolidchemicals + Cn[self.layers[0].nchemicals:]
                else:
                    Cn = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn[self.deplayer[0].nchemicals:]
                for i in range(points):
                    CnL = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + CnL
                    Cn  = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn
                    for component in self.components:
                        Fis[component.name].insert(0, 0)
                        FisL[component.name].insert(0, 0)
                    for component in self.deplayer[0].init_components:
                        Fis[component.name][0] = component.fraction
                        FisL[component.name][0] = component.fraction
                    self.z.insert(0, round(self.z[0]-self.delzdep, 8))
                updatecheck = 1

        elif self.depcheck == 2:

            self.toppoints = 2
            self.ptot   = [0] + [p +  self.toppoints for p  in self.p]
            self.cptot  = [0] + [cp + self.toppoints * self.deplayer[0].nchemicals + self.deplayer[0].nlowersolidchemicals  for cp in self.cp]
            self.cptot[1] = self.cptot[1] - self.deplayer[0].nlowersolidchemicals
            self.cptotL = [0] + [cp + self.toppoints * self.deplayer[0].nchemicals + self.deplayer[0].nlowersolidchemicals for cp in self.cpL]
            try: Cn = list(Cn.copy())
            except: pass
            for component in self.components:
                Fis[component.name][0] = 0
            for component in self.deplayer[0].init_components:
                Fis[component.name][0] = component.fraction
            Cn_solid = []
            for solidchemical in self.deplayer[0].solidchemicals:
                Cn_solid.append(self.ICs[self.deplayer[0].name][solidchemical.name].uniform)
            if self.z[0] == 0:
                Lowersolidchemicals = []
                for soldichemical in self.deplayer[0].lowersolidchemicals:
                    Lowersolidchemicals.append(self.ICs[self.deplayer[0].name][soldichemical.name].uniform)
                CnL = CnL[:self.layers[0].nchemicals] + Lowersolidchemicals + CnL[self.layers[0].nchemicals:]
                Cn = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn[self.nchemicals:self.layers[0].nchemicals] + Lowersolidchemicals + Cn[self.layers[0].nchemicals:]
            else:
                Cn = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn[self.deplayer[0].nchemicals:]
            for j in range(self.toppoints):
                CnL = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + CnL
                Cn  = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn
                for component in self.components:
                    Fis[component.name].insert(0, 0)
                    FisL[component.name].insert(0, 0)
                for component in self.deplayer[0].init_components:
                    Fis[component.name][0]  = component.fraction
                    FisL[component.name][0] = component.fraction
                self.z.insert(0, - round(self.depsize * (j+1), 8))
            self.layertot = self.deplayer + self.layers

            updatecheck = 1

        elif self.depcheck == 1:

            updatecheck = 1
            points = int(round((self.Vdep*t-abs(self.z[0]))/self.depsize, 8))

            self.toppoints = self.toppoints + points
            self.ptot     = [ptot +  points for ptot in self.ptot]
            self.cptot    = [cptot + self.deplayer[0].nchemicals * points for cptot in self.cptot]
            self.ptot[0]  = 0
            self.cptot[0] = 0

            try: Cn = list(Cn.copy())
            except: pass
            Cn_solid = []
            for solidchemical in self.deplayer[0].solidchemicals:
                Cn_solid.append(self.ICs[self.deplayer[0].name][solidchemical.name].uniform)
            CnL = CnL
            Cn = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn[self.deplayer[0].nchemicals:]
            for i in range(points):
                CnL = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + CnL
                Cn  = [self.BCs[chemical.name].Cw for chemical in (self.chemicals)] + Cn_solid + Cn
                for component in self.components:
                    Fis[component.name].insert(0, 0)
                    FisL[component.name].insert(0, 0)
                for component in self.deplayer[0].init_components:
                    Fis[component.name][0]  = component.fraction
                    FisL[component.name][0] = component.fraction
                self.z.insert(0, round(self.z[0] - self.depsize, 8))

        if self.depcheck >= 0 and round(self.Vdep*t-abs(self.z[self.toppoints]), 8) > round(self.delzdep*1.5, 8) and round(self.Vdep*t, 8) > round(self.delzdep*2.5, 8):

            updatecheck = 1
            points = int(round(self.toppoints/self.depgrid, 8))
            toppoints_new = self.toppoints - points * self.depgrid

            if self.z[self.toppoints] == 0:
                self.ptot   = [0] + [p +  points * (1-self.depgrid) for p  in self.ptot]
                self.cptot  = [0] + [cp + points * (1-self.depgrid) * self.deplayer[0].nchemicals for cp in self.cptot]
                self.cptotL = [0] + [cp + points * (1-self.depgrid) * self.deplayer[0].nchemicals for cp in self.cptotL]
                self.ptot[1]  = toppoints_new
                self.cptot[1] = toppoints_new* self.deplayer[0].nchemicals
                self.cptotL[1] = toppoints_new* self.deplayer[0].nchemicals
            else:
                self.ptot      = [ptot +  points * (1-self.depgrid) for ptot in self.ptot]
                self.cptot     = [cptot + points * (1-self.depgrid) * self.deplayer[0].nchemicals for cptot in self.cptot]
                self.cptotL    = [cptot + points * (1-self.depgrid) * self.deplayer[0].nchemicals for cptot in self.cptotL]
                self.ptot[0]   = 0
                self.cptot[0]  = 0
                self.cptotL[0] = 0
                self.ptot[1]   = toppoints_new
                self.cptot[1]  = toppoints_new* self.deplayer[0].nchemicals
                self.cptotL[1] = toppoints_new* self.deplayer[0].nchemicals

            try: Cn = list(Cn.copy())
            except: pass

            z_top    = []
            Cn_top   = []
            Fis_top  = {}
            CnL_top  = []
            FisL_top = {}

            z_mid    = []
            Cn_mid   = []
            Fis_mid  = {}
            CnL_mid  = []
            FisL_mid = {}

            for toppoint in range(toppoints_new):
                for i in range(self.deplayer[0].nchemicals):
                    Cn_top.append(Cn[self.deplayer[0].nchemicals*toppoint+i])
                    CnL_top.append(CnL[self.deplayer[0].nchemicals*toppoint+i])
                z_top.append(self.z[toppoint])
            for component in self.components:
                Fis_top[component.name]  = []
                FisL_top[component.name] = []
                for toppoint in range(toppoints_new):
                    Fis_top[component.name].append(Fis[component.name][toppoint])
                    FisL_top[component.name].append(FisL[component.name][toppoint])

            for point in range(points):
                midpoint = toppoints_new + point * self.depgrid
                for i in range(self.deplayer[0].nchemicals):
                    Cn_mid.append(Cn[self.deplayer[0].nchemicals*midpoint+i])
                    CnL_mid.append(CnL[self.deplayer[0].nchemicals*midpoint+i])
            for component in self.components:
                Fis_mid[component.name] = []
                FisL_mid[component.name] = []
                for point in range(points):
                    midpoint = toppoints_new + point * self.depgrid
                    Fis_mid[component.name].append(Fis[component.name][midpoint])
                    FisL_mid[component.name].append(FisL[component.name][midpoint])

            for j in range (self.toppoints):
                self.z.remove(self.z[0])
            for point in range(points):
                self.z.insert(0, round(self.z[0]-self.delzdep, 8))
            for point in range(toppoints_new):
                self.z.insert(0, round(self.z[0]-self.depsize, 8))
            self.layertot = self.deplayer + self.deplayer + self.layers

            Cn  = Cn[self.toppoints* self.deplayer[0].nchemicals:]
            CnL = CnL[self.toppoints* self.deplayer[0].nchemicals:]
            for component in self.components:
                Fis[component.name]  = Fis[component.name][self.toppoints:]
                FisL[component.name] = FisL[component.name][self.toppoints:]

            Cn  = Cn_top  + Cn_mid + Cn
            CnL = CnL_top  + CnL_mid + CnL
            for component in self.components:
                Fis[component.name]  = Fis_top[component.name]  + Fis_mid[component.name] + Fis[component.name]
                FisL[component.name] = FisL_top[component.name] + FisL_mid[component.name] + FisL[component.name]

            self.toppoints = toppoints_new

        if updatecheck == 1:

            if self.bio == 1:
                self.update_bioturbation()
                self.make_grid_Dbiops()
            else:
                self.ptotbio     = [ptot               for ptot    in self.ptot]
                self.cptotbio    = [cptot              for cptot   in self.cptot]
                self.cptotLbio   = [cptotL             for cptotL  in self.cptotL]
                self.layertotbio = [layer.bio_copy()   for layer   in self.layertot]

            self.make_matrix_parameter_vectors(Cn, CnL, Fis, FisL)
            self.make_transport_parameter_vectors()
            self.make_reaction_parameter_vectors(Cn, CnL, Fis, FisL)

        return Cn, CnL, Fis, FisL, updatecheck #, FisM



    def update_consolidation(self, time, U):
        """Updates the Darcy velocity and then the equations if there is
        consolidation in the underlying sediment."""

        self.U_plus_1 = U + consolidation(time, self.Vcon0, self.kcon)

    def update_tidal(self, time, U):
        """Updates the Darcy velocity and then the equations if there is
        consolidation in the underlying sediment."""

        self.U_plus_1 = U + tidal(time, self.Vtidal, self.ptidal)

        if self.topBCtype == 'Finite mixed water column':
            Q = self.taucoefs['Q'] - self.taucoefs['Qevap'] + self.taucoefs['V']/self.taucoefs['h']*self.U

            if self.dep == 1:
                if self.lengthunit == 'cm':         kdep  = self.Vdep/100
                elif self.lengthunit == u'\u03BCm': kdep  = self.Vdep/100/1000
                else:                               kdep  = self.Vdep
                epsilon = self.matrices[self.deplayer[0].type_index].e
                matrixa  = self.matrices[self.deplayer[0].type_index]
                dep_rho = matrixa.rho
                dep_fraction = []
                for component in matrixa.components:
                    dep_fraction.append(component.mfraction)
            else:
                kdep = 0
                epsilon = 0
                dep_rho = 0
                dep_fraction = []

            for n in range(len(self.chemicals)):
                chemical = self.chemicals[n]
                if self.taucoefs['Evap'] == 'None':  kevap = 0
                else:                                kevap = self.BCs[chemical.name].kevap
                if self.taucoefs['Decay'] == 'None': kdecay = 0
                else:                                kdecay = self.BCs[chemical.name].kdecay

                K = []
                if self.dep == 1:
                    matrixa  = self.matrices[self.deplayer[0].type_index]
                    for component in matrixa.components:
                        if self.sorptions[component.name][chemical.name].isotherm == 'Linear--Kd specified': K.append(self.sorptions[component.name][chemical.name].K)
                        elif self.sorptions[component.name][chemical.name].isotherm == 'Linear--Kocfoc': K.append(10**self.sorptions[component.name][chemical.name].Koc)

                self.BCs[chemical.name].tau = tauwater(Q, self.taucoefs['DOC'], chemical.Kdoc, kdep, epsilon, dep_rho, dep_fraction, K, self.taucoefs['V'], self.taucoefs['h'], kevap, kdecay)

    def update_Cw(self, C):

        for n in range(len(self.chemicals)):
            chemical = self.chemicals[n]
            h = (self.taucoefs['h']*self.h_factor)/(1+self.taucoefs['DOC']/10**6*10**chemical.Kdoc)
            if self.Us_plus_1[n] >= 0:
                self.BCs[chemical.name].Cw_plus_1 = (self.BCs[chemical.name].Cw+(self.BCs[chemical.name].k+self.Us_plus_1[n])*self.delt/h*C[n])/(1. + self.BCs[chemical.name].tau*self.delt+self.BCs[chemical.name].k/h*self.delt)
            else:
                self.BCs[chemical.name].Cw_plus_1 = (self.BCs[chemical.name].Cw+(self.BCs[chemical.name].k)*self.delt/h*C[n])/(1. + self.BCs[chemical.name].tau*self.delt+(self.BCs[chemical.name].k-self.Us_plus_1[n])/h*self.delt)



    def update_nonlinear(self, Cn, CnL, Fis, FisL):
        """Updates the retardation factors, delt, and governing equations for
        nonlinear sorption."""


        Cn  = transpose(array([i for i in Cn]))
        CnL = transpose(array([i for i in CnL]))

        if self.sorp == 1:
            self.make_grid_Rs(Cn, CnL, Fis, FisL)
            self.make_grid_Ds()

        if self.reac == 1:
            self.make_grid_rates_plus_1(Cn, CnL, Fis, FisL)
            self.make_grid_rates_diff(Cn, CnL, Fis, FisL)

        self.A = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, self.cptot[-1] + self.layertot[-1].nchemicals]))
        self.B = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, self.cptot[-1] + self.layertot[-1].nchemicals]))
        self.a = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, 1]))
        self.b = matrix(zeros([self.cptot[-1] + self.layertot[-1].nchemicals, 1]))

        self.make_boundary_equations()
        self.make_governing_equations()
        self.make_Newton_Raphson_equations(Cn, CnL, Fis, FisL)

    def non_linear_solver(self, Cn, CnL, Fis, FisL):

        convergence_check = 0
        B_old = self.B.copy()
        b_old = self.b.copy()
        RMSE_old = 100000

        #self.update_nonlinear(Cn, CnL, Fis, FisL)
        Cn_init = [C for C in CnL]

        Newton_fail_check = 0

        count = 0
        Cn_old = [C for C in Cn_init]

        while convergence_check == 0 and Newton_fail_check == 0:

            self.update_nonlinear(Cn_old, Cn_old, Fis, FisL)
            Cn_new = Cn_old - array(transpose(linalg.solve(self.NR, self.A * transpose(matrix(Cn_old)) + self.a - B_old * transpose(matrix(Cn)) - b_old)))[0]

            SE              = []
            AG              = 0

            for i in range(len(Cn_new)):
                SE.append(((Cn_new[i]-Cn_old[i])**2)**0.5)
                AG = AG + abs(Cn_new[i])/2+abs(Cn_old[i])/2
                if Cn_new[i] < 0:  negative_check = 1

            RMSE = sum(SE)/AG

            if RMSE_old <= RMSE:
                Newton_fail_check = 1

            if RMSE < self.nlerror:
                convergence_check = 1
            else:
                Cn_old = Cn_new

            RMSE_old = RMSE

        if Newton_fail_check == 1:

            Cn_old = [C for C in Cn_init]

            self.update_nonlinear(Cn_old, Cn_old, Fis, FisL)

            while convergence_check == 0:

                Cn_new = array((transpose(linalg.solve(self.A, B_old * transpose(matrix(Cn)) + b_old - self.a))))[0]

                SE              = []
                AG              = 0

                for i in range(len(Cn_new)):
                    SE.append(((Cn_new[i]-Cn_old[i])**2)**0.5)
                    AG = AG + abs(Cn_new[i])/2+abs(Cn_old[i])/2
                    if Cn_new[i] < 0:  negative_check = 1

                RMSE = sum(SE)/AG

                if RMSE < self.nlerror:
                    convergence_check = 1
                else:
                    Cn_old = Cn_new
                self.update_nonlinear(Cn_old, Cn_old, Fis, FisL)

        return Cn_new

    def get_Fis_plus_1(self, Fis, FisL):
        """Uses the matrices to solve the system.  Returns the concentrations
        at the next time step in a one-dimensional row array."""
        Fis_plus_1  = {}
        FisL_plus_1 = {}

        for component in self.components:
            Fis_plus_1[component.name]  = list(array((transpose(linalg.solve(self.Acomp, self.Bcomp * transpose(matrix(Fis[component.name]))))))[0])

            if self.pbio < self.ptot[-1]:
                FisL_plus_1[component.name]    = Fis_plus_1[component.name][:self.pbio] + FisL[component.name][self.pbio:]
            else:
                FisL_plus_1[component.name]    = Fis_plus_1[component.name][:]

        return Fis_plus_1, FisL_plus_1

    def get_Cn_plus_1(self, Cn):
        """Uses the matrices to solve the system.  Returns the concentrations
        at the next time step in a one-dimensional row array."""

        return array((transpose(linalg.solve(self.A, self.B * transpose(matrix(Cn)) + self.b - self.a))))[0]



    def get_Cws(self):

        Cws         = zeros(self.nchemicals)
        Cws_plus_1  = zeros(self.nchemicals)

        for n in range (self.nchemicals):
            Cws[n]        = self.BCs[self.chemicals[n].name].Cw * self.chemicals[n].MW
            Cws_plus_1[n] = self.BCs[self.chemicals[n].name].Cw_plus_1 * self.chemicals[n].MW

        return Cws, Cws_plus_1

    def get_DM(self, DM, C, Fis, DM_depth):

        DM_plus_1 = DM.copy()

        for n in range (self.nchemicals):
            chemical = self.chemicals[n]
            for i in range(len(self.z)):
                if self.z < DM_depth:
                    for component in self.components:
                        DM_plus_1[n] = DM_plus_1[n] + self.sorptions[component.name][chemical.name].get_q(C[i*self.layertot[0].chemicals + n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho * (self.z[i+1]-self.z[i])*self.flux_factor *chemical.MW
                    DM_plus_1[n] = DM_plus_1[n] + C[i*self.layertot[0].chemicals + n] * self.es_plus_1[i]*(1+self.layertot[0].doc/(10**6)*10**chemical.Kdoc) * (self.z[i+1]-self.z[i])*self.flux_factor *chemical.MW

        DM_depth = self.z[0]

        return DM_plus_1, DM_depth


    def get_FM(self, FM, CL):

        FM_plus_1 = FM.copy()

        for n in range (self.nchemicals):
            chemical = self.chemicals[n]
            if self.topBCtype =='Mass transfer' or self.topBCtype == 'Finite mixed water column':
                FM_plus_1[0,n] = FM_plus_1[0,n] + get_surface_flux(self.Us_plus_1[n], self.BCs[chemical.name].k, CL[n], self.BCs[chemical.name].Cw ) * self.flux_factor *chemical.MW * self.delt
            else:
                if self.bio == 1 and self.Dbiop > 0:
                    KDbiops_plus_1  = [self.DbiopsL_plus_1[0]*self.KsL_plus_1[n], self.DbiopsL_plus_1[0]*self.Ks_plus_1[self.layertot[0].nchemicals+n]]
                else:
                    KDbiops_plus_1 = [0,0]
                FM_plus_1[0,n] = FM_plus_1[0,n] + get_top_flux(self.Us_plus_1[n], self.Ds_plus_1[n], KDbiops_plus_1, array([CL[n], CL[self.layertot[0].nchemicals+n]]), self.z[0:2])*chemical.MW * self.flux_factor * self.delt

            i   = self.ptot[-1]
            if self.botBCtype =='Flux-matching':
                FM_plus_1[1,n] = FM_plus_1[1,n] + self.U_plus_1 * self.BCs[chemical.name].Cb * chemical.MW * self.flux_factor * self.delt
            else:
                if self.bio == 1 and self.Dbiop > 0:
                    KDbiops_plus_1  = [self.Dbiops_plus_1[i]*self.Ks_plus_1[self.cptot[-1]-self.layertot[-1].nchemicals+n], self.Dbiops_plus_1[i]*self.Ks_plus_1[self.cptot[-1]+n]]
                else:
                    KDbiops_plus_1 = [0,0]
                FM_plus_1[1,n] = FM_plus_1[1,n] + get_boundary_flux(self.Us_plus_1[self.cptot[-1]+n], self.Ds_plus_1[self.cptot[-1]+n] , KDbiops_plus_1, array([CL[self.cptot[-1]-self.layertot[-1].nchemicals+n], CL[self.cptot[-1]+n]]), self.z[-2:])*chemical.MW * self.flux_factor * self.delt

        return FM_plus_1

    def get_fluxes(self, C, CL, O = None, OL = None, flag = None, Cn_top = None, O_top = None):
        """Calculates the fluxes in the domain "z" using the concentration at
        each grid point "C," the Darcy velocity "U," and the diffusion
        coefficient at each grid point "D." Converts units to ug/m2/yr."""

        U       = self.U_plus_1
        Ds      = self.Ds_plus_1

        ptot    = []
        cptot   = []
        cptotL  = []
        layers  = []
        for i in self.ptot:   ptot.append(i)
        for i in self.cptot:  cptot.append(i)
        for i in self.cptotL: cptotL.append(i)

        for layer in self.layertot:
            layers.append(layer)

        i    = 0
        if self.bio == 1:

            if ptot.count(self.pbio) == 0:
                ptot.insert(self.layerbio + 1, self.pbio)
                cptot.insert(self.layerbio + 1, self.cpbio)
                cptotL.insert(self.layerbio + 1, self.cpbioL)
                layers.insert(self.layerbio + 1, self.layertot[self.layerbio])

            if self.dep == 1 and self.toppoints > 1:
                F = zeros(((len(self.z)-self.toppoints + 1), self.nchemicals))

                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for jj in range(len(ptot) - 2):
                        j = jj + 1
                        for i in range(ptot[j] + 1, ptot[j+1]):
                            ii = i - ptot[j]
                            iii = i - self.toppoints + 1
                            if ii == 1: cp = cptot[j]
                            else:       cp = cptotL[j]
                            if self.Dbiop > 0:
                                KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cp+(ii-1)*layers[j].nchemicals+n],
                                                   self.DbiopsL_plus_1[i]*self.Ks_plus_1[cptotL[j]+ii*layers[j].nchemicals+n],
                                                   self.DbiopsL_plus_1[i]*self.Ks_plus_1[cptotL[j]+(ii+1)*layers[j].nchemicals+n]]
                            else:
                                KDbiops_plus_1 = [0,0,0]
                            F[iii,n] =  get_point_flux(U* (1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptotL[j]+ii*layers[j].nchemicals+n], KDbiops_plus_1, C[iii-1:iii+2,n], self.z[i-1:i+2]) * self.flux_factor

                        if j < len(ptot) - 2:
                            i = ptot[j+1]
                            iii = i - self.toppoints + 1
                            if self.Dbiop > 0:
                                KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j+1]+n],
                                                   self.DbiopsL_plus_1[i]*self.Ks_plus_1[cptotL[j+1]+layers[j+1].nchemicals+n]]
                            else:
                                KDbiops_plus_1 = [0,0]
                            F[iii,n] = get_top_flux(U* (1+layers[j+1].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[cptot[j+1]+n], KDbiops_plus_1, CL[iii:iii+2,n], self.z[i:i+2]) * self.flux_factor


                    if self.topBCtype =='Mass transfer' or self.topBCtype == 'Finite mixed water column':
                        F[0,n] = get_surface_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.BCs[chemical.name].k, Cn_top[0,n], self.BCs[chemical.name].Cw ) * self.flux_factor
                    else:
                        KDbiops_plus_1  = [self.Dbiops_plus_1[0]*self.KsL_plus_1[n],
                                           self.Dbiops_plus_1[0]*self.Ks_plus_1[layers[0].nchemicals+n]]
                        F[0,n] = get_top_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[n], KDbiops_plus_1, Cn_top[0:2, n], self.z[0:2]) * self.flux_factor

                    i   = ptot[-1]
                    KDbiops_plus_1  = [self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[-1]-layers[-1].nchemicals+n],
                                       self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[-1]+n]]
                    F[-1,n] = get_boundary_flux(U* (1+layers[-1].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[-1]+n] , KDbiops_plus_1, C[-2:,n], self.z[-2:]) * self.flux_factor

                if self.nsolidchemicals > 0:
                    for jj in range(len(ptot) - 2):
                        j = jj + 1
                        num = layers[j].nchemicals
                        layer = layers[j]
                        for solidchemical in layers[j].solidchemicals:
                            n  = self.solidchemical_list.index(solidchemical.name)
                            nn = self.chemical_list.index(solidchemical.chemical_name)
                            na = layer.solidchemical_list.index(solidchemical.name)+self.nchemicals
                            for i in range(ptot[j] + 1, ptot[j+1]-1):
                                ii = i - ptot[j]
                                iii = i - self.toppoints + 1
                                if ii == 1: cp = cptot[j]
                                else:       cp = cptotL[j]
                                if i <= self.pbio:
                                    KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cp+(ii-1)*num+na],
                                                       self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+ii*num+na],
                                                       self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+(ii+1)*num+na]]
                                    F[iii,nn] =  F[iii,nn] + get_point_flux(0, 0, KDbiops_plus_1, O[iii-1:iii+2,n], self.z[i-1:i+2]) * self.flux_factor

                            i = ptot[j+1]-1
                            iii = i - self.toppoints + 1
                            if i <= self.pbio:
                                na = layers[j].solidchemical_list.index(solidchemical.name) + self.nchemicals
                                nb = layers[j].lowersolidchemicalnums[layers[j].solidchemical_list.index(solidchemical.name)]
                                KDbiops_plus_1  = [ self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j+1]-2*num+na],
                                                    self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]-num+na],
                                                    self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]+nb]]
                                F[iii,nn] =  F[iii,nn] + get_point_flux(0, 0, KDbiops_plus_1, O[iii-1:iii+2,n], self.z[i-1:i+2]) * self.flux_factor

                            i = ptot[j]
                            iii = i - self.toppoints + 1
                            if i <= self.pbio and j < len(ptot) - 2:
                                nb = layers[j].solidchemical_list.index(solidchemical.name)+self.nchemicals
                                KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j]+nb],
                                                   self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptotL[j]+layers[j].nchemicals+nb]]
                                F[iii,nn] =  F[iii,nn] + get_top_flux(0, 0, KDbiops_plus_1, OL[iii:iii+2,n], self.z[i:i+2]) * self.flux_factor

                    for solidchemical in layers[0].solidchemicals:
                        n  = self.solidchemical_list.index(solidchemical.name)
                        nn = self.chemical_list.index(solidchemical.chemical_name)
                        nb = layers[0].solidchemical_list.index(solidchemical.name)+self.nchemicals
                        KDbiops_plus_1  = [self.Dbiops_plus_1[0]*self.Ks_plus_1[nb],
                                           self.Dbiops_plus_1[0]*self.Ks_plus_1[layers[0].nchemicals+nb]]
                        F[0,nn] = F[0,nn] + get_top_flux(0,0, KDbiops_plus_1, O[0:2,nn], self.z[0:2]) * self.flux_factor

                    i = ptot[-1]
                    if i <= self.pbio:
                        for solidchemical in layers[-1].solidchemicals:
                            nn = self.chemical_list.index(solidchemical.chemical_name)
                            nb = layers[-1].solidchemical_list.index(solidchemical.name)+self.nchemicals
                            KDbiops_plus_1  = [self.Dbiops_plus_1[-1]*self.Ks_plus_1[cptot[-1]-layers[-1].nchemicals+nb],
                                               self.Dbiops_plus_1[-1]*self.Ks_plus_1[cptot[-1]+nb]]
                            F[-1,nn] = get_boundary_flux(0, 0, KDbiops_plus_1, C[-2:,nn], self.z[-2:]) * self.flux_factor

            else:
                F = zeros((len(self.z), self.nchemicals))
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for j in range(len(ptot) - 1):
                        for i in range(ptot[j] + 1, ptot[j+1]-1):
                            ii = i - ptot[j]
                            if ii == 1: cp = cptot[j]
                            else:       cp = cptotL[j]
                            KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cp+(ii-1)*layers[j].nchemicals+n],
                                               self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+ii*layers[j].nchemicals+n],
                                               self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+(ii+1)*layers[j].nchemicals+n]]
                            F[i,n] =  get_point_flux(U*(1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptotL[j]+ii*layers[j].nchemicals+n], KDbiops_plus_1, CL[i-1:i+2,n], self.z[i-1:i+2]) * self.flux_factor

                        i = ptot[j+1]-1
                        KDbiops_plus_1  = [ self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j+1]-2*layers[j].nchemicals+n],
                                            self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]-layers[j].nchemicals+n],
                                            self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]+n]]
                        F[i,n] =  get_point_flux(U*(1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[j+1]-layers[j].nchemicals+n], KDbiops_plus_1, C[i-1:i+2,n], self.z[i-1:i+2]) * self.flux_factor

                        if j < len(ptot) - 2:
                            i = ptot[j+1]
                            KDbiopsL_plus_1  = [self.DbiopsL_plus_1[i] * self.KsL_plus_1[cptot[j+1]+n],
                                                self.DbiopsL_plus_1[i] * self.KsL_plus_1[cptotL[j+1]+layers[j+1].nchemicals+n]]
                            F[i,n] = get_top_flux(U* (1+layers[j+1].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[cptot[j]+n], KDbiopsL_plus_1, CL[i:i+2,n], self.z[i:i+2]) * self.flux_factor

                    if self.topBCtype =='Mass transfer' or self.topBCtype == 'Finite mixed water column':
                        F[0,n] = get_surface_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.BCs[chemical.name].k, C[0,n], self.BCs[chemical.name].Cw ) * self.flux_factor
                    else:
                        KDbiops_plus_1  = [self.DbiopsL_plus_1[0]*self.KsL_plus_1[n],
                                           self.DbiopsL_plus_1[0] *self.Ks_plus_1[layers[0].nchemicals+n]]
                        F[0,n] = get_top_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[n], KDbiops_plus_1, C[0:2,n], self.z[0:2]) * self.flux_factor

                    i = ptot[-1]
                    KDbiops_plus_1  = [self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[-1]-layers[-1].nchemicals+n],
                                       self.Dbiops_plus_1[i] *self.Ks_plus_1[cptot[-1]+n]]
                    F[-1,n] = get_boundary_flux(U* (1+layers[-1].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[-1]+n] ,KDbiops_plus_1,  C[-2:,n], self.z[-2:]) * self.flux_factor
                if self.nsolidchemicals > 0:
                    for j in range(len(ptot) - 1):
                        layer = layers[j]
                        for solidchemical in layer.solidchemicals:
                            n  = self.solidchemical_list.index(solidchemical.name)
                            nn = self.chemical_list.index(solidchemical.chemical_name)
                            na = layer.solidchemical_list.index(solidchemical.name)+self.nchemicals
                            nb = layers[j].lowersolidchemicalnums[layers[j].solidchemical_list.index(solidchemical.name)]
                            num = layer.nchemicals
                            for i in range(ptot[j] + 1, ptot[j+1] - 1):
                                if i <= self.pbio:
                                    ii = i - ptot[j]
                                    if ii == 1: cp = cptot[j]
                                    else:       cp = cptotL[j]
                                    KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cp+(ii-1)*num+na],
                                                       self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+ii*num+na],
                                                       self.Dbiops_plus_1[i]*self.Ks_plus_1[cptotL[j]+(ii+1)*num+na]]
                                    F[i,nn] =  F[i,nn] + get_point_flux(0, 0, KDbiops_plus_1, OL[i-1:i+2,n], self.z[i-1:i+2]) * self.flux_factor

                            i = ptot[j+1]-1
                            if i < self.pbio:
                                KDbiops_plus_1  = [ self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j+1]-2*num+na],
                                                    self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]-num+na],
                                                    self.Dbiops_plus_1[i]*self.Ks_plus_1[cptot[j+1]+nb]]
                                F[i,nn] =  F[i,nn] + get_point_flux(0, 0, KDbiops_plus_1, O[i-1:i+2,n], self.z[i-1:i+2]) * self.flux_factor

                            i = ptot[j]
                            if i < self.pbio:
                                KDbiops_plus_1  = [self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptot[j]+na],
                                                   self.DbiopsL_plus_1[i]*self.KsL_plus_1[cptotL[j]+num+na]]
                                F[i,nn] =  F[i,nn] + get_top_flux(0, 0, KDbiops_plus_1, OL[i:i+2,n], self.z[i:i+2]) * self.flux_factor

                    if ptot[-1] <= self.pbio:
                        layer = layers[-1]
                        for solidchemical in layer.solidchemicals:
                            nn = self.chemical_list.index(solidchemical.chemical_name)
                            na = layer.solidchemical_list.index(solidchemical.name)
                            num = layer.nchemicals
                            KDbiops_plus_1  = [self.Dbiops_plus_1[-1]*self.Ks_plus_1[cptot[-1]-num+na],
                                               self.Dbiops_plus_1[-1]*self.Ks_plus_1[cptot[-1]+na]]
                            F[-1,nn] = get_boundary_flux(0, 0, KDbiops_plus_1, C[-2:,nn], self.z[-2:]) * self.flux_factor

        else:
            if self.dep == 1 and self.toppoints > 1:
                F = zeros(((len(self.z)-self.toppoints + 1), self.nchemicals))

                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for jj in range(len(ptot) - 2):
                        j = jj + 1
                        for i in range(ptot[j] + 1, ptot[j+1]):
                            ii = i - ptot[j]
                            iii = i - self.toppoints + 1
                            KDbiops_plus_1 =  [0,0,0]
                            F[iii,n] =  get_point_flux(U* (1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptotL[j]+ii*layers[j].nchemicals+n], KDbiops_plus_1, C[iii-1:iii+2,n], self.z[i-1:i+2]) * self.flux_factor

                        if j < len(ptot) - 2:
                            i = ptot[j+1]
                            iii = i - self.toppoints + 1
                            KDbiops_plus_1 =  [0,0]
                            F[iii,n] = get_top_flux(U* (1+layers[j+1].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[cptot[j+1]+n], KDbiops_plus_1, C[iii:iii+2,n], self.z[i:i+2]) * self.flux_factor

                    if self.topBCtype =='Mass transfer' or self.topBCtype == 'Finite mixed water column':
                        F[0,n] = get_surface_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.BCs[chemical.name].k, Cn_top[0,n], self.BCs[chemical.name].Cw ) * self.flux_factor
                    else:
                        KDbiops_plus_1 =  [0,0]
                        F[0,n] = get_top_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[n], KDbiops_plus_1, Cn_top[0:2, n], self.z[0:2]) * self.flux_factor

                    #if self.bio == 1 and self.pbio <> ptot[-1]:
                    KDbiops_plus_1 =  [0,0]
                    F[-1,n] = get_boundary_flux(U* (1+layers[-1].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[-1]+n] , KDbiops_plus_1, C[-2:,n], self.z[-2:]) * self.flux_factor
            else:
                F = zeros((len(self.z), self.nchemicals))

                for n in range (self.nchemicals):
                    #top boundary
                    chemical = self.chemicals[n]

                    for j in range(len(ptot) - 1):
                        #loop through the interior points of the layers
                        for i in range(ptot[j] + 1, ptot[j+1]):
                            ii = i - ptot[j]
                            KDbiops_plus_1 =  [0,0,0]
                            F[i,n] =  get_point_flux(U*(1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptotL[j]+ii*layers[j].nchemicals+n], KDbiops_plus_1, C[i-1:i+2,n], self.z[i-1:i+2]) * self.flux_factor

                        i = ptot[j+1]
                        KDbiops_plus_1 =  [0,0]
                        F[i,n] = get_boundary_flux(U* (1+layers[j].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[j+1]+n], KDbiops_plus_1, C[i-1:i+1,n], self.z[i-1:i+1]) * self.flux_factor

                    if self.topBCtype =='Mass transfer' or self.topBCtype == 'Finite mixed water column':
                        F[0,n] = get_surface_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.BCs[chemical.name].k, C[0,n], self.BCs[chemical.name].Cw ) * self.flux_factor
                    else:
                        KDbiops_plus_1 =  [0,0]
                        F[0,n] = get_top_flux(U* (1+layers[0].doc/(10**6)*10**chemical.Kdoc), self.DsL_plus_1[n], KDbiops_plus_1, C[0:2,n], self.z[0:2]) * self.flux_factor

                    #if self.bio == 1 and self.pbio <> ptot[-1]:
                    KDbiops_plus_1 =  [0,0]
                    F[-1,n] = get_boundary_flux(U* (1+layers[-1].doc/(10**6)*10**chemical.Kdoc), self.Ds_plus_1[cptot[-1]+n] ,KDbiops_plus_1,  C[-2:,n], self.z[-2:]) * self.flux_factor

        return F

    def get_qs(self, C, CL, Fis, FisL, O = None, OL = None):
        """Calculates the solid-phase concentration "q" in the domain "z"
        using the pore water concentration at each grid point "C." """

        if self.dep == 1 and self.toppoints > 1:
            q = zeros(((len(self.z)-self.toppoints + 1), self.nchemicals))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for jj in range(len(self.ptotbio) - 2):
                    j = jj + 1
                    for i in range(self.ptotbio[j], self.ptotbio[j+1]+1):
                        iii = i - self.toppoints + 1
                        q[iii, n] = 0
                        matrix_rho = 0
                        for component in self.components:
                            q[iii,n] = q[iii,n] + self.sorptions[component.name][chemical.name].get_q(C[iii,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho
                            matrix_rho = matrix_rho + (Fis[component.name][i]+FisL[component.name][i])/2 * component.rho
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for solidchemical in self.solidchemicals:
                                    if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                        q[iii,n] = q[iii,n] + O[iii,self.solidchemicals.index(solidchemical)]* Fis[component.name][i] * component.rho
                        q[iii,n] = q[iii,n]/matrix_rho

            qL = q.copy()
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for jj in range(len(self.ptotbio) - 2):
                    j = jj + 1
                    i = self.ptotbio[j+1]
                    iii = i - self.toppoints + 1
                    qL[iii, n] = 0
                    matrix_rho = 0
                    for component in self.components:
                        qL[iii,n] = qL[iii,n] + self.sorptions[component.name][chemical.name].get_q(CL[iii,n], self.Cmax[chemical.name]) * FisL[component.name][i] * component.rho
                        matrix_rho = matrix_rho + FisL[component.name][i] * component.rho
                        if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            for solidchemical in self.solidchemicals:
                                if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                    qL[iii,n] = qL[iii,n] + OL[iii,self.solidchemicals.index(solidchemical)]* FisL[component.name][i] * component.rho
                    qL[iii,n] = qL[iii,n]/matrix_rho

        else:
            q = zeros((len(self.z), self.nchemicals))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptotbio) - 1):
                    for i in range(self.ptotbio[j], self.ptotbio[j+1]+1):
                        q[i, n] = 0
                        matrix_rho = 0
                        for component in self.components:
                            q[i,n] = q[i,n] + self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho
                            matrix_rho = matrix_rho + Fis[component.name][i] * component.rho
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for solidchemical in self.solidchemicals:
                                    if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                        q[i,n] = q[i,n] + O[i,self.solidchemicals.index(solidchemical)]* Fis[component.name][i] * component.rho
                        q[i,n] = q[i,n]/matrix_rho

            qL = q.copy()

            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptotbio) - 1):
                    i = self.ptotbio[j+1]
                    qL[i, n] = 0
                    matrix_rho = 0
                    for component in self.components:
                        qL[i,n] = qL[i,n] + self.sorptions[component.name][chemical.name].get_q(CL[i,n], self.Cmax[chemical.name]) * FisL[component.name][i] * component.rho
                        matrix_rho = matrix_rho + FisL[component.name][i] * component.rho
                        if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            for solidchemical in self.solidchemicals:
                                if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                    qL[i,n] = qL[i,n] + OL[i,self.solidchemicals.index(solidchemical)]* FisL[component.name][i] * component.rho
                    qL[i,n] = qL[i,n]/matrix_rho

        return q, qL

    def get_qms(self, C, CL, Fis, FisL,O = None, OL = None):
        """Calculates the solid-phase concentration "q" in the domain "z"
        using the pore water concentration at each grid point "C." """

        if self.dep == 1 and self.toppoints > 1:
            qm = zeros(((len(self.z)-self.toppoints + 1), self.nchemicals, len(self.components)))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for jj in range(len(self.ptotbio) - 2):
                    j = jj +1
                    for i in range(self.ptotbio[j], self.ptotbio[j+1]+1):
                        iii = i - self.toppoints + 1
                        for component in self.components:
                            m = self.component_list.index(component.name)
                            if Fis[component.name][i] > 0.0001:
                                qm[iii,n,m] = self.sorptions[component.name][chemical.name].get_q(C[iii,n], self.Cmax[chemical.name])
                                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                    for solidchemical in self.solidchemicals:
                                        if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                            qm[iii,n,m] = qm[iii,n,m] + O[iii,self.solidchemicals.index(solidchemical)]



            qmL = qm.copy()
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for jj in range(len(self.ptotbio) - 2):
                    j = jj +1
                    i = self.ptotbio[j+1]
                    iii = i - self.toppoints + 1
                    for component in self.components:
                        m = self.component_list.index(component.name)
                        qmL[iii,n,m] = 0
                        if Fis[component.name][i] > 0.0001:
                            qmL[iii,n,m] = self.sorptions[component.name][chemical.name].get_q(CL[iii,n], self.Cmax[chemical.name])
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for solidchemical in self.solidchemicals:
                                    if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                        qmL[iii,n,m] = qmL[iii,n,m] + OL[iii,self.solidchemicals.index(solidchemical)]

        else:
            qm = zeros((len(self.z), self.nchemicals, len(self.components)))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptotbio) - 1):
                    for i in range(self.ptotbio[j], self.ptotbio[j+1]+1):
                        for component in self.components:
                            m = self.component_list.index(component.name)
                            if Fis[component.name][i] > 0.0001:
                                qm[i,n,m] = self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name])
                                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                    for solidchemical in self.solidchemicals:
                                        if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                            qm[i,n,m] = qm[i,n,m] + O[i,self.solidchemicals.index(solidchemical)]

            qmL = qm.copy()
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptotbio) - 2):
                    i = self.ptotbio[j+1]
                    for component in self.components:
                        m = self.component_list.index(component.name)
                        qmL[i,n,m] = 0
                        if FisL[component.name][i] > 0.0001:
                            qmL[i,n,m] = self.sorptions[component.name][chemical.name].get_q(CL[i,n], self.Cmax[chemical.name])
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for solidchemical in self.solidchemicals:
                                    if solidchemical.chemical_name == chemical.name and solidchemical.component_name == component.name:
                                        qmL[i,n,m] = qmL[i,n,m] + OL[i,self.solidchemicals.index(solidchemical)]

        return qm, qmL


    def get_Ws(self,C):
        """Calculates the total (aqueous + solid-phase) concentration "W" in
        the domain "z" using the pore water concentration at each grid point
        "C." """

        if self.dep == 1 and self.toppoints > 1:
            W = zeros(((len(self.z)-self.toppoints + 1), self.nchemicals))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for jj in range(len(self.ptot) - 2):
                    j = jj + 1
                    layer  = self.layertot[j]
                    for i in range(self.ptot[j], self.ptot[j+1]+1):
                        iii = i - self.toppoints + 1
                        W[iii, n] = C[iii, n] *(1+layer.doc/(10**6)*10**chemical.Kdoc)

        else:
            W = zeros((len(self.z), self.nchemicals))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptot) - 1):
                    layer  = self.layertot[j]
                    for i in range(self.ptot[j], self.ptot[j+1]+1):
                        W[i, n] = C[i, n] *(1+layer.doc/(10**6)*10**chemical.Kdoc)

        return W

    def get_Ms(self,C, CL, Fis, FisL, O = None, OL = None):
        """Calculates the total (aqueous + solid-phase) concentration "W" in
        the domain "z" using the pore water concentration at each grid point
        "C." """

        M = zeros((len(self.layers+self.deplayer), self.nchemicals))

        if self.dep == 1 and len(self.layertot) < len(self.layers+self.deplayer):
            if self.toppoints > 1:
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for jj in range(len(self.ptot) - 2):
                        j = jj + 2
                        layer  = self.layertot[j]
                        M[j, n] = 0
                        for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                            iii = i - self.toppoints + 1
                            if i == self.pbio:
                                for component in self.components:
                                    M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(CL[iii,n], self.Cmax[chemical.name]) * FisL[component.name][i] * component.rho * (self.z[i+1]-self.z[i])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j,n] = M[j,n] + OL[iii,nn]* FisL[component.name][i] * component.rho* (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                                for component in self.components:
                                    M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(C[iii,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho * (self.z[i]-self.z[i-1])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j,n] = M[j,n] + O[iii,nn]* Fis[component.name][i] * component.rho*(self.z[i]-self.z[i-1])/2*self.flux_factor
                                M[j, n] = (M[j, n] + CL[iii,n] * self.esL_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i+1]-self.z[i])  /2*self.flux_factor
                                                      + C[iii,n]  * self.es_plus_1[i] *(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i]  -self.z[i-1])/2*self.flux_factor)
                            else:
                                for component in self.components:
                                    M[j, n] = M[j, n] + (self.sorptions[component.name][chemical.name].get_q(C[iii,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho) * (self.z[i+1]-self.z[i-1])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j,n] = M[j,n] + O[iii,nn]* Fis[component.name][i] * component.rho * (self.z[i+1]-self.z[i-1])/2*self.flux_factor
                                M[j, n] = M[j, n] + C[iii,n] * self.es_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc)*(self.z[i+1]-self.z[i-1])/2*self.flux_factor

                        for component in self.components:
                            M[j, n] = M[j, n] + (self.sorptions[component.name][chemical.name].get_q(CL[self.ptot[j]- self.toppoints + 1,n], self.Cmax[chemical.name]) * FisL[component.name][self.ptot[j]] * component.rho) * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for nn in range(len(self.solidchemicals)):
                                    if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                        M[j,n] = M[j,n] + OL[self.ptot[j]- self.toppoints + 1,nn]* FisL[component.name][self.ptot[j]] * component.rho* (self.z[i+1]-self.z[i-1])/2*self.flux_factor

                            M[j, n] = M[j, n] + (self.sorptions[component.name][chemical.name].get_q(C[self.ptot[j+1]- self.toppoints + 1,n], self.Cmax[chemical.name]) * Fis[component.name][self.ptot[j+1]] * component.rho) * (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for nn in range(len(self.solidchemicals)):
                                    if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                        M[j,n] = M[j,n] + O[self.ptot[j+1]- self.toppoints + 1,nn]* Fis[component.name][self.ptot[j+1]] * component.rho*(self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor

                        M[j, n] = (M[j, n] + CL[self.ptot[j]- self.toppoints + 1,n] * self.esL_plus_1[self.ptot[j]]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                                           + C[self.ptot[j+1]- self.toppoints + 1,n] * self.es_plus_1[self.ptot[j+1]]*(1+layer.doc/(10**6)*10**chemical.Kdoc)* (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor)

            else:
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for j in range(len(self.ptot) - 1):
                        layer  = self.layertot[j]
                        M[j+1, n] = 0
                        for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                            if i == self.pbio:
                                for component in self.components:
                                    M[j+1, n] = M[j+1, n] + self.sorptions[component.name][chemical.name].get_q(CL[i,n], self.Cmax[chemical.name]) * FisL[component.name][i] * component.rho * (self.z[i+1]-self.z[i])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j+1,n] = M[j+1,n] + OL[i,nn]* FisL[component.name][i] * component.rho* (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                                for component in self.components:
                                    M[j+1, n] = M[j+1, n] + self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho * (self.z[i]-self.z[i-1])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j+1,n] = M[j+1,n] + O[i,nn]* Fis[component.name][i] * component.rho*(self.z[i]-self.z[i-1])/2*self.flux_factor
                                M[j+1, n] = (M[j+1, n] + CL[i,n] * self.esL_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i+1]-self.z[i])  /2*self.flux_factor
                                                       + C[i,n]  * self.es_plus_1[i] *(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i]  -self.z[i-1])/2*self.flux_factor)
                            else:
                                for component in self.components:
                                    M[j+1, n] = M[j+1, n] + (self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho) * (self.z[i+1]-self.z[i-1])/2*self.flux_factor
                                    if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                        for nn in range(len(self.solidchemicals)):
                                            if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                                M[j+1,n] = M[j+1,n] + O[i,nn]* Fis[component.name][i] * component.rho * (self.z[i+1]-self.z[i-1])/2*self.flux_factor
                                M[j+1, n] = M[j+1, n] + C[i,n] * self.es_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc)*(self.z[i+1]-self.z[i-1])/2*self.flux_factor

                        for component in self.components:
                            M[j+1, n] = M[j+1, n] + (self.sorptions[component.name][chemical.name].get_q(CL[self.ptot[j]- self.toppoints + 1,n], self.Cmax[chemical.name]) * FisL[component.name][self.ptot[j]] * component.rho) * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for nn in range(len(self.solidchemicals)):
                                    if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                        M[j+1,n] = M[j+1,n] + OL[self.ptot[j]- self.toppoints + 1, nn]* FisL[component.name][self.ptot[j]] * component.rho* (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor

                        for component in self.components:
                            M[j+1, n] = M[j+1, n] + (self.sorptions[component.name][chemical.name].get_q(C[self.ptot[j+1]- self.toppoints + 1,n], self.Cmax[chemical.name]) * Fis[component.name][self.ptot[j+1]] * component.rho) * (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor
                            if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                for nn in range(len(self.solidchemicals)):
                                    if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                        M[j+1,n] = M[j+1,n] + O[self.ptot[j+1]- self.toppoints + 1,nn]* Fis[component.name][self.ptot[j+1]] * component.rho*(self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor


                        M[j+1, n] = (M[j+1, n] + CL[self.ptot[j]- self.toppoints + 1,n] * self.esL_plus_1[self.ptot[j]]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                                                + C[self.ptot[j+1]- self.toppoints + 1,n] * self.es_plus_1[self.ptot[j+1]]*(1+layer.doc/(10**6)*10**chemical.Kdoc)* (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor)

        else:
            M = zeros((len(self.layertot), self.nchemicals))
            for n in range (self.nchemicals):
                chemical = self.chemicals[n]
                for j in range(len(self.ptot) - 1):
                    layer  = self.layertot[j]
                    M[j, n] = 0
                    for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                        if self.bio == 1 and i == self.pbio:
                            for component in self.components:
                                M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(CL[i,n], self.Cmax[chemical.name]) * FisL[component.name][i] * component.rho * (self.z[i+1]-self.z[i])/2*self.flux_factor
                                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                    for nn in range(len(self.solidchemicals)):
                                        if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                            M[j,n] = M[j,n] + OL[i,nn]* FisL[component.name][i] * component.rho* (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                            for component in self.components:
                                M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho * (self.z[i]-self.z[i-1])/2*self.flux_factor
                                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                    for nn in range(len(self.solidchemicals)):
                                        if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                            M[j,n] = M[j,n] + O[i,nn]* Fis[component.name][i] * component.rho*(self.z[i]-self.z[i-1])/2*self.flux_factor
                            M[j, n] = (M[j, n] + CL[i,n] * self.esL_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i+1]-self.z[i])  /2*self.flux_factor
                                               + C[i,n]  * self.es_plus_1[i] *(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[i]  -self.z[i-1])/2*self.flux_factor)
                        else:
                            for component in self.components:
                                M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(C[i,n], self.Cmax[chemical.name]) * Fis[component.name][i] * component.rho*(self.z[i+1]-self.z[i-1])/2*self.flux_factor
                                if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                                    for nn in range(len(self.solidchemicals)):
                                        if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                            M[j,n] = M[j,n] + O[i,nn]* Fis[component.name][i] * component.rho*(self.z[i+1]-self.z[i-1])/2*self.flux_factor
                            M[j, n] = M[j, n] + C[i,n] * self.es_plus_1[i]*(1+layer.doc/(10**6)*10**chemical.Kdoc)*(self.z[i+1]-self.z[i-1])/2*self.flux_factor

                    for component in self.components:
                        M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(CL[self.ptot[j],n], self.Cmax[chemical.name]) * FisL[component.name][self.ptot[j]] * component.rho * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                        if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            for nn in range(len(self.solidchemicals)):
                                if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                    M[j,n] = M[j,n] + OL[self.ptot[j],nn]* FisL[component.name][self.ptot[j]] * component.rho* (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor

                    for component in self.components:
                        M[j, n] = M[j, n] + self.sorptions[component.name][chemical.name].get_q(C[self.ptot[j+1],n], self.Cmax[chemical.name]) * Fis[component.name][self.ptot[j+1]] * component.rho * (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor
                        if self.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            for nn in range(len(self.solidchemicals)):
                                if self.solidchemicals[nn].chemical_name == chemical.name and self.solidchemicals[nn].component_name == component.name:
                                    M[j,n] = M[j,n] + O[self.ptot[j+1],nn]* Fis[component.name][self.ptot[j+1]] * component.rho*(self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor

                    M[j, n] = (M[j, n] + CL[self.ptot[j],n] * self.esL_plus_1[self.ptot[j]]*(1+layer.doc/(10**6)*10**chemical.Kdoc) * (self.z[self.ptot[j]+1]-self.z[self.ptot[j]])/2*self.flux_factor
                                       + C[self.ptot[j+1],n] * self.es_plus_1[self.ptot[j+1]]*(1+layer.doc/(10**6)*10**chemical.Kdoc)* (self.z[self.ptot[j+1]]-self.z[self.ptot[j+1]-1])/2*self.flux_factor)

        return M


    def get_RM(self, RM, CL, CL_plus_1):
        """Calculates the total (aqueous + solid-phase) concentration "W" in
        the domain "z" using the pore water concentration at each grid point
        "C." """

        if len(self.layertot) > len(RM):
            RM_plus_1 = zeros((len(self.layertot), self.nchemicals))
            for i in range(len(RM)):
                for n in range(self.nchemicals):
                    RM_plus_1[i+1,n] = RM_plus_1[i,n]

        RM_plus_1 = RM.copy()

        elams           = self.elams.copy()
        rates           = self.rates.copy()
        rates_plus_1    = self.rates_plus_1.copy()

        for j in range(len(self.ptotbio) - 1):
            layer  = self.layertotbio[j]
            num    = self.layertotbio[j].nchemicals
            for nn in range(self.nchemicals, num):
                solidchemical = layer.chemicals[nn]
                n = solidchemical.kineticchemnum
                nnn = layer.lowersolidchemicalnums[nn - self.nchemicals]
                for ii in range(self.ptotbio[j+1] - self.ptotbio[j]-1):
                    i = ii + 1
                    elams [self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+nn] = 0
                    elams [self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+n] = elams [self.cptotLbio[j]+i*num+n, self.cptotLbio[j]+i*num+n] + elams [self.cptotLbio[j]+i*num+nn, self.cptotLbio[j]+i*num+n]
                    rates_plus_1[self.cptotLbio[j]+i*num+n] = rates_plus_1[self.cptotLbio[j]+i*num+n] + rates_plus_1[self.cptotLbio[j]+i*num+nn]
                    rates[self.cptotLbio[j]+i*num+n]        = rates[self.cptotLbio[j]+i*num+n]        + rates[self.cptotLbio[j]+i*num+nn]
                elams [self.cptotbio[j]+n, self.cptotbio[j]+nn]    = 0
                elams [self.cptotbio[j]+n, self.cptotbio[j]+n]     = elams [self.cptotbio[j]+n, self.cptotbio[j]+n] + elams [self.cptotbio[j]+nn, self.cptotbio[j]+n]
                rates_plus_1[self.cptotbio[j]+n]                   = rates_plus_1[self.cptotbio[j]+n]               + rates_plus_1[self.cptotbio[j]+nn]
                rates[self.cptotbio[j]+n]                          = rates[self.cptotbio[j]+n]                      + rates[self.cptotbio[j]+nn]
                if j == len(self.ptotbio) - 2:
                    elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+nnn]   = 0
                    elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+n] = elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+n] + elams [self.cptotbio[j+1]+nnn, self.cptotbio[j+1]+n]
                    rates[self.cptotbio[j+1]+n]                        = rates[self.cptotbio[j+1]+n]          + rates[self.cptotbio[j+1]+nnn]
                    rates_plus_1[self.cptotbio[j+1]+n]                 = rates_plus_1[self.cptotbio[j+1]+n]   + rates_plus_1[self.cptotbio[j+1]+nnn]
                elif nnn >= self.layertotbio[j+1].nchemicals:
                    elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+nnn]   = 0
                    elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+n] = elams [self.cptotbio[j+1]+n, self.cptotbio[j+1]+n] + elams [self.cptotbio[j+1]+nnn, self.cptotbio[j+1]+n]
                    rates[self.cptotbio[j+1]+n]                        = rates[self.cptotbio[j+1]+n]          + rates[self.cptotbio[j+1]+nnn]
                    rates_plus_1[self.cptotbio[j+1]+n]                 = rates_plus_1[self.cptotbio[j+1]+n]   + rates_plus_1[self.cptotbio[j+1]+nnn]

        if len(self.layertot) < len(self.layers+self.deplayer):
            elam_matrix_plus_1 = elams * transpose(matrix(CL_plus_1))
            elam_array_plus_1 = array([elam_matrix_plus_1[n,0] for n in range(len(elam_matrix_plus_1))])
            reacs_imp  = (elam_array_plus_1 + rates_plus_1) * self.delt
            if self.timeoption == 'Implicit':
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for j in range(len(self.ptot) - 1):
                        layer  = self.layertot[j]
                        for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                            ii = i - self.ptot[j]
                            RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_imp[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                        RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_imp[self.cptot[j] + n]   * self.flux_factor * chemical.MW
                        RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_imp[self.cptot[j+1] + n] * self.flux_factor * chemical.MW
            else:
                elam_matrix = elams * transpose(matrix(CL))
                elam_array  = array([elam_matrix[n,0] for n in range(len(elam_matrix))])
                reacs_CN  = (elam_array  + elam_array_plus_1 + rates + rates_plus_1)/2 * self.delt
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    for j in range(len(self.ptot) - 1):
                        layer  = self.layertot[j]
                        for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                            ii = i - self.ptot[j]
                            RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_CN[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                        RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_CN[self.cptot[j] + n] * self.flux_factor * chemical.MW
                        RM_plus_1[j+1, n] = RM_plus_1[j+1, n] + reacs_CN[self.cptot[j+1] + n] * self.flux_factor * chemical.MW

        else:
            elam_matrix_plus_1 = elams * transpose(matrix(CL_plus_1))
            elam_array_plus_1 = array([elam_matrix_plus_1[n,0] for n in range(len(elam_matrix_plus_1))])
            reacs_imp  = (elam_array_plus_1 + rates_plus_1) * self.delt
            if self.timeoption == 'Implicit':
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    if self.dep == 1 and self.toppoints > 1:
                        for jj in range(len(self.ptot) - 2):
                            j = jj + 1
                            layer  = self.layertot[j]
                            for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                                ii = i - self.ptot[j]
                                RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_imp[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                            RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_imp[self.cptot[j] + n]   * self.flux_factor * chemical.MW
                            RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_imp[self.cptot[j+1] + n] * self.flux_factor * chemical.MW
                    else:
                        for j in range(len(self.ptot) - 1):
                            layer  = self.layertot[j]
                            for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                                ii = i - self.ptot[j]
                                RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_imp[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                            RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_imp[self.cptot[j] + n]   * self.flux_factor * chemical.MW
                            RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_imp[self.cptot[j+1] + n] * self.flux_factor * chemical.MW
            else:
                elam_matrix = elams * transpose(matrix(CL))
                elam_array  = array([elam_matrix[n,0] for n in range(len(elam_matrix))])
                reacs_CN  = (elam_array  + elam_array_plus_1 + rates + rates_plus_1)/2 * self.delt
                for n in range (self.nchemicals):
                    chemical = self.chemicals[n]
                    if self.dep == 1 and self.toppoints > 1:
                        for jj in range(len(self.ptot) - 2):
                            j = jj + 1
                            layer  = self.layertot[j]
                            for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                                ii = i - self.ptot[j]
                                RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_CN[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                            RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_CN[self.cptot[j] + n]   * self.flux_factor * chemical.MW
                            RM_plus_1[jj, n] = RM_plus_1[jj, n] + reacs_CN[self.cptot[j+1] + n] * self.flux_factor * chemical.MW
                    else:
                        for j in range(len(self.ptot) - 1):
                            layer  = self.layertot[j]
                            for i in range(self.ptot[j] + 1, self.ptot[j+1]):
                                ii = i - self.ptot[j]
                                RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_CN[self.cptot[j]+ii*layer.nchemicals + n] * (self.z[i+1]-self.z[i-1])/2*self.flux_factor*chemical.MW
                            RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_CN[self.cptot[j] + n]   * self.flux_factor * chemical.MW
                            RM_plus_1[j, n] = RM_plus_1[j, n] + reacs_CN[self.cptot[j+1] + n] * self.flux_factor * chemical.MW

        return RM_plus_1

class Output:
    """Stores the output information from a simulation."""

    def __init__(self, parameters):
        """Constructor method.  Essentially this just lumps the output data
        together for portability.  Variable definitions:

        C     -- Porewater concentrations at 250 times and all grid points
        F     -- Fluxes at 250 times and all grid points
        q     -- Solid concentrations at 250 times and all grid points
        W     -- Total concentrations at 250 times and all grid points
        z     -- The grid
        t     -- Times for the profiles for the other variables
        tCFqw -- A list of tuples containing the flux, pore water concentration,
                 solid and total concentration at a given time
                 at the depth of interest
        tplot -- A list of the times to get for plotting
        Cplot -- The porewater concentrations at the times in tplot
        Fplot -- The fluxes at the times in tplot
        qplot -- The solid concentrations at the times in tplot
        Wplot -- The total concentrations at the times in tplot
        """
        self.p                  = parameters.p
        self.z                  = parameters.zdep + parameters.z
        self.deplayer           = parameters.deplayer
        self.layers             = parameters.layers
        self.layertot           = parameters.layertot
        self.chemicals          = parameters.chemicals
        self.nchemicals         = parameters.nchemicals
        self.outputsteps        = parameters.outputsteps

        self.solidchemicals     = parameters.solidchemicals
        self.nsolidchemicals    = parameters.nsolidchemicals
        self.solidchemical_list = parameters.solidchemical_list

        self.times      = [parameters.tstart + round(i * (parameters.tfinal_ori-parameters.tstart)/self.outputsteps, 10) for i in range(self.outputsteps+1)]
        self.sizeflag = 0

        try:
            self.O   = zeros((len(self.times), len(self.z), parameters.nsolidchemicals))
            self.OL  = zeros((len(self.times), len(self.z), parameters.nsolidchemicals))
            self.C   = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.CL  = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.F   = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.q   = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.qL  = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.W   = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.WL  = zeros((len(self.times), len(self.z), parameters.nchemicals))
            self.qm  = zeros((len(self.times), len(self.z), parameters.nchemicals, len(parameters.components)))
            self.qmL = zeros((len(self.times), len(self.z), parameters.nchemicals, len(parameters.components)))
            self.Cw  = zeros((len(self.times), parameters.nchemicals))
            self.Fi  = zeros((len(self.times), len(self.z), parameters.ncomponents))
            self.FiL = zeros((len(self.times), len(self.z), parameters.ncomponents))
            self.FM  = zeros((len(self.times), 2, parameters.nchemicals))
            self.DM  = zeros((len(self.times), parameters.nchemicals))
            self.M   = zeros((len(self.times), len(parameters.layers+parameters.deplayer), parameters.nchemicals))
            self.RM  = zeros((len(self.times), len(parameters.layers+parameters.deplayer), parameters.nchemicals))

        except :    self.sizeflag = 1
        self.n     = 0     #index of next time to collect data

        self.parameters = parameters

    def copy(self):

        output = Output(self.parameters)
        output.z = self.z

        return output


    def converter (self, parameters, Cn, CnL, Fis, FisL, FM, RM, DM):

        """Stores the concentrations and fluxes from a given time."""
        self.ptotbio       = parameters.ptotbio
        self.cptotbio      = parameters.cptotbio
        self.cptotLbio     = parameters.cptotLbio
        self.layertotbio   = parameters.layertotbio

        z_t             = parameters.z

        if parameters.dep == 1 and parameters.toppoints > 1:


            toppoints           = parameters.toppoints
            pdep                = len(self.z) - (len(parameters.z) - toppoints + 1)
            CnL_2D              = zeros((len(self.z), parameters.nchemicals))
            Solid_CnL_2D        = zeros((len(self.z), self.nsolidchemicals))
            Fis_2D              = zeros((len(self.z), parameters.ncomponents))

            F_2D                = zeros((len(self.z), parameters.nchemicals))
            q_2D                = zeros((len(self.z), parameters.nchemicals))
            qL_2D               = zeros((len(self.z), parameters.nchemicals))
            qm_2D               = zeros((len(self.z), parameters.nchemicals, len(parameters.components)))
            qmL_2D              = zeros((len(self.z), parameters.nchemicals, len(parameters.components)))
            W_2D                = zeros((len(self.z), parameters.nchemicals))
            WL_2D               = zeros((len(self.z), parameters.nchemicals))

            Cn_top_2D           = zeros((toppoints + 1, parameters.nchemicals))
            Solid_Cn_top_2D     = zeros((toppoints + 1, self.nsolidchemicals))

            for jj in range(len(self.layertotbio) - 1):
                j = jj + 1
                num = self.layertotbio[j].nchemicals
                for i in range(1, self.ptotbio[j+1] - self.ptotbio[j]):
                    for n in range(self.nchemicals):
                        CnL_2D[pdep + self.ptotbio[j] - toppoints + 1 + i,n] = CnL[self.cptotLbio[j] + i*num + n] * self.chemicals[n].MW
                for n in range(self.nchemicals):
                    CnL_2D[pdep + self.ptotbio[j] - toppoints + 1, n]        = CnL[self.cptotbio[j] + n] * self.chemicals[n].MW
            for n in range(self.nchemicals):
                CnL_2D[pdep + self.ptotbio[-1] - toppoints + 1,n]        = CnL[self.cptotbio[-1] + n] * self.chemicals[n].MW
                CnL_2D[pdep,n]                                           = CnL[n] * self.chemicals[n].MW

            for n in range(self.nchemicals):
                for i in range(toppoints + 1):
                    Cn_top_2D[i,n]        = CnL[self.layertotbio[0].nchemicals*i+n] * self.chemicals[n].MW

            Cn_2D         = CnL_2D.copy()

            for jj in range(len(self.layertotbio) - 1):
                j = jj + 1
                for n in range(self.nchemicals):
                    Cn_2D[pdep + self.ptotbio[j] - toppoints + 1, n]  = Cn[self.cptotbio[j] + n] * self.chemicals[n].MW


            for jj in range(len(self.layertotbio) - 1):
                j = jj + 1
                for i in range((self.ptotbio[j+1]-self.ptotbio[j])):
                    for m in range(parameters.ncomponents):
                        Fis_2D[pdep + self.ptotbio[j] - toppoints + 1 + i,m]  = Fis[parameters.components[m].name][self.ptotbio[j] + i]

                for m in range(parameters.ncomponents):
                    Fis_2D[pdep + self.ptotbio[-1]- toppoints + 1,m]         = Fis[parameters.components[m].name][self.ptotbio[-1]]
                    Fis_2D[pdep,m]                                           = Fis[parameters.components[m].name][0]

            FisL_2D = Fis_2D.copy()
            for jj in range(len(self.layertotbio) - 1):
                j = jj + 1
                for m in range(parameters.ncomponents):
                    FisL_2D[pdep + self.ptotbio[j] - toppoints + 1,m]  = FisL[parameters.components[m].name][self.ptotbio[j]]

            if self.nsolidchemicals > 0:
                for jj in range(len(self.layertotbio) - 1):
                    j = jj + 1
                    num = self.layertotbio[j].nsolidchemicals
                    for i in range(1, self.ptotbio[j+1]-self.ptotbio[j]):
                        for n in range(num):
                            nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                            Solid_CnL_2D[pdep + self.ptotbio[j] - toppoints + 1 + i,nn] = CnL[self.cptotLbio[j] + i*(self.nchemicals + num) + self.nchemicals + n] * self.solidchemicals[nn].MW

                    for n in range(num):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                        Solid_CnL_2D[pdep + self.ptotbio[j] - toppoints + 1,nn] = CnL[self.cptotbio[j] + self.nchemicals + n] * self.solidchemicals[nn].MW

                for n in range(self.layertotbio[-1].nsolidchemicals):
                    nn = parameters.solidchemical_list.index(self.layertotbio[-1].solidchemicals[n].name)
                    Solid_CnL_2D[pdep + self.ptotbio[-1] - toppoints + 1, nn] = CnL[self.cptotbio[-1] + self.nchemicals + n] * self.solidchemicals[nn].MW
                    Solid_CnL_2D[pdep,nn]                                     = CnL[self.nchemicals + n] * self.solidchemicals[nn].MW

                for n in range(self.nsolidchemicals):
                    for i in range(toppoints + 1):
                        Solid_Cn_top_2D[i,n]           = Cn[self.layertotbio[0].nchemicals*i+self.nchemicals+n] * self.chemicals[n].MW

                Solid_Cn_2D         = Solid_CnL_2D.copy()

                for jj in range(len(self.layertotbio)-1):
                    j = jj + 1
                    for n in range(self.layertotbio[j].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                        Solid_Cn_2D[pdep + self.ptotbio[j] - toppoints + 1,nn] = Cn[self.cptotbio[j] + n] * self.solidchemicals[nn].MW

                    for n in range(self.layertotbio[j-1].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j-1].solidchemicals[n].name)
                        na = self.layertotbio[j-1].lowersolidchemicalnums[n]
                        Solid_Cn_2D[pdep + self.ptotbio[j] - toppoints + 1,nn] = Cn[self.cptotbio[j] + na] * self.solidchemicals[nn].MW

                F_2D[pdep:, :]                          = parameters.get_fluxes(Cn_2D[pdep:, :], CnL_2D[pdep:, :], O = Solid_Cn_2D[pdep:, :], OL = Solid_CnL_2D[pdep:, :], Cn_top = Cn_top_2D, O_top = Solid_Cn_top_2D )
                q_2D[pdep:, :], qL_2D[pdep:, :]         = parameters.get_qs(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
                W_2D[pdep:, :]                          = parameters.get_Ws(Cn_2D[pdep:, :], Solid_Cn_2D[pdep:, :])
                WL_2D[pdep:, :]                         = parameters.get_Ws(CnL_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
                qm_2D[pdep:, :, :], qmL_2D[pdep:, :, :] = parameters.get_qms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
                M_2D                                    = parameters.get_Ms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
            else:
                F_2D[pdep:, :]                          = parameters.get_fluxes(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Cn_top = Cn_top_2D)
                q_2D[pdep:, :],qL_2D[pdep:, :]          = parameters.get_qs(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL)
                W_2D[pdep:, :]                          = parameters.get_Ws(Cn_2D[pdep:, :])
                WL_2D[pdep:, :]                         = parameters.get_Ws(CnL_2D[pdep:, :])
                qm_2D[pdep:, :, :], qmL_2D[pdep:, :, :] = parameters.get_qms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL)
                M_2D                                    = parameters.get_Ms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL)

            if parameters.topBCtype == 'Finite mixed water column' or parameters.topBCtype == 'Mass transfer':
                Cw, Cw_plus_1= parameters.get_Cws()
            else:
                Cw = zeros(parameters.nchemicals)

        else:

            pdep            = len(self.z) - len(parameters.z)


            CnL_2D          = zeros((len(self.z), parameters.nchemicals))
            Solid_CnL_2D    = zeros((len(self.z), self.nsolidchemicals))
            F_2D            = zeros((len(self.z), parameters.nchemicals))
            q_2D            = zeros((len(self.z), parameters.nchemicals))
            qL_2D           = zeros((len(self.z), parameters.nchemicals))
            qm_2D           = zeros((len(self.z), parameters.nchemicals, len(parameters.components)))
            qmL_2D          = zeros((len(self.z), parameters.nchemicals, len(parameters.components)))
            W_2D            = zeros((len(self.z), parameters.nchemicals))
            WL_2D           = zeros((len(self.z), parameters.nchemicals))

            for j in range(len(self.layertotbio)):
                num = self.layertotbio[j].nchemicals
                for ii in range((self.ptotbio[j+1] - self.ptotbio[j] - 1)):
                    i = ii + 1
                    for n in range(self.nchemicals):
                        CnL_2D[pdep + self.ptotbio[j] + i,n] = CnL[self.cptotLbio[j] + i*num + n] * self.chemicals[n].MW
                for n in range(self.nchemicals):
                    CnL_2D[pdep + self.ptotbio[j], n]            = CnL[self.cptotbio[j] + n]          * self.chemicals[n].MW

            for n in range(self.nchemicals):
                CnL_2D[pdep + self.ptotbio[0],n]             = CnL[self.cptotbio[0] + n] * self.chemicals[n].MW
                CnL_2D[pdep + self.ptotbio[-1],n]            = CnL[self.cptotbio[-1] + n] * self.chemicals[n].MW

            Cn_2D         = CnL_2D.copy()
            for jj in range(len(self.layertotbio)-1):
                j = jj+1
                for n in range(self.nchemicals):
                    Cn_2D[pdep + self.ptotbio[j], n]  = Cn[self.cptotbio[j] + n] * self.chemicals[n].MW

            Fis_2D           = zeros((len(self.z), parameters.ncomponents))
            for j in range(len(self.layertotbio)):
                for i in range((self.ptotbio[j+1]-self.ptotbio[j])):
                    for m in range(parameters.ncomponents):
                        Fis_2D[pdep + self.ptotbio[j] + i,m]        = Fis[parameters.components[m].name][self.ptotbio[j] + i]
                for m in range(parameters.ncomponents):
                    Fis_2D[pdep + self.ptotbio[-1],m]               = Fis[parameters.components[m].name][self.ptotbio[-1]]

            FisL_2D = Fis_2D.copy()
            for j in range(len(self.layertotbio)):
                for m in range(parameters.ncomponents):
                    FisL_2D[pdep + self.ptotbio[j],m]  = FisL[parameters.components[m].name][self.ptotbio[j]]

            W_2D[pdep:, :]                          = parameters.get_Ws(Cn_2D[pdep:, :])
            WL_2D[pdep:, :]                         = parameters.get_Ws(CnL_2D[pdep:, :])

            if self.nsolidchemicals > 0:
                for j in range(len(self.layertotbio)-1):
                    num = self.layertotbio[j].nsolidchemicals
                    for i in range((self.ptotbio[j+1]-self.ptotbio[j])):
                        for n in range(num):
                            nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                            Solid_CnL_2D[pdep + self.ptotbio[j] + i,nn]        = CnL[self.cptotLbio[j] + i*(self.nchemicals + num) + self.nchemicals + n] * self.solidchemicals[nn].MW

                    for n in range(self.layertotbio[j+1].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j+1].solidchemicals[n].name)
                        Solid_CnL_2D[pdep + self.ptotbio[j+1],nn]        = CnL[self.cptotbio[j+1] + self.nchemicals + n] * self.solidchemicals[nn].MW

                    for n in range(self.layertotbio[j].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                        na = self.layertotbio[j].lowersolidchemicalnums[self.layertotbio[j].solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)]
                        Solid_CnL_2D[pdep + self.ptotbio[j+1],nn]        = CnL[self.cptotbio[j+1] + na] * self.solidchemicals[nn].MW

                for n in range(self.layertotbio[0].nsolidchemicals):
                    nn = parameters.solidchemical_list.index(self.layertotbio[0].solidchemicals[n].name)
                    Solid_CnL_2D[pdep + self.ptotbio[0],nn]        = CnL[self.nchemicals + n] * self.solidchemicals[nn].MW

                for n in range(self.layertotbio[-1].nsolidchemicals):
                    nn = parameters.solidchemical_list.index(self.layertotbio[-1].solidchemicals[n].name)
                    Solid_CnL_2D[pdep + self.ptotbio[-1],nn]        = CnL[self.cptotbio[-1] + self.nchemicals + n] * self.solidchemicals[nn].MW

                Solid_Cn_2D         = Solid_CnL_2D.copy()

                for jj in range(len(self.layertotbio)-1):
                    j = jj + 1
                    for n in range(self.layertotbio[j].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                        Solid_Cn_2D[pdep + self.ptotbio[j],nn] = Cn[self.cptotbio[j] + n] * self.solidchemicals[nn].MW

                    for n in range(self.layertotbio[j-1].nsolidchemicals):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j-1].solidchemicals[n].name)
                        na = self.layertotbio[j-1].lowersolidchemicalnums[n]
                        Solid_Cn_2D[pdep + self.ptotbio[j],nn] = Cn[self.cptotbio[j] + na] * self.solidchemicals[nn].MW

                F_2D[pdep:, :]                          = parameters.get_fluxes(Cn_2D[pdep:, :], CnL_2D[pdep:, :], O = Solid_Cn_2D[pdep:, :],  OL = Solid_CnL_2D[pdep:, :])
                q_2D[pdep:, :], qL_2D[pdep:, :]         = parameters.get_qs(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
                qm_2D[pdep:, :, :], qmL_2D[pdep:, :, :] = parameters.get_qms(Cn_2D[pdep:, :],CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])
                M_2D                                    = parameters.get_Ms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL, Solid_Cn_2D[pdep:, :], Solid_CnL_2D[pdep:, :])

            else:
                F_2D[pdep:, :]                          = parameters.get_fluxes(Cn_2D[pdep:, :], CnL_2D[pdep:, :])
                q_2D[pdep:, :], qL_2D[pdep:, :]         = parameters.get_qs(Cn_2D[pdep:, :],  CnL_2D[pdep:, :], Fis, FisL)
                qm_2D[pdep:, :, :], qmL_2D[pdep:, :, :] = parameters.get_qms(Cn_2D[pdep:, :], CnL_2D[pdep:, :], Fis, FisL)
                M_2D                                    = parameters.get_Ms(Cn_2D[pdep:, :],  CnL_2D[pdep:, :], Fis, FisL)

            if parameters.topBCtype == 'Finite mixed water column' or parameters.topBCtype == 'Mass transfer':
                Cw, Cw_plus_1 = parameters.get_Cws()
            else:
                Cw = zeros(parameters.nchemicals)

        return [Cn_2D, CnL_2D, Fis_2D, FisL_2D, F_2D, q_2D, qL_2D, qm_2D, qmL_2D, W_2D, WL_2D, Cw, M_2D, FM, RM, DM]

    def store(self, t, t_plus_1, parameters, variables, variables_plus_1):

        delt                  = t_plus_1 - t

        Cn_2D                 = variables[0]
        CnL_2D                = variables[1]
        Fis_2D                = variables[2]
        FisL_2D               = variables[3]
        F_2D                  = variables[4]
        q_2D                  = variables[5]
        qL_2D                 = variables[6]
        qm_2D                 = variables[7]
        qmL_2D                = variables[8]
        W_2D                  = variables[9]
        WL_2D                 = variables[10]
        Cw_2D                 = variables[11]
        M_2D                  = variables[12]
        FM                    = variables[13]
        RM                    = variables[14]
        DM                    = variables[15]

        Cn_plus_1_2D          = variables_plus_1[0]
        CnL_plus_1_2D         = variables_plus_1[1]
        Fis_plus_1_2D         = variables_plus_1[2]
        FisL_plus_1_2D        = variables_plus_1[3]
        F_plus_1_2D           = variables_plus_1[4]
        q_plus_1_2D           = variables_plus_1[5]
        qL_plus_1_2D          = variables_plus_1[6]
        qm_plus_1_2D          = variables_plus_1[7]
        qmL_plus_1_2D         = variables_plus_1[8]
        W_plus_1_2D           = variables_plus_1[9]
        WL_plus_1_2D          = variables_plus_1[10]
        Cw_plus_1_2D          = variables_plus_1[11]
        M_plus_1_2D           = variables_plus_1[12]
        FM_plus_1_2D          = variables_plus_1[13]
        RM_plus_1_2D          = variables_plus_1[14]
        DM_plus_1_2D          = variables_plus_1[15]

        while self.n < len(self.times) and round(self.times[self.n], 8) <= round(t_plus_1, 8):

            time = self.times[self.n]
            self.C[self.n, :, :]         = time_interpolate(time, t, delt, Cn_2D,   Cn_plus_1_2D)
            self.CL[self.n, :, :]        = time_interpolate(time, t, delt, CnL_2D,  CnL_plus_1_2D)
            self.Fi[self.n, :, :]        = time_interpolate(time, t, delt, Fis_2D,  Fis_plus_1_2D)
            self.FiL[self.n, :, :]       = time_interpolate(time, t, delt, FisL_2D, FisL_plus_1_2D)
            self.F[self.n, :, :]         = time_interpolate(time, t, delt, F_2D,    F_plus_1_2D)
            self.q[self.n, :, :]         = time_interpolate(time, t, delt, q_2D,    q_plus_1_2D)
            self.qL[self.n, :, :]        = time_interpolate(time, t, delt, qL_2D,   qL_plus_1_2D)
            self.qm[self.n, :, :, :]     = time_interpolate(time, t, delt, qm_2D,   qm_plus_1_2D)
            self.qmL[self.n, :, :, :]    = time_interpolate(time, t, delt, qmL_2D,  qmL_plus_1_2D)
            self.W[self.n, :, :]         = time_interpolate(time, t, delt, W_2D,    W_plus_1_2D)
            self.WL[self.n, :, :]        = time_interpolate(time, t, delt, WL_2D,   WL_plus_1_2D)
            self.FM[self.n, :, :]        = time_interpolate(time, t, delt, FM,      FM_plus_1_2D)
            self.DM[self.n, :]           = time_interpolate(time, t, delt, DM,      DM_plus_1_2D)
            self.M[self.n, :, :]         = time_interpolate(time, t, delt, M_2D,    M_plus_1_2D)
            self.RM[self.n, :, :]        = time_interpolate(time, t, delt, RM,      RM_plus_1_2D)

            if parameters.topBCtype == 'Finite mixed water column' or parameters.topBCtype == 'Mass transfer':
                for n in range(parameters.nchemicals):
                    self.Cw[self.n, :] = time_interpolate(self.times[self.n], t, delt, Cw_2D, Cw_plus_1_2D)

            self.n = self.n + 1



    def store_no_dep(self, Cn, CnL, CnL_plus_1, t, t_plus_1, parameters, Fis, Fis_plus_1, FisL, FisL_plus_1, FM, FM_plus_1, RM, RM_plus_1, DM, DM_plus_1):
        """Stores the concentrations and fluxes from a given time."""

        self.ptot           = parameters.ptot
        self.cptot          = parameters.cptot
        self.cptotL         = parameters.cptotL
        self.layertot       = parameters.layertot

        self.ptotbio        = parameters.ptotbio
        self.cptotbio       = parameters.cptotbio
        self.layertotbio    = parameters.layertotbio

        delt            = t_plus_1 - t
        pdep            = len(self.z) - len(parameters.z)

        z_t             = len(parameters.z)
        CnL_2D          = zeros((z_t, parameters.nchemicals))
        CnL_plus_1_2D   = zeros((z_t, parameters.nchemicals))

        for j in range(len(self.layertotbio)):
            num = self.layertotbio[j].nchemicals
            for ii in range((self.ptotbio[j+1] - self.ptotbio[j]-1)):
                i = ii + 1
                for n in range(self.nchemicals):
                    CnL_2D[self.ptotbio[j] + i,n]          = CnL[self.cptotLbio[j] + i*num + n] * self.chemicals[n].MW
                    CnL_plus_1_2D[self.ptotbio[j] + i,n]   = CnL_plus_1[self.cptotLbio[j] + i*num + n] * self.chemicals[n].MW
            for n in range(self.nchemicals):
                CnL_2D[self.ptotbio[j],n]                  = CnL[self.cptotbio[j] + n] * self.chemicals[n].MW
                CnL_plus_1_2D[self.ptotbio[j],n]           = CnL_plus_1[self.cptotbio[j] + n] * self.chemicals[n].MW

        for n in range(self.nchemicals):
            CnL_2D[self.ptotbio[-1],n]                 = CnL[self.cptotbio[-1] + n] * self.chemicals[n].MW
            CnL_plus_1_2D[self.ptotbio[-1],n]          = CnL_plus_1[self.cptotbio[-1] + n] * self.chemicals[n].MW

        Cn_2D         = CnL_2D.copy()
        Cn_plus_1_2D  = CnL_plus_1_2D.copy()

        for j in range(len(self.layertotbio)-1):
            for n in range(self.nchemicals):
                Cn_2D[self.ptotbio[j+1],n]  = Cn[self.cptotbio[j+1] + n] * self.chemicals[n].MW

        if self.nsolidchemicals > 0:
            Solid_CnL_2D        = zeros((z_t, self.nsolidchemicals))
            Solid_CnL_plus_1_2D = zeros((z_t, self.nsolidchemicals))

            for j in range(len(self.layertotbio)):
                num = self.layertotbio[j].nsolidchemicals
                for ii in range((self.ptotbio[j+1]-self.ptotbio[j]-1)):
                    i = ii + 1
                    for n in range(num):
                        nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                        Solid_CnL_2D[self.ptotbio[j] + i,nn]        = CnL[self.cptotLbio[j] + i*(self.nchemicals + num) + self.nchemicals + n] * self.solidchemicals[nn].MW
                        Solid_CnL_plus_1_2D[self.ptotbio[j] + i,nn] = CnL_plus_1[self.cptotLbio[j] + i*(self.nchemicals + num) + self.nchemicals + n] * self.solidchemicals[nn].MW

                for n in range(self.layertotbio[j].nsolidchemicals):
                    nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                    Solid_CnL_2D[self.ptotbio[j],nn]        = CnL[self.cptotbio[j] + self.nchemicals + n] * self.solidchemicals[nn].MW
                    Solid_CnL_plus_1_2D[self.ptotbio[j],nn] = CnL_plus_1[self.cptotbio[j]  + self.nchemicals + n] * self.solidchemicals[nn].MW

            for n in range(self.layertotbio[-1].nsolidchemicals):
                nn = parameters.solidchemical_list.index(self.layertotbio[-1].solidchemicals[n].name)
                Solid_CnL_2D[self.ptotbio[-1],nn]        = CnL[self.cptotbio[-1] + self.nchemicals + n] * self.solidchemicals[nn].MW
                Solid_CnL_plus_1_2D[self.ptotbio[-1],nn] = CnL_plus_1[self.cptotbio[-1]  + self.nchemicals + n] * self.solidchemicals[nn].MW

            Solid_Cn_2D         = Solid_CnL_2D.copy()
            Solid_Cn_plus_1_2D  = Solid_CnL_plus_1_2D.copy()

            for j in range(len(self.layertotbio)-1):
                num = self.layertotbio[j+1].nsolidchemicals
                for n in range(self.layertotbio[j+1].nsolidchemicals):
                    if self.layertotbio[j].solidchemical_list.count(self.layertotbio[j+1].solidchemicals[n].name)> 0:
                        nn = parameters.solidchemical_list.index(self.layertotbio[j+1].solidchemicals[n].name)
                        Solid_Cn_2D[self.ptotbio[j+1],nn]           = Cn[self.cptotbio[j+1] + n] * self.solidchemicals[nn].MW

            for j in range(len(self.layertotbio)-1):
                num = self.layertotbio[j].nsolidchemicals
                for n in range(self.layertotbio[j].nsolidchemicals):
                    nn = parameters.solidchemical_list.index(self.layertotbio[j].solidchemicals[n].name)
                    na = self.layertotbio[j].lowersolidchemicalnums[n]
                    Solid_Cn_2D[self.ptotbio[j+1],nn]           = Cn[self.cptotbio[j+1] + na] * self.solidchemicals[nn].MW
                    Solid_Cn_plus_1_2D[self.ptotbio[j+1],nn]    = CnL_plus_1[self.cptotbio[j+1] + na] * self.solidchemicals[nn].MW

        Fis_2D           = zeros((z_t, parameters.ncomponents))
        Fis_plus_1_2D    = zeros((z_t, parameters.ncomponents))
        for j in range(len(self.layertotbio)):
            for i in range((self.ptotbio[j+1]-self.ptotbio[j])):
                for m in range(parameters.ncomponents):
                    Fis_2D[self.ptotbio[j] + i,m]        = Fis[parameters.components[m].name][self.ptotbio[j] + i]
                    Fis_plus_1_2D[self.ptotbio[j] + i,m] = Fis_plus_1[parameters.components[m].name][self.ptotbio[j] + i]
            for m in range(parameters.ncomponents):
                Fis_2D[self.ptotbio[-1],m]               = Fis[parameters.components[m].name][self.ptotbio[-1]]
                Fis_plus_1_2D[self.ptotbio[-1],m]        = Fis_plus_1[parameters.components[m].name][self.ptotbio[-1]]


        FisL_2D         = Fis_2D.copy()
        FisL_plus_1_2D  = Fis_plus_1_2D.copy()

        for j in range(len(self.layertotbio)-1):
            i = self.ptotbio[j+1]
            for m in range(parameters.ncomponents):
                FisL_2D[i,m]        = FisL[parameters.components[m].name][i]
                FisL_plus_1_2D[i,m] = FisL_plus_1[parameters.components[m].name][i]

        while self.n < len(self.times) and round(self.times[self.n], 8) <= round(t_plus_1, 8):
            self.C[self.n, pdep:, :]            = time_interpolate(self.times[self.n], t, delt, Cn_2D,    Cn_plus_1_2D)
            self.CL[self.n, pdep:, :]           = time_interpolate(self.times[self.n], t, delt, CnL_2D,   CnL_plus_1_2D)
            self.Fi[self.n, pdep:, :]           = time_interpolate(self.times[self.n], t, delt, Fis_2D,   Fis_plus_1_2D)
            self.FiL[self.n, pdep:, :]          = time_interpolate(self.times[self.n], t, delt, FisL_2D,  FisL_plus_1_2D)
            self.W[self.n, pdep:, :]            = parameters.get_Ws( self.C[self.n, pdep:])
            self.WL[self.n, pdep:, :]           = parameters.get_Ws( self.CL[self.n, pdep:])
            self.RM[self.n, :, :]           = time_interpolate(self.times[self.n], t, delt, RM, RM_plus_1)
            self.FM[self.n, :, :]           = time_interpolate(self.times[self.n], t, delt, FM, FM_plus_1)
            self.DM[self.n, :]              = time_interpolate(self.times[self.n], t, delt, DM, DM_plus_1)
            if self.nsolidchemicals > 0:
                self.O[self.n, pdep:, :]        = time_interpolate(self.times[self.n], t, delt, Solid_Cn_2D,  Solid_Cn_plus_1_2D)
                self.OL[self.n, pdep:, :]       = time_interpolate(self.times[self.n], t, delt, Solid_CnL_2D, Solid_CnL_plus_1_2D)
                self.F[self.n, pdep:, :]        = time_interpolate(self.times[self.n], t, delt, parameters.get_fluxes(Cn_2D, CnL_2D, O = Solid_Cn_2D, OL = Solid_CnL_2D, flag = 1), parameters.get_fluxes(Cn_plus_1_2D, CnL_plus_1_2D, O = Solid_CnL_plus_1_2D, OL = Solid_CnL_plus_1_2D))
                q,qL                            = parameters.get_qs(Cn_2D, CnL_2D, Fis, FisL, Solid_Cn_2D, Solid_CnL_2D)
                q_plus_1, qL_plus_1             = parameters.get_qs(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1, Solid_Cn_plus_1_2D, Solid_CnL_plus_1_2D)
                self.q[self.n, pdep:, :]        = time_interpolate(self.times[self.n], t, delt, q,  q_plus_1)
                self.qL[self.n, pdep:, :]       = time_interpolate(self.times[self.n], t, delt, qL, qL_plus_1)
                self.M[self.n, :, :]            = time_interpolate(self.times[self.n], t, delt, parameters.get_Ms(Cn_2D, CnL_2D, Fis, FisL, O = Solid_Cn_2D, OL = Solid_CnL_2D), parameters.get_Ms(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1, O = Solid_CnL_2D, OL = Solid_CnL_plus_1_2D))

                qm, qmL                         = parameters.get_qms(Cn_2D, CnL_2D, Fis, FisL, Solid_Cn_2D, Solid_CnL_2D)
                qm_plus_1, qmL_plus_1           = parameters.get_qms(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1, Solid_Cn_plus_1_2D, Solid_CnL_plus_1_2D)
                self.qm[self.n, pdep:, :, :]    = time_interpolate(self.times[self.n], t, delt, qm,  qm_plus_1)
                self.qmL[self.n, pdep:, :, :]   = time_interpolate(self.times[self.n], t, delt, qmL, qmL_plus_1)

            else:
                self.F[self.n, pdep:, :]        = time_interpolate(self.times[self.n], t, delt, parameters.get_fluxes(Cn_2D, CnL_2D, flag = 1), parameters.get_fluxes(Cn_plus_1_2D, CnL_plus_1_2D))
                q,qL                            = parameters.get_qs(Cn_2D, CnL_2D, Fis, FisL)
                q_plus_1, qL_plus_1             = parameters.get_qs(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1)
                #qL = parameters.get_qs(Cn_2D, CnL_2D, Fis, FisL)[1]
                self.q[self.n, pdep:, :]        = time_interpolate(self.times[self.n], t, delt, q, q_plus_1)
                self.qL[self.n, pdep:, :]       = time_interpolate(self.times[self.n], t, delt, qL, qL_plus_1)
                self.M[self.n, :, :]            = time_interpolate(self.times[self.n], t, delt, parameters.get_Ms(Cn_2D, CnL_2D, Fis, FisL), parameters.get_Ms(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1))
                qm, qmL                         = parameters.get_qms(Cn_2D, CnL_2D, Fis, FisL)
                qm_plus_1, qmL_plus_1           = parameters.get_qms(Cn_plus_1_2D, CnL_plus_1_2D, Fis_plus_1, FisL_plus_1)
                self.qm[self.n, pdep:, :, :]    = time_interpolate(self.times[self.n], t, delt, qm,  qm_plus_1)
                self.qmL[self.n, pdep:, :, :]   = time_interpolate(self.times[self.n], t, delt, qmL, qmL_plus_1)

            if parameters.topBCtype == 'Finite mixed water column' or parameters.topBCtype == 'Mass transfer':
                Cw, Cw_plus_1 = parameters.get_Cws()
                self.Cw[self.n, :] = time_interpolate(self.times[self.n], t, delt, Cw, Cw_plus_1)

            self.n = self.n + 1

def first_deriv_2pt_fwd(x):
    """Returns the finite difference coefficients for the first derivative for
    a two-point forward finite difference equation with uneven grid spacing.
    The variable "x" is an array with the values of the independent variable
    at points xi and xi+1."""

    return array([-1. / (x[1] - x[0]), 1. / (x[1] - x[0])])


def first_deriv_3pt_fwd(x):
    """Returns the finite difference coefficients for the first derivative for
    a three-point forward finite difference equation with uneven grid spacing.
    The variable "x" is an array with the values of the independent variable
    at points xi, xi+1, xi+2."""

    y = []
    y.append((2. * x[0] - x[1] - x[2]) / (x[1] - x[0]) / (x[2] - x[0]))
    y.append((x[2] - x[0]) / float(x[1] - x[0]) / (x[2] - x[1]))
    y.append((x[1] - x[0]) / float(x[2] - x[0]) / (x[1] - x[2]))
    return array(y)

def first_deriv_3pt_bwd(x):
    """Returns the finite difference coefficients for the first derivative for
    a three-point backward finite difference equation with uneven grid spacing.
    The variable "x" is an array with the values of the independent variable
    at points xi-2, xi-1, and xi."""

    y = []
    y.append((x[2] - x[1]) / float(x[2] - x[0]) / (x[1] - x[0]))
    y.append((x[0] - x[2]) / float(x[2] - x[1]) / (x[1] - x[0]))
    y.append((2 * x[2] - x[1] - x[0]) / float(x[2] - x[1]) / (x[2] - x[0]))
    return array(y)

def first_deriv_3pt_cen(x):
    """Returns the finite difference coefficients for the first derivative for
    a three-point forward finite difference equation with uneven grid spacing.
    The variable "x" is an array with the values of the independent variable
    at points xi-1, xi, and xi+1."""

    y = []
    y.append((x[1] - x[2]) / float(x[1] - x[0]) / (x[2] - x[0]))
    y.append((x[0] + x[2] - 2. * x[1]) / float(x[1] - x[0]) / (x[2] - x[1]))
    y.append((x[1] - x[0]) / float(x[2] - x[0]) / (x[2] - x[1]))
    return array(y)

def first_deriv_4pt_upw(x):
    """Returns the finite difference coefficients for the first derivative for
    a four-point finite difference equation with uneven grid spacing.  The
    algorithm uses a linear com+bination of the 3-point centered scheme and the
    3-point forward scheme (2/3 to 1/3 ratio, respectively).  The scheme
    provides good stability for hyperbolic problems with better accuracy than
    a simple forward difference."""

    y      = zeros(4)
    z      = zeros(4)
    y[0:3] = first_deriv_3pt_cen(x[0:3])
    z[1:4] = first_deriv_3pt_fwd(x[1:4])
    return 2. / 3 * y + 1. / 3 * z

def second_deriv_3pt_cen(x):
    """Returns the finite difference coefficients for the second derivative for
    a three-point finite difference equation with uneven grid spacing. The
    variable "x" is an array with the values of the independent variable at
    points xi-1, xi, and xi+1."""

    y = []
    y.append(2.  / (x[2] - x[0]) / (x[1] - x[0]))
    y.append(-2. / (x[1] - x[0]) / (x[2] - x[1]))
    y.append(2.  / (x[2] - x[0]) / (x[2] - x[1]))
    return array(y)

def get_delzmax(U, D, elam):
    """Calculates the maximum grid spacing for each layer to prevent
    oscillations in the discretized solution (mesh Peclet number = 2).
    The function uses the effective diffusion for a layer "D" and the
    Darcy velocity "U." """

    if U != 0:    delzmax = 2. * D / abs(U)
    else:         delzmax = 0

    if elam != 0: delzmax_e = 2 * (D / elam)**0.5
    else:         delzmax_e = 0

    if U != 0 and elam != 0: delzmax = min(delzmax,delzmax_e)
    else:                    delzmax = max(delzmax,delzmax_e)

    return delzmax

def get_max_time_step(delz, D, R):
    """Calculates the maximum time step size from the Courant-Friedrichs-Lewy
    condition (R * delz**2 / 2D).  "R" is the retardation factor, "delz" is
    the grid spacing, and "D" the effective diffusion coefficient."""

    return R * delz**2 / 2. / D

def top_mass_transfer_boundary(D, KDbiop, kbl, Cw, z):
    """Returns the finite difference equation for the top mass transfer
    boundary condition using the benthic boundary layer mass transfer
    coefficient "kbl," the diffusion coefficient in the top layer "D," the
    overlying water concentration "Cw," and the top three grid points
    (in array "z")."""

    kbl = kbl #convert the unit of k from cm/hr to cm/yr

    LHS = -D   * first_deriv_3pt_fwd(z) + array([kbl, 0 , 0]) - array(KDbiop) * first_deriv_3pt_fwd(z[0:3])
    RHS = kbl * Cw

    return LHS, RHS

def top_mass_transfer_boundary_kinetic(D, U, KDbiop, kbl, Cw, z):
    """Returns the finite difference equation for the top mass transfer
    boundary condition using the benthic boundary layer mass transfer
    coefficient "kbl," the diffusion coefficient in the top layer "D," the
    overlying water concentration "Cw," and the top three grid points
    (in array "z")."""

    kbl = kbl #convert the unit of k from cm/hr to cm/yr
    if U >= 0:
        LHS = -D   * first_deriv_2pt_fwd(z) + array([kbl, 0]) - array(KDbiop) * first_deriv_2pt_fwd(z) - [-U, U]
        RHS = kbl * Cw
    else:
        LHS = -D   * first_deriv_2pt_fwd(z) + array([kbl, 0]) - array(KDbiop) * first_deriv_2pt_fwd(z) - [U, 0]
        RHS = kbl * Cw - U *Cw

    return LHS, RHS

def top_CSTR_boundary(D, KDbiop, U, kbl, tau, h, z):
    """Returns the finite difference equation for the top mass transfer
    boundary condition using the benthic boundary layer mass transfer
    coefficient "kbl," the diffusion coefficient in the top layer "D," the
    overlying water concentration "Cw," and the top three grid points
    (in array "z")."""

    kbl = kbl #convert the unit of k from cm/hr to cm/yr
    k_Cw  = (U+kbl)/h/(1/tau+kbl/h)

    LHS = -D   * first_deriv_3pt_fwd(z) - array(KDbiop) * first_deriv_3pt_fwd(z[0:3]) + array([kbl*(1-k_Cw), 0 , 0])
    RHS = 0
    return LHS, RHS

def top_CSTR_boundary_kinetic(D, KDbiop, U, kbl, tau, h, z):
    """Returns the finite difference equation for the top mass transfer
    boundary condition using the benthic boundary layer mass transfer
    coefficient "kbl," the diffusion coefficient in the top layer "D," the
    overlying water concentration "Cw," and the top three grid points
    (in array "z")."""

    kbl = kbl #convert the unit of k from cm/hr to cm/yr
    #k_Cw  = (U+kbl)/h/(1/tau+kbl/h)

    if U >= 0:
        k_Cw  = (U+kbl)/h/(1/tau+kbl/h)
        LHS = -D   * first_deriv_2pt_fwd(z) - array(KDbiop) * first_deriv_2pt_fwd(z) + array([kbl*(1-k_Cw), 0]) - [-U, U]
    else:
        k_Cw  = (kbl)/h/(1/tau+kbl/h-U/h)
        LHS = -D   * first_deriv_2pt_fwd(z) - array(KDbiop) * first_deriv_2pt_fwd(z) + array([kbl*(1-k_Cw), 0]) - [U * (1-k_Cw), 0]

    RHS = 0

    return LHS, RHS

def interface_boundary_kinetic(D1, U1, KDbiop1, D2, U2, KDbiop2, z):
    """Returns the finite difference equation for the interface between two
    layers.  The diffusion coefficient in the first layer is "D1," the diffusion
    coefficient in the second layer is "D2", and "z" is an array containing the
    grid points two above the interface, at the interface, and two below the
    interface (total length five)."""

    equation      = zeros(3)
    if U1 >= 0:
        equation[0:2] = D1 * first_deriv_2pt_fwd(z[0:2]) + array([0, U1]) + array(KDbiop1) * first_deriv_2pt_fwd(z[0:2])
        equation[1:3] = equation[1:3] - D2 * first_deriv_2pt_fwd(z[1:3]) - array([0, U2]) - array(KDbiop2) * first_deriv_2pt_fwd(z[1:3])
    else:
        equation[0:2] = D1 * first_deriv_2pt_fwd(z[0:2]) + array([U1, 0]) + array(KDbiop1) * first_deriv_2pt_fwd(z[0:2])
        equation[1:3] = equation[1:3] - D2 * first_deriv_2pt_fwd(z[1:3]) - array([U2, 0]) - array(KDbiop2) * first_deriv_2pt_fwd(z[1:3])

    return equation

def interface_boundary(D1, U1, KDbiop1, D2, U2, KDbiop2, z):
    """Returns the finite difference equation for the interface between two
    layers.  The diffusion coefficient in the first layer is "D1," the diffusion
    coefficient in the second layer is "D2", and "z" is an array containing the
    grid points two above the interface, at the interface, and two below the
    interface (total length five)."""

    equation      = zeros(5)
    equation[0:3] = D1 * first_deriv_3pt_bwd(z[0:3]) + array([0, 0, U1]) + array(KDbiop1) * first_deriv_3pt_bwd(z[0:3])
    equation[2:5] = equation[2:5] - D2 * first_deriv_3pt_fwd(z[2:5]) - array([U2, 0, 0]) - array(KDbiop2) * first_deriv_3pt_fwd(z[2:5])

    return equation

def comp_interface_boundary(D1, U1, KDbiop1, D2, U2, KDbiop2, z):
    """Returns the finite difference equation for the interface between two
    layers.  The diffusion coefficient in the first layer is "D1," the diffusion
    coefficient in the second layer is "D2", and "z" is an array containing the
    grid points two above the interface, at the interface, and two below the
    interface (total length five)."""

    equation      = zeros(5)
    equation[0:3] = -D1 * first_deriv_3pt_bwd(z[0:3]) - array([0, 0, U1]) - array(KDbiop1) * first_deriv_3pt_bwd(z[0:3])
    equation[2:5] = equation[2:5] + D2 * first_deriv_3pt_fwd(z[2:5]) + array([U2, 0, 0]) + array(KDbiop2) * first_deriv_3pt_fwd(z[2:5])

    #equation      = zeros(5)
    #equation[1:3] = -D1 * first_deriv_2pt_fwd(z[1:3]) - array([0, U1]) - array(KDbiop1) * first_deriv_2pt_fwd(z[1:3])
    #equation[2:4] = equation[2:4] + D2 * first_deriv_2pt_fwd(z[2:4]) + array([U2, 0]) + array(KDbiop2) * first_deriv_2pt_fwd(z[2:4])

    return equation

def flux_bottom_boundary(U, D, Dbiops, C0, z):
    """Returns the finite difference equation for a flux-matching bottom
    boundary condition using the Darcy velocity "U," the diffusion coefficient
    in the bottom layer "D," the underlying sediment concentration "C0," and
    the bottom three grid points "x." """

    y    = D * first_deriv_3pt_bwd(z) + array(Dbiops) * first_deriv_3pt_bwd(z)
    y[2] = y[2] + U
    return array(y), U * C0

def flux_bottom_boundary_kinetic(U, D, Dbiops, C0, z):
    """Returns the finite difference equation for a flux-matching bottom
    boundary condition using the Darcy velocity "U," the diffusion coefficient
    in the bottom layer "D," the underlying sediment concentration "C0," and
    the bottom three grid points "x." """

    y    = D * first_deriv_2pt_fwd(z) + array(Dbiops) * first_deriv_2pt_fwd(z)
    if U >= 0:   y[1] = y[1] + U
    else:        y[0] = y[0] + U

    return array(y), U * C0


def get_4pt_adr_fde_imp(R, R_plus_1, U_plus_1, D_plus_1, K_plus_1, delt, z):
    """Returns the LHS of the finite difference equation at a point for the
    advection-diffusion-reaction equation with sorption using one point
    downwind and two points upwind.  "U" is the Darcy velocity, "D" is the
    diffusion coefficient, "e" is the porosity, "lam" is the first-order
    decay rate constant, "delt" is the time step size, and "z" is an array
    of the four points used in the discretization."""

    a      = zeros(4)
    a[0:3] = second_deriv_3pt_cen(z)
    b      = zeros(4)
    if U_plus_1 >= 0:    b[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                b[0:2] = first_deriv_2pt_fwd(z[0:2])
    c      = array([0, 1, 0, 0])

    return - array(D_plus_1 + [0]) * a - array(K_plus_1 + [0]) * a - U_plus_1 * b + (R_plus_1/delt)*c, (R/delt)*c

def get_3pt_adr_fde_imp(R, R_plus_1, U_plus_1, D_plus_1, K_plus_1, delt, z):
    """Returns the LHS of the finite difference equation at a point for the
    advection-diffusion-reaction equation with sorption using one point
    downwind and one point upwind.  "U" is the Darcy velocity, "D" is the
    diffusion coefficient, "e" is the porosity, "lam" is the first-order
    decay rate constant, "delt" is the time step size, and "z" is an array
    of the four points used in the discretization."""

    a      = second_deriv_3pt_cen(z)
    b      = zeros(3)

    #
    if U_plus_1 >= 0:    b[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                b[0:2] = first_deriv_2pt_fwd(z[0:2])
    c      = array([0, 1, 0])

    return -array(D_plus_1) * a - array(K_plus_1) * a  - U_plus_1 * b + (R_plus_1/delt)*c, (R/delt)*c


def get_4pt_adr_fde_CN(R, R_plus_1, U, U_plus_1, D, D_plus_1, K, K_plus_1, delt, z):
    """Returns the LHS of the finite difference equation at a point for the 
    advection-diffusion-reaction equation with sorption using one point 
    downwind and two points upwind.  "U" is the Darcy velocity, "D" is the 
    diffusion coefficient, "e" is the porosity, "lam" is the first-order 
    decay rate constant, "delt" is the time step size, and "z" is an array 
    of the four points used in the discretization."""

    a               = zeros(4)
    a[0:3]          = second_deriv_3pt_cen(z)
    b      = zeros(4)
    if U >= 0:           b[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                b[0:2] = first_deriv_2pt_fwd(z[0:2])
    b_plus_1 = zeros(4)
    if U_plus_1 >= 0:    b_plus_1[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                b_plus_1[0:2] = first_deriv_2pt_fwd(z[0:2])
    c               = array([0, 1, 0, 0])

    return (-array(D_plus_1 + [0])*a - array(K_plus_1 + [0])*a - U_plus_1 * b + (2*R_plus_1/delt)*c)/2, (array(D + [0])*a + array(K + [0])*a + U*b + (2*R/delt)*c)/2

def get_3pt_adr_fde_CN(R, R_plus_1, U, U_plus_1, D, D_plus_1, K, K_plus_1, delt, z):
    """Returns the LHS of the finite difference equation at a point for the 
    advection-diffusion-reaction equation with sorption using one point 
    downwind and one point upwind.  "U" is the Darcy velocity, "D" is the 
    diffusion coefficient, "e" is the porosity, "lam" is the first-order 
    decay rate constant, "delt" is the time step size, and "z" is an array 
    of the four points used in the discretization."""

    a           = second_deriv_3pt_cen(z)
    b           = zeros(3)
    b_plus_1    = zeros(3)
    if U >= 0:              b[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                   b[0:2] = first_deriv_2pt_fwd(z[0:2])
    if U_plus_1 >= 0:       b_plus_1[1:3] = first_deriv_2pt_fwd(z[1:3])
    else:                   b_plus_1[0:2] = first_deriv_2pt_fwd(z[0:2])
    c      = array([0, 1, 0])

    return (-array(D_plus_1) * a - array(K_plus_1) * a- U_plus_1*b_plus_1 + (2*R_plus_1/delt)*c)/2, (array(D)*a + array(K)*a + U*b + (2*R/delt)*c)/2

def time_interpolate(tint, t, delt, Cn, Cn_plus_1):
    """Returns the interpolated concentrations at time "tint" using the
    concentrations "Cn" at time "t" and concentrations "Cn_plus_1" at time 
    "t + delt." """

    return (Cn + (tint - t) / delt * (Cn_plus_1 - Cn))

def get_surface_flux(U, kbl, C, Cw):
    """Returns the surface flux using the benthic boundary layer mass transfer 
    coefficient "kbl," the concentration at the cap-water interface "C," and
    the overlying water concentration "Cw." """

    if U >= 0:  F = (kbl * (C - Cw))+ U * C
    else:       F = (kbl * (C - Cw))+ U * Cw

    return F

def get_top_flux(U, D, KD, C, z):

    return U * C[0] + D * sum(first_deriv_2pt_fwd(z) *  C) + sum(array(KD) * first_deriv_2pt_fwd(z) * C)

def get_boundary_flux(U, D, KD, C, z):
    """Returns the flux at a boundary between two layers using the properties
    of the overlying layer.  Uses the three-point backwards difference, the
    overlying effective diffusion coefficient, the concentrations at the three
    points above in array "C," and the grid depths in array "z." """

    return U * C[1] + D * sum(first_deriv_2pt_fwd(z) * C) + sum(array(KD) * first_deriv_2pt_fwd(z) * C)

def get_point_flux(U, D, KD, C, z):
    """Returns the flux at a grid point that is not a boundary.  Uses the three-
    point centered difference, the diffusion coefficient at the point "D," the
    concentrations at the points before, at and after "C," and the grid points
    "z." """

    return U * C[1] + D * sum(first_deriv_3pt_cen(z) *  C) + sum(array(KD) * first_deriv_3pt_cen(z) * C)

def get_2_point_upwind_flux(U, D, C, z):
    """Returns the flux at a grid point that is not a boundary.  Uses the three-
    point centered difference, the diffusion coefficient at the point "D," the
    concentrations at the points before, at and after "C," and the grid points
    "z." """

    return U * C[0] + D * sum(first_deriv_2pt_fwd(z) *  C)

def get_2_point_downwind_flux(U, D, C, z):
    """Returns the flux at a grid point that is not a boundary.  Uses the three-
    point centered difference, the diffusion coefficient at the point "D," the
    concentrations at the points before, at and after "C," and the grid points
    "z." """

    return U * C[1] + D * sum(first_deriv_2pt_fwd(z) *  C)

def get_nonpoint_flux(U, D, C1, C2, z1, z2):
    """Returns the flux at a point that is not part of the grid.  Uses the two-
    point centered difference between points "z1" and "z2" and the 
    concentrations "C1" and "C2" at those points and the diffusion coefficient
    at the point "D." """
    
    return U * (C1 + C2) / 2 + D * (C2 - C1) / (z2 - z1)
