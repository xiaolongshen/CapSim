# -*- coding: utf-8 -*-

#!/usr/bin/env python
#These subroutines are used in post-processing for CapSim.

import matplotlib.pyplot as plt, tkMessageBox as tkmb, math, sys, os, codecs, csv, _winreg as wreg
import matplotlib.patches  as patches
import matplotlib.gridspec as gridspec
import tkFileDialog as tkfd

if sys.path[0][-3:] == 'zip': 
    os.chdir(sys.path[0][:-12])
    path = sys.path[0][:-12]
else: path = sys.path[0]

CapSimReg = wreg.ConnectRegistry(None, wreg.HKEY_CURRENT_USER)
CapSimKey = wreg.OpenKey(CapSimReg, r'Software\CapSim')
Filepath =wreg.QueryValueEx(CapSimKey, 'FilePath')[0]
CapSimKey.Close()

from numpy               import zeros, transpose, interp, array, isnan
from capsim_object_types import System, CapSimWindow
from postprocess         import FigureEditor, GraphEditor
from datetime            import datetime
from Tkinter             import Tk, Toplevel, Canvas, Frame, Label, Entry, Text, Button, Scrollbar, OptionMenu, StringVar, DoubleVar, IntVar, FLAT, RAISED
from PIL                 import Image, ImageTk
from textwrap            import fill

class Graphprocess:
    """Makes a window to display the simulation plots."""

    def __init__(self, master, system):
        """Constructor method."""

        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.system    = system
        rgb            = self.master.frame.winfo_rgb(self.master.frame.cget('bg'))
        self.color     = '#%02x%02x%02x' %rgb
        self.mplrgb    = (rgb[0] / 65536., rgb[1] / 65536., rgb[2] / 65536.) 
        self.tframe    = Frame(master.tframe, bg = self.color)
        self.frame     = Frame(master.frame,  bg = self.color)
        self.bframe    = Frame(master.bframe, bg = self.color)
        self.top       = None
        self.flag      = 0
        self.filename  = 'plot'

        self.lengthunit = system.lengthunit
        self.concunit   = system.concunit
        self.timeunit   = system.timeunit
        self.diffunit   = system.diffunit

        self.lengthunits = [u'\u03BCm', 'cm', 'm']
        self.concunits   = [u'\u03BCg/L', 'mg/L', 'g/L', u'\u03BCmol/L', 'mmol/L', 'mol/L']
        self.timeunits   = ['s', 'min', 'hr', 'day', 'yr']
        self.diffunits   = [u'cm\u00B2/s', u'cm\u00B2/yr']

        self.colors      = ['b', 'g', 'r', 'c', 'm', 'y', '#ff7f0e', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', 'k']
        self.colornames  = ['Blue', 'Green', 'Red', 'Cyan', 'Magenta', 'Yellow', 'Orange', 'Purple', 'Brown', 'Pink', 'Grey', 'Black']

        self.styles      = ['-', '--', '-.', ':', '.', 'o', '^', 's', '*', 'x']
        self.stylenames  = ['Solid', 'Dashed', 'Dash-dot', 'Dotted', 'Point', 'Circle', 'Triangle', 'Square', 'Star', 'Cross']

        self.types        = ['Spatial profile', 'Time profile']
        self.variables    = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Water concentration', 'Material fraction', 'Mass']

        self.type        = self.types[0]
        self.variable    = self.variables[0]

        self.sketch       = {}
        self.sketch['Sketches'] = ['Hide sketch']
        self.sketch['Sketch']   = 'Hide sketch'
        self.sketch['Label']    = 1
        self.sketch['Time']     = 0
        self.sketch['Depth']    = 0

        self.profiles    = []

        self.output_type = []
        self.chemicals  = []
        self.components = []
        self.solids     = []

        self.z          = []
        self.times      = []
        self.dep        = []
        self.Vdep       = []
        self.C          = []
        self.F          = []
        self.Fi         = []
        self.q          = []
        self.W          = []
        self.qm         = []
        self.Cw         = []
        self.M          = []
        self.FDM        = []
        self.RM         = []

        self.C_smax      = []
        self.F_smax      = []
        self.Fi_smax     = []
        self.q_smax      = []
        self.W_smax      = []
        self.qm_smax     = []
        self.Cw_smax     = []
        self.M_smax      = []
        self.FDM_smax    = []
        self.RM_smax     = []

        self.C_smin      = []
        self.F_smin      = []
        self.Fi_smin     = []
        self.q_smin      = []
        self.W_smin      = []
        self.qm_smin     = []
        self.Cw_smin     = []
        self.M_smin      = []
        self.FDM_smin    = []
        self.RM_smin     = []

        self.C_max      = []
        self.F_max      = []
        self.Fi_max     = []
        self.q_max      = []
        self.W_max      = []
        self.qm_max     = []
        self.Cw_max     = []
        self.M_max      = []
        self.FDM_max    = []
        self.RM_max     = []

        self.C_min      = []
        self.F_min      = []
        self.Fi_min     = []
        self.q_min      = []
        self.W_min      = []
        self.qm_min     = []
        self.Cw_min     = []
        self.M_min      = []
        self.FDM_min    = []
        self.RM_min     = []

        self.solids     = []
        self.layer_type = []
        self.layer_list = []
        self.layer_h    = []
        self.hbio       = []

        self.timeplotdatas      = []
        self.spatialplotdatas   = []

        self.figuresize = [600, 450]

        if master.tk.winfo_screenheight()/master.tk.winfo_screenwidth() > 3/4:
            if master.tk.winfo_screenheight()*3/5 < 450:
                self.figuresize[1] = int(master.tk.winfo_screenheight()*3/5)
                self.figuresize[0] = int(self.figuresize[1]/3*4)
        else:
            if master.tk.winfo_screenwidth()*4/5 < 600:
                self.figuresize[0] = int(master.tk.winfo_screenwidth()*4/5)
                self.figuresize[1] = int(self.figuresize[0]*3/4)

        self.axislimit  = [0, 1.2, 0, 1]
        self.legendposition = 0

        self.figuretexts = {}
        self.figuretexts['Flag']    = 'Default'
        self.figuretexts['title']   = 'Porewater Concentration Profiles'
        self.figuretexts['xlabel']  = 'Concentration ('+self.concunit+')'
        self.figuretexts['ylabel' ] = 'Depth ('+self.lengthunit+')'

        self.fontsize = {}
        self.fontsize['Flag']   = 'Default'
        self.fontsize['Title']  = 12
        self.fontsize['Label']  = [10, 10]
        self.fontsize['Axis' ]  = [8, 8]
        self.fontsize['Legend'] = 9

        self.MCplots = {}
        self.MCplots['Flag']        = 'Default'
        self.MCplots['AVG']         = 1.
        self.MCplots['SD']          = 0.2
        self.MCplots['Max/Min' ]    = 0.1

        self.graphpath = Filepath + '\output\\'
        self.make_graphs()

        self.gs = gridspec.GridSpec(1, 2, width_ratios=[1, 5])
        self.gs.update(left = 0.01, right = 0.95)


    def make_widgets(self):
        """Makes the widgets for the window."""

        self.master.frame.config(bg      = self.color)
        self.master.mainbutton.config(bg = self.color)
        self.master.exitbutton.config(bg = self.color)

        self.graph          = Image.open(Filepath + r'/output/%s.png' % self.filename)
        self.image          = ImageTk.PhotoImage(self.graph)
        self.graphlabel     = Label(self.frame, image = self.image, bg = self.color)

        self.blank1         = Label(self.bframe, text = ' ')
        self.blank2         = Label(self.bframe, text = ' ')
        self.blank3         = Label(self.bframe, text = ' ')

        self.importbutton = Button(self.bframe, text = 'Import profiles',   bg = self.color, command = self.importdata,  width = 20, activebackground = 'white', highlightbackground = 'black')
        self.editbutton   = Button(self.bframe, text = 'Edit plot',         bg = self.color, command = self.editplot,    width = 20, activebackground = 'white', highlightbackground = 'black')
        self.editgbutton  = Button(self.bframe, text = 'Edit profiles',     bg = self.color, command = self.editprofile, width = 20, activebackground = 'white', highlightbackground = 'black')
        self.figurebutton = Button(self.bframe, text = 'Edit figure',       bg = self.color, command = self.editfigure,  width = 20, activebackground = 'white', highlightbackground = 'black')
        self.savebutton   = Button(self.bframe, text = 'Save figure',       bg = self.color, command = self.savegraph,   width = 20, activebackground = 'white', highlightbackground = 'black')

        self.graphlabel.grid(   row = 0, column = 1, padx = 40)

        self.blank1.grid(       row = 1, column = 0, columnspan = 7, padx = 40)
        self.importbutton.grid( row = 2, column = 0, columnspan = 7, padx = 40)

        self.focusbutton = self.importbutton

        self.master.geometry()
        self.master.center()

    def importdata(self, event = None):

        self.dataname = tkfd.askopenfilename(initialdir = Filepath +r'\output',  filetypes = [('Output files','*.csv')])

        if self.dataname != '':

            content       = []
            file          = open(self.dataname, 'r')
            file_content  = csv.reader(file)
            for row in file_content:    content.append(row)
            rows     = len(content)
            row = 0

            self.chemicals.append([])
            self.components.append([])
            self.solids.append([])
            self.layer_type.append([])
            self.layer_h .append([])
            self.layer_list.append([])
            self.solids[-1].append('Total solid')
            self.z.append([])
            self.times.append([])
            self.dep.append(0)
            self.Vdep.append(0)
            self.hbio.append(0)

            if content[0][0] == 'CapSim':

                self.output_type.append(content[2][1])

                property_check = 0
                while row < rows and property_check == 0:
                    if len(content[row]) > 0:
                        if content[row][0] == 'System Units':                   row_unit        = row
                        if content[row][0] == 'Chemicals':                      row_chem        = row
                        if content[row][0] == 'Matrices':                       row_matrix      = row
                        if content[row][0] == 'Matrix Components':              row_comp        = row
                        if content[row][0] == 'Layers':                         row_layer       = row
                        if content[row][0] == 'System Properties':              row_system      = row
                        if content[row][0] == 'Auxiliary Conditions':           row_BCIC        = row
                        if content[row][0] == 'Output Results':
                            row_output     = row
                            property_check = 1
                    row = row + 1

                nchemicals  = int(content[row_chem][1])
                ncomponents = int(content[row_comp][1])
                nlayers     = int(content[row_layer][1])
                topBCtype   = content[row_layer+2][2]

                for i in range(nlayers):
                    if content[row_layer+2+i][1] == 'Deposition':   self.layer_list[-1].append('Deposition Layer')
                    else:                                           self.layer_list[-1].append(content[row_layer+2+i][1])
                    self.layer_type[-1].append(content[row_layer+2+i][2])
                    self.layer_h[-1].append(float(content[row_layer+2+i][3]))

                if self.layer_list[-1][0] == 'Deposition Layer':
                    self.dep[-1] = 1
                    self.Vdep[-1] = self.layer_h[-1][0]

                if content[row_system + 5][2] == 'Uniform':
                    self.hbio[-1] = float(content[row_system + 6][2])

                row_solids     = []
                row_conc       = []
                row_flux       = []
                row_solid      = []
                row_whole      = []
                row_water      = []
                row_mass_layer = []
                row_mass_IOD   = []
                row_mass_RG    = []
                row_Fis        = []

                for n in range(nchemicals):
                    chemical_check = 0
                    row_solids.append([])
                    while row < rows and chemical_check == 0:
                        if len(content[row]) > 0:
                            if content[row][0] == 'Porewater concentrations':       row_conc.append(row)
                            if content[row][0] == 'Fluxes':                         row_flux.append(row)
                            if content[row][0] == 'Solid concentrations':
                                if content[row][1] == 'Total solid':                row_solid.append(row)
                                else:
                                    row_solids[-1].append(row)
                            if content[row][0] == 'Total concentrations':            row_whole.append(row)
                            if content[row][0] == 'Overlying water concentrations':  row_water.append(row)
                            if content[row][0] == 'Mass in layers':                  row_mass_layer.append(row)
                            if content[row][0] == 'Mass in top/bottom/deposition': row_mass_IOD.append(row)
                            if content[row][0] == 'Mass by reaction/generation':
                                row_mass_RG.append(row)
                                chemical_check = 1
                        row = row + 1

                for m in range(ncomponents):
                    chemical_check = 0
                    while row < rows and chemical_check == 0:
                        if len(content[row]) > 0:
                            if content[row][0] == 'Component':
                                row_Fis.append(row)
                                chemical_check = 1
                        row = row + 1

                lengthunit = content[row_unit+1][2]
                concunit   = content[row_unit+2][2]
                timeunit   = content[row_unit+3][2]

                lengthunits = ['um^2', 'cm^2', 'm^2']
                concunits   = ['ug/L', 'mg/L', 'g/L', 'umol/L', 'mmol/L', 'mol/L']
                timeunits   = ['s', 'min', 'hr', 'day', 'yr']

                length_converter = 1.
                time_converter   = 1.

                if self.lengthunit == self.lengthunits[0]:  length_converter = length_converter * 10000
                if self.lengthunit == self.lengthunits[1]:  length_converter = length_converter
                if self.lengthunit == self.lengthunits[2]:  length_converter = length_converter / 100

                if lengthunit == lengthunits[0]:            length_converter = length_converter / 10000
                if lengthunit == lengthunits[1]:            length_converter = length_converter
                if lengthunit == lengthunits[2]:            length_converter = length_converter * 100

                if self.timeunit == self.timeunits[0]:      time_converter = time_converter*365.25*24*3600
                if self.timeunit == self.timeunits[1]:      time_converter = time_converter*24*3600
                if self.timeunit == self.timeunits[2]:      time_converter = time_converter*3600
                if self.timeunit == self.timeunits[3]:      time_converter = time_converter*60
                if self.timeunit == self.timeunits[4]:      time_converter = time_converter

                if timeunit == timeunits[0]:                time_converter = time_converter/365.25/24/3600
                if timeunit == timeunits[1]:                time_converter = time_converter/24/3600
                if timeunit == timeunits[2]:                time_converter = time_converter/3600
                if timeunit == timeunits[3]:                time_converter = time_converter/60
                if timeunit == timeunits[4]:                time_converter = time_converter

                row = row_conc[0] + 2

                num_t = int(content[row_output+2][1])
                for t in range(num_t):
                    self.times[-1].append(float(content[row][t+1]) * time_converter)

                num_z = int(content[row_output+1][1])
                for z in range(num_z):
                    self.z[-1].append(float(content[row+z+1][0]) * length_converter)

                self.C.append(zeros([num_t, num_z, nchemicals]))
                self.F.append(zeros([num_t, num_z, nchemicals]))
                self.q.append(zeros([num_t, num_z, nchemicals, len(row_solids[-1])+1]))
                self.W.append(zeros([num_t, num_z, nchemicals]))
                self.Cw.append(zeros([num_t, nchemicals]))
                self.M.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                self.RM.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                self.FDM.append(zeros([3, num_t, nchemicals]))
                self.Fi.append(zeros([num_t, num_z, ncomponents]))

                if self.output_type[-1] == 'Monte Carlo Simulation':

                    self.Fi_smax.append(zeros([num_t, num_z, ncomponents]))
                    self.C_smax.append(zeros([num_t, num_z, nchemicals]))
                    self.F_smax.append(zeros([num_t, num_z, nchemicals]))
                    self.q_smax.append(zeros([num_t, num_z, nchemicals, len(row_solids[-1])+1]))
                    self.W_smax.append(zeros([num_t, num_z, nchemicals]))
                    self.Cw_smax.append(zeros([num_t, nchemicals]))
                    self.M_smax.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.RM_smax.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.FDM_smax.append(zeros([3, num_t, nchemicals]))

                    self.Fi_smin.append(zeros([num_t, num_z, ncomponents]))
                    self.C_smin.append(zeros([num_t, num_z, nchemicals]))
                    self.F_smin.append(zeros([num_t, num_z, nchemicals]))
                    self.q_smin.append(zeros([num_t, num_z, nchemicals, len(row_solids[-1])+1]))
                    self.W_smin.append(zeros([num_t, num_z, nchemicals]))
                    self.Cw_smin.append(zeros([num_t, nchemicals]))
                    self.M_smin.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.RM_smin.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.FDM_smin.append(zeros([3, num_t, nchemicals]))

                    self.Fi_max.append(zeros([num_t, num_z, ncomponents]))
                    self.C_max.append(zeros([num_t, num_z, nchemicals]))
                    self.F_max.append(zeros([num_t, num_z, nchemicals]))
                    self.q_max.append(zeros([num_t, num_z, nchemicals, len(row_solids[-1])+1]))
                    self.W_max.append(zeros([num_t, num_z, nchemicals]))
                    self.Cw_max.append(zeros([num_t, nchemicals]))
                    self.M_max.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.RM_max.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.FDM_max.append(zeros([3, num_t, nchemicals]))

                    self.Fi_min.append(zeros([num_t, num_z, ncomponents]))
                    self.C_min.append(zeros([num_t, num_z, nchemicals]))
                    self.F_min.append(zeros([num_t, num_z, nchemicals]))
                    self.q_min.append(zeros([num_t, num_z, nchemicals, len(row_solids[-1])+1]))
                    self.W_min.append(zeros([num_t, num_z, nchemicals]))
                    self.Cw_min.append(zeros([num_t, nchemicals]))
                    self.M_min.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.RM_min.append(zeros([len(self.layer_h[-1]), num_t, nchemicals]))
                    self.FDM_min.append(zeros([3, num_t, nchemicals]))
                    self.Fi_min.append(zeros([num_t, num_z, ncomponents]))
                else:
                    self.C_smax.append([])
                    self.F_smax.append([])
                    self.q_smax.append([])
                    self.W_smax.append([])
                    self.qm_smax.append([])
                    self.Cw_smax.append([])
                    self.M_smax.append([])
                    self.RM_smax.append([])
                    self.FDM_smax.append([])
                    self.Fi_smax.append([])

                    self.C_smin.append([])
                    self.F_smin.append([])
                    self.q_smin.append([])
                    self.W_smin.append([])
                    self.qm_smin.append([])
                    self.Cw_smin.append([])
                    self.M_smin.append([])
                    self.RM_smin.append([])
                    self.FDM_smin.append([])
                    self.Fi_smin.append([])

                    self.C_max.append([])
                    self.F_max.append([])
                    self.q_max.append([])
                    self.W_max.append([])
                    self.qm_max.append([])
                    self.Cw_max.append([])
                    self.M_max.append([])
                    self.RM_max.append([])
                    self.FDM_max.append([])
                    self.Fi_max.append([])

                    self.C_min.append([])
                    self.F_min.append([])
                    self.q_min.append([])
                    self.W_min.append([])
                    self.qm_min.append([])
                    self.Cw_min.append([])
                    self.M_min.append([])
                    self.RM_min.append([])
                    self.FDM_min.append([])
                    self.Fi_min.append([])

                for n in range(nchemicals):
                    self.chemicals[-1].append(content[row_chem+2+n][1])
                    MW = float(content[row_chem+2+n][3])
                    conc_converter   = 1.

                    if self.concunit == self.concunits[0]: conc_converter = conc_converter
                    if self.concunit == self.concunits[1]: conc_converter = conc_converter/1000
                    if self.concunit == self.concunits[2]: conc_converter = conc_converter/1000/1000
                    if self.concunit == self.concunits[3]: conc_converter = conc_converter/MW
                    if self.concunit == self.concunits[4]: conc_converter = conc_converter/1000/MW
                    if self.concunit == self.concunits[5]: conc_converter = conc_converter/1000/1000/MW

                    if concunit == concunits[0]:    conc_converter = conc_converter
                    if concunit == concunits[1]:    conc_converter = conc_converter * 1000
                    if concunit == concunits[2]:    conc_converter = conc_converter * 1000 *1000
                    if concunit == concunits[3]:    conc_converter = conc_converter * MW
                    if concunit == concunits[4]:    conc_converter = conc_converter * MW* 1000
                    if concunit == concunits[5]:    conc_converter = conc_converter * MW* 1000 *1000

                    flux_converter  = conc_converter

                    if self.lengthunit == self.lengthunits[0]:  flux_converter = flux_converter/10000/10000
                    if self.lengthunit == self.lengthunits[1]:  flux_converter = flux_converter
                    if self.lengthunit == self.lengthunits[2]:  flux_converter = flux_converter*100*100

                    if lengthunit == lengthunits[0]:            flux_converter = flux_converter * 10000 * 10000
                    if lengthunit == lengthunits[1]:            flux_converter = flux_converter
                    if lengthunit == lengthunits[2]:            flux_converter = flux_converter /100/100

                    if self.timeunit == self.timeunits[0]:      flux_converter = flux_converter/365.25/24/3600
                    if self.timeunit == self.timeunits[1]:      flux_converter = flux_converter/24/3600
                    if self.timeunit == self.timeunits[2]:      flux_converter = flux_converter/3600
                    if self.timeunit == self.timeunits[3]:      flux_converter = flux_converter/60
                    if self.timeunit == self.timeunits[4]:      flux_converter = flux_converter

                    if timeunit == timeunits[0]:                flux_converter = flux_converter*365.25*24*3600
                    if timeunit == timeunits[1]:                flux_converter = flux_converter*24*3600
                    if timeunit == timeunits[2]:                flux_converter = flux_converter*3600
                    if timeunit == timeunits[3]:                flux_converter = flux_converter*60
                    if timeunit == timeunits[4]:                flux_converter = flux_converter

                    for z in range(num_z):
                        for t in range(num_t):
                            self.C[-1][t,z,n]   = float(content[row_conc[n] + 3 + z][1 + t]) * conc_converter
                            self.F[-1][t,z,n]   = float(content[row_flux[n] + 3 + z][1 + t]) * flux_converter
                            self.q[-1][t,z,n,0] = float(content[row_solid[n] + 3 + z][1 + t])* conc_converter
                            self.W[-1][t,z,n]   = float(content[row_whole[n] + 3 + z][1 + t])* conc_converter
                            for m in range(len(row_solids[-1])):
                                self.q[-1][t,z,n,m+1] = float(content[row_solids[n][m] + 3 + z][1 + t])* conc_converter

                    if row_water != 0:
                        for t in range(num_t):
                            self.Cw[-1][t,n]      = float(content[row_water[n] + 3][1 + t]) * conc_converter

                    if row_mass_layer != 0:
                        for layer in range(len(self.layer_h[-1])):
                            for t in range(num_t):
                                self.M[-1][layer, t, n]      = float(content[row_mass_layer[n] + 3 + layer][1 + t]) * flux_converter

                    if row_mass_IOD != 0:
                        for t in range(num_t):
                            self.FDM[-1][0, t, n]     = float(content[row_mass_IOD[n] + 3][1 + t]) * flux_converter
                            self.FDM[-1][1, t, n]     = float(content[row_mass_IOD[n] + 4][1 + t]) * flux_converter
                            self.FDM[-1][2, t, n]     = float(content[row_mass_IOD[n] + 5][1 + t]) * flux_converter

                    if row_mass_RG != 0:
                        for layer in range(len(self.layer_h[-1])):
                            for t in range(num_t):
                                self.RM[-1][layer, t, n]     = float(content[row_mass_RG[n] + 3 + layer][1 + t]) * flux_converter

                    if self.output_type[-1] == 'Monte Carlo Simulation':

                        for z in range(num_z):
                            for t in range(num_t):
                                self.C_smax[-1][t,z,n]   = float(content[row_conc[n] + 3 + z + 4 + num_z][1 + t]) * conc_converter + self.C[-1][t,z,n]
                                self.F_smax[-1][t,z,n]   = float(content[row_flux[n] + 3 + z + 4 + num_z][1 + t]) * flux_converter + self.F[-1][t,z,n]
                                self.q_smax[-1][t,z,n,0] = float(content[row_solid[n] + 3 + z+ 4 + num_z][1 + t]) * conc_converter + self.q[-1][t,z,n, 0]
                                self.W_smax[-1][t,z,n]   = float(content[row_whole[n] + 3 + z+ 4 + num_z][1 + t]) * conc_converter + self.W[-1][t,z,n]
                                for m in range(len(row_solids)):
                                    self.q_smax[-1][t,z,n,m+1] = float(content[row_solids[n][m] + 3 + z + 4 + num_z][1 + t])* conc_converter + self.q[-1][t,z,n,m+1]
                                self.C_smin[-1][t,z,n]   = self.C[-1][t,z,n] - float(content[row_conc[n] + 3 + z + 4 + num_z][1 + t]) * conc_converter
                                self.F_smin[-1][t,z,n]   = self.F[-1][t,z,n] - float(content[row_flux[n] + 3 + z + 4 + num_z][1 + t]) * flux_converter
                                self.q_smin[-1][t,z,n,0] = self.q[-1][t,z,n, 0] - float(content[row_solid[n] + 3 + z+ 4 + num_z][1 + t]) * conc_converter
                                self.W_smin[-1][t,z,n]   = self.W[-1][t,z,n] - float(content[row_whole[n] + 3 + z+ 4 + num_z][1 + t]) * conc_converter
                                for m in range(len(row_solids)):
                                    self.q_smin[-1][t,z,n,m+1] = self.q[-1][t,z,n,m+1] - float(content[row_solids[n][m] + 3 + z + 4 + num_z][1 + t])* conc_converter
                                self.C_max[-1][t,z,n] = float(content[row_conc[n] + 3 + z + (4 + num_z)*2][1 + t]) * conc_converter
                                self.F_max[-1][t,z,n] = float(content[row_flux[n] + 3 + z + (4 + num_z)*2][1 + t]) * flux_converter
                                self.q_max[-1][t,z,n,0] = float(content[row_solid[n] + 3 + z+ (4 + num_z)*2][1 + t]) * conc_converter
                                self.W_max[-1][t,z,n] = float(content[row_whole[n] + 3 + z+ (4 + num_z)*2][1 + t]) * conc_converter
                                for m in range(len(row_solids)):
                                    self.q_max[-1][t,z,n,m+1] = float(content[row_solids[n][m] + 3 + z + (4 + num_z)*2][1 + t])* conc_converter
                                self.C_min[-1][t,z,n] = float(content[row_conc[n] + 3 + z +(4 + num_z)*3][1 + t]) * conc_converter
                                self.F_min[-1][t,z,n] = float(content[row_flux[n] + 3 + z + (4 + num_z)*3][1 + t]) * flux_converter
                                self.q_min[-1][t,z,n,0] = float(content[row_solid[n] + 3 + z+ (4 + num_z)*3][1 + t]) * conc_converter
                                self.W_min[-1][t,z,n] = float(content[row_whole[n] + 3 + z+ (4 + num_z)*3][1 + t]) * conc_converter
                                for m in range(len(row_solids)):
                                    self.q_min[-1][t,z,n,m+1] = float(content[row_solids[n][m] + 3 + z + (4 + num_z)*3][1 + t])* conc_converter

                        if row_water != 0 and topBCtype == 'Finite mixed water column':
                            for t in range(num_t):
                                self.Cw_smax[-1][t,n]     = self.Cw[-1][t,n] + float(content[row_water[n] + 3 + 5][1 + t]) * conc_converter
                                self.Cw_smin[-1][t,n]     = self.Cw[-1][t,n] - float(content[row_water[n] + 3 + 5][1 + t]) * conc_converter
                                self.Cw_max[-1][t,n]      = float(content[row_water[n] + 3 + 10][1 + t]) * conc_converter
                                self.Cw_min[-1][t,n]      = float(content[row_water[n] + 3 + 15][1 + t]) * conc_converter

                        if row_mass_layer != 0:
                            for layer in range(len(self.layer_h[-1])):
                                for t in range(num_t):
                                    self.M_smax[-1][layer, t, n]     = self.M[-1][layer, t, n] + float(content[row_mass_layer[n] + 3 + layer + 4 + nlayers][1 + t]) * flux_converter
                                    self.M_smin[-1][layer, t, n]     = self.M[-1][layer, t, n] - float(content[row_mass_layer[n] + 3 + layer + 4 + nlayers][1 + t]) * flux_converter
                                    self.M_max[-1][layer, t, n]      = float(content[row_mass_layer[n] + 3 + layer + (4 + nlayers)*2][1 + t]) * flux_converter
                                    self.M_min[-1][layer, t, n]      = float(content[row_mass_layer[n] + 3 + layer + (4 + nlayers)*3][1 + t]) * flux_converter

                        if row_mass_IOD != 0:
                            for t in range(num_t):
                                self.FDM_smax[-1][0, t, n]    = self.FDM[-1][0, t, n] + float(content[row_mass_IOD[n] + 3 + 7][1 + t]) * flux_converter
                                self.FDM_smax[-1][1, t, n]    = self.FDM[-1][1, t, n] + float(content[row_mass_IOD[n] + 4 + 7][1 + t]) * flux_converter
                                self.FDM_smax[-1][2, t, n]    = self.FDM[-1][2, t, n] + float(content[row_mass_IOD[n] + 5 + 7][1 + t]) * flux_converter
                                self.FDM_smin[-1][0, t, n]    = self.FDM[-1][0, t, n] - float(content[row_mass_IOD[n] + 3 + 7][1 + t]) * flux_converter
                                self.FDM_smin[-1][1, t, n]    = self.FDM[-1][1, t, n] - float(content[row_mass_IOD[n] + 4 + 7][1 + t]) * flux_converter
                                self.FDM_smin[-1][2, t, n]    = self.FDM[-1][2, t, n] - float(content[row_mass_IOD[n] + 5 + 7][1 + t]) * flux_converter
                                self.FDM_max[-1][0, t, n]     = float(content[row_mass_IOD[n] + 3 + 14][1 + t]) * flux_converter
                                self.FDM_max[-1][1, t, n]     = float(content[row_mass_IOD[n] + 4 + 14][1 + t]) * flux_converter
                                self.FDM_max[-1][2, t, n]     = float(content[row_mass_IOD[n] + 5 + 14][1 + t]) * flux_converter
                                self.FDM_min[-1][0, t, n]     = float(content[row_mass_IOD[n] + 3 + 21][1 + t]) * flux_converter
                                self.FDM_min[-1][1, t, n]     = float(content[row_mass_IOD[n] + 4 + 21][1 + t]) * flux_converter
                                self.FDM_min[-1][2, t, n]     = float(content[row_mass_IOD[n] + 5 + 21][1 + t]) * flux_converter

                        if row_mass_RG != 0:
                            for layer in range(len(self.layer_h[-1])):
                                for t in range(num_t):
                                    self.RM_smax[-1][layer, t, n]    = self.RM[-1][layer, t, n] + float(content[row_mass_RG[n] + 3 + layer + 4 + nlayers][1 + t]) * flux_converter
                                    self.RM_smin[-1][layer, t, n]    = self.RM[-1][layer, t, n] - float(content[row_mass_RG[n] + 3 + layer + 4 + nlayers][1 + t]) * flux_converter
                                    self.RM_max[-1][layer, t, n]     = float(content[row_mass_RG[n] + 3 + layer + (4 + nlayers)*2][1 + t]) * flux_converter
                                    self.RM_min[-1][layer, t, n]     = float(content[row_mass_RG[n] + 3 + layer + (4 + nlayers)*3][1 + t]) * flux_converter


                for m in range(ncomponents):
                    self.components[-1].append(content[row_Fis[m]][1])
                    self.solids[-1].append(content[row_Fis[m]][1])
                    for z in range(num_z):
                        for t in range(num_t):
                            self.Fi[-1][t,z,m] = float(content[row_Fis[m] + 4 + z][1 + t])
                    if self.output_type[-1] == 'Monte Carlo Simulation':
                        for z in range(num_z):
                            for t in range(num_t):
                                self.Fi_smax[-1][t,z,m] = self.Fi[-1][t,z,m] + float(content[row_Fis[m] + 4 + z + 4 + num_z][1 + t])
                                self.Fi_smin[-1][t,z,m] = self.Fi[-1][t,z,m] - float(content[row_Fis[m] + 4 + z + 4 + num_z][1 + t])
                                self.Fi_max[-1][t,z,m] = float(content[row_Fis[m] + 4 + z + (4 + num_z)*2 ][1 + t])
                                self.Fi_min[-1][t,z,m] = float(content[row_Fis[m] + 4 + z + (4 + num_z)*3 ][1 + t])
            else:

                self.output_type.append('Non-CapSim data')

                self.C.append([])
                self.F.append([])
                self.q.append([])
                self.W.append([])
                self.Cw.append([])
                self.M.append([])
                self.RM.append([])
                self.FDM.append([])
                self.Fi.append([])

                self.C_smax.append([])
                self.F_smax.append([])
                self.q_smax.append([])
                self.W_smax.append([])
                self.qm_smax.append([])
                self.Cw_smax.append([])
                self.M_smax.append([])
                self.RM_smax.append([])
                self.FDM_smax.append([])
                self.Fi_smax.append([])

                self.C_smin.append([])
                self.F_smin.append([])
                self.q_smin.append([])
                self.W_smin.append([])
                self.qm_smin.append([])
                self.Cw_smin.append([])
                self.M_smin.append([])
                self.RM_smin.append([])
                self.FDM_smin.append([])
                self.Fi_smin.append([])

                self.C_max.append([])
                self.F_max.append([])
                self.q_max.append([])
                self.W_max.append([])
                self.qm_max.append([])
                self.Cw_max.append([])
                self.M_max.append([])
                self.RM_max.append([])
                self.FDM_max.append([])
                self.Fi_max.append([])

                self.C_min.append([])
                self.F_min.append([])
                self.q_min.append([])
                self.W_min.append([])
                self.qm_min.append([])
                self.Cw_min.append([])
                self.M_min.append([])
                self.RM_min.append([])
                self.FDM_min.append([])
                self.Fi_min.append([])

                row_plots      = []
                while row < rows:
                    if len(content[row]) > 0:
                        if content[row][0] == 'Profile name': row_plots.append(row)
                    row = row + 1

                for row_plot in row_plots:
                    self.chemicals[-1].append(content[row_plot][1])
                    self.components[-1].append('N/A')
                    self.solids[-1].append('N/A')
                    self.times[-1].append([])
                    self.z[-1].append([])
                    self.C[-1].append([])
                    self.C_smax[-1].append([])
                    if content[row_plot+1][0].count('Independent variable')>0:
                        for i in range(len(content[row_plot+1])-1):
                            self.times[-1][-1].append(float(content[row_plot+1][i+1]))
                            self.z[-1][-1].append(float(content[row_plot+1][i+1]))
                    if content[row_plot+2][0].count('Plot variable')>0:
                        for i in range(len(content[row_plot+1])-1):
                            self.C[-1][-1].append(float(content[row_plot+2][i+1]))
                    if row_plot+3 < rows and len(content[row_plot+3])>0 and content[row_plot+3][0].count('standard deviation') > 0:
                        for i in range(len(content[row_plot+1])-1):
                            self.C_smax[-1][-1].append(float(content[row_plot+3][i+1]))

            profilename = self.dataname[:-4]
            p = len(profilename) - 1
            while profilename[p] <> '/':
                p = p - 1

            try: self.profiles.remove('None profile')
            except:pass
            self.profiles.append(profilename[p+1:])

            self.update_widgets()
            self.updategraph()


    def show_error(self, event = None):

        tkmb.showerror(title = 'Run Error', message = 'Unable to process the file')


    def update_widgets(self):

        try:
            self.editbutton.grid_forget()
            self.figurebutton.grid_forget()
            self.savebutton.grid_forget()
        except: pass

        if len(self.profiles) > 1:
            self.editbutton.grid(   row = 3, column = 0, columnspan = 7, padx = 40)
            self.editgbutton.grid(  row = 4, column = 0, columnspan = 7, padx = 40)
            self.figurebutton.grid( row = 5, column = 0, columnspan = 7, padx = 40)
            self.savebutton.grid(   row = 6, column = 0, columnspan = 7)

        elif len(self.profiles) == 1:

            self.editbutton.grid(   row = 3, column = 0, columnspan = 7, padx = 40)
            self.editgbutton.grid(  row = 4, column = 0, columnspan = 7, padx = 40)
            self.figurebutton.grid( row = 5, column = 0, columnspan = 7, padx = 40)
            self.savebutton.grid(   row = 6, column = 0, columnspan = 7)

            if self.output_type[0] == 'Monte Carlo Simulation': spatialplotnumber = 2
            elif self.output_type[0] == 'Non-CapSim data'     : spatialplotnumber = 1
            else:                                               spatialplotnumber = 6

            for i in range(spatialplotnumber):
                if self.output_type[0] == 'Non-CapSim data':
                    self.spatialplotdatas.append(PlotData(i, self.timeunit, self.lengthunit))
                    self.spatialplotdatas[-1].name      = self.profiles[0]
                    self.spatialplotdatas[-1].chemical  = self.chemicals[0][0]
                    self.spatialplotdatas[-1].value     = 0
                    self.spatialplotdatas[-1].type      = self.spatialplotdatas[-1].types[0]
                    self.spatialplotdatas[-1].component = self.components[0][0]
                    self.spatialplotdatas[-1].solid     = self.solids[0][0]
                    self.spatialplotdatas[-1].layer     = 'All layers'
                    self.spatialplotdatas[-1].style     = self.stylenames[i+4]
                    self.spatialplotdatas[-1].color     = self.colornames[i]
                    self.spatialplotdatas[-1].size      = 10
                    self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name + ' ' + self.spatialplotdatas[-1].chemical
                else:
                    self.spatialplotdatas.append(PlotData(i, self.timeunit, self.lengthunit))
                    self.spatialplotdatas[-1].name      = self.profiles[0]
                    self.spatialplotdatas[-1].chemical  = self.chemicals[0][0]
                    self.spatialplotdatas[-1].value     = round(self.times[0][-1]/(spatialplotnumber-1)*i, 2)
                    self.spatialplotdatas[-1].type      = self.spatialplotdatas[-1].types[0]
                    self.spatialplotdatas[-1].component = self.components[0][0]
                    self.spatialplotdatas[-1].solid     = 'Total solid'
                    self.spatialplotdatas[-1].layer     = 'All layers'
                    self.spatialplotdatas[-1].style     = 'Solid'
                    self.spatialplotdatas[-1].color     = self.colornames[i]
                    self.spatialplotdatas[-1].size      = 0.5
                    self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name + ' ' + self.spatialplotdatas[-1].chemical +' '+str(self.spatialplotdatas[-1].value) + ' ' + self.timeunit

            if self.output_type[0] == 'Non-CapSim data':
                self.timeplotdatas               = [PlotData(0, self.timeunit, self.lengthunit)]
                self.timeplotdatas[-1].name      = self.profiles[0]
                self.timeplotdatas[-1].chemical  = self.chemicals[0][0]
                self.timeplotdatas[-1].value     = 0
                self.timeplotdatas[-1].type      = self.timeplotdatas[-1].types[0]
                self.timeplotdatas[-1].component = self.components[0][0]
                self.timeplotdatas[-1].solid     = self.solids[0][0]
                self.timeplotdatas[-1].layer     = 'All layers'
                self.timeplotdatas[-1].style     = self.stylenames[4]
                self.timeplotdatas[-1].color     = self.colornames[0]
                self.timeplotdatas[-1].size      = 0.5
                self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name + ' ' + self.spatialplotdatas[-1].chemical
            else:
                self.timeplotdatas               = [PlotData(0, self.timeunit, self.lengthunit)]
                self.timeplotdatas[-1].name      = self.profiles[0]
                self.timeplotdatas[-1].chemical  = self.chemicals[0][0]
                self.timeplotdatas[-1].value     = 0
                self.timeplotdatas[-1].type      = self.timeplotdatas[-1].types[0]
                self.timeplotdatas[-1].component = self.components[0][0]
                self.timeplotdatas[-1].solid     = 'Total solid'
                self.timeplotdatas[-1].layer     = 'All layers'
                self.timeplotdatas[-1].style     = 'Solid'
                self.timeplotdatas[-1].color     = self.colornames[0]
                self.timeplotdatas[-1].size      = 0.5
                if self.dep[-1] == 1:
                    self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name + ' ' + self.spatialplotdatas[-1].chemical + ' ' +  str(self.timeplotdatas[-1].value) + ' ' + self.lengthunit + ' from ' + self.timeplotdatas[-1].type
                else:
                    self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name + ' ' + self.spatialplotdatas[-1].chemical + ' ' +  str(self.timeplotdatas[-1].value) + ' ' + self.lengthunit

            self.plot = []
            self.plot_smax = []
            self.plot_smin = []
            self.plot_max  = []
            self.plot_min  = []
            data_smax = self.C_smax
            data_smin = self.C_smin
            data_max  = self.C_max
            data_min  = self.C_min
            for n in range(len(self.spatialplotdatas)):
                filenum = 0
                chemnum = 0
                if self.output_type[0] == 'Non-CapSim data':
                    self.plot.append(self.C[filenum][chemnum])
                    self.plot_smax.append(self.C_smax[filenum][chemnum])
                    self.plot_smin.append([])
                    self.plot_max.append([])
                    self.plot_min.append([])
                else:
                    i = 0
                    while self.spatialplotdatas[n].value > self.times[0][i+1]:  i = i + 1
                    self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i + 1] - self.times[filenum][i], self.C[filenum][i,:,chemnum], self.C[filenum][i+1,:,chemnum]))
                    if self.output_type[0] == 'Monte Carlo Simulation':
                        self.plot_smax.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_smax[filenum][i,:,chemnum], data_smax[filenum][i+1,:,chemnum]))
                        self.plot_smin.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_smin[filenum][i,:,chemnum], data_smin[filenum][i+1,:,chemnum]))
                        self.plot_max.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_max[filenum][i,:,chemnum], data_max[filenum][i+1,:,chemnum]))
                        self.plot_min.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_min[filenum][i,:,chemnum], data_min[filenum][i+1,:,chemnum]))
                    else:
                        self.plot_smax.append([])
                        self.plot_smin.append([])
                        self.plot_max.append([])
                        self.plot_min.append([])

            if self.output_type[0] == 'Non-CapSim data':
                self.axislimit  = [0, 1.2 * max(self.plot[0]), min(self.z[0][0])*0.8, max(self.z[0][0])*1.2]
            elif self.output_type[0] == 'Monte Carlo Simulation':
                self.axislimit  = [0, 1.2 * max([self.plot_max[n].max() for n in range(len(self.spatialplotdatas))]), min(self.z[0]), max(self.z[0])]
            else:
                self.axislimit  = [0, 1.2 * max([self.plot[n].max() for n in range(len(self.spatialplotdatas))]), min(self.z[0]), max(self.z[0])]
            self.legendposition = 0
        self.sketch['Sketches'] = ['Hide sketch']
        for i in self.profiles: self.sketch['Sketches'].append(i)
        self.sketch['Sketch'] = self.sketch['Sketches'][0]

        self.master.geometry()
        self.master.center()


    def editplot(self, event = None):

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 2)
            self.top.make_window(PlotEditor(self.top, self.system, self.type, self.variable, self.profiles, self.chemicals, self.components, self.solids, self.dep, self.spatialplotdatas, self.timeplotdatas, self.z, self.times, self.layer_list, self.output_type))
            self.top.tk.mainloop()

            if self.top.window.cancelflag == 0:
                self.type     = self.top.window.type.get()
                self.variable = self.top.window.variable.get()

                for spatialplotdata in self.top.window.spatialplotdatas:  spatialplotdata.get_plotdata()
                for timeplotdata    in self.top.window.timeplotdatas:     timeplotdata.get_plotdata()

                self.spatialplotdatas = [spatialplotdata.copy() for spatialplotdata in self.top.window.spatialplotdatas]
                self.timeplotdatas    = [timeplotdata.copy()    for timeplotdata    in self.top.window.timeplotdatas]
                if self.type == self.types[0]:
                    axislimit3_temp = 0
                    axislimit2_temp = 0
                    for n in range(len(self.spatialplotdatas)):
                        filenum = self.profiles.index(self.spatialplotdatas[n].name)
                        chemnum = self.chemicals[filenum].index(self.spatialplotdatas[n].chemical)
                        if self.output_type[self.profiles.index(self.spatialplotdatas[n].name)] !=  'Non-CapSim data':
                            axislimit3_temp = max(axislimit3_temp, max(self.z[filenum]))
                            axislimit2_temp = min(axislimit2_temp, min(self.z[filenum]))
                        else:
                            axislimit3_temp = max(axislimit3_temp, max(self.z[filenum][chemnum]))
                            axislimit2_temp = min(axislimit2_temp, min(self.z[filenum][chemnum]))
                    self.axislimit[3] = axislimit3_temp
                    self.axislimit[2] = axislimit2_temp
                    self.plot      = []
                    self.plot_smax = []
                    self.plot_smin = []
                    self.plot_max  = []
                    self.plot_min  = []

                    if self.variable == 'Material fraction':
                        self.figuretexts['title']  = 'Solid material fraction Profiles'
                        self.figuretexts['xlabel'] = 'Volumetric fraction'
                    elif self.variable == 'Solid concentration':
                        self.figuretexts['title']  = 'Solid Concentration Profiles'
                        self.figuretexts['xlabel'] = 'Solid concentration ('+self.concunit[:-1]+'kg)'
                    elif self.variable == 'Concentration':
                        self.figuretexts['title']  = 'Porewater Concentration Profiles'
                        self.figuretexts['xlabel'] = 'Porewater concentration ('+self.concunit+')'
                    elif self.variable == 'Flux':
                        self.figuretexts['title']  = 'Flux Profiles'
                        self.figuretexts['xlabel'] = 'Flux (' + self.concunit[:-1] + self.lengthunit + u'\u00B2' + '/' + self.timeunit+')'
                    elif self.variable == 'Total concentration':
                        self.figuretexts['title']  = 'Porespace Concentration Profiles'
                        self.figuretexts['xlabel'] = 'Porespace Concentration ('+self.concunit+')'
                    self.figuretexts['ylabel'] = 'Depth ('+self.lengthunit+')'

                    for n in range(len(self.spatialplotdatas)):
                        filenum = self.profiles.index(self.spatialplotdatas[n].name)
                        chemnum = self.chemicals[filenum].index(self.spatialplotdatas[n].chemical)
                        if self.output_type[self.profiles.index(self.spatialplotdatas[n].name)] == 'Non-CapSim data':
                            self.plot.append(self.C[filenum][chemnum])
                            self.plot_smax.append(self.C_smax[filenum][chemnum])
                            self.plot_smin.append([])
                            self.plot_max.append([])
                            self.plot_min.append([])
                        else:
                            compnum  = self.components[filenum].index(self.spatialplotdatas[n].component)
                            solidnum = self.solids[filenum].index(self.spatialplotdatas[n].solid)
                            i = 0
                            while self.spatialplotdatas[n].value > self.times[filenum][i+1]: i = i + 1
                            if self.variable == 'Material fraction':
                                self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.Fi[filenum][i,:,compnum], self.Fi[filenum][i+1, :, compnum]))
                            elif self.variable == 'Solid concentration':
                                self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.q[filenum][i,:,chemnum,solidnum], self.q[filenum][i+1,:,chemnum,solidnum]))
                            else:
                                if self.variable == 'Concentration':       data = self.C
                                if self.variable == 'Flux':                data = self.F
                                if self.variable == 'Total concentration': data = self.W
                                self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data[filenum][i,:,chemnum], data[filenum][i+1,:,chemnum]))

                            if self.output_type[filenum] == 'Monte Carlo Simulation':
                                if self.variable == 'Material fraction':
                                    self.plot_smax.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.Fi_smax[filenum][i,:,compnum], self.Fi_smax[filenum][i+1, :, compnum]))
                                    self.plot_smin.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.Fi_smin[filenum][i,:,compnum], self.Fi_smin[filenum][i+1, :, compnum]))
                                    self.plot_max.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.Fi_max[filenum][i,:,compnum], self.Fi_max[filenum][i+1, :, compnum]))
                                    self.plot_min.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.Fi_min[filenum][i,:,compnum], self.Fi_min[filenum][i+1, :, compnum]))
                                elif self.variable == 'Solid concentration':
                                    self.plot_smax.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.q_smax[filenum][i,:,chemnum,solidnum], self.q_smax[filenum][i+1,:,chemnum,solidnum]))
                                    self.plot_smin.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.q_smin[filenum][i,:,chemnum,solidnum], self.q_smin[filenum][i+1,:,chemnum,solidnum]))
                                    self.plot_max.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.q_max[filenum][i,:,chemnum,solidnum], self.q_max[filenum][i+1,:,chemnum,solidnum]))
                                    self.plot_min.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], self.q_min[filenum][i,:,chemnum,solidnum], self.q_min[filenum][i+1,:,chemnum,solidnum]))
                                else:
                                    if self.variable == 'Concentration':
                                        data_smax = self.C_smax
                                        data_smin = self.C_smin
                                        data_max  = self.C_max
                                        data_min  = self.C_min
                                    if self.variable == 'Flux':
                                        data_smax = self.F_smax
                                        data_smin = self.F_smin
                                        data_max  = self.F_max
                                        data_min  = self.F_min
                                    if self.variable == 'Total concentration':
                                        data_smax = self.W_smax
                                        data_smin = self.W_smin
                                        data_max  = self.W_max
                                        data_min  = self.W_min

                                    self.plot_smax.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_smax[filenum][i,:,chemnum], data_smax[filenum][i+1,:,chemnum]))
                                    self.plot_smin.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_smin[filenum][i,:,chemnum], data_smin[filenum][i+1,:,chemnum]))
                                    self.plot_max.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_max[filenum][i,:,chemnum], data_max[filenum][i+1,:,chemnum]))
                                    self.plot_min.append(time_interpolate(self.spatialplotdatas[n].value, self.times[filenum][i], self.times[filenum][i+1] - self.times[filenum][i], data_min[filenum][i,:,chemnum], data_min[filenum][i+1,:,chemnum]))
                            else:
                                self.plot_smax.append([])
                                self.plot_smin.append([])
                                self.plot_max.append([])
                                self.plot_min.append([])

                    if self.variable == 'Flux':
                        Fmax = 0
                        Fmin = 0
                        for n in range(len(self.spatialplotdatas)):
                            filenum = self.profiles.index(self.spatialplotdatas[n].name)
                            if self.output_type[filenum] == 'Monte Carlo Simulation':
                                max_temp = [self.plot[n].max(), self.plot_max[n].max(), self.plot_min[n].max()]
                                min_temp = [self.plot[n].min(), self.plot_max[n].min(), self.plot_min[n].min()]
                            else:
                                max_temp = max(self.plot[n])
                                min_temp = min(self.plot[n])

                            if max_temp < 0:     Fmax = max(Fmax, 0.8 * max_temp)
                            else:                Fmax = max(Fmax, 1.2 * max_temp)

                            if min_temp < 0:     Fmin = min(Fmin, 1.2 * min_temp)
                            else:                Fmin = min(Fmin, 0.8 * min_temp)


                        self.axislimit[0]  = Fmin
                        self.axislimit[1]  = Fmax
                    else:
                        self.axislimit[0]  = 0
                        Cmax = 0
                        for n in range(len(self.spatialplotdatas)):
                            filenum = self.profiles.index(self.spatialplotdatas[n].name)
                            if self.output_type[filenum] == 'Monte Carlo Simulation':
                                if Cmax < max(self.plot_max[n]): Cmax = max(self.plot_max[n])
                            else:
                                if Cmax < max(self.plot[n]): Cmax = max(self.plot[n])
                        self.axislimit[1]  = 1.2  * Cmax

                else:
                    self.axislimit[0] = 0
                    axislimit1_temp = 0
                    for n in range(len(self.timeplotdatas)):
                        filenum  = self.profiles.index(self.timeplotdatas[n].name)
                        chemnum  = self.chemicals[filenum].index(self.timeplotdatas[n].chemical)
                        if self.output_type[self.profiles.index(self.timeplotdatas[n].name)] != 'Non-CapSim data':
                            axislimit1_temp = max(axislimit1_temp, max(self.times[filenum]))
                        else:
                            axislimit1_temp = max(axislimit1_temp, max(self.times[filenum][chemnum])*1.2)
                    self.axislimit[1] = axislimit1_temp
                    self.plotinterest         = []
                    self.plotinterest_smax    = []
                    self.plotinterest_smin    = []
                    self.plotinterest_max     = []
                    self.plotinterest_min     = []

                    if self.variable == 'Material fraction':
                        self.figuretexts['title']  = 'Solid material fraction Time Profiles'
                        self.figuretexts['ylabel'] = 'Volumetric fraction ('+self.concunit+')'
                    elif self.variable == 'Solid concentration':
                        self.figuretexts['title']  = 'Solid Concentration Time Profiles'
                        self.figuretexts['ylabel'] = 'Solid concentration ('+self.concunit[:-1]+'kg)'
                    elif self.variable == 'Water concentration':
                        self.figuretexts['title']  = 'Overlying Water Concentration Time Profiles'
                        self.figuretexts['ylabel'] = 'Overlying water concentration ('+self.concunit+')'
                    elif self.variable == 'Mass':
                        self.figuretexts['title']  = 'Accumulated Mass Time Profiles'
                        self.figuretexts['ylabel'] = 'Mass(' + self.concunit[:-1] + 'cm'+u'\u00B2'+')'
                    elif self.variable == 'Concentration':
                        self.figuretexts['title']  = 'Porewater Concentration Time Profiles'
                        self.figuretexts['ylabel'] = 'Porewater concentration ('+self.concunit+')'
                    elif self.variable == 'Flux':
                        self.figuretexts['title']  = 'Flux Time Profiles'
                        self.figuretexts['ylabel'] = 'Flux (' + self.concunit[:-1]+ self.lengthunit + u'\u00B2' + '/' + self.timeunit+')'
                    elif self.variable == 'Total concentration':
                        self.figuretexts['title']  = 'Porespace Concentration Time Profiles'
                        self.figuretexts['ylabel'] = 'Porespace Concentration ('+self.concunit+')'
                    self.figuretexts['xlabel'] = 'Time (' + self.timeunit + ')'

                    plotmax = 0
                    for n in range(len(self.timeplotdatas)):
                        filenum  = self.profiles.index(self.timeplotdatas[n].name)
                        chemnum  = self.chemicals[filenum].index(self.timeplotdatas[n].chemical)
                        if self.output_type[self.profiles.index(self.timeplotdatas[n].name)] == 'Non-CapSim data':
                            self.plotinterest.append(self.C[filenum][chemnum])
                            self.plotinterest_smax.append(self.C_smax[filenum][chemnum])
                            self.plotinterest_smin.append([])
                            self.plotinterest_max.append([])
                            self.plotinterest_min.append([])
                        else:
                            self.plotinterest.append([])
                            compnum  = self.components[filenum].index(self.timeplotdatas[n].component)
                            solidnum = self.solids[filenum].index(self.timeplotdatas[n].solid)
                            if self.variable == 'Material fraction':
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.Fi[filenum][i,:, compnum]))
                                    else:
                                        self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.Fi[filenum][i,:,compnum]))
                            elif self.variable == 'Solid concentration':
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.q[filenum][i,:,chemnum,solidnum]))
                                    else:
                                        self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.q[filenum][i,:,chemnum,compnum]))
                            elif self.variable == 'Water concentration':
                                for i in range (len(self.times[filenum])):
                                    self.plotinterest[-1].append(self.Cw[filenum][i, chemnum])

                            elif self.variable == 'Mass':
                                for i in range (len(self.times[filenum])):
                                    if self.timeplotdatas[n].layer == 'All layers':
                                        self.plotinterest[-1].append(0)
                                        for j in range(len(self.layer_h[filenum])):
                                            self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.M[filenum][j,i,chemnum]
                                    elif self.timeplotdatas[n].layer == 'Top accumulation':
                                        self.plotinterest[-1].append(self.FDM[filenum][0,i,chemnum])
                                    elif self.timeplotdatas[n].layer == 'Bottom accumulation':
                                        self.plotinterest[-1].append(self.FDM[filenum][1,i,chemnum])
                                    elif self.timeplotdatas[n].layer == 'Deposition':
                                        self.plotinterest[-1].append(self.FDM[filenum][2,i,chemnum])
                                    elif self.timeplotdatas[n].layer == 'Reaction/Generation':
                                        self.plotinterest[-1].append(0)
                                        for j in range(len(self.layer_h[filenum])):
                                            self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.RM[filenum][j,i,chemnum]
                                    elif self.timeplotdatas[n].layer == 'Universal':
                                        self.plotinterest[-1].append(0)
                                        for j in range(len(self.layer_h[filenum])):
                                            self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.M[filenum][j,i,chemnum]
                                            self.plotinterest[-1][-1] = self.plotinterest[-1][-1] - self.RM[filenum][j,i,chemnum]
                                        self.plotinterest[-1][-1] = self.plotinterest[-1][-1]+self.FDM[filenum][0,i,chemnum]-self.FDM[filenum][1,i,chemnum]-self.FDM[filenum][2,i,chemnum]
                                    else:
                                        self.plotinterest[-1].append(self.M[filenum][self.layer_list[filenum].index(self.timeplotdatas[n].layer),i,chemnum])
                            else:
                                if self.variable == 'Concentration':       data = self.C
                                if self.variable == 'Flux':                data = self.F
                                if self.variable == 'Total concentration': data = self.W
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], data[filenum][i,:,chemnum]))
                                    else:
                                        self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], data[filenum][i,:,chemnum]))

                        self.plotinterest_smax.append([])
                        self.plotinterest_smin.append([])
                        self.plotinterest_max.append([])
                        self.plotinterest_min.append([])
                        if self.output_type[filenum] == 'Monte Carlo Simulation':
                            if self.variable == 'Material fraction':
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest_smax[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.Fi_smax[filenum][i,:, compnum]))
                                        self.plotinterest_smin[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.Fi_smin[filenum][i,:, compnum]))
                                        self.plotinterest_max[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.Fi_max[filenum][i,:, compnum]))
                                        self.plotinterest_min[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.Fi_min[filenum][i,:, compnum]))
                                    else:
                                        self.plotinterest_smax[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.Fi_smax[filenum][i,:,compnum]))
                                        self.plotinterest_smin[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.Fi_smin[filenum][i,:,compnum]))
                                        self.plotinterest_max[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.Fi_max[filenum][i,:,compnum]))
                                        self.plotinterest_min[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.Fi_min[filenum][i,:,compnum]))

                            elif self.variable == 'Solid concentration':
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest_smax[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.q_smax[filenum][i,:,chemnum,solidnum]))
                                        self.plotinterest_smin[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.q_smin[filenum][i,:,chemnum,solidnum]))
                                        self.plotinterest_max[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.q_max[filenum][i,:,chemnum,solidnum]))
                                        self.plotinterest_min[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], self.q_min[filenum][i,:,chemnum,solidnum]))
                                    else:
                                        self.plotinterest_smax[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.q_smax[filenum][i,:,chemnum,compnum]))
                                        self.plotinterest_smin[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.q_smin[filenum][i,:,chemnum,compnum]))
                                        self.plotinterest_max[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.q_max[filenum][i,:,chemnum,compnum]))
                                        self.plotinterest_min[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], self.q_min[filenum][i,:,chemnum,compnum]))

                            elif self.variable == 'Water concentration':
                                for i in range (len(self.times[filenum])):
                                    self.plotinterest_smax[-1].append(self.Cw_smax[filenum][i, chemnum])
                                    self.plotinterest_smin[-1].append(self.Cw_smin[filenum][i, chemnum])
                                    self.plotinterest_max[-1].append(self.Cw_max[filenum][i, chemnum])
                                    self.plotinterest_min[-1].append(self.Cw_min[filenum][i, chemnum])

                            elif self.variable == 'Mass':
                                for i in range (len(self.times[filenum])):
                                    if self.timeplotdatas[n].layer == 'All layers':
                                        self.plotinterest_smax[-1].append(0)
                                        self.plotinterest_smin[-1].append(0)
                                        self.plotinterest_max[-1].append(0)
                                        self.plotinterest_min[-1].append(0)
                                    elif self.timeplotdatas[n].layer == 'Top accumulation':
                                        self.plotinterest_smax[-1][-1] = self.FDM_smax[filenum][0,i,chemnum]
                                        self.plotinterest_smin[-1][-1] = self.FDM_smin[filenum][0,i,chemnum]
                                        self.plotinterest_max[-1][-1]  = self.FDM_max[filenum][0,i,chemnum]
                                        self.plotinterest_min[-1][-1]  = self.FDM_min[filenum][0,i,chemnum]
                                    elif self.timeplotdatas[n].layer == 'Bottom accumulation':
                                        self.plotinterest_smax[-1][-1] = self.FDM_smax[filenum][1,i,chemnum]
                                        self.plotinterest_smin[-1][-1] = self.FDM_smin[filenum][1,i,chemnum]
                                        self.plotinterest_max[-1][-1]  = self.FDM_max[filenum][1,i,chemnum]
                                        self.plotinterest_min[-1][-1]  = self.FDM_min[filenum][1,i,chemnum]
                                    elif self.timeplotdatas[n].layer == 'Deposition':
                                        self.plotinterest_smax[-1][-1] = self.FDM_smax[filenum][2,i,chemnum]
                                        self.plotinterest_smin[-1][-1] = self.FDM_smin[filenum][2,i,chemnum]
                                        self.plotinterest_max[-1][-1]  = self.FDM_max[filenum][2,i,chemnum]
                                        self.plotinterest_min[-1][-1]  = self.FDM_min[filenum][2,i,chemnum]
                                    elif self.timeplotdatas[n].layer == 'Reaction/Generation':
                                        self.plotinterest_smax[-1].append(0)
                                        self.plotinterest_smin[-1].append(0)
                                        self.plotinterest_max[-1].append(0)
                                        self.plotinterest_min[-1].append(0)
                                    elif self.timeplotdatas[n].layer == 'Universal':
                                        self.plotinterest_smax[-1].append(0)
                                        self.plotinterest_smin[-1].append(0)
                                        self.plotinterest_max[-1].append(0)
                                        self.plotinterest_min[-1].append(0)
                                    else:
                                        self.plotinterest_smax[-1].append(self.M_smax[filenum][self.layer_list[filenum].index(self.timeplotdatas[n].layer),i,chemnum])
                                        self.plotinterest_smin[-1].append(self.M_smin[filenum][self.layer_list[filenum].index(self.timeplotdatas[n].layer),i,chemnum])
                                        self.plotinterest_max[-1].append(self.M_max[filenum][self.layer_list[filenum].index(self.timeplotdatas[n].layer),i,chemnum])
                                        self.plotinterest_min[-1].append(self.M_min[filenum][self.layer_list[filenum].index(self.timeplotdatas[n].layer),i,chemnum])
                            else:
                                if self.variable == 'Concentration':
                                    data_smax = self.C_smax
                                    data_smin = self.C_smin
                                    data_max = self.C_max
                                    data_min = self.C_min
                                if self.variable == 'Flux':
                                    data_smax = self.F_smax
                                    data_smin = self.F_smin
                                    data_max = self.F_max
                                    data_min = self.F_min
                                if self.variable == 'Total concentration':
                                    data_smax = self.W_smax
                                    data_smin = self.W_smin
                                    data_max = self.W_max
                                    data_min = self.W_min
                                for i in range (len(self.times[filenum])):
                                    if self.dep[filenum] == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                        self.plotinterest_smax[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], data_smax[filenum][i,:,chemnum]))
                                        self.plotinterest_smin[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], data_smin[filenum][i,:,chemnum]))
                                        self.plotinterest_max[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], data_max[filenum][i,:,chemnum]))
                                        self.plotinterest_min[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep[filenum]*self.times[filenum][i], (self.z[filenum][0]-self.z[filenum][1]), self.z[filenum], data_min[filenum][i,:,chemnum]))
                                    else:
                                        self.plotinterest_smax[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], data_smax[filenum][i,:,chemnum]))
                                        self.plotinterest_smin[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], data_smin[filenum][i,:,chemnum]))
                                        self.plotinterest_max[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], data_max[filenum][i,:,chemnum]))
                                        self.plotinterest_min[-1].append(interp(self.timeplotdatas[n].value, self.z[filenum], data_min[filenum][i,:,chemnum]))

                    if self.variable == 'Flux':
                        Fmax = 0
                        Fmin = 0
                        for n in range(len(self.timeplotdatas)):
                            filenum  = self.profiles.index(self.timeplotdatas[n].name)
                            if self.output_type[filenum] == 'Monte Carlo Simulation':
                                max_temp = max([max(self.plotinterest[n]), max(self.plotinterest_max[n]), max(self.plotinterest_min[n])])
                                min_temp = max([max(self.plotinterest[n]), max(self.plotinterest_max[n]), max(self.plotinterest_min[n])])
                            else:
                                max_temp = max(self.plotinterest[n])
                                min_temp = max(self.plotinterest[n])

                            if max_temp < 0:     Fmax = max(Fmax, 0.8 * max_temp)
                            else:                Fmax = max(Fmax, 1.2 * max_temp)

                            if min_temp < 0:     Fmin = min(Fmin, 1.2 * min_temp)
                            else:                Fmin = min(Fmin, 0.8 * min_temp)

                        self.axislimit[2]  = Fmin
                        self.axislimit[3]  = Fmax

                    else:
                        self.axislimit[0]  = 0
                        for n in range(len(self.timeplotdatas)):
                            filenum  = self.profiles.index(self.timeplotdatas[n].name)
                            if self.output_type[filenum] == 'Monte Carlo Simulation':
                                if max(self.plotinterest_max[n]) > plotmax: plotmax = max(self.plotinterest_max[n])
                            else:
                                if max(self.plotinterest[n]) > plotmax:   plotmax = max(self.plotinterest[n])

                            self.axislimit[2]  = 0
                            self.axislimit[3]  = plotmax*1.2

                self.updategraph()

            if self.top is not None:
                self.top.destroy()
                self.top = None

        elif self.top is not None:
            tkmb.showerror(title = self.system.version, message = 'Please close the existing parameter input window first.')
            self.top.tk.focus()

    def editprofile(self, event = None):

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 2)
            if self.type == self.types[0]:  plotdatas = self.spatialplotdatas
            else:                           plotdatas = self.timeplotdatas

            self.top.make_window(GraphEditor(self.top, self.system, self.type, self.variable, plotdatas))
            self.top.tk.mainloop()

            if self.top.window.cancelflag == 0:
                for plotdata in self.top.window.plotdatas:  plotdata.get_graphicdata()
                if self.type == self.types[0]:  self.spatialplotdatas = [plotdata.copy() for plotdata in self.top.window.plotdatas]
                else:                           self.timeplotdatas    = [plotdata.copy() for plotdata in self.top.window.plotdatas]
                self.updategraph()

            if self.top is not None:
                self.top.destroy()
                self.top = None

        elif self.top is not None:
            tkmb.showerror(title = self.system.version, message = 'Please close the existing parameter input window first.')
            self.top.tk.focus()

    def editfigure(self, event = None):

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 2)
            self.top.make_window(FigureEditor(self.top, self.system, self.type, self.variable, self.sketch, self.axislimit, self.figuresize, self.legendposition, self.fontsize, self.figuretexts, self.MCplots))

            self.top.tk.mainloop()

            if self.top.window.cancelflag == 0:

                self.sketch['Sketch']       = self.top.window.sketch.get()
                self.sketch['Label']        = self.top.window.sketchlabel.get()
                self.sketch['Time']         = self.top.window.sketchline.get()
                self.sketch['Depth']        = self.top.window.locationline.get()

                self.legendposition         = self.top.window.legendpositions.index(self.top.window.legendposition.get())
                self.axislimit[0]           = self.top.window.xaxismin.get()
                self.axislimit[1]           = self.top.window.xaxismax.get()
                self.axislimit[2]           = self.top.window.yaxismin.get()
                self.axislimit[3]           = self.top.window.yaxismax.get()

                self.figuresize[0]          = self.top.window.figurewidth.get()
                self.figuresize[1]          = self.top.window.figureheight.get()

                self.figuretexts['Flag']    = self.top.window.textflag.get()
                self.figuretexts['title']   = self.top.window.title.get()
                self.figuretexts['xlabel']  = self.top.window.xlabel.get()
                self.figuretexts['ylabel' ] = self.top.window.ylabel.get()

                self.fontsize['Flag']       = self.top.window.fontsize.get()
                self.fontsize['Title']      = self.top.window.titlefontsize.get()
                self.fontsize['Axis'][0]    = self.top.window.xaxisfontsize.get()
                self.fontsize['Axis'][1]    = self.top.window.yaxisfontsize.get()
                self.fontsize['Label'][0]   = self.top.window.xlabelfontsize.get()
                self.fontsize['Label'][1]   = self.top.window.ylabelfontsize.get()
                self.fontsize['Legend']     = self.top.window.legendfontsize.get()

                self.MCplots['Flag']        = self.top.window.MCflag.get()
                self.MCplots['AVG']         = self.top.window.MCavg.get()
                self.MCplots['SD']          = self.top.window.MCsd.get()
                self.MCplots['Max/Min' ]    = self.top.window.MCmaxmin.get()

                self.updategraph()

            if self.top is not None:
                self.top.destroy()
                self.top = None

        elif self.top is not None:
            tkmb.showerror(title = self.system.version, message = 'Please close the existing parameter input window first.')
            self.top.tk.focus()

    def updategraph(self, event = None):

        try:
            self.graphlabel.grid_forget()
        except: pass

        self.make_graphs()
        
        self.graph      = Image.open(Filepath + r'/output/%s.png' % self.filename)
        self.image      = ImageTk.PhotoImage(self.graph)
        self.graphlabel = Label(self.frame, image = self.image, bg = self.color)

        self.graphlabel.grid( row = 0, column = 0, columnspan = 7, padx = 40)
            
    def savegraph(self, event = None):

        self.graphname = tkfd.asksaveasfilename(initialdir=self.graphpath,  title='Please define the directory to save the graph', filetypes = [('PNG files', '*.png')])
        self.outputgraphs.savefig(self.graphname)

        #update the default directory
        try:
            i = 1
            while self.graphname[-i] != '/':
                i = i + 1
            self.graphpath = self.graphname[0:-i+1]
        except: pass
               

    def make_graphs(self):

        """Makes graphs for the postprocessing window."""
        self.family         = 'Calibri'

        plt.rcParams['mathtext.default']      = 'regular'
        plt.rcParams['axes.facecolor']        = self.mplrgb
        plt.rcParams['savefig.edgecolor']     = self.mplrgb
        plt.rcParams['savefig.facecolor']     = self.mplrgb
        plt.rcParams['axes.formatter.limits'] = [-4,4]
        plt.rcParams['xtick.labelsize']       = self.fontsize['Axis'][0]
        plt.rcParams['ytick.labelsize']       = self.fontsize['Axis'][1]

        if self.sketch['Sketch'] == self.sketch['Sketches'][0]:
            self.outputgraphs = plt.figure(figsize=(self.figuresize[0]/100, self.figuresize[1]/100))
        else:
            self.outputgraphs = plt.figure(figsize=(self.figuresize[0]/100, self.figuresize[1]/100))
            self.make_sketch_image()

        if len(self.profiles) >= 1:
            if self.type == self.types[0]:
                self.make_profiles(self.plot, self.plot_smax, self.plot_smin, self.plot_max, self.plot_min)
            elif self.type == self.types[1]:
                self.make_time_profiles(self.plotinterest, self.plotinterest_smax, self.plotinterest_smin, self.plotinterest_max, self.plotinterest_min)
        else:
            self.make_blank_graphs()

        if 0.15 * self.figuresize[0] < 90:   leftspace    = 90./self.figuresize[0]
        else:                                leftspace    = 0.15
        if 0.05 * self.figuresize[0] < 30:   rightspace   = 1. - 30./self.figuresize[0]
        else:                                rightspace   = 0.95
        if 0.35 * self.figuresize[0] < 210:  widthspace   = 210./self.figuresize[0]
        else:                                widthspace   = 0.35
        if 0.1  * self.figuresize[1] < 45:   bottomspace  = 45./self.figuresize[1]
        else:                                bottomspace  = 0.1
        if 0.07 * self.figuresize[1] < 31.5: topspace     = 1. - 45./self.figuresize[1]
        else:                                topspace     = 0.93
        if 0.33 * self.figuresize[0] < 121.5:heightspace  = 121.5/self.figuresize[0]
        else:                                heightspace  = 0.33

        plt.subplots_adjust(left = leftspace, right = rightspace, bottom =bottomspace, top = topspace,  wspace = widthspace, hspace = heightspace)
        self.outputgraphs.savefig(Filepath + r'/output/%s.png' % self.filename)

    def make_blank_graphs(self):
        """Makes plots of concentration profiles over time."""

        blankprofiles  = self.outputgraphs.add_subplot(111)

        blankprofiles.text(0.35,0.5,'Please import profiles',fontsize=15)

        blankprofiles.set_xlabel('Concentration ('+self.concunit+')', fontsize = self.fontsize['Label'][0], family = self.family)
        blankprofiles.set_ylabel('Depth ('+self.lengthunit+')', fontsize = self.fontsize['Label'][1], family =self.family)

    def make_sketch_image(self):
        """Makes plots of concentration profiles over time."""

        Sketchprofiles  = self.outputgraphs.add_subplot(self.gs[0])

        pnum    = self.profiles.index(self.sketch['Sketch'])
        Htotal  = self.z[pnum][-1] - self.z[pnum][0]

        sfont = {9:6, 10:7, 11:7.5, 12:8, 13:9}

        H = 0
        color_step = round(0.5/len(self.layer_h[pnum]),2)
        Blocks = []

        for i in range(len(self.layer_h[pnum])):
            if self.dep[pnum] == 0 or i > 0:
                Blocks.append(patches.Rectangle((0.0, H), 0.99, self.layer_h[pnum][i], facecolor = str( 0.75 - color_step * i), edgecolor = 'k', linestyle = 'solid', linewidth = '1'))
                Sketchprofiles.add_patch(Blocks[-1])
                if self.sketch['Label'] == 1:
                    rx, ry = Blocks[-1].get_xy()
                    cx = rx + Blocks[-1].get_width()/2.0
                    cy = ry + Blocks[-1].get_height()/2.0
                    if Blocks[-1].get_height()/Htotal*self.figuresize[1] >= 14:
                        Sketchprofiles.annotate(fill(self.layer_type[pnum][i], 10), (cx, cy), color='w', fontsize = 10, ha='center', va='center')
                    elif Blocks[-1].get_height()/Htotal*self.figuresize[1] >= 9:
                        Sketchprofiles.annotate(fill(self.layer_type[pnum][i], 10), (cx, cy), color='w', fontsize = sfont[int(Blocks[-1].get_height()/Htotal*self.figuresize[1])], ha='center', va='center')
                H = H + self.layer_h[pnum][i]

        if self.dep[pnum] == 1:
            Blockdep = patches.Polygon(array([[0.995, self.z[pnum][0]], [0.0, 0.0], [0.995, 0.0]]), facecolor = '0.75', edgecolor = 'k', linestyle = 'solid', linewidth = '1')
            Sketchprofiles.add_patch(Blockdep)
            if self.sketch['Label'] == 1:
                rx, ry = 0.0, 0.0
                cx = rx + 0.995/2.0
                cy = ry + self.z[pnum][0]/2.0
                if abs(self.z[pnum][0])/Htotal*self.figuresize[1] >= 14:
                    Sketchprofiles.annotate(fill(self.layer_type[pnum][0], 10), (cx, cy), color='k', fontsize = 10, ha='center', va='center')
                elif abs(self.z[pnum][0])/Htotal*self.figuresize[1] >= 9:
                    Sketchprofiles.annotate(fill(self.layer_type[pnum][0], 10), (cx, cy), color='k', fontsize = sfont[int(abs(self.z[pnum][0])/Htotal*self.figuresize[1])], ha='center', va='center')

        if self.dep[pnum] == 1 or self.sketch['Time'] == 1:
            Sketchprofiles.annotate('0', (0, H),                                         color='k', fontsize = 10, ha='center', va='top')
            Sketchprofiles.annotate(fill(str(self.times[pnum][-1]) + self.timeunit, 6), (1.0, H), color='k', fontsize = 10, ha='center', va='top')

        if self.type == self.types[0]:
            if self.sketch['Time'] == 1:
                if self.dep[pnum] == 1:
                    for plotdata in self.spatialplotdatas:
                        Sketchprofiles.plot( [plotdata.value/self.times[pnum][-1], plotdata.value/self.times[pnum][-1]],[0-plotdata.value*self.Vdep[pnum], self.z[pnum][-1]])
                else:
                    for plotdata in self.spatialplotdatas:
                        Sketchprofiles.plot( [plotdata.value/self.times[pnum][-1], plotdata.value/self.times[pnum][-1]], [0, self.z[pnum][-1]])

        if self.type == self.types[1]:
            if self.sketch['Depth'] == 1:
                for plotdata in self.timeplotdatas:
                    if plotdata.type == plotdata.types[0]:
                        Sketchprofiles.plot([0.0, 1.0], [plotdata.value, plotdata.value])
                    else:
                        Sketchprofiles.plot([0.0, 1.0], [plotdata.value, plotdata.value + self.z[pnum][0]])

        if self.hbio[pnum] <> 0 :
            Sketchprofiles.plot([0.0, 1.0], [self.hbio[pnum], self.hbio[pnum] + self.z[pnum][0]], linestyle = ':', color = 'k')

        Sketchprofiles.axis('off')
        if self.type == self.types[0]:
            Sketchprofiles.set_ylim([self.axislimit[3], self.axislimit[2]])
        else:
            Sketchprofiles.set_ylim([max(self.z[pnum]), min(self.z[pnum])])


    def make_profiles(self, plots, plots_smax, plots_smin, plots_max, plots_min):
        """Makes plots of concentration profiles over time."""

        if self.sketch['Sketch'] == self.sketch['Sketches'][0]: Graph  = self.outputgraphs.add_subplot(111)
        else:                                                   Graph  = self.outputgraphs.add_subplot(self.gs[1])

        profiles  = []
        tlegend   = []

        for n in range(len(self.spatialplotdatas)):
            filenum = self.profiles.index(self.spatialplotdatas[n].name)
            color     = self.colors[self.colornames.index(self.spatialplotdatas[n].color)]
            style     = self.styles[self.stylenames.index(self.spatialplotdatas[n].style)]
            linewidth = self.spatialplotdatas[n].size

            if self.output_type[filenum] == 'Non-CapSim data':
                chemnum = self.chemicals[filenum].index(self.spatialplotdatas[n].chemical)
                if self.stylenames.index(self.spatialplotdatas[n].style) < 4:
                    p, = Graph.plot(plots[n], self.z[filenum][chemnum], color = color, linestyle = style, linewidth = linewidth)
                    if len(plots_smax[n]) >0: Graph.errorbar(plots[n], self.z[filenum][chemnum], xerr = plots_smax[n], fmt = style, ecolor = color, elinewidth = 0.5)
                else:
                    p, = Graph.plot(plots[n], self.z[filenum][chemnum], color = color, linestyle = 'None', marker = style,  markersize = linewidth)
                    if len(plots_smax[n]) >0: Graph.errorbar(plots[n], self.z[filenum][chemnum], xerr = plots_smax[n], fmt = style, ecolor = color, elinewidth = 0.5)
            else:
                if self.dep[filenum] == 1:
                    time = self.spatialplotdatas[n].value
                    dep = -self.Vdep[filenum] * time
                    i = 0
                    while round(self.z[filenum][i] + 0.5*(self.z[filenum][1]-self.z[filenum][0]), 8) < round(dep, 8)  : i = i + 1
                    if round(-dep, 8) < 1.5 * (self.z[filenum][1]-self.z[filenum][0]) and round(-dep, 8) >= 0.5 * (self.z[filenum][1]-self.z[filenum][0]):
                        Cs = plots[n][i+2:]
                        zs = self.z[filenum][i+2:]
                        if self.output_type[filenum] == 'Monte Carlo Simulation':
                            Cs_smax = plots_smax[n][i+2:]
                            Cs_smin = plots_smin[n][i+2:]
                            Cs_max  = plots_max[n][i+2:]
                            Cs_min  = plots_min[n][i+2:]
                    else:
                        Cs = plots[n][i+1:]
                        zs = self.z[filenum][i+1:]
                        if self.output_type[filenum] == 'Monte Carlo Simulation':
                            Cs_smax = plots_smax[n][i+1:]
                            Cs_smin = plots_smin[n][i+1:]
                            Cs_max  = plots_max[n][i+1:]
                            Cs_min  = plots_min[n][i+1:]
                    zs[0] = -self.Vdep[filenum] * time

                else:
                    Cs = plots[n][:]
                    zs = self.z[filenum][:]
                    if self.output_type[filenum] == 'Monte Carlo Simulation':
                        Cs_smax = plots_smax[n][:]
                        Cs_smin = plots_smin[n][:]
                        Cs_max  = plots_max[n][:]
                        Cs_min  = plots_min[n][:]

                if self.output_type[filenum] == 'Monte Carlo Simulation':
                    Graph.fill_betweenx(zs, Cs_smax, Cs_smin, alpha = self.MCplots['SD'], facecolor = color, edgecolor = color)
                    Graph.fill_betweenx(zs, Cs_max, Cs_min,   alpha = self.MCplots['Max/Min'], facecolor = color, edgecolor = color)
                    if self.stylenames.index(self.spatialplotdatas[n].style) < 4:
                        p, = Graph.plot(Cs, zs, color = color, alpha = self.MCplots['AVG'], linestyle = style, linewidth = linewidth)
                    else:
                        p, = Graph.plot(Cs, zs, color = color, alpha = self.MCplots['AVG'], linestyle = 'None', marker = style,  markersize = linewidth)
                else:
                    if self.stylenames.index(self.spatialplotdatas[n].style) < 4:
                        p, = Graph.plot(Cs, zs, color = color, linestyle = style, linewidth = linewidth)
                    else:
                        p, = Graph.plot(Cs, zs, color = color, linestyle = 'None', marker = style,  markersize = linewidth)

            if self.spatialplotdatas[n].legend != '':
                profiles.append(p)
                tlegend.append(self.spatialplotdatas[n].legend)

        Graph.set_xlim([self.axislimit[0], self.axislimit[1]])
        Graph.set_ylim([self.axislimit[3], self.axislimit[2]])
        Graph.set_xlabel(self.figuretexts['xlabel'], fontsize = self.fontsize['Label'][0], family = self.family)
        Graph.set_ylabel(self.figuretexts['ylabel' ], fontsize = self.fontsize['Label'][0], family =self.family)
        Graph.set_title(self.figuretexts['title'], fontsize = self.fontsize['Title'], family = self.family)
        if len(tlegend) > 0 and self.legendposition < 11:
            Graphleg = Graph.legend(profiles, tlegend, loc = self.legendposition, labelspacing = 0, borderpad = 0.3, handlelength = 1.5, handletextpad = 0.1, fancybox = 0)
            Graphleg.legendPatch.set_fc(plt.rcParams['axes.facecolor'])
            for caption in Graphleg.get_texts():
                caption.set_fontsize(self.fontsize['Legend'])
                caption.set_fontname(self.family)

    def make_time_profiles(self, plots, plots_smax, plots_smin, plots_max, plots_min):
        """Makes a graph of the pore water concentration at the depth of
        interest."""

        if self.sketch['Sketch'] == self.sketch['Sketches'][0]: Time_graph  = self.outputgraphs.add_subplot(111)
        else:                                                   Time_graph  = self.outputgraphs.add_subplot(self.gs[1])

        zlegend         = []
        profiles        = []
        max_value = 0

        for n in range(len(self.timeplotdatas)):
            filenum = self.profiles.index(self.timeplotdatas[n].name)
            color     = self.colors[self.colornames.index(self.timeplotdatas[n].color)]
            style     = self.styles[self.stylenames.index(self.timeplotdatas[n].style)]
            linewidth = self.timeplotdatas[n].size

            if self.output_type[filenum] == 'Non-CapSim data':
                chemnum = self.chemicals[filenum].index(self.timeplotdatas[n].chemical)
                if self.stylenames.index(self.timeplotdatas[n].style) < 4:
                    p, = Time_graph.plot(self.times[filenum][chemnum], plots[n], color = color, linestyle = style, linewidth = linewidth)
                    if len(plots_smax[n]) > 0: Time_graph.errorbar(self.times[filenum][chemnum], plots[n], yerr = plots_smax[n], fmt = style, ecolor = color, elinewidth = 0.5)
                else:
                    p, = Time_graph.plot(self.times[filenum][chemnum], plots[n], color = color, linestyle = 'None', marker = style,  markersize = linewidth)
                    if len(plots_smax[n]) > 0: Time_graph.errorbar(self.times[filenum][chemnum], plots[n], yerr = plots_smax[n], fmt = style, ecolor = color, elinewidth = 0.5)
            elif self.output_type[filenum] == 'Monte Carlo Simulation':
                Time_graph.fill_between(self.times[filenum], plots_smax[n], plots_smin[n], alpha = self.MCplots['SD'],      facecolor = color, edgecolor = color)
                Time_graph.fill_between(self.times[filenum], plots_max[n], plots_min[n],   alpha = self.MCplots['Max/Min'], facecolor = color, edgecolor = color)
                if self.stylenames.index(self.timeplotdatas[n].style) < 4:
                    p, = Time_graph.plot(self.times[filenum], plots[n], alpha = self.MCplots['AVG'], color = color, linestyle = style, linewidth = linewidth)
                else:
                    p, = Time_graph.plot(self.times[filenum], plots[n], alpha = self.MCplots['AVG'], color = color, linestyle = 'None', marker = style,  markersize = linewidth)
            else:
                if self.stylenames.index(self.timeplotdatas[n].style) < 4:
                    p, = Time_graph.plot(self.times[filenum], plots[n],  color = color, linestyle = style, linewidth = linewidth)
                else:
                    p, = Time_graph.plot(self.times[filenum], plots[n],  color = color, linestyle = 'None', marker = style,  markersize = linewidth)

            if self.timeplotdatas[n].legend != '':
                profiles.append(p)
                zlegend.append(self.timeplotdatas[n].legend)


            for i in range (len(self.times[filenum])):
                if abs(plots[n][i]) > max_value:
                    max_value = abs(plots[n][i])

        Time_graph.set_xlabel(self.figuretexts['xlabel'], fontsize = self.fontsize['Label'][0], family = self.family)
        Time_graph.set_ylabel(self.figuretexts['ylabel'], fontsize = self.fontsize['Label'][0], family =self.family)
        Time_graph.set_title(self.figuretexts['title'],   fontsize = self.fontsize['Title'], family = self.family)

        if max_value < 10**-3 or max_value >= 10**4:
            Time_graph.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1e'))

        Time_graph.set_xlim([self.axislimit[0], self.axislimit[1]])
        Time_graph.set_ylim([self.axislimit[2], self.axislimit[3]])
        if len(zlegend) > 0 and self.legendposition < 11:
            Time_graphleg = Time_graph.legend(profiles, zlegend, loc = self.legendposition, labelspacing = 0, borderpad = 0.3, handlelength = 1.5, handletextpad = 0.1, fancybox = 0)
            Time_graphleg.legendPatch.set_fc(plt.rcParams['axes.facecolor'])
            for caption in Time_graphleg.get_texts():
                caption.set_fontsize(self.fontsize['Legend'])
                caption.set_fontname(self.family)

class PlotEditor:

    def __init__(self, master, system, type, variable, profiles, chemicals, components, solids, dep, spatialplotdatas, timeplotdatas, outputz, outputt, layer_list, output_type):
        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.tframe    = Frame(master.tframe)
        self.frame     = Frame(master.frame)
        self.bframe    = Frame(master.bframe)
        self.top       = None
        self.flag      = 0
        self.system    = system

        rgb            = self.master.frame.winfo_rgb(self.master.frame.cget('bg'))
        self.color     = '#%02x%02x%02x' %rgb
        self.mplrgb    = (rgb[0] / 65536., rgb[1] / 65536., rgb[2] / 65536.)

        self.types          = ['Spatial profile', 'Time profile']
        self.variables      = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Material fraction']
        self.extravariables = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Water concentration', 'Material fraction', 'Mass']

        self.layers  = []
        layer_extend = ['All layers', 'Top accumulation', 'Bottom accumulation', 'Deposition', 'Reaction/Generation', 'Universal']
        for layers in layer_list:
            self.layers.append([])
            for layer_name in layer_extend:
                self.layers[-1].append(layer_name)
            for layer_name in layers:
                self.layers[-1].append(layer_name)

        self.type      = StringVar(value = type)
        self.variable  = StringVar(value = variable)

        self.profiles    = profiles
        self.chemicals   = chemicals
        self.components  = components
        self.solids      = solids
        self.dep         = dep
        self.output_type = output_type

        self.timeunit   = system.timeunit
        self.lengthunit = system.lengthunit

        self.cancelflag =  0

        self.spatialplotdatas   = [plotdata.copy() for plotdata in spatialplotdatas]
        self.timeplotdatas      = [plotdata.copy() for plotdata in timeplotdatas]

        for plotdata in (self.spatialplotdatas + self.timeplotdatas):
            plotdata.name     = StringVar(value = plotdata.name)
            plotdata.chemical = StringVar(value = plotdata.chemical)
            plotdata.value    = DoubleVar(value = plotdata.value)
            plotdata.type     = StringVar(value = plotdata.type)
            plotdata.component= StringVar(value = plotdata.component)
            plotdata.solid    = StringVar(value = plotdata.solid)
            plotdata.layer    = StringVar(value = plotdata.layer)

        self.outputz = outputz
        self.outputt = outputt

    def make_widgets(self):

        self.instructions   = Label(self.tframe, text = ' Please provides the following output information                    ')

        self.startcolumn    = Label(self.tframe, text = ' ', width = 2)
        self.delcolumn      = Label(self.tframe, text = ' ', width = 8)
        self.numcolumn      = Label(self.tframe, text = ' ', width = 6)
        self.namecolumn     = Label(self.tframe, text = ' ', width = 20)
        self.chemcolumn     = Label(self.tframe, text = ' ', width = 20)
        self.valuecolumn    = Label(self.tframe, text = ' ', width = 12)
        self.matrixcolumn   = Label(self.tframe, text = ' ', width = 20)
        self.depthcolumn    = Label(self.tframe, text = ' ', width = 21)
        self.endcolumn      = Label(self.tframe, text = ' ', width = 2)

        self.typelabel      = Label(self.tframe, bg = self.color, text = 'Plot type:')
        self.varilabel      = Label(self.tframe, bg = self.color, text = 'Plot variable:')
        self.sizelabel      = Label(self.tframe, bg = self.color, text = 'Plot size:')
        self.sizeunit       = Label(self.tframe, bg = self.color, text = 'Pixels')

        self.typewidget     = OptionMenu(self.tframe, self.type,     *self.types,     command = self.updateplots)
        self.variwidget     = OptionMenu(self.tframe, self.variable, *self.variables, command = self.updateplots)
        self.extravariwidget= OptionMenu(self.tframe, self.variable, *self.extravariables, command = self.updateplots)

        self.typewidget.config(width = 16)
        self.variwidget.config(width = 16)
        self.extravariwidget.config(width = 16)

        self.namelabel      = Label(self.tframe, text = 'Profile Name')
        self.chemlabel      = Label(self.tframe, text = 'Chemical')
        self.timelabel      = Label(self.tframe, text = 'Plot Time(' + self.timeunit + ')')
        self.depthlabel     = Label(self.tframe, text = 'Plot Depth (' + self.lengthunit + ')')
        self.ztlabel        = Label(self.tframe, bg = self.color, text = 'Depth from:')
        self.matrixlabel    = Label(self.tframe, bg = self.color, text = 'Matrix:')

        self.blank3         = Label (self.tframe, text = ' ')

        self.botstartcolumn = Label(self.frame, text = ' ', width = 2)
        self.botdelcolumn   = Label(self.frame, text = ' ', width = 8)
        self.botnumcolumn   = Label(self.frame, text = ' ', width = 6)
        self.botnamecolumn  = Label(self.frame, text = ' ', width = 20)
        self.botchemcolumn  = Label(self.frame, text = ' ', width = 20)
        self.botvaluecolumn = Label(self.frame, text = ' ', width = 12)
        self.botmatrixcolumn= Label(self.frame, text = ' ', width = 20)
        self.botdepthcolumn = Label(self.frame, text = ' ', width = 21)
        self.botendcolumn   = Label(self.frame, text = ' ', width = 2)

        self.blank1         = Label (self.bframe, text = ' ')
        self.blank2         = Label (self.bframe, text = ' ')

        self.addwidget      = Button(self.bframe, text = 'Add plots', command = self.addplotdata, width = 20)
        self.okbutton       = Button(self.bframe, text = 'OK', width = 20, command = self.OK)
        self.cancelbutton   = Button(self.bframe, text = 'Cancel', width = 20, command = self.Cancel)

        self.instructions.grid( row = 0, column = 0, columnspan = 6, padx = 8, sticky = 'W')

        self.startcolumn.grid(  row = 1, column = 0)
        self.delcolumn.grid(    row = 1, column = 1)
        self.numcolumn.grid(    row = 1, column = 2)
        self.namecolumn.grid(   row = 1, column = 3)
        self.chemcolumn.grid(   row = 1, column = 4)
        self.valuecolumn.grid(  row = 1, column = 5)
        self.matrixcolumn.grid( row = 1, column = 6)
        self.endcolumn.grid(    row = 1, column = 8)

        self.typelabel.grid(    row = 2, column = 1, sticky = 'E',  pady = 1, columnspan = 2)
        self.typewidget.grid(   row = 2, column = 3, sticky = 'W',  pady = 1)
        self.varilabel.grid(    row = 2, column = 4, sticky = 'E',  pady = 1)

        self.blank3.grid(       row = 3, column = 3)
        self.namelabel.grid(    row = 4, column = 3)
        self.chemlabel.grid(    row = 4, column = 4)
        self.matrixlabel.grid(  row = 4, column = 6)

        self.botstartcolumn.grid(  row = 0, column = 0)
        self.botdelcolumn.grid(    row = 0, column = 1)
        self.botnumcolumn.grid(    row = 0, column = 2)
        self.botnamecolumn.grid(   row = 0, column = 3)
        self.botchemcolumn.grid(   row = 0, column = 4)
        self.botvaluecolumn.grid(  row = 0, column = 5)
        self.botmatrixcolumn.grid( row = 0, column = 6)
        self.botendcolumn.grid(    row = 0, column = 8)

        self.updateplots()

    def updateplots(self, event = None):

        try:
            self.timelabel.grid_forget()
            self.depthlabel.grid_forget()
            self.ztlabel.grid_forget()
            self.depthcolumn.grid_forget()
            self.botdepthcolumn.grid_forget()
            self.ztlabel.grid_forget()
        except: pass

        for plotdata in self.spatialplotdatas:
            try: plotdata.remove_propertieswidgets()
            except:pass

        for plotdata in self.timeplotdatas:
            try: plotdata.remove_propertieswidgets()
            except:pass

        try:
            self.variwidget.grid_forget()
        except: pass
        try:
            self.extravariwidget.grid_forget()
        except: pass


        if self.type.get() == 'Spatial profile':
            if self.variable.get() == 'Water concentration' or self.variable.get() == 'Mass': self.variable.set('Concentration')
            self.variwidget.grid(   row = 2, column = 5,  pady = 1, sticky = 'W',  columnspan = 2)
        else:
            self.extravariwidget.grid(row = 2, column = 5,  pady = 1, sticky = 'W',  columnspan = 2)

        row = 2

        if self.type.get() == 'Spatial profile':
            self.timelabel.grid(  row = 4, column = 5)
            for plotdata in self.spatialplotdatas:
                plotdata.number = self.spatialplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.profiles, self.chemicals, self.output_type, self.type.get(), self.variable.get(), self.components, self.solids)
                row = row + 1

        elif self.dep.count(1) >= 1 and self.variable.get() != 'Mass' and self.variable.get() != 'Water concentration':
            self.depthcolumn.grid(      row = 1, column = 7)
            self.depthlabel.grid(       row = 4, column = 5)
            self.ztlabel.grid(          row = 4, column = 7)
            self.botdepthcolumn.grid(   row = 0, column = 7)
            for plotdata in self.timeplotdatas:
                plotdata.number = self.timeplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.profiles, self.chemicals, self.output_type, self.type.get(), self.variable.get(), self.components, self.solids, depositionflag = 1, layer_list = self.layers)
                row = row + 1
        else:
            self.depthlabel.grid(       row = 4, column = 5)
            for plotdata in self.timeplotdatas:
                plotdata.number = self.timeplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.profiles, self.chemicals, self.output_type, self.type.get(), self.variable.get(), self.components, self.solids, layer_list = self.layers)
                row = row + 1

        self.blank1.grid(row = row)
        row = row + 1
        self.blank2.grid(row = row)
        row = row + 1
        self.addwidget.grid(row = row, columnspan = 11)
        row = row + 1
        self.okbutton.grid(row = row, columnspan = 11)
        row = row + 1
        self.cancelbutton.grid(row = row, columnspan = 11)
        row = row + 1

        self.focusbutton = self.okbutton
        self.okbutton.bind('<Return>',   self.OK)

        self.master.geometry()
        self.master.center()

    def addplotdata(self):

        if self.output_type[0] == 'Non-CapSim data':
            if self.type.get() == 'Spatial profile':
                self.spatialplotdatas.append(PlotData(len(self.spatialplotdatas)+ 1, self.timeunit, self.lengthunit))

                self.spatialplotdatas[-1].name      = StringVar(value = self.profiles[0])
                self.spatialplotdatas[-1].chemical  = StringVar(value = self.chemicals[0][0])
                self.spatialplotdatas[-1].value     = DoubleVar(value = 0)
                self.spatialplotdatas[-1].component = StringVar(value = self.components[0][0])
                self.spatialplotdatas[-1].solid     = StringVar(value = 'Total solid')
                self.spatialplotdatas[-1].type      = StringVar(value = 'Initial benthic surface')
                self.spatialplotdatas[-1].layer     = StringVar(value = 'All layers')
                self.spatialplotdatas[-1].style     = 'Point'
                self.spatialplotdatas[-1].color     = self.spatialplotdatas[-1].colornames[self.spatialplotdatas[-1].number - 1]
                self.spatialplotdatas[-1].size      = 10
                self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name.get() + ' ' + self.spatialplotdatas[-1].chemical.get()

            else:
                self.timeplotdatas.append(PlotData(len(self.timeplotdatas)+ 1, self.timeunit, self.lengthunit))

                self.timeplotdatas[-1].name         = StringVar(value = self.profiles[0])
                self.timeplotdatas[-1].chemical     = StringVar(value = self.chemicals[0][0])
                self.timeplotdatas[-1].value        = DoubleVar(value = 0)
                self.timeplotdatas[-1].type         = StringVar(value = 'Initial benthic surface')
                self.timeplotdatas[-1].component    = StringVar(value = self.components[0][0])
                self.timeplotdatas[-1].solid        = StringVar(value = 'Total solid')
                self.timeplotdatas[-1].layer        = StringVar(value = 'All layers')
                self.timeplotdatas[-1].style        = 'Point'
                self.timeplotdatas[-1].color        = self.timeplotdatas[-1].colornames[self.timeplotdatas[-1].number - 1]
                self.timeplotdatas[-1].size         = 10
                self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name.get() + ' ' + self.timeplotdatas[-1].chemical.get()
        else:
            if self.type.get() == 'Spatial profile':
                self.spatialplotdatas.append(PlotData(len(self.spatialplotdatas)+ 1, self.timeunit, self.lengthunit))

                self.spatialplotdatas[-1].name      = StringVar(value = self.profiles[0])
                self.spatialplotdatas[-1].chemical  = StringVar(value = self.chemicals[0][0])
                self.spatialplotdatas[-1].value     = DoubleVar(value = 0)
                self.spatialplotdatas[-1].component = StringVar(value = self.components[0][0])
                self.spatialplotdatas[-1].solid     = StringVar(value = 'Total solid')
                self.spatialplotdatas[-1].type      = StringVar(value = 'Initial benthic surface')
                self.spatialplotdatas[-1].layer     = StringVar(value = 'All layers')
                self.spatialplotdatas[-1].style     = 'Solid'
                self.spatialplotdatas[-1].color     = self.spatialplotdatas[-1].colornames[self.spatialplotdatas[-1].number - 1]
                self.spatialplotdatas[-1].size      = 0.5
                self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name.get() + ' ' + self.spatialplotdatas[-1].chemical.get() + ' ' +  str(self.spatialplotdatas[-1].value.get()) + ' ' + self.timeunit

            else:
                self.timeplotdatas.append(PlotData(len(self.timeplotdatas)+ 1, self.timeunit, self.lengthunit))

                self.timeplotdatas[-1].name         = StringVar(value = self.profiles[0])
                self.timeplotdatas[-1].chemical     = StringVar(value = self.chemicals[0][0])
                self.timeplotdatas[-1].value        = DoubleVar(value = 0)
                self.timeplotdatas[-1].type         = StringVar(value = 'Initial benthic surface')
                self.timeplotdatas[-1].component    = StringVar(value = self.components[0][0])
                self.timeplotdatas[-1].solid        = StringVar(value = 'Total solid')
                self.timeplotdatas[-1].layer        = StringVar(value = 'All layers')
                self.timeplotdatas[-1].style        = 'Solid'
                self.timeplotdatas[-1].color        = self.timeplotdatas[-1].colornames[self.timeplotdatas[-1].number - 1]
                self.timeplotdatas[-1].size         = 0.5
                if self.dep[-1] == 1:
                    self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name.get() + ' ' + self.timeplotdatas[-1].chemical.get() +' ' +  str(self.timeplotdatas[-1].value.get()) + ' ' + self.lengthunit + ' from ' + self.timeplotdatas[-1].type.get()
                else:
                    self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name.get() + ' ' + self.timeplotdatas[-1].chemical.get() +' ' +  str(self.timeplotdatas[-1].value.get()) + ' ' + self.lengthunit

        self.updateplots()

    def delplotdata(self, number):

        if self.type.get() == 'Spatial profile':
            self.spatialplotdatas[number - 1].remove_propertieswidgets()
            self.spatialplotdatas.remove(self.spatialplotdatas[number - 1])
        else:
            self.timeplotdatas[number - 1].remove_propertieswidgets()
            self.timeplotdatas.remove(self.timeplotdatas[number - 1])

        self.updateplots()

    def OK(self, event = None):
        """Finish and move on.  Checks that the number chemicals are less than the
        total number of chemicals in database."""

        error_flag = 0

        if self.type.get() == 'Spatial profile':
            for plotdata in self.spatialplotdatas:
                if self.output_type[self.profiles.index(plotdata.name.get())] != 'Non-CapSim data':
                    if plotdata.value.get() < self.outputt[self.profiles.index(plotdata.name.get())][0]:  error_flag = 1
                    if plotdata.value.get() > self.outputt[self.profiles.index(plotdata.name.get())][-1]: error_flag = 1
                plotdata.update_legend()
        else:
            for plotdata in self.timeplotdatas:
                if self.output_type[self.profiles.index(plotdata.name.get())] != 'Non-CapSim data':
                    if plotdata.type == 'Initial benthic surface':
                        if plotdata.value.get() < self.outputz[self.profiles.index(plotdata.name.get())][0]:  error_flag = 1
                        if plotdata.value.get() > self.outputz[self.profiles.index(plotdata.name.get())][-1]: error_flag = 1
                    else:
                        if plotdata.value.get() < self.outputz[self.profiles.index(plotdata.name.get())][0] - self.outputz[self.profiles.index(plotdata.name.get())][0]: error_flag = 1
                        if plotdata.value.get() > self.outputz[self.profiles.index(plotdata.name.get())][-1]- self.outputz[self.profiles.index(plotdata.name.get())][0]: error_flag = 1
                plotdata.update_legend()
        if self.master.window.top is not None: self.master.open_toplevel()
        elif error_flag == 1: self.warning()
        else: self.master.tk.quit()

    def warning(self):

        tkmb.showerror(title = self.version, message = 'The input depth/time is out of range, please correct')
        self.focusbutton = None
        self.master.tk.lift()

    def Cancel(self):

        self.cancelflag = 1

        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()

class PlotData:

    def __init__(self, number, timeunit, lengthunit):

        self.number     = number
        self.types      = ['Initial benthic surface', 'New benthic surface']
        self.colornames = ['Blue', 'Green', 'Red', 'Cyan', 'Magenta', 'Yellow', 'Orange', 'Purple', 'Brown', 'Pink', 'Grey', 'Black']
        self.stylenames = ['Solid', 'Dashed', 'Dash-dot', 'Dotted', 'Point', 'Circle', 'Triangle', 'Square', 'Star', 'Cross']
        self.timeunit   = timeunit
        self.lengthunit = lengthunit

    def copy(self):

        plotdata = PlotData(self.number, self.timeunit, self.lengthunit)

        plotdata.name      = self.name
        plotdata.value     = self.value
        plotdata.type      = self.type
        plotdata.chemical  = self.chemical
        plotdata.component = self.component
        plotdata.solid     = self.solid
        plotdata.layer     = self.layer
        plotdata.number    = self.number
        plotdata.color     = self.color
        plotdata.style     = self.style
        plotdata.size      = self.size
        plotdata.legend    = self.legend

        return plotdata

    def propertieswidgets(self, frame, row, master, profiles, chemical_list, output_type, type, variable, component_list, solid_list, depositionflag = None, layer_list = None):

        self.master = master

        self.profiles       = profiles
        self.chemical_list  = chemical_list
        self.component_list = component_list
        self.solid_list     = solid_list
        self.layer_list     = layer_list
        self.variable       = variable
        self.output_type    = output_type
        self.plot_type      = type
        self.depositionflag = depositionflag

        self.row        = row
        self.frame      = frame

        self.delwidget      = Button(frame, width = 5,  justify = 'center', text = 'Delete', command = self.del_plotdata)
        self.numwidget      = Label (frame, width = 4,  justify = 'center', text = self.number)
        self.namewidget     = OptionMenu(frame, self.name, *profiles,       command = self.click_name)
        self.valuewidget    = Entry (frame, width = 8,  justify = 'center', textvariable = self.value)
        self.valuelabel     = Label (frame, width = 8,  justify = 'center', text = 'N/A')

        self.delwidget.grid(      row  = row, column = 1, padx = 2 ,pady = 1)

        self.numwidget.grid(      row  = row, column = 2, padx = 2 ,pady = 1)
        self.namewidget.grid(     row  = row, column = 3, padx = 2 ,pady = 1, sticky = 'WE')

        self.click_name()


    def click_name(self, default = None):

        try: self.componentwidget.grid_forget()
        except: pass

        try:
            self.chemwidget.grid_forget()
        except: pass

        try: self.typewidget.grid_forget()
        except: pass

        try:
            self.valuelabel.grid_forget()
            self.valuewidget.grid_forget()
        except: pass

        if default <> None: self.update_legend()

        filenum = self.profiles.index(self.name.get())

        if self.variable == 'Water concentration' or self.variable == 'Mass' or self.output_type[filenum] == 'Non-CapSim data':
            self.valuelabel.grid(    row  = self.row, column = 5, padx = 2, pady = 1)
        else:
            self.valuewidget.grid(    row  = self.row, column = 5, padx = 2, pady = 1)

        if self.output_type[filenum] == 'Non-CapSim data':
            self.componentwidget = Label (self.frame, justify = 'center', text = 'N/A')
            self.chemwidget = OptionMenu(self.frame, self.chemical, *self.chemical_list[filenum], command = self.update_legend)
            self.typewidget = Label (self.frame, justify = 'center', text = 'N/A')
        else:
            if self.variable == 'Solid concentration':
                if self.component_list[filenum].count(self.component.get()) == 0:
                    self.component.set(self.component_list[filenum][0])
                self.componentwidget = OptionMenu(self.frame, self.component, *self.component_list[filenum], command = self.update_legend)
            elif self.variable == 'Mass':
                if self.layer_list[filenum].count(self.layer.get()) == 0:
                    self.layer.set(self.layer_list[filenum][0])
                self.componentwidget = OptionMenu(self.frame, self.layer, *self.layer_list[filenum], command = self.update_legend)
            elif self.variable == 'Concentration':
                self.componentwidget = Label (self.frame, width = 16,  justify = 'center', text = 'Porewater')
            elif self.variable == 'Total concentration':
                self.componentwidget = Label (self.frame, width = 16,  justify = 'center', text = 'Porewater/DOC')
            elif self.variable == 'Flux':
                self.componentwidget = Label (self.frame, width = 16,  justify = 'center', text = 'Porewater/DOC/Solid')
            elif self.variable == 'Water concentration':
                self.componentwidget = Label (self.frame, width = 16,  justify = 'center', text = 'Overlying water column')
            elif self.variable == 'Material fraction':
                self.componentwidget = OptionMenu(self.frame, self.component, *self.component_list[filenum], command = self.update_legend)

            if self.variable == 'Material fraction':
                self.chemwidget = Label (self.frame, width = 16,  justify = 'center', text = 'N/A')
            else:
                self.chemwidget = OptionMenu(self.frame, self.chemical, *self.chemical_list[filenum], command = self.update_legend)

            if self.depositionflag == 1 and self.variable != 'Mass':
                self.typewidget     = OptionMenu(self.frame, self.type, *self.types, command = self.update_legend)
                self.typewidget.grid(     row  = self.row, column = 7, padx = 2, pady = 1, sticky = 'WE')


        if self.chemical_list[filenum].count(self.chemical.get()) == 0:
            self.chemical.set(self.chemical_list[filenum][0])

        if self.component_list[filenum].count(self.component.get()) == 0:
            self.component.set(self.component_list[filenum][0])

        self.chemwidget.grid(     row  = self.row, column = 4, padx = 2 ,pady = 1, sticky = 'WE')
        self.componentwidget.grid(row  = self.row, column = 6, padx = 2 ,pady = 1, sticky = 'WE')

    def update_legend(self, default = None):

        filenum = self.profiles.index(self.name.get())

        if self.output_type[filenum] == 'Non-CapSim data':
            self.legend = self.name.get() + ' ' + self.chemical.get()
        elif self.variable == 'Material fraction':
            if self.plot_type == 'Spatial profile':
                self.legend = self.name.get() + ' ' + self.component.get() + ' at ' +  str(self.value.get()) + ' ' + self.timeunit
            else:
                self.legend = self.name.get() + ' ' + self.component.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit
        else:
            if self.plot_type == 'Spatial profile':
                if self.variable == 'Solid concentration':
                    self.legend = self.name.get() + ' ' + self.chemical.get() + ' in ' +  self.component.get() + ' at ' + str(self.value.get()) + ' ' + self.timeunit
                else:
                    self.legend = self.name.get() + ' ' + self.chemical.get() + ' at ' +  str(self.value.get()) + ' ' + self.timeunit
            else:
                if self.depositionflag == 1:
                    if self.variable == 'Solid concentration':
                        self.legend = self.name.get() + ' ' + self.chemical.get()+ ' in '+ self.component.get()  + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit + ' from ' + self.type.get()
                    else:
                        self.legend = self.name.get() + ' ' + self.chemical.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit + ' from ' + self.type.get()
                else:
                    if self.variable == 'Solid concentration':
                        self.legend = self.name.get() + ' ' + self.chemical.get()+ ' in '+ self.component.get()  + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit
                    else:
                        self.legend = self.name.get() + ' ' + self.chemical.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit

    def remove_propertieswidgets(self):

        self.delwidget.grid_forget()
        self.namewidget.grid_forget()
        self.chemwidget.grid_forget()
        self.numwidget.grid_forget()

        try: self.typewidget.grid_forget()
        except:pass
        try: self.componentwidget.grid_forget()
        except:pass
        try: self.valuewidget.grid_forget()
        except:pass
        try: self.valuelabel.grid_forget()
        except:pass

        self.master          = 0
        self.frame           = 0
        self.delwidget       = 0
        self.namewidget      = 0
        self.chemwidget      = 0
        self.numwidget       = 0
        self.valuewidget     = 0
        self.typewidget      = 0
        self.componentwidget = 0

    def graphicwidgets(self, frame, row, master):

        self.master = master
        self.frame  = frame
        self.row        = row

        self.numwidget       = Label (frame, width = 4,  justify = 'center', text = str(self.number))
        self.stylewidget     = OptionMenu(frame, self.style, *self.stylenames)
        self.colorwidget     = OptionMenu(frame, self.color, *self.colornames)
        self.sizewidget      = Entry(frame, width = 10,   justify = 'center', textvariable = self.size)
        self.legendwidget    = Entry(frame, width = 20,  justify = 'center', textvariable = self.legend)

        self.numwidget.grid(        row  = row, column = 2, padx = 2 ,pady = 1)
        self.stylewidget.grid(      row  = row, column = 3, padx = 2 ,pady = 1, sticky = 'WE')
        self.colorwidget.grid(      row  = row, column = 4, padx = 2 ,pady = 1, sticky = 'WE')
        self.sizewidget.grid(       row  = row, column = 5, padx = 2 ,pady = 1)
        self.legendwidget.grid(     row  = row, column = 6, padx = 2 ,pady = 1)

    def get_plotdata(self):

        self.name       = self.name.get()
        self.chemical   = self.chemical.get()
        self.value      = self.value.get()
        self.type       = self.type.get()
        self.component  = self.component.get()
        self.solid      = self.solid.get()
        self.layer      = self.layer.get()

    def del_plotdata(self):

        self.master.window.delplotdata(self.number)

    def get_graphicdata(self):

        self.style     = self.style.get()
        self.color     = self.color.get()
        self.size      = self.size.get()
        self.legend    = self.legend.get()


def graphprocess_data(system):
    """Shows the results of the simulation."""

    root = CapSimWindow(buttons = 2)
    root.make_window(Graphprocess(root, system))
    root.tk.config(background = root.window.color)
    root.buttonframe.config(background = root.window.color)
    root.blank.config(background = root.window.color)

    root.mainloop()
    main = root.main.get()
    root.destroy()
    
    return main

def time_interpolate(tint, t, delt, Cn, Cn_plus_1):
    """Returns the interpolated concentrations at time "tint" using the
    concentrations "Cn" at time "t" and concentrations "Cn_plus_1" at time 
    "t + delt." """

    return (Cn + (tint - t) / delt * (Cn_plus_1 - Cn))

def interp_dep(zpoint, dep, delz, z, Cn):
    """Returns the interpolated concentrations at time "tint" using the
    concentrations "Cn" at time "t" and concentrations "Cn_plus_1" at time
    "t + delt." """

    jj = 0

    for j in range(len(z)-1):
        if round(z[j], 8) == 0.0: jo = j

    if round(dep, 8) > round(0.5*delz, 8):
        ans = interp(zpoint+dep, [dep]+z[jo+1:], Cn[jo:])

    else:
        for j in range(len(z)-1):
            if round(dep, 8) >= round(z[j] + 0.5 * delz,8) and round(dep,8) < round(z[j+1] + 0.5 * delz, 8):
                jj = j + 1
        ans = interp(zpoint+dep, [dep]+z[jj+1:], Cn[jj:])

    return ans
