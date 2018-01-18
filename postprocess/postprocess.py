# -*- coding: utf-8 -*-

#!/usr/bin/env python
#These subroutines are used in post-processing for CapSim.

import matplotlib.pyplot as plt, tkMessageBox as tkmb, math, sys, os, codecs, _winreg as wreg
import tkFileDialog as tkfd
import matplotlib.patches  as patches
import matplotlib.gridspec as gridspec
import matplotlib.tri      as tri
from capsim_functions import text_converter

if sys.path[0][-3:] == 'zip':
    os.chdir(sys.path[0][:-12])
    path = sys.path[0][:-12]
else: path = sys.path[0]

CapSimReg = wreg.ConnectRegistry(None, wreg.HKEY_CURRENT_USER)
CapSimKey = wreg.OpenKey(CapSimReg, r'Software\CapSim')
Filepath =wreg.QueryValueEx(CapSimKey, 'FilePath')[0]
CapSimKey.Close()

from numpy               import zeros, transpose, interp, array, isnan, shape, mean, std, nanmax, nanmin
from capsim_object_types import System, CapSimWindow
from capsim_functions    import round_to_n
from datetime            import datetime
from Tkinter             import Tk, Toplevel, Canvas, Frame, Label, Entry, Text, Button, Scrollbar, OptionMenu, StringVar, DoubleVar, IntVar, FLAT, RAISED, Checkbutton
from PIL                 import Image, ImageTk
from textwrap            import fill

class Postprocess:
    """Makes a window to display the simulation plots."""

    def __init__(self, master, system, output, batch_type = None):
        """Constructor method."""

        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.system    = system
        self.output    = output
        rgb            = self.master.frame.winfo_rgb(self.master.frame.cget('bg'))
        self.color     = '#%02x%02x%02x' %rgb
        self.mplrgb    = (rgb[0] / 65536., rgb[1] / 65536., rgb[2] / 65536.) 
        self.tframe    = Frame(master.tframe, bg = self.color)
        self.frame     = Frame(master.frame,  bg = self.color)
        self.bframe    = Frame(master.bframe, bg = self.color)
        self.top       = None
        self.flag      = 0
        self.filename  = 'plot_output'

        self.lengthunit = system.lengthunit
        self.concunit   = system.concunit
        self.timeunit   = system.timeunit
        self.diffunit   = system.diffunit

        self.colors      = ['b', 'g', 'r', 'c', 'm', 'y', '#ff7f0e', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', 'k']
        self.colornames  = ['Blue', 'Green', 'Red', 'Cyan', 'Magenta', 'Yellow', 'Orange', 'Purple', 'Brown', 'Pink', 'Grey', 'Black']

        self.styles      = ['-', '--', '-.', ':', '.', 'o', '^', 's', '*', 'x']
        self.stylenames  = ['Solid', 'Dashed', 'Dash-dot', 'Dotted', 'Point', 'Circle', 'Triangle', 'Square', 'Star', 'Cross']

        self.lengthunits = [u'\u03BCm', 'cm', 'm']
        self.concunits   = [u'\u03BCg/L', 'mg/L', 'g/L', u'\u03BCmol/L', 'mmol/L', 'mol/L']
        self.timeunits   = ['s', 'min', 'hr', 'day', 'yr']
        self.diffunits   = [u'cm\u00B2/s', u'cm\u00B2/yr']

        self.types        = ['Spatial profile', 'Time profile']
        self.variables    = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Water concentration', 'Material fraction', 'Mass']

        self.sketch       = {}
        self.sketch['Sketches'] = ['Show sketch', 'Hide sketch']
        self.sketch['Sketch']   = 'Show sketch'
        self.sketch['Label']    = 1
        self.sketch['Time']     = 0
        self.sketch['Depth']    = 0

        self.layer_list = [layer.name for layer in system.layers]
        self.chemicals  = [chemical.name for chemical in system.chemicals]
        self.components = [component for component in system.component_list]
        self.solids     = [component for component in system.component_list]
        self.solids.insert(0, 'Total solid')
        self.nchemicals = system.nchemicals

        if system.dep == 'Deposition':
            self.dep = 1
            self.Vdep = system.Vdep
        else:
            self.dep = 0

        self.times = self.output.times
        self.z     = self.output.z
        self.O     = self.output.O
        self.OL    = self.output.OL
        self.C     = self.output.C
        self.CL    = self.output.CL
        self.F     = self.output.F
        self.W     = self.output.W
        self.WL    = self.output.WL
        self.Cw    = self.output.Cw
        self.FM    = self.output.FM
        self.DM    = self.output.DM
        self.M     = self.output.M
        self.RM    = self.output.RM
        self.Fi    = self.output.Fi
        self.FiL   = self.output.FiL

        self.q     = zeros([len(self.times), len(self.z), self.nchemicals, len(self.solids)])
        self.qL    = zeros([len(self.times), len(self.z), self.nchemicals, len(self.solids)])
        self.q[:,:,:,0]  = self.output.q
        self.qL[:,:,:,0] = self.output.qL
        self.q[:,:,:,1:] = self.output.qm
        self.qL[:,:,:,1:] = self.output.qmL

        self.layer_h   =[0]
        for layer in system.layers:
            if layer.name !='Deposition':
                self.layer_h.append(self.layer_h[-1]+layer.h)

        self.type        = self.types[0]
        self.variable    = self.variables[0]

        self.spatialplotdatas   = []
        for i in range(6):
            self.spatialplotdatas.append(PlotData(i, self.timeunit, self.lengthunit))
            self.spatialplotdatas[-1].name      = self.chemicals[0]
            self.spatialplotdatas[-1].value     = round_to_n(self.times[-1]/5*i, 2)
            self.spatialplotdatas[-1].type      = self.spatialplotdatas[-1].types[0]
            self.spatialplotdatas[-1].component = self.components[0]
            self.spatialplotdatas[-1].solid     = 'Total solid'
            self.spatialplotdatas[-1].layer     = 'All layers'
            self.spatialplotdatas[-1].style     = 'Solid'
            self.spatialplotdatas[-1].color     = self.colornames[i]
            self.spatialplotdatas[-1].size      = 0.5
            self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name + ' ' +  str(self.spatialplotdatas[-1].value) + ' ' + self.timeunit

            if self.spatialplotdatas[-1].value > self.times[-1]: self.spatialplotdatas[-1].value = self.times[-1]

        self.timeplotdatas               = [PlotData(0, self.timeunit, self.lengthunit)]
        self.timeplotdatas[-1].name      = self.chemicals[0]
        self.timeplotdatas[-1].value     = 0
        self.timeplotdatas[-1].type      = self.timeplotdatas[-1].types[0]
        self.timeplotdatas[-1].component = self.components[0]
        self.timeplotdatas[-1].solid     = 'Total solid'
        self.timeplotdatas[-1].layer     = 'All layers'
        self.timeplotdatas[-1].style     = 'Solid'
        self.timeplotdatas[-1].color     = self.colornames[0]
        self.timeplotdatas[-1].size      = 0.5
        if self.system.dep == 'Deposition':
            self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name + ' ' +  str(self.timeplotdatas[-1].value) + ' ' + self.lengthunit + ' from ' + self.timeplotdatas[-1].type
        else:
            self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name + ' ' +  str(self.timeplotdatas[-1].value) + ' ' + self.lengthunit
        self.plot  = []
        self.Lplot = []

        for n in range(len(self.spatialplotdatas)):
            chemnum = self.chemicals.index(self.spatialplotdatas[n].name)
            i = 0
            if self.spatialplotdatas[n].value >= self.times[-1]:
                i =  len(self.times) - 2
            else:
                while self.spatialplotdatas[n].value > self.times[i+1]: i = i + 1
            self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i + 1] - self.times[i], self.C[i, :, chemnum], self.C[i + 1, :, chemnum]))
            self.Lplot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i + 1] - self.times[i], self.CL[i, :, chemnum], self.CL[i + 1, :, chemnum]))

        Cmax = 0
        for n in range(6):
            if max(self.plot[n]) > Cmax:
                Cmax = max(self.plot[n])

        self.figuresize = [600, 450]

        if master.tk.winfo_screenheight()/master.tk.winfo_screenwidth() > 3/4:
            if master.tk.winfo_screenheight()*3/5 < 450:
                self.figuresize[1] = int(master.tk.winfo_screenheight()*3/5)
                self.figuresize[0] = int(self.figuresize[1]/3*4)
        else:
            if master.tk.winfo_screenwidth()*4/5 < 600:
                self.figuresize[0] = int(master.tk.winfo_screenwidth()*4/5)
                self.figuresize[1] = int(self.figuresize[0]*3/4)

        self.axislimit      = [0, 1.2 * Cmax, min(self.output.z), max(self.output.z)]
        self.linewidth      = 0.5
        self.layerlinewidth = 0.5

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

        self.graphpath = Filepath + '\output\\'

        self.gs = gridspec.GridSpec(1, 2, width_ratios=[1, 5])
        self.gs.update(left = 0.01, right = 0.95)
        self.batch_type = batch_type

        if batch_type == None:  self.make_graphs()


    def make_widgets(self):
        """Makes the widgets for the window."""

        self.master.frame.config(bg      = self.color)
        self.master.mainbutton.config(bg = self.color)
        self.master.exitbutton.config(bg = self.color)

        self.graph      = Image.open(Filepath + r'/output/%s.png' % self.filename)
        self.image      = ImageTk.PhotoImage(self.graph)
        self.graphlabel = Label(self.frame, image = self.image,       bg = self.color)

        self.blank1     = Label(self.bframe, text = ' ')
        self.blank2     = Label(self.bframe, text = ' ')
        self.blank3     = Label(self.bframe, text = ' ')

        self.editbutton   = Button(self.bframe, text = 'Edit plot',     bg = self.color, command = self.editplot,    width = 20, activebackground = 'white', highlightbackground = 'black')
        self.editgbutton  = Button(self.bframe, text = 'Edit profiles', bg = self.color, command = self.editprofile, width = 20, activebackground = 'white', highlightbackground = 'black')
        self.figurebutton = Button(self.bframe, text = 'Edit figure',   bg = self.color, command = self.editfigure,  width = 20, activebackground = 'white', highlightbackground = 'black')
        self.savebutton   = Button(self.bframe, text = 'Save figure',   bg = self.color, command = self.savegraph,   width = 20, activebackground = 'white', highlightbackground = 'black')
        self.exportbutton = Button(self.bframe, text = 'Export result', bg = self.color, command = self.export_data, width = 20, activebackground = 'white', highlightbackground = 'black')
        self.modifybutton = Button(self.bframe, text = 'Modify System', bg = self.color, command = self.modify,      width = 20, activebackground = 'white', highlightbackground = 'black')

        self.graphlabel.grid(   row = 0, column = 1, padx = 40)

        self.blank1.grid(       row = 1, column = 0, columnspan = 7, padx = 40)
        self.editbutton.grid(   row = 2, column = 0, columnspan = 7, padx = 40)
        self.editgbutton.grid(  row = 3, column = 0, columnspan = 7, padx = 40)
        self.figurebutton.grid( row = 4, column = 0, columnspan = 7, padx = 40)
        self.savebutton.grid(   row = 5, column = 0, columnspan = 7, padx = 40)
        self.exportbutton.grid( row = 6, column = 0, columnspan = 7, padx = 40)
        self.savebutton.grid(   row = 7, column = 0, columnspan = 7)

        if self.master.master is None: self.modifybutton.grid(row = 8, columnspan = 7, column = 0)

        self.modifybutton.bind('<Return>', self.modify)

        self.focusbutton = self.editbutton

        self.master.geometry()
        self.master.center()

    def editplot(self, event = None):

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 2)
            self.top.make_window(PlotEditor(self.top, self.system, self.type, self.variable, self.chemicals, self.components, self.solids, self.dep, self.spatialplotdatas, self.timeplotdatas, self.z, self.times, self.layer_list))

            self.top.tk.mainloop()

            if self.top.window.cancelflag == 0:
                self.type     = self.top.window.type.get()
                self.variable = self.top.window.variable.get()

                for spatialplotdata in self.top.window.spatialplotdatas:  spatialplotdata.get_plotdata()
                for timeplotdata in self.top.window.timeplotdatas:        timeplotdata.get_plotdata()

                self.spatialplotdatas = [spatialplotdata.copy() for spatialplotdata in self.top.window.spatialplotdatas]
                self.timeplotdatas    = [timeplotdata.copy()    for timeplotdata    in self.top.window.timeplotdatas]


                if self.type == self.types[0]:
                    self.axislimit[2]  = min(self.output.z)
                    self.axislimit[3]  = max(self.output.z)
                    self.plot  = []
                    self.Lplot = []

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
                        chemnum  = self.chemicals.index(self.spatialplotdatas[n].name)
                        compnum  = self.components.index(self.spatialplotdatas[n].component)
                        solidnum = self.solids.index(self.spatialplotdatas[n].solid)
                        i = 0
                        while self.spatialplotdatas[n].value > self.times[i+1]: i = i + 1

                        if self.variable == 'Material fraction':
                            self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], self.Fi[i,:,compnum], self.Fi[i+1, :, compnum]))
                            self.Lplot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], self.FiL[i,:,compnum], self.FiL[i+1, :, compnum]))
                        elif self.variable == 'Solid concentration':
                            self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], self.q[i,:,chemnum,solidnum], self.q[i+1,:,chemnum,solidnum]))
                            self.Lplot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], self.qL[i,:,chemnum,solidnum], self.qL[i+1,:,chemnum,solidnum]))
                        else:
                            if self.variable == 'Concentration':
                                data  = self.C
                                dataL = self.CL
                            if self.variable == 'Flux':
                                data  = self.F
                                dataL = self.F
                            if self.variable == 'Total concentration':
                                data  = self.W
                                dataL = self.WL
                            self.plot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], data[i,:,chemnum], data[i+1,:,chemnum]))
                            self.Lplot.append(time_interpolate(self.spatialplotdatas[n].value, self.times[i], self.times[i+1] - self.times[i], dataL[i,:,chemnum], dataL[i+1,:,chemnum]))

                        #if self.variable == 'Material fraction':
                        #    self.spatialplotdatas[n].legend = self.spatialplotdatas[n].component + ' ' +  str(self.spatialplotdatas[n].value) + ' ' + self.timeunit
                        #else:
                        #    self.spatialplotdatas[n].legend = self.spatialplotdatas[n].name + ' ' +  str(self.spatialplotdatas[n].value) + ' ' + self.timeunit

                    if self.variable == 'Flux':
                        Fmax = 0
                        Fmin = 0
                        for n in range(len(self.spatialplotdatas)):
                            if self.plot[n].max() < 0:     Fmax = max(Fmax, 0.8 * self.plot[n].max())
                            else:                          Fmax = max(Fmax, 1.2 * self.plot[n].max())

                            if self.plot[n].min() < 0:     Fmin = min(Fmin, 1.2 * self.plot[n].min())
                            else:                          Fmin = min(Fmin, 0.8 * self.plot[n].min())

                        self.axislimit[0]  = Fmin
                        self.axislimit[1]  = Fmax
                    else:
                        self.axislimit[0]  = 0
                        self.axislimit[1]  = 1.2  * max([self.plot[n].max() for n in range(len(self.spatialplotdatas))])

                else:
                    self.axislimit[0]  = 0
                    self.axislimit[1]  = self.times[-1]
                    self.plotinterest  = []
                    plotmax = 0

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

                    for n in range(len(self.timeplotdatas)):
                        chemnum  = self.chemicals.index(self.timeplotdatas[n].name)
                        compnum  = self.components.index(self.timeplotdatas[n].component)
                        solidnum = self.solids.index(self.timeplotdatas[n].solid)
                        self.plotinterest.append([])

                        if self.variable == 'Material fraction':
                            for i in range (len(self.times)):
                                if self.dep == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                    self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep*self.times[i], (self.z[0]-self.z[1]), self.z, self.Fi[i,:, compnum]))
                                else:
                                    self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z, self.Fi[i,:,compnum]))

                        elif self.variable == 'Solid concentration':
                            for i in range (len(self.times)):
                                if self.dep == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                    self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep*self.times[i], (self.z[0]-self.z[1]), self.z, self.q[i,:,chemnum,solidnum]))
                                else:
                                    self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z, self.q[i,:,chemnum,compnum]))

                        elif self.variable == 'Water concentration':
                            for i in range (len(self.times)):
                                self.plotinterest[-1].append(self.Cw[i, chemnum])

                        elif self.variable == 'Mass':
                            for i in range (len(self.times)):
                                if self.timeplotdatas[n].layer == 'All layers':
                                    self.plotinterest[-1].append(0)
                                    for j in range(len(self.output.layers + self.output.deplayer)):
                                        self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.M[i,j,chemnum]
                                elif self.timeplotdatas[n].layer == 'Top accumulation':
                                    self.plotinterest[-1].append(self.FM[i,0,chemnum])
                                elif self.timeplotdatas[n].layer == 'Bottom accumulation':
                                    self.plotinterest[-1].append(self.FM[i,1,chemnum])
                                elif self.timeplotdatas[n].layer == 'Deposition':
                                    self.plotinterest[-1].append(self.DM[i,chemnum])
                                elif self.timeplotdatas[n].layer == 'Reaction/Generation':
                                    self.plotinterest[-1].append(0)
                                    for j in range(len(self.output.layers + self.output.deplayer)):
                                        self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.RM[i,j,chemnum]
                                elif self.timeplotdatas[n].layer == 'Universal':
                                    self.plotinterest[-1].append(0)
                                    for j in range(len(self.output.layers + self.output.deplayer)):
                                        self.plotinterest[-1][-1] = self.plotinterest[-1][-1] + self.M[i,j,chemnum]
                                        self.plotinterest[-1][-1] = self.plotinterest[-1][-1] - self.RM[i,j,chemnum]
                                    self.plotinterest[-1][-1] = self.plotinterest[-1][-1]+self.FM[i,0,chemnum]-self.FM[i,1,chemnum]-self.DM[i,chemnum]
                                else:
                                    self.plotinterest[-1].append(self.M[i, self.layer_list.index(self.timeplotdatas[n].layer),chemnum])
                        else:
                            if self.variable == 'Concentration':       data = self.C
                            if self.variable == 'Flux':                data = self.F
                            if self.variable == 'Total concentration': data = self.W

                            for i in range (len(self.times)):
                                if self.dep == 1 and self.timeplotdatas[n].type == self.timeplotdatas[n].types[1]:
                                    self.plotinterest[-1].append(interp_dep(self.timeplotdatas[n].value, -self.Vdep*self.times[i], (self.z[0]-self.z[1]), self.z, data[i,:,chemnum]))
                                else:
                                    self.plotinterest[-1].append(interp(self.timeplotdatas[n].value, self.z, data[i,:,chemnum]))

                        #if self.variable == 'Material fraction':
                        #    self.timeplotdatas[n].legend = self.timeplotdatas[n].component + ' ' +  str(self.timeplotdatas[n].value) + ' ' + self.lengthunit
                        #else:
                        #    self.timeplotdatas[n].legend = self.timeplotdatas[n].name + ' ' +  str(self.timeplotdatas[n].value) + ' ' + self.lengthunit

                        if max(self.plotinterest[-1]) > plotmax: plotmax = max(self.plotinterest[-1])

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
            self.top.make_window(FigureEditor(self.top, self.system, self.type, self.variable, self.sketch, self.axislimit, self.figuresize, self.legendposition, self.fontsize, self.figuretexts))

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

                self.updategraph()

            if self.top is not None:
                self.top.destroy()
                self.top = None

        elif self.top is not None:
            tkmb.showerror(title = self.system.version, message = 'Please close the existing parameter input window first.')
            self.top.tk.focus()

    def updategraph(self, event = None):

        try: self.graphlabel.grid_forget()
        except: pass

        self.make_graphs()

        self.graph      = Image.open(Filepath + r'/output/%s.png' % self.filename)
        self.image      = ImageTk.PhotoImage(self.graph)
        self.graphlabel = Label(self.frame, image = self.image, bg = self.color)

        self.graphlabel.grid(  row = 0, column = 1, padx = 40, sticky = 'WE')

        self.master.geometry()
        self.master.center()

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


    def export_data(self, event = None):

        if self.top is None:
            self.top = CapSimWindow(master = self.master, buttons = 2)
            self.top.make_window(ResultExporter(self.top, self.system, self.output))
            self.top.mainloop()

            if self.top.window.cancelflag == 0:
                filename = self.top.window.filepath.get() + self.top.window.name.get()
                if self.top.window.type.get() == self.top.window.types[0]:
                    self.make_csv_file(filename,  'Regular')

            if self.top is not None:
                self.top.destroy()
                self.top = None

        elif self.top is not None:
            tkmb.showerror(title = self.system.version, message = 'Please close the existing result export window first.')
            self.top.tk.focus()

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
            self.make_sketch_image()
        else:
            self.outputgraphs = plt.figure(figsize=(self.figuresize[0]/100, self.figuresize[1]/100))

        if self.type == self.types[0]:
            self.make_profiles(self.plot, self.Lplot)
        elif self.type == self.types[1]:
            self.make_time_profiles(self.plotinterest)

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

    def make_sketch_image(self):
        """Makes plots of concentration profiles over time."""

        Sketchprofiles  = self.outputgraphs.add_subplot(self.gs[0])

        Htotal = self.output.z[-1] - self.output.z[0]

        sfont = {9:6, 10:7, 11:7.5, 12:8, 13:9}

        H = 0
        color_step = round(0.5/len(self.system.layers),2)
        Blocks = []
        for layer in self.system.layers:
            if layer.number <> 0:
                Blocks.append(patches.Rectangle((0.0, H), 0.99, layer.h, facecolor = str( 0.75 - color_step * layer.number), edgecolor = 'k', linestyle = 'solid', linewidth = '1'))
                Sketchprofiles.add_patch(Blocks[-1])
                if self.sketch['Label'] == 1:
                    rx, ry = Blocks[-1].get_xy()
                    cx = rx + Blocks[-1].get_width()/2.0
                    cy = ry + Blocks[-1].get_height()/2.0
                    if Blocks[-1].get_height()/Htotal*self.figuresize[1] >= 14:
                        Sketchprofiles.annotate(fill(layer.type, 10), (cx, cy), color='w', fontsize = 10, ha='center', va='center')
                    elif Blocks[-1].get_height()/Htotal*self.figuresize[1] >= 9:
                        Sketchprofiles.annotate(fill(layer.type, 10), (cx, cy), color='w', fontsize = sfont[int(Blocks[-1].get_height()/Htotal*self.figuresize[1])], ha='center', va='center')
                H = H + layer.h

        if self.system.dep == 'Deposition':
            Blockdep = patches.Polygon(array([[0.995, self.output.z[0]], [0.0, 0.0], [0.995, 0.0]]), facecolor = '0.75', edgecolor = 'k', linestyle = 'solid', linewidth = '1')
            Sketchprofiles.add_patch(Blockdep)
            if self.sketch['Label'] == 1:
                rx, ry = 0.0, 0.0
                cx = rx + 0.995/2.0
                cy = ry + self.output.z[0]/2.0
                if abs(self.output.z[0])/Htotal*self.figuresize[1] >= 14:
                    Sketchprofiles.annotate(fill(self.system.layers[0].type, 10), (cx, cy), color='k', fontsize = 10, ha='center', va='center')
                elif abs(self.output.z[0])/Htotal*self.figuresize[1] >= 9:
                    Sketchprofiles.annotate(fill(self.system.layers[0].type, 10), (cx, cy), color='k', fontsize = sfont[int(abs(self.output.z[0])/Htotal*self.figuresize[1])], ha='center', va='center')

        if self.system.dep == 'Deposition' and self.sketch['Time'] == 1:
            Sketchprofiles.annotate('0', (0, H),                                                color='k', fontsize = 10, ha='center', va='top')
            Sketchprofiles.annotate(fill(str(self.output.times[-1]) + self.system.timeunit, 6), (1.0, H),color='k', fontsize = 10, ha='center', va='top')

        if self.type == self.types[0]:
            if self.sketch['Time'] == 1:
                if self.system.dep == 'Deposition':
                    for plotdata in self.spatialplotdatas:
                        Sketchprofiles.plot( [plotdata.value/self.output.times[-1], plotdata.value/self.output.times[-1]],[0-plotdata.value*self.system.Vdep, self.output.z[-1]])
                else:
                    for plotdata in self.spatialplotdatas:
                        Sketchprofiles.plot( [plotdata.value/self.output.times[-1], plotdata.value/self.output.times[-1]], [0, self.output.z[-1]])

        if self.type == self.types[1]:
            if self.sketch['Depth'] == 1:
                for plotdata in self.timeplotdatas:
                    if plotdata.type == plotdata.types[0]:
                        Sketchprofiles.plot([0.0, 1.0], [plotdata.value, plotdata.value])
                    else:
                        Sketchprofiles.plot([0.0, 1.0], [plotdata.value, plotdata.value + self.output.z[0]])

        if self.system.bio <> 'None':
            Sketchprofiles.plot([0.0, 1.0], [self.system.hbio, self.system.hbio + self.output.z[0]], linestyle = ':', color = 'k')

        Sketchprofiles.axis('off')
        Sketchprofiles.set_ylim([max(self.output.z), min(self.output.z)])

    def make_profiles(self, plots, Lplots):
        """Makes plots of concentration profiles over time."""

        if self.sketch['Sketch'] == self.sketch['Sketches'][0]: Graph  = self.outputgraphs.add_subplot(self.gs[1])
        else:                                                   Graph  = self.outputgraphs.add_subplot(111)

        profiles  = []
        tlegend   = []

        for n in range(len(self.spatialplotdatas)):
            if self.dep == 1:
                time = self.spatialplotdatas[n].value
                dep = -self.Vdep * time
                i = 0
                while round(self.z[i] + 0.5*(self.z[1]-self.z[0]), 8) < round(dep, 8)  : i = i + 1
                if round(-dep, 8) < 1.5 * (self.z[1]-self.z[0]) and round(-dep, 8) >= 0.5 * (self.z[1]-self.z[0]):
                    Cs  = plots[n][i+2:]
                    CLs = Lplots[n][i+2:]
                    zs  = self.z[i+2:]
                else:
                    Cs = plots[n][i+1:]
                    CLs = Lplots[n][i+1:]
                    zs = self.z[i+1:]
                zs[0] = -self.Vdep * time
            else:
                Cs  = plots[n][:]
                CLs = Lplots[n][:]
                zs  = self.z[:]

            Cs_plot = []
            zs_plot = []

            for i in range(len(zs)):
                if Cs[i] <> CLs[i]:
                    zs_plot.append(zs[i])
                    zs_plot.append(zs[i])
                    Cs_plot.append(Cs[i])
                    Cs_plot.append(CLs[i])
                else:
                    zs_plot.append(zs[i])
                    Cs_plot.append(Cs[i])

            color     = self.colors[self.colornames.index(self.spatialplotdatas[n].color)]
            style     = self.styles[self.stylenames.index(self.spatialplotdatas[n].style)]
            linewidth = self.spatialplotdatas[n].size
            if self.stylenames.index(self.spatialplotdatas[n].style) < 4:
                p, = Graph.plot(Cs_plot, zs_plot, color = color, linestyle = style, linewidth = linewidth)
            else:
                p, = Graph.plot(Cs_plot, zs_plot, color = color, linestyle = 'None', marker = style,  markersize = linewidth)
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

    def make_time_profiles(self, plots):
        """Makes a graph of the pore water concentration at the depth of
        interest."""

        if self.sketch['Sketch'] == self.sketch['Sketches'][0]: Time_graph  = self.outputgraphs.add_subplot(self.gs[1])
        else:                                                   Time_graph  = self.outputgraphs.add_subplot(111)

        zlegend         = []
        profiles        = []
        max_value = 0

        for n in range(len(self.timeplotdatas)):
            color     = self.colors[self.colornames.index(self.timeplotdatas[n].color)]
            style     = self.styles[self.stylenames.index(self.timeplotdatas[n].style)]
            linewidth = self.timeplotdatas[n].size

            if self.stylenames.index(self.timeplotdatas[n].style) < 4:
                p, = Time_graph.plot(self.times, plots[n], color = color, linestyle = style, linewidth = linewidth)
            else:
                p, = Time_graph.plot(self.times, plots[n], color = color, linestyle = 'None', marker = style,  markersize = linewidth)

            if self.timeplotdatas[n].legend != '':
                profiles.append(p)
                zlegend.append(self.timeplotdatas[n].legend)

            for i in range (len(self.times)):
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

    def modify(self, event = None):
        """Takes the user back to the Summary window for the simulation."""

        self.frame.quit()

    def make_csv_file(self, filename, output_type):
        """This subroutine is used to create the output file for CapSim."""

        lengthunits = ['um', 'cm', 'm']
        concunits   = ['ug/L', 'mg/L', 'g/L', 'umol/L', 'mmol/L', 'mol/L']
        diffunits   = ['cm^2/s', 'cm^2/yr']

        timeunit    = self.timeunit

        if self.lengthunits.count(self.lengthunit) > 0:   lengthunit     = lengthunits[self.lengthunits.index(self.lengthunit)]
        else:                                             lengthunit     = self.lengthunit
        if self.concunits.count(self.concunit) > 0:       concunit       = concunits[self.concunits.index(self.concunit)]
        else:                                             concunit       = self.concunit
        if self.diffunits.count(self.diffunit) > 0:       diffunit       = diffunits[self.system.diffunits.index(self.system.diffunit)]
        else:                                             diffunit       = self.diffunit

        fluxunit    = concunit[:-2]+ ',/,'+ lengthunit + '^2' + ',/,' + timeunit

        file = open(filename + '.csv', 'w')

        file.write('CapSim,'+self.system.version + '\n\n')

        file.write('Output type,'+output_type + '\n\n')

        file.write('System Units' + '\n')
        file.write(',Length: ,'         + lengthunit    + '\n')
        file.write(',Concentration: ,'  + concunit      + '\n')
        file.write(',Time: ,'           + timeunit      + '\n')
        file.write(',Diffusivity: ,'    + diffunit      + '\n')
        file.write('\n')

        file.write('Chemicals,'+str(self.system.nchemicals)  + '\n')
        file.write(', Name, Formula, MW, Temperature, Dw, Koc, Kdoc, Reference\n')
        for chemical in self.system.chemicals:
            file.write(',' + text_converter(chemical.name) + ',')
            file.write(text_converter(chemical.formula) + ',')
            file.write(str(chemical.MW) + ',')
            file.write(str(chemical.temp) + ',')
            file.write(str(chemical.Dw) + ',')
            file.write(str(chemical.Koc) + ',')
            file.write(str(chemical.Kdoc) + ',')
            file.write(str(chemical.Ref) + ',')
            file.write('\n')
        file.write('\n')

        file.write('Matrices,' + str(self.system.nmatrices)+ '\n')
        file.write(', Name, Components, Porosity, Bulk density(L/kg), foc, Mixture\n')
        for matrix in self.system.matrices:
            file.write(',' + matrix.name            + ',')
            file.write(str(len(matrix.components))  + ',')
            file.write(str(matrix.e)                + ',')
            file.write(str(matrix.rho)              + ',')
            file.write(str(matrix.foc)              + ',')
            file.write(matrix.model                 + ',')
            file.write('\n')
        file.write('\n')

        file.write('Matrix Components,' + str(len(self.system.components))+'\n')
        file.write(',Matrix, Component, Weight fraction, Porosity, Bulk density(L/kg), foc\n')
        for matrix in self.system.matrices:
            for component in matrix.components:
                file.write( ',' + str(matrix.name)  + ',')
                file.write(component.name           + ',')
                file.write(str(component.mfraction) + ',')
                file.write(str(component.e)         + ',')
                file.write(str(component.rho)       + ',')
                file.write(str(component.foc)       + ',')
                file.write('\n')
        file.write('\n')

        file.write('Sorptions,'+ '\n')
        file.write(',Component, Chemical, Isotherm, Kinetic, Kd/Koc/Kf/qmax, foc/N/b, ksorp, kdesorp, equilibrium time\n')
        for component in self.system.components:
            for chemical in self.system.chemicals:
                file.write(',' + component.name + ',')
                file.write(text_converter(chemical.name) + ',')
                file.write(self.system.sorptions[component.name][chemical.name].isotherm         + ',')
                file.write(self.system.sorptions[component.name][chemical.name].kinetic          + ',')
                if self.system.sorptions[component.name][chemical.name].isotherm == self.system.sorptions[component.name][chemical.name].isotherms[0]:
                    file.write(str(self.system.sorptions[component.name][chemical.name].K)      + ',')
                    file.write(',')
                elif self.system.sorptions[component.name][chemical.name].isotherm == self.system.sorptions[component.name][chemical.name].isotherms[1]:
                    file.write(str(self.system.sorptions[component.name][chemical.name].Koc)     + ',')
                    file.write(str(self.system.sorptions[component.name][chemical.name].foc)     + ',')
                elif self.system.sorptions[component.name][chemical.name].isotherm == self.system.sorptions[component.name][chemical.name].isotherms[2]:
                    file.write(str(self.system.sorptions[component.name][chemical.name].Kf)      + ',')
                    file.write(str(self.system.sorptions[component.name][chemical.name].N)       + ',')
                else:
                    file.write(str(self.system.sorptions[component.name][chemical.name].qmax)    + ',')
                    file.write(str(self.system.sorptions[component.name][chemical.name].b)       + ',')

                file.write(str(self.system.sorptions[component.name][chemical.name].ksorp)   + ',')
                file.write(str(self.system.sorptions[component.name][chemical.name].kdesorp) + ',')
                file.write(str(self.system.sorptions[component.name][chemical.name].thalf)   + ',')
                file.write('\n')
        file.write('\n')

        file.write('Layers,'+str(self.system.nlayers)  + '\n')
        file.write(', Name, Matrix, Thickness, Tortuosity, Hydraulic Dispersivity, DOC concentration\n')
        for layer in self.system.layers:
            file.write(',' + layer.name + ',')
            file.write(layer.type + ',')
            file.write(str(layer.h) + ',')
            file.write(layer.tort + ',')
            file.write(str(layer.alpha) + ',')
            file.write(str(layer.doc) + ',')
            file.write('\n')
        file.write('\n')

        file.write('Reactions,'+str(len(self.system.reactions))  + '\n')
        for reaction in self.system.reactions:
            file.write(',' + str(reaction.number) + ',')
            file.write(reaction.name  + ',')
            file.write(reaction.model + ',')
            file.write(text_converter(reaction.equation) + ',')
            file.write('\n')
            file.write(',,Reactants,'+str(len(reaction.reactants)) + ',' * (len(reaction.reactants))  + 'Products,' + str(len(reaction.products))+ '\n')

            file.write(',Name,')
            for reactant in reaction.reactants:
                file.write(text_converter(reactant.name) + ',')
            file.write('-->,')
            for product in reaction.products:
                file.write(text_converter(product.name) + ',')
            file.write('\n')

            file.write(',Formula,')
            for reactant in reaction.reactants:
                file.write(text_converter(reactant.formula) + ',')
            file.write('-->,')
            for product in reaction.products:
                file.write(text_converter(product.formula) + ',')
            file.write('\n')

            file.write(',MW,')
            for reactant in reaction.reactants:
                file.write(str(reactant.MW) + ',')
            file.write('-->,')
            for product in reaction.products:
                file.write(str(product.MW) + ',')
            file.write('\n')

            file.write(',Coef,')
            for reactant in reaction.reactants:
                file.write(str(reactant.coef) + ',')
            file.write('-->,')
            for product in reaction.products:
                file.write(str(product.coef) + ',')
            file.write('\n')

            file.write(',Index,')
            for reactant in reaction.reactants:
                file.write(str(reactant.index) + ',')
            file.write('-->,')
            file.write('\n')
        file.write('\n')


        file.write('Reaction Coefficients,'+str(len(self.system.reactions))  + '\n')
        file.write(',Layer, Reaction, Coefficient\n')
        for layer in self.system.layers:
            for reaction in self.system.reactions:
                file.write(',' + layer.name + ',')
                file.write(reaction.name + ',')
                file.write(str(self.system.coefficients[layer.name][reaction.name].lam) + ',')
                file.write('\n')
        file.write('\n')

        file.write('System Properties\n')
        file.write(',Advection,'            +self.system.adv+'\n')
        file.write(',Darcy Velocity,'       +str(self.system.Vdar)+'\n')
        file.write(',Tidal Max Velocity,'   +str(self.system.Vtidal)+'\n')
        file.write(',Tidal Period,'         +str(self.system.ptidal)+'\n')
        file.write(',Bioturbation,'         +self.system.bio+'\n')
        file.write(',Bioturbation Depth,'   +str(self.system.hbio)+'\n')
        file.write(',Bioturbation Gaussian model coefficient,'   +str(self.system.sigma)+'\n')
        file.write(',Particle biodiffusion coefficient,'        +str(self.system.Dbiop)+'\n')
        file.write(',Porewater biodiffusion coeffcient,'       +str(self.system.Dbiopw)+'\n')
        file.write(',Consolidation,'                +self.system.con+'\n')
        file.write(',Consolidation thickness,'      +str(self.system.hcon)+'\n')
        file.write(',Time to reach 90% Consolidation,'       +str(self.system.t90)+'\n')
        file.write('\n')

        file.write('Auxiliary Conditions\n')
        file.write(',Layer, Type, Parameter,')
        for chemical in self.system.chemicals:
            file.write(text_converter(chemical.name) + ',')
        solidcheck = 0
        for component in self.system.components:
            for chemical in self.system.chemicals:
                if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                    file.write(component.name +'_'+text_converter(chemical.name) + ',')
        file.write('\n')

        file.write(', Benthic,'+ self.system.topBCtype + ',')
        if self.system.topBCtype == 'Fixed Concentration':
            file.write('Surface concentration,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].Co) + ',')
        elif self.system.topBCtype == 'Mass transfer':
            file.write('Mass transfer coefficient,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].k) + ',')
            file.write('\n')
            file.write(',,,Water concentration,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].Cw) + ',')
        else:
            file.write('Mass transfer coefficient,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].k) + ',')
            file.write('\n')
            file.write(',,,Water column retention time,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].tau) + ',')

        file.write('\n\n')

        for layer in self.system.layers:
            file.write(','+layer.name+','+ layer.ICtype + ',')
            if layer.ICtype == 'Uniform':
                file.write('Initial concentration,')
                for chemical in self.system.chemicals:
                    file.write(str(self.system.ICs[layer.name][chemical.name].uniform) + ',')
                for component in self.system.matrices[layer.type_index].components:
                    for chemical in self.system.chemicals:
                        if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            file.write(str(self.system.SolidICs[layer.name][component.name][chemical.name].uniform) + ',')
            else:
                file.write('Top concentration,')
                for chemical in self.system.chemicals:
                    file.write(str(self.system.ICs[layer.name][chemical.name].top) + ',')
                for component in self.system.matrices[layer.type_index].components:
                    for chemical in self.system.chemicals:
                        if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            file.write(str(self.system.SolidICs[layer.name][component.name][chemical.name].uniform) + ',')
                file.write('\n')
                file.write(',,,Bottom concentration,')
                for chemical in self.system.chemicals:
                    file.write(str(self.system.ICs[layer.name][chemical.name].bot) + ',')
                for component in self.system.matrices[layer.type_index].components:
                    for chemical in self.system.chemicals:
                        if self.system.sorptions[component.name][chemical.name].kinetic == 'Transient':
                            file.write(str(self.system.SolidICs[layer.name][component.name][chemical.name].uniform) + ',')
            file.write('\n')
        file.write('\n')

        file.write(', Underlying,'+ self.system.botBCtype + ',')
        if self.system.botBCtype == 'Fixed Concentration' or self.system.botBCtype == 'Flux-matching':
            file.write('Concentration,')
            for chemical in self.system.chemicals:
                file.write(str(self.system.BCs[chemical.name].Cb) + ',')
        else:
            file.write('No flux')
        file.write('\n\n')

        file.write('Solver Options, 1\n')
        file.write(',Simulation time,'+str(self.system.tfinal)+'\n')
        file.write(',Output time steps,'+str(self.system.outputsteps)+'\n')
        file.write(',Discrete time option,'+str(self.system.timeoption)+'\n')
        file.write(',Discrete,'+self.system.discrete + '\n')
        file.write(',Total number of grid points,'+str(self.system.ptotal)+'\n')
        file.write(',Time step,'+str(self.system.delt)+'\n')
        file.write(',Grid option,'+self.system.ptype+'\n')
        file.write(',Time step option,'+self.system.tvariable+'\n')
        file.write(',,')
        for i in range(self.system.nlayers):
            file.write(self.system.layers[i].name + ',')
        file.write('\n')
        file.write(',Layer Grid size,')
        for i in range(self.system.nlayers):
            file.write(str(self.system.delz[i]) +',')
        file.write('\n')
        file.write(',Layer Number of Grids,')
        for i in range(self.system.nlayers):
            file.write( str(self.system.players[i]) + ',')
        file.write('\n')
        file.write(',Time step in a tidal circle,'+str(self.system.tidalsteps)+'\n')
        file.write(',nonlinear options,'+self.system.nonlinear+'\n')
        file.write(',nonlinear error(%),'+str(self.system.nlerror)+'\n')
        file.write(',Steps in deposition grid,'+str(self.system.depgrid)+'\n')
        file.write(',Oscillation output options,'+str(self.system.averageoption)+'\n')
        file.write(',Mass tracking option,'+str(self.system.massoption)+'\n')
        file.write('\n\n')

        file.write('Overlying Water Column Properties,')
        if self.system.topBCtype == 'Finite mixing water column':   file.write('Involved\n')
        else:                                                       file.write('Not applicable\n')
        file.write(',Inflow rate (m^3/s),'+str(self.system.taucoefs['Q'])+'\n')
        file.write(',Water body volume (m^3),'+str(self.system.taucoefs['V'])+'\n')
        file.write(',Water body depth (m),'+str(self.system.taucoefs['h'])+'\n')
        file.write(',Water evaporation rate (m^3/s),'+str(self.system.taucoefs['Qevap'])+'\n')
        file.write(',Water body DOC concentration (mg/L),'+str(self.system.taucoefs['DOC']) + '\n')
        file.write(',,')
        for chemical in self.system.chemicals:
            file.write(chemical.name + ',')
        file.write('\n')
        file.write(',Contaminant decay ('+timeunit+'),')
        for chemical in self.system.chemicals:
            file.write(str(self.system.BCs[chemical.name].kdecay) + ',')
        file.write('\n')
        file.write(',Contaminant evaporation ('+ lengthunit+ '/'+ timeunit + '),')
        for chemical in self.system.chemicals:
            file.write(str(self.system.BCs[chemical.name].kevap)+ ',')
        file.write('\n')

        file.write('\n')

        Double_z  = []
        for target_name in self.chemicals:
            n = self.chemicals.index(target_name)
            for j in range(len(self.output.times)):
                for i in range(len(self.output.z)):
                    if self.output.C[j, i, n] <> self.output.CL[j, i, n] and Double_z.count(self.output.z[i]) == 0:
                        Double_z.append(self.output.z[i])
                    if self.output.W[j, i, n] <> self.output.WL[j, i, n] and Double_z.count(self.output.z[i]) == 0:
                        Double_z.append(self.output.z[i])
                    if self.output.q[j, i, n] <> self.output.qL[j, i, n] and Double_z.count(self.output.z[i]) == 0:
                        Double_z.append(self.output.z[i])
                    for m in range(len(self.system.components)):
                        if self.output.qm[j, i, n, m] <> self.output.qmL[j, i, n, m] and Double_z.count(self.output.z[i]) == 0:
                            Double_z.append(self.output.z[i])

        for target_name in self.system.component_list:
            n = self.system.component_list.index(target_name)
            for j in range(len(self.output.times)):
                for i in range(len(self.output.z)):
                    if self.output.Fi[j, i, n] <> self.output.Fi[j, i, n] and Double_z.count(self.output.z[i]) == 0:
                        Double_z.append(self.output.z[i])

        file.write('\n')
        file.write('Output Results\n')
        file.write('Output Grid,'+str(len(self.output.z)+len(Double_z))+'\n')
        file.write('Output Time steps,'+str(len(self.output.times))+'\n')
        file.write('\n')

        for target_name in self.chemicals:
            n = self.chemicals.index(target_name)

            file.write('Chemical,' + target_name +'\n')

            title    = []
            data     = []
            dataL    = []

            title.append('Porewater concentrations,' + concunit)
            data.append(self.output.C)
            dataL.append(self.output.CL)
            if output_type == 'Monte Carlo Simulation':
                title.append('Porewater concentrations standard deviation,' + concunit)
                title.append('Porewater concentrations maximum,' + concunit)
                title.append('Porewater concentrations minimum,' + concunit)
                data.append(self.output.C_std)
                data.append(self.output.C_max)
                data.append(self.output.C_min)
                dataL.append(self.output.CL_std)
                dataL.append(self.output.CL_max)
                dataL.append(self.output.CL_min)
            title.append('Fluxes,' + fluxunit)
            data.append(self.output.F)
            dataL.append(self.output.F)
            if output_type == 'Monte Carlo Simulation':
                title.append('Fluxes standard deviation,' + fluxunit)
                title.append('Fluxes maximum,' + fluxunit)
                title.append('Fluxes minimum,' + fluxunit)
                data.append(self.output.F_std)
                data.append(self.output.F_max)
                data.append(self.output.F_min)
                dataL.append(self.output.F_std)
                dataL.append(self.output.F_max)
                dataL.append(self.output.F_min)
            title.append('Total concentrations,' + concunit)
            data.append(self.output.W)
            dataL.append(self.output.WL)
            if output_type == 'Monte Carlo Simulation':
                title.append('Total concentrations standard deviation,' + concunit)
                title.append('Total concentrations maximum,' + concunit)
                title.append('Total concentrations minimum,' + concunit)
                data.append(self.output.W_std)
                data.append(self.output.W_max)
                data.append(self.output.W_min)
                dataL.append(self.output.WL_std)
                dataL.append(self.output.WL_max)
                dataL.append(self.output.WL_min)

            title.append('Solid concentrations,' + 'Total solid,' + concunit[:-1] + 'kg')
            data.append(self.output.q)
            dataL.append(self.output.qL)
            if output_type == 'Monte Carlo Simulation':
                title.append('Solid concentrations standard deviation,' + 'Total solid,' + concunit[:-1] + 'kg')
                title.append('Solid concentrations maximum,' + 'Total solid,' + concunit[:-1] + 'kg')
                title.append('Solid concentrations minimum,' + 'Total solid,' + concunit[:-1] + 'kg')
                data.append(self.output.q_std)
                data.append(self.output.q_max)
                data.append(self.output.q_min)
                dataL.append(self.output.qL_std)
                dataL.append(self.output.qL_max)
                dataL.append(self.output.qL_min)
            for component in self.system.components:
                title.append('Solid concentrations,' + component.name + ',' + concunit[:-1] + 'kg' )
                data.append(self.output.qm[:,:,:,self.system.components.index(component)])
                dataL.append(self.output.qmL[:,:,:,self.system.components.index(component)])
                if output_type == 'Monte Carlo Simulation':
                    title.append('Solid concentrations standard deviation,' + component.name + ',' + concunit[:-1] + 'kg' )
                    title.append('Solid concentrations maximum,' + component.name + ',' + concunit[:-1] + 'kg' )
                    title.append('Solid concentrations minimum,' + component.name + ',' + concunit[:-1] + 'kg' )
                    data.append(self.output.qm_std[:,:,:,self.system.components.index(component)])
                    data.append(self.output.qm_max[:,:,:,self.system.components.index(component)])
                    data.append(self.output.qm_min[:,:,:,self.system.components.index(component)])
                    dataL.append(self.output.qmL_std[:,:,:,self.system.components.index(component)])
                    dataL.append(self.output.qmL_max[:,:,:,self.system.components.index(component)])
                    dataL.append(self.output.qmL_min[:,:,:,self.system.components.index(component)])

            for ii in range(len(title)):
                file.write(title[ii] + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.z)):
                    if Double_z.count(self.output.z[i]) > 0:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % data[ii][j, i, n])
                        file.write('\n')
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % dataL[ii][j, i, n])
                    else:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % data[ii][j, i, n])
                    file.write('\n')
                file.write('\n')

            if self.system.topBCtype == 'Finite mixed water column':
                file.write('Overlying water concentrations,' + concunit + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write(',')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.Cw[j, n])
                file.write('\n')

                if output_type == 'Monte Carlo Simulation':
                    file.write('Overlying water concentrations standard deviation,' + concunit + '\n,Times\n Depths,')
                    for t in self.output.times: file.write('%.3e,' % t)
                    file.write('\n')
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.Cw_std[j, n])
                    file.write('\n')

                    file.write('Overlying water concentrations maximum,' + concunit + '\n,Times\n Depths,')
                    for t in self.output.times: file.write('%.3e,' % t)
                    file.write('\n')
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.Cw_max[j, n])
                    file.write('\n')

                    file.write('Overlying water concentrations minimum,' + concunit + '\n,Times\n Depths,')
                    for t in self.output.times: file.write('%.3e,' % t)
                    file.write('\n')
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.Cw_min[j, n])
                    file.write('\n')

            elif self.system.topBCtype == 'Mass transfer':
                file.write('Overlying water concentrations,' + concunit + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write(',')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.system.BCs[target_name].Cw)
                file.write('\n')

            else:
                file.write('Overlying water concentrations,' + concunit + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write(',')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.system.BCs[target_name].Co)
                file.write('\n')

            file.write('\n')

            file.write('Mass in layers,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
            for t in self.output.times: file.write('%.3e,' % t)
            file.write('\n')
            for i in range(len(self.output.deplayer+self.output.layers)):
                file.write((self.output.deplayer+self.output.layers)[i].name)
                file.write(',')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.M[j, i, n])
                file.write('\n')
            file.write('\n')

            if output_type == 'Monte Carlo Simulation':
                file.write('Mass in layers standard deviation,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.M_std[j, i, n])
                    file.write('\n')
                file.write('\n')

                file.write('Mass in layers maximum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.M_max[j, i, n])
                    file.write('\n')
                file.write('\n')

                file.write('Mass in layers minimum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.M_min[j, i, n])
                    file.write('\n')
                file.write('\n')

            file.write('Mass in top/bottom/deposition,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
            for t in self.output.times: file.write('%.3e,' % t)
            file.write('\n')
            file.write('Top,')
            for j in range(len(self.output.times)):
                file.write('%.3e,' % self.output.FM[j, 0, n])
            file.write('\n')
            file.write('Bottom,')
            for j in range(len(self.output.times)):
                file.write('%.3e,' % self.output.FM[j, 1, n])
            file.write('\n')
            file.write('Deposition,')
            for j in range(len(self.output.times)):
                file.write('%.3e,' % self.output.DM[j, n])
            file.write('\n')
            file.write('\n')

            if output_type == 'Monte Carlo Simulation':
                file.write('Mass in top/bottom/deposition standard deviation,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write('Top,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_std[j, 0, n])
                file.write('\n')
                file.write('Bottom,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_std[j, 1, n])
                file.write('\n')
                file.write('Deposition,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.DM_std[j, n])
                file.write('\n')
                file.write('\n')

                file.write('Mass in top/bottom/deposition maximum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write('Top,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_max[j, 0, n])
                file.write('\n')
                file.write('Bottom,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_max[j, 1, n])
                file.write('\n')
                file.write('Deposition,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.DM_max[j, n])
                file.write('\n')
                file.write('\n')

                file.write('Mass in top/bottom/deposition minimum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                file.write('Top,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_min[j, 0, n])
                file.write('\n')
                file.write('Bottom,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.FM_min[j, 1, n])
                file.write('\n')
                file.write('Deposition,')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.DM_min[j, n])
                file.write('\n')
                file.write('\n')

            file.write('Mass by reaction/generation,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
            for t in self.output.times: file.write('%.3e,' % t)
            file.write('\n')
            for i in range(len(self.output.deplayer+self.output.layers)):
                file.write((self.output.deplayer+self.output.layers)[i].name)
                file.write(',')
                for j in range(len(self.output.times)):
                    file.write('%.3e,' % self.output.RM[j, i, n])
                file.write('\n')
            file.write('\n')

            if output_type == 'Monte Carlo Simulation':
                file.write('Mass by reaction/generation standard deviation,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.RM_std[j, i, n])
                    file.write('\n')
                file.write('\n')
                file.write('Mass by reaction/generation maximum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.RM_max[j, i, n])
                    file.write('\n')
                file.write('\n')
                file.write('Mass by reaction/generation minimum,' + concunit[:-1] + 'cm^2' + '\n,Times\n,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.deplayer+self.output.layers)):
                    file.write((self.output.deplayer+self.output.layers)[i].name)
                    file.write(',')
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.RM_min[j, i, n])
                    file.write('\n')
                file.write('\n')

            file.write('\n')

        file.write('Solid material fractions\n')

        for target_name in self.system.component_list:
            n = self.system.component_list.index(target_name)
            file.write('Component,' + target_name + '\n')
            file.write('Volumetric fraction,' + '\n,Times\n Depths,')
            for t in self.output.times: file.write('%.3e,' % t)
            file.write('\n')
            for i in range(len(self.output.z)):
                if Double_z.count(self.output.z[i]) > 0:
                    file.write('%.3e,' % self.output.z[i])
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.Fi[j, i, n])
                    file.write('\n')
                    file.write('%.3e,' % self.output.z[i])
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.FiL[j, i, n])
                else:
                    file.write('%.3e,' % self.output.z[i])
                    for j in range(len(self.output.times)):
                        file.write('%.3e,' % self.output.Fi[j, i, n])
                file.write('\n')
            file.write('\n')

            if output_type == 'Monte Carlo Simulation':
                file.write('Volumetric fraction standard deviation,' + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.z)):
                    if Double_z.count(self.output.z[i]) > 0:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_std[j, i, n])
                        file.write('\n')
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.FiL_std[j, i, n])
                    else:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_std[j, i, n])
                    file.write('\n')
                file.write('\n')
                file.write('Volumetric fraction maximum,' + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.z)):
                    if Double_z.count(self.output.z[i]) > 0:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_max[j, i, n])
                        file.write('\n')
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.FiL_max[j, i, n])
                    else:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_max[j, i, n])
                    file.write('\n')
                file.write('\n')
                file.write('Volumetric fraction minimum,' + '\n,Times\n Depths,')
                for t in self.output.times: file.write('%.3e,' % t)
                file.write('\n')
                for i in range(len(self.output.z)):
                    if Double_z.count(self.output.z[i]) > 0:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_min[j, i, n])
                        file.write('\n')
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.FiL_min[j, i, n])
                    else:
                        file.write('%.3e,' % self.output.z[i])
                        for j in range(len(self.output.times)):
                            file.write('%.3e,' % self.output.Fi_min[j, i, n])
                    file.write('\n')
                file.write('\n')

class PlotEditor:

    def __init__(self, master, system, type, variable, chemicals, components, solids, dep, spatialplotdatas, timeplotdatas, outputz, outputt, layer_list):

        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.tframe    = Frame(master.tframe)
        self.frame     = Frame(master.frame)
        self.bframe    = Frame(master.bframe)
        self.top       = None
        self.flag      = 0
        self.system    = system

        self.layers = ['All layers', 'Top accumulation', 'Bottom accumulation', 'Deposition', 'Reaction/Generation', 'Universal']
        for layer in (layer_list):
            if layer == 'Deposition':
                self.layers.append('Deposition Layer')
            else:
                self.layers.append(layer)

        rgb            = self.master.frame.winfo_rgb(self.master.frame.cget('bg'))
        self.color     = '#%02x%02x%02x' %rgb
        self.mplrgb    = (rgb[0] / 65536., rgb[1] / 65536., rgb[2] / 65536.)

        self.types          = ['Spatial profile', 'Time profile']
        if system.biomix == 1:
            self.variables =      ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Material fraction']
            self.extravariables = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Mass', 'Water concentration', 'Material fraction']
            self.timevariables  = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Mass', 'Material fraction']
        else:
            self.variables =      ['Concentration', 'Flux', 'Solid concentration', 'Total concentration']
            self.extravariables = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Mass', 'Water concentration']
            self.timevariables  = ['Concentration', 'Flux', 'Solid concentration', 'Total concentration', 'Mass']

        self.type      = StringVar(value = type)
        self.variable  = StringVar(value = variable)

        self.chemicals  = chemicals
        self.components = components
        self.solids     = solids

        self.dep        = dep

        self.timeunit   = system.timeunit
        self.lengthunit = system.lengthunit

        self.cancelflag =  0

        self.spatialplotdatas   = [plotdata.copy() for plotdata in spatialplotdatas]
        self.timeplotdatas      = [plotdata.copy() for plotdata in timeplotdatas]

        for plotdata in (self.spatialplotdatas+self.timeplotdatas):
            plotdata.name     = StringVar(value = plotdata.name)
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
        self.valuecolumn    = Label(self.tframe, text = ' ', width = 12)
        self.matrixcolumn   = Label(self.tframe, text = ' ', width = 20)
        self.depthcolumn    = Label(self.tframe, text = ' ', width = 21)
        self.endcolumn      = Label(self.tframe, text = ' ', width = 2)

        self.typelabel      = Label(self.tframe, bg = self.color, text = 'Type:')
        self.varilabel      = Label(self.tframe, bg = self.color, text = 'Variable:')

        self.typewidget     = OptionMenu(self.tframe, self.type,     *self.types,           command = self.updateplots)
        self.variwidget     = OptionMenu(self.tframe, self.variable, *self.variables,       command = self.updateplots)
        self.extravariwidget= OptionMenu(self.tframe, self.variable, *self.extravariables,  command = self.updateplots)
        self.timevariwidget = OptionMenu(self.tframe,  self.variable, *self.timevariables,  command = self.updateplots)

        self.numlabel       = Label(self.tframe, text = '#')
        self.namelabel      = Label(self.tframe, text = 'Chemical Name')
        self.timelabel      = Label(self.tframe, text = 'Plot Time(' + self.timeunit + ')')
        self.depthlabel     = Label(self.tframe, text = 'Plot Depth (' + self.lengthunit + ')')
        self.ztlabel        = Label(self.tframe, bg = self.color, text = 'Depth from:')
        self.matrixlabel    = Label(self.tframe, bg = self.color, text = 'Phase:')

        self.blank3         = Label (self.tframe, text = ' ')

        self.botstartcolumn     = Label(self.frame, text = ' ', width = 2)
        self.botdelcolumn       = Label(self.frame, text = ' ', width = 8)
        self.botnumcolumn       = Label(self.frame, text = ' ', width = 6)
        self.botnamecolumn      = Label(self.frame, text = ' ', width = 20)
        self.botvaluecolumn     = Label(self.frame, text = ' ', width = 12)
        self.botmatrixcolumn    = Label(self.frame, text = ' ', width = 20)
        self.botdepthcolumn     = Label(self.frame, text = ' ', width = 21)
        self.botendcolumn       = Label(self.frame, text = ' ', width = 2)

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
        self.valuecolumn.grid(  row = 1, column = 4)
        self.matrixcolumn.grid( row = 1, column = 5)
        self.endcolumn.grid(    row = 1, column = 7)

        self.typelabel.grid(    row = 2, column = 1, sticky = 'E',  pady = 1, columnspan = 2)
        self.typewidget.grid(   row = 2, column = 3, sticky = 'WE', pady = 1)
        self.varilabel.grid(    row = 2, column = 4, sticky = 'E',  pady = 1)

        self.blank3.grid(       row = 3, column = 3)
        self.numlabel.grid(     row = 4, column = 2)
        self.namelabel.grid(    row = 4, column = 3)
        self.matrixlabel.grid(  row = 4, column = 5)

        self.botstartcolumn.grid(  row = 0, column = 0)
        self.botdelcolumn.grid(    row = 0, column = 1)
        self.botnumcolumn.grid(    row = 0, column = 2)
        self.botnamecolumn.grid(   row = 0, column = 3)
        self.botvaluecolumn.grid(  row = 0, column = 4)
        self.botmatrixcolumn.grid( row = 0, column = 5)
        self.botendcolumn.grid(    row = 0, column = 7)

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

        try:
            self.timevariwidget.grid_forget()
        except: pass

        if self.type.get() == 'Spatial profile':
            if self.variable.get() == 'Water concentration' or self.variable.get() == 'Mass': self.variable.set('Concentration')
            self.variwidget.grid(   row = 2, column = 5, sticky = 'WE', pady = 1)
        else:
            if self.system.topBCtype == 'Finite mixed water column':
                self.extravariwidget.grid(row = 2, column = 5, sticky = 'WE', pady = 1)
            else:
                self.timevariwidget.grid( row = 2, column = 5, sticky = 'WE', pady = 1)

        row = 2

        if self.type.get() == 'Spatial profile':
            self.timelabel.grid(  row = 4, column = 4)
            for plotdata in self.spatialplotdatas:
                plotdata.number = self.spatialplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.chemicals, self.type.get(), self.variable.get(), component_list = self.components, solid_list = self.solids)
                row = row + 1

        elif self.system.dep == 'Deposition' and self.variable.get() != 'Mass' and self.variable.get() != 'Water concentration':
            self.depthcolumn.grid(      row = 1, column = 6)
            self.depthlabel.grid(       row = 4, column = 4)
            self.ztlabel.grid(          row = 4, column = 6)
            self.botdepthcolumn.grid(   row = 0, column = 6)
            for plotdata in self.timeplotdatas:
                plotdata.number = self.timeplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.chemicals, self.type.get(), self.variable.get(), component_list = self.components, solid_list = self.solids, depositionflag = 1, layer_list = self.layers)
                row = row + 1
        else:
            self.depthlabel.grid(       row = 4, column = 4)
            for plotdata in self.timeplotdatas:
                plotdata.number = self.timeplotdatas.index(plotdata) + 1
                plotdata.propertieswidgets(self.frame, row, self.master, self.chemicals, self.type.get(), self.variable.get(), component_list = self.components, solid_list = self.solids, layer_list = self.layers)
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

        if self.type.get() == 'Spatial profile':
            self.spatialplotdatas.append(PlotData(len(self.spatialplotdatas)+ 1, self.timeunit, self.lengthunit))

            self.spatialplotdatas[-1].name      = StringVar(value = self.chemicals[0])
            self.spatialplotdatas[-1].value     = DoubleVar(value = 0)
            self.spatialplotdatas[-1].type      = StringVar(value = 'Initial benthic surface')
            self.spatialplotdatas[-1].component = StringVar(value = self.components[0])
            self.spatialplotdatas[-1].layer     = StringVar(value = 'All layers')
            self.spatialplotdatas[-1].solid     = StringVar(value = 'Total solid')
            self.spatialplotdatas[-1].style     = 'Solid'
            self.spatialplotdatas[-1].color     = self.spatialplotdatas[-1].colornames[self.spatialplotdatas[-1].number - 1]
            self.spatialplotdatas[-1].size      = 0.5
            self.spatialplotdatas[-1].legend    = self.spatialplotdatas[-1].name.get() + ' ' +  str(self.spatialplotdatas[-1].value.get()) + ' ' + self.timeunit

        else:
            self.timeplotdatas.append(PlotData(len(self.timeplotdatas)+ 1, self.timeunit, self.lengthunit))

            self.timeplotdatas[-1].name         = StringVar(value = self.chemicals[0])
            self.timeplotdatas[-1].value        = DoubleVar(value = 0)
            self.timeplotdatas[-1].type         = StringVar(value = 'Initial benthic surface')
            self.timeplotdatas[-1].component    = StringVar(value = self.components[0])
            self.timeplotdatas[-1].solid        = StringVar(value = 'Total solid')
            self.timeplotdatas[-1].layer        = StringVar(value = 'All layers')
            self.timeplotdatas[-1].style        = 'Solid'
            self.timeplotdatas[-1].color        = self.timeplotdatas[-1].colornames[self.timeplotdatas[-1].number - 1]
            self.timeplotdatas[-1].size         = 0.5
            if self.system.dep == 'Deposition':
                self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name.get() + ' ' +  str(self.timeplotdatas[-1].value.get()) + ' ' + self.lengthunit + ' from ' + self.timeplotdatas[-1].type.get()
            else:
                self.timeplotdatas[-1].legend = self.timeplotdatas[-1].name.get() + ' ' +  str(self.timeplotdatas[-1].value.get()) + ' ' + self.lengthunit

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
            try:
                for plotdata in self.spatialplotdatas:
                    if plotdata.value.get() < self.outputt[0]:  error_flag = 1
                    if plotdata.value.get() > self.outputt[-1]: error_flag = 1
            except:
                for plotdata in self.spatialplotdatas:
                    if plotdata.value < self.outputt[0]:  error_flag = 1
                    if plotdata.value > self.outputt[-1]: error_flag = 1
            for plotdata in self.spatialplotdatas:
                plotdata.update_legend()
        else:
            try:
                for plotdata in self.timeplotdatas:
                    if plotdata.type == 'Initial benthic surface':
                        if plotdata.value.get() < self.outputz[0]:  error_flag = 1
                        if plotdata.value.get() > self.outputz[-1]: error_flag = 1
                    else:
                        if plotdata.value.get() < self.outputz[0] - self.outputz[0]: error_flag = 1
                        if plotdata.value.get() > self.outputz[-1]- self.outputz[0]: error_flag = 1
            except:
                for plotdata in self.timeplotdatas:
                    if plotdata.type == 'Initial benthic surface':
                        if plotdata.value < self.outputz[0]:  error_flag = 1
                        if plotdata.value > self.outputz[-1]: error_flag = 1
                    else:
                        if plotdata.value < self.outputz[0] - self.outputz[0]: error_flag = 1
                        if plotdata.value > self.outputz[-1]- self.outputz[0]: error_flag = 1
            for plotdata in self.timeplotdatas:
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


class GraphEditor:

    def __init__(self, master, system, type, variable, plotdatas):

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
        self.type           = type
        self.variable       = variable

        self.cancelflag =  0

        self.plotdatas   = [plotdata.copy() for plotdata in plotdatas]

        for plotdata in self.plotdatas:

            plotdata.legend   = StringVar(value = plotdata.legend)
            plotdata.color    = StringVar(value = plotdata.color)
            plotdata.style    = StringVar(value = plotdata.style)
            plotdata.size     = DoubleVar(value = plotdata.size)

    def make_widgets(self):

        self.instructions    = Label(self.tframe, text = ' Please provides the following output information                    ')

        self.startcolumn     = Label(self.tframe, text = ' ', width = 2)
        self.numcolumn       = Label(self.tframe, text = ' ', width = 6)
        self.stylecolumn     = Label(self.tframe, text = ' ', width = 14)
        self.colorcolumn     = Label(self.tframe, text = ' ', width = 14)
        self.widthcolumn     = Label(self.tframe, text = ' ', width = 12)
        self.legendcolumn    = Label(self.tframe, text = ' ', width = 20)
        self.endcolumn       = Label(self.tframe, text = ' ', width = 2)

        self.numlabel        = Label(self.tframe, text = '#')
        self.stylelabel      = Label(self.tframe, text = 'Style')
        self.colorlabel      = Label(self.tframe, text = 'Color')
        self.widthlabel      = Label(self.tframe, text = 'Weight/Size')
        self.legendlabel     = Label(self.tframe, text = 'Legend')

        self.blank3          = Label (self.tframe, text = ' ')

        self.botstartcolumn  = Label(self.frame, text = ' ', width = 2)
        self.botnumcolumn    = Label(self.frame, text = ' ', width = 6)
        self.botstylecolumn  = Label(self.frame, text = ' ', width = 14)
        self.botcolorcolumn  = Label(self.frame, text = ' ', width = 14)
        self.botwidthcolumn  = Label(self.frame, text = ' ', width = 12)
        self.botlegendcolumn = Label(self.frame, text = ' ', width = 20)
        self.botendcolumn    = Label(self.frame, text = ' ', width = 2)

        self.blank1         = Label (self.bframe, text = ' ')
        self.blank2         = Label (self.bframe, text = ' ')

        self.okbutton       = Button(self.bframe, text = 'OK', width = 20, command = self.OK)
        self.cancelbutton   = Button(self.bframe, text = 'Cancel', width = 20, command = self.Cancel)

        self.instructions.grid( row = 0, column = 0, columnspan = 6, padx = 8, sticky = 'W')

        self.startcolumn.grid(  row = 1, column = 0)
        self.numcolumn.grid(    row = 1, column = 2)
        self.stylecolumn.grid(  row = 1, column = 3)
        self.colorcolumn.grid(  row = 1, column = 4)
        self.widthcolumn.grid(  row = 1, column = 5)
        self.legendcolumn.grid( row = 1, column = 6)
        self.endcolumn.grid(    row = 1, column = 7)

        self.blank3.grid(       row = 3, column = 2)
        self.numlabel.grid(     row = 4, column = 2)
        self.stylelabel.grid(   row = 4, column = 3)
        self.colorlabel.grid(   row = 4, column = 4)
        self.widthlabel.grid(   row = 4, column = 5)
        self.legendlabel.grid(  row = 4, column = 6)

        self.botstartcolumn.grid(  row = 0, column = 0)
        self.botnumcolumn.grid(    row = 0, column = 2)
        self.botstylecolumn.grid(  row = 0, column = 3)
        self.botcolorcolumn.grid(  row = 0, column = 4)
        self.botwidthcolumn.grid(  row = 0, column = 5)
        self.botlegendcolumn.grid( row = 0, column = 6)
        self.botendcolumn.grid(    row = 0, column = 7)

        self.updateplots()


    def updateplots(self, event = None):

        row = 2

        for plotdata in self.plotdatas:
            plotdata.graphicwidgets(self.frame, row, self.master)
            row = row + 1

        self.blank1.grid(row = row)
        row = row + 1
        self.blank2.grid(row = row)
        row = row + 1
        self.okbutton.grid(row = row, columnspan = 11)
        row = row + 1
        self.cancelbutton.grid(row = row, columnspan = 11)
        row = row + 1

        self.focusbutton = self.okbutton
        self.okbutton.bind('<Return>',   self.OK)

        self.master.geometry()
        self.master.center()

    def OK(self, event = None):
        """Finish and move on.  Checks that the number chemicals are less than the
        total number of chemicals in database."""

        error_flag = 0

        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()

    def Cancel(self):

        self.cancelflag = 1

        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()


class PlotData:

    def __init__(self, number, timeunit, lengthunit):

        self.number       = number
        self.types        = ['Initial benthic surface', 'New benthic surface']
        self.colornames   = ['Blue', 'Green', 'Red', 'Cyan', 'Magenta', 'Yellow', 'Orange', 'Purple', 'Brown', 'Pink', 'Grey', 'Black']
        self.stylenames   = ['Solid', 'Dashed', 'Dash-dot', 'Dotted', 'Point', 'Circle', 'Triangle', 'Square', 'Star', 'Cross']
        self.timeunit     = timeunit
        self.lengthunit   = lengthunit

    def copy(self):

        plotdata = PlotData(self.number, self.timeunit, self.lengthunit)

        plotdata.name         = self.name
        plotdata.value        = self.value
        plotdata.component    = self.component
        plotdata.solid        = self.solid
        plotdata.type         = self.type
        plotdata.layer        = self.layer
        plotdata.number       = self.number
        plotdata.color        = self.color
        plotdata.style        = self.style
        plotdata.size         = self.size
        plotdata.legend       = self.legend

        return plotdata

    def propertieswidgets(self, frame, row, master, chemical_list, type, variable, component_list = None, solid_list = None, depositionflag = None, layer_list = None):

        self.master = master
        self.frame  = frame

        self.chemical_list  = chemical_list
        self.row            = row
        self.variable       = variable
        self.plot_type      = type
        self.depositionflag = depositionflag

        self.delwidget      = Button(frame, width = 5,  justify = 'center', text = 'Delete', command = self.del_plotdata)
        self.valuewidget    = Entry (frame, width = 8,  justify = 'center', textvariable = self.value)
        self.valuelabel     = Label (frame, width = 8,  justify = 'center', text = 'N/A')
        self.numwidget      = Label (frame, width = 4,  justify = 'center', text = self.number)

        self.delwidget.grid(        row  = row, column = 1, padx = 2 ,pady = 1)
        self.numwidget.grid(         row  = row, column = 2, padx = 2 ,pady = 1)

        if variable == 'Water concentration' or variable == 'Mass':
            self.valuelabel.grid(   row  = row, column = 4, padx = 2, pady = 1)
        else:
            self.valuewidget.grid(  row  = row, column = 4, padx = 2, pady = 1)

        if variable == 'Material fraction':
            self.namewidget      = OptionMenu(frame, self.component, *component_list[1:], command = self.update_legend)
        else:
            self.namewidget      = OptionMenu(frame, self.name, *chemical_list, command = self.update_legend)

        self.namewidget.grid(     row  = row, column = 3, padx = 2 ,pady = 1, sticky = 'WE')

        if variable == 'Solid concentration':
            self.componentwidget = OptionMenu(frame, self.solid, *solid_list, command = self.update_legend)
        elif variable == 'Mass':
            self.componentwidget = OptionMenu(frame, self.layer, *layer_list, command = self.update_legend)
        elif variable == 'Concentration':
            self.componentwidget = Label (frame, width = 16,  justify = 'center', text = 'Porewater')
        elif variable == 'Total concentration':
            self.componentwidget = Label (frame, width = 16,  justify = 'center', text = 'Porewater/DOC')
        elif variable == 'Flux':
            self.componentwidget = Label (frame, width = 16,  justify = 'center', text = 'Porewater/DOC/Solid')
        elif variable == 'Water concentration':
            self.componentwidget = Label (frame, width = 16,  justify = 'center', text = 'Overlying water column')
        elif variable == 'Material fraction':
            self.componentwidget = Label (frame, width = 16,  justify = 'center', text = 'N/A')

        self.componentwidget.grid(row  = row, column = 5, padx = 2 ,pady = 1, sticky = 'WE')

        if depositionflag == 1:
            self.typewidget     = OptionMenu(frame, self.type, *self.types, command = self.update_legend)
            self.typewidget.grid(     row  = row, column = 6, padx = 2, pady = 1, sticky = 'WE')

    def update_legend(self, default = None):

        if self.variable == 'Material fraction':
            if self.plot_type == 'Spatial profile':
                self.legend = self.name.get() + ' ' + self.component.get() + ' at ' +  str(self.value.get()) + ' ' + self.timeunit
            else:
                self.legend = self.name.get() + ' ' + self.component.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit
        else:
            if self.plot_type == 'Spatial profile':
                if self.variable == 'Solid concentration':
                    self.legend = self.name.get() + ' in ' +  self.component.get() + ' at ' + str(self.value.get()) + ' ' + self.timeunit
                else:
                    self.legend = self.name.get() + ' at ' +  str(self.value.get()) + ' ' + self.timeunit
            else:
                if self.depositionflag == 1:
                    if self.variable == 'Solid concentration':
                        self.legend = self.name.get() + ' in '+ self.component.get()  + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit + ' from ' + self.type.get()
                    else:
                        self.legend = self.name.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit + ' from ' + self.type.get()
                else:
                    if self.variable == 'Solid concentration':
                        self.legend = self.name.get() + ' in '+ self.component.get()  + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit
                    else:
                        self.legend = self.name.get() + ' at ' +  str(self.value.get()) + ' ' + self.lengthunit


    def remove_propertieswidgets(self):

        self.delwidget.grid_forget()
        self.namewidget.grid_forget()
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

        self.name          = self.name.get()
        self.value         = self.value.get()
        self.type          = self.type.get()
        self.component     = self.component.get()
        self.solid         = self.solid.get()
        self.layer         = self.layer.get()

    def del_plotdata(self):

        self.master.window.delplotdata(self.number)

    def get_graphicdata(self):

        self.style     = self.style.get()
        self.color     = self.color.get()
        self.size      = self.size.get()
        self.legend    = self.legend.get()


class FigureEditor:

    def __init__(self, master, system, type, variable, sketch, axislimit, figuresize, legendposition, fontsize, figuretexts, MCplots = None):

        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.tframe    = Frame(master.tframe)
        self.frame     = Frame(master.frame)
        self.bframe    = Frame(master.bframe)
        self.top       = None
        self.flag      = 0
        self.system    = system
        self.type      = type
        self.variable  = variable

        rgb            = self.master.frame.winfo_rgb(self.master.frame.cget('bg'))
        self.color     = '#%02x%02x%02x' %rgb
        self.mplrgb    = (rgb[0] / 65536., rgb[1] / 65536., rgb[2] / 65536.)

        self.timeunit   = system.timeunit
        self.lengthunit = system.lengthunit
        self.concunit   = system.concunit

        self.xaxismin = DoubleVar(value = '%.4g' % axislimit[0])
        self.xaxismax = DoubleVar(value = '%.4g' % axislimit[1])
        self.yaxismin = DoubleVar(value = '%.4g' % axislimit[2])
        self.yaxismax = DoubleVar(value = '%.4g' % axislimit[3])

        self.figurewidth  = IntVar(value = figuresize[0])
        self.figureheight = IntVar(value = figuresize[1])

        self.titlefontsize  = IntVar(value = fontsize['Title'])
        self.xaxisfontsize  = IntVar(value = fontsize['Axis'][0])
        self.yaxisfontsize  = IntVar(value = fontsize['Axis'][1])
        self.xlabelfontsize = IntVar(value = fontsize['Label'][0])
        self.ylabelfontsize = IntVar(value = fontsize['Label'][1])
        self.legendfontsize = IntVar(value = fontsize['Legend'])

        self.legendpositions = ['Default', 'Upper right', 'Upper left', 'Lower left', 'Lower right', 'Right', 'Center left', 'Center right', 'Lower center', 'Upper center', 'Center', 'None']
        self.legendposition = StringVar(value = self.legendpositions[legendposition])

        self.sketches       = sketch['Sketches']
        self.plotsizes      = ['Default', 'User-defined']
        self.fontsizes      = ['Default', 'User-defined']
        self.textflags      = ['Default', 'User-defined']
        self.MCflags        = ['Default', 'User-defined']

        self.sketch         = StringVar(value = sketch['Sketch'])
        self.sketchlabel    = IntVar(value = sketch['Label'])
        self.sketchline     = IntVar(value = sketch['Time'])
        self.locationline   = IntVar(value = sketch['Depth'])
        self.fontsize       = StringVar(value = fontsize['Flag'])

        self.textflag       = StringVar(value = figuretexts['Flag'])
        self.title          = StringVar(value = figuretexts['title'])
        self.xlabel         = StringVar(value = figuretexts['xlabel'])
        self.ylabel         = StringVar(value = figuretexts['ylabel'])

        if MCplots <> None:
            self.MCflag       = StringVar(value = MCplots['Flag'])
            self.MCavg        = DoubleVar(value = MCplots['AVG'])
            self.MCsd         = DoubleVar(value = MCplots['SD'])
            self.MCmaxmin     = DoubleVar(value = MCplots['Max/Min'])

        self.MCplots = MCplots
        self.cancelflag =  0

    def make_widgets(self):

        self.instructions   = Label(self.frame, text = ' Please provides the following output information                    ')

        self.startcolumn    = Label(self.frame, text = ' ', width = 2)
        self.delcolumn      = Label(self.frame, text = ' ', width = 12)
        self.name1column    = Label(self.frame, text = ' ', width = 6)
        self.name2column    = Label(self.frame, text = ' ', width = 2)
        self.name3column    = Label(self.frame, text = ' ', width = 6)
        self.unitcolumn     = Label(self.frame, text = ' ', width = 12)
        self.endcolumn      = Label(self.frame, text = ' ', width = 2)

        self.sizelabel      = Label(self.frame, bg = self.color, text = 'Plot size:')
        self.widthwidget    = Entry(self.frame, textvariable = self.figurewidth ,  width = 8, justify = 'center')
        self.xwidget        = Label(self.frame, bg = self.color, text = 'x', width = 2, justify = 'center')
        self.heightwidget   = Entry(self.frame, textvariable = self.figureheight , width = 8, justify = 'center')
        self.sizeunit       = Label(self.frame, bg = self.color, text = 'Pixels')

        self.xaxiswidget1    = Entry(self.frame, textvariable = self.xaxismin,  width = 8, justify = 'center')
        self.xaxiswidget2    = Label(self.frame, text = 'to'         ,  width = 2, justify = 'center')
        self.xaxiswidget3    = Entry(self.frame, textvariable = self.xaxismax,  width = 8, justify = 'center')

        self.yaxiswidget1    = Entry(self.frame, textvariable = self.yaxismin,  width = 8, justify = 'center')
        self.yaxiswidget2    = Label(self.frame, text = 'to'         ,  width = 2, justify = 'center')
        self.yaxiswidget3    = Entry(self.frame, textvariable = self.yaxismax,  width = 8, justify = 'center')

        if self.type == 'Spatial profile':
            self.yaxislabel = Label(self.frame, bg = self.color, text = 'Depth:')
            self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.lengthunit)
            if self.variable == 'Concentration':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Concentration:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Flux':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Flux:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit[:-1] + self.system.lengthunit + '^2' + '/' + self.system.timeunit)
            if self.variable == 'Solid concentration':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Solid concentration:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit[:-1] + 'kg')
            if self.variable == 'Total concentration':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Total concentration:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Water concentration':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Water column concentration:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Material fraction':
                self.xaxislabel = Label(self.frame, bg = self.color, text = 'Material fraction:')
                self.xaxisunit  = Label(self.frame, bg = self.color, text = '')
        else:
            self.xaxislabel = Label(self.frame, bg = self.color, text = 'Time:')
            self.xaxisunit  = Label(self.frame, bg = self.color, text = self.system.timeunit)
            if self.variable == 'Concentration':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Concentration:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Flux':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Flux:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit[:-1]+ '/'+ self.system.lengthunit + '^2' + '/' + self.system.timeunit)
            if self.variable == 'Solid concentration':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Solid concentration:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit[:-1] + 'kg')
            if self.variable == 'Total concentration':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Total concentration:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Water concentration':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Water column concentration:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit)
            if self.variable == 'Material fraction':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Material fraction:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = '')
            if self.variable == 'Mass':
                self.yaxislabel = Label(self.frame, bg = self.color, text = 'Mass:')
                self.yaxisunit  = Label(self.frame, bg = self.color, text = self.system.concunit[:-1])

        self.legendlabel    = Label(self.frame, bg = self.color, text = 'Legend:')
        self.legendwidget   = OptionMenu(self.frame, self.legendposition,  *self.legendpositions)

        self.sketchslabel   = Label(self.frame, bg = self.color, text = 'Column sketch:')
        self.sketchwidget   = OptionMenu(self.frame, self.sketch, *self.sketches, command = self.updatewidgets)

        self.typelabel      = Label(self.frame, bg = self.color, text = 'Materials:')
        self.typewidget     = Checkbutton(self.frame, variable = self.sketchlabel)

        if self.type == 'Spatial profile':
            self.linelabel  = Label(self.frame, bg = self.color, text = 'Time:')
            self.linewidget = Checkbutton(self.frame, variable = self.sketchline)
        else:
            self.linelabel  = Label(self.frame, bg = self.color, text = 'Depth:')
            self.linewidget = Checkbutton(self.frame, variable = self.locationline)

        self.plottextlabel  = Label(self.frame, bg = self.color, text = 'Graph texts:')
        self.plottextwidget = OptionMenu(self.frame, self.textflag,  *self.textflags, command = self.updatewidgets)

        self.titlelabel     = Label(self.frame, bg = self.color, text = 'Title:')
        self.titleentry     = Entry(self.frame, textvariable = self.title ,  width = 28, justify = 'center')

        self.xlabellabel     = Label(self.frame, bg = self.color, text = 'X-label:')
        self.xlabelentry     = Entry(self.frame, textvariable = self.xlabel ,  width = 28, justify = 'center')

        self.ylabellabel     = Label(self.frame, bg = self.color, text = 'Y-label:')
        self.ylabelentry     = Entry(self.frame, textvariable = self.ylabel ,  width = 28, justify = 'center')

        self.fontsizelabel  = Label(self.frame, bg = self.color, text = 'Fontsize:')
        self.fontsizewidget = OptionMenu(self.frame, self.fontsize,  *self.fontsizes, command = self.updatewidgets)

        self.titlefontlabel  = Label(self.frame, bg = self.color, text = 'Title:')
        self.titlefontentry  = Entry(self.frame, textvariable = self.titlefontsize ,  width = 8, justify = 'center')
        self.titlefontunit   = Label(self.frame, bg = self.color, text = 'pt')

        self.xaxisfontlabel  = Label(self.frame, bg = self.color, text = 'X-axis:')
        self.xaxisfontentry  = Entry(self.frame, textvariable = self.xaxisfontsize ,   width = 8, justify = 'center')
        self.xaxisfontunit   = Label(self.frame, bg = self.color, text = 'pt')
        self.yaxisfontlabel  = Label(self.frame, bg = self.color, text = 'Y-axis:')
        self.yaxisfontentry  = Entry(self.frame, textvariable = self.yaxisfontsize ,   width = 8, justify = 'center')
        self.yaxisfontunit   = Label(self.frame, bg = self.color, text = 'pt')

        self.xlabelfontlabel = Label(self.frame, bg = self.color, text = 'X-Label:')
        self.xlabelfontentry = Entry(self.frame, textvariable = self.xlabelfontsize ,  width = 8, justify = 'center')
        self.xlabelfontunit  = Label(self.frame, bg = self.color, text = 'pt')
        self.ylabelfontlabel = Label(self.frame, bg = self.color, text = 'Y-Label:')
        self.ylabelfontentry = Entry(self.frame, textvariable = self.ylabelfontsize ,  width = 8, justify = 'center')
        self.ylabelfontunit  = Label(self.frame, bg = self.color, text = 'pt')

        self.legendfontlabel = Label(self.frame, bg = self.color, text = 'Legend:')
        self.legendfontentry = Entry(self.frame, textvariable = self.legendfontsize , width = 8, justify = 'center')
        self.legendfontunit  = Label(self.frame, bg = self.color, text = 'pt')

        if self.MCplots <> None:

            self.MCflaglabel     = Label(self.frame, bg = self.color, text = 'Batch file plot:')
            self.MCflagwidget    = OptionMenu(self.frame, self.MCflag, *self.MCflags, command = self.updatewidgets)

            self.MCavglabel      = Label(self.frame, bg = self.color, text = 'Average:')
            self.MCavgwidget     = Entry(self.frame, textvariable = self.MCavg,   width = 8, justify = 'center')

            self.MCsdlabel       = Label(self.frame, bg = self.color, text = 'StDev:')
            self.MCsdwidget      = Entry(self.frame, textvariable = self.MCsd,   width = 8, justify = 'center')

            self.MCmaxminlabel   = Label(self.frame, bg = self.color, text = 'Max/Min:')
            self.MCmaxminwidget  = Entry(self.frame, textvariable = self.MCmaxmin,   width = 8, justify = 'center')

        self.blank1         = Label (self.bframe, text = ' ')
        self.blank2         = Label (self.bframe, text = ' ')

        self.okbutton       = Button(self.bframe, text = 'OK', width = 20, command = self.OK)
        self.cancelbutton   = Button(self.bframe, text = 'Cancel', width = 20, command = self.Cancel)

        self.instructions.grid( row = 0, column = 0, columnspan = 6, padx = 8, sticky = 'W')

        self.startcolumn.grid(  row = 1, column = 0)
        self.delcolumn.grid(    row = 1, column = 1)
        self.name1column.grid(  row = 1, column = 2)
        self.name2column.grid(  row = 1, column = 3)
        self.name3column.grid(  row = 1, column = 4)
        self.endcolumn.grid(    row = 1, column = 5)

        self.sizelabel.grid(    row = 2, column = 1, sticky = 'E', pady = 4)
        self.widthwidget.grid(  row = 2, column = 2, sticky = 'E', pady = 2)
        self.xwidget.grid(      row = 2, column = 3)
        self.heightwidget.grid( row = 2, column = 4, sticky = 'W', pady = 2)
        self.sizeunit.grid(     row = 2, column = 5, sticky = 'W', pady = 2)

        self.xaxislabel.grid(   row = 3, column = 1, sticky = 'E', pady = 4)
        self.xaxiswidget1.grid( row = 3, column = 2, sticky = 'E', pady = 2)
        self.xaxiswidget2.grid( row = 3, column = 3, pady = 2)
        self.xaxiswidget3.grid( row = 3, column = 4, sticky = 'W', pady = 2)
        self.xaxisunit.grid(    row = 3, column = 5, sticky = 'W', pady = 2)

        self.yaxislabel.grid(   row = 4, column = 1, sticky = 'E', pady = 4)
        self.yaxiswidget1.grid( row = 4, column = 2, sticky = 'E', pady = 2)
        self.yaxiswidget2.grid( row = 4, column = 3)
        self.yaxiswidget3.grid( row = 4, column = 4, sticky = 'W', pady = 2)
        self.yaxisunit.grid(    row = 4, column = 5, sticky = 'W', pady = 2)

        self.legendlabel.grid(  row = 5, column = 1, sticky = 'E', pady = 4)
        self.legendwidget.grid( row = 5, column = 2, sticky = 'WE', columnspan = 3, pady = 1)

        self.sketchslabel.grid( row = 6, column = 1, sticky = 'E', pady = 4)
        self.sketchwidget.grid( row = 6, column = 2, sticky = 'WE', columnspan = 3, pady = 1)

        row = 0
        self.blank1.grid(row = row)
        row = row + 1
        self.blank2.grid(row = row)
        row = row + 1
        self.okbutton.grid(row = row, columnspan = 11)
        row = row + 1
        self.cancelbutton.grid(row = row, columnspan = 11)
        row = row + 1

        self.okbutton.bind('<Return>',   self.OK)
        self.updatewidgets()

    def updatewidgets(self, default = None):

        try:
            self.typelabel.grid_forget()
            self.typewidget.grid_forget()
            self.linelabel.grid_forget()
            self.linewidget.grid_forget()
        except: pass

        try:
            self.plottextlabel.grid_forget()
            self.plottextwidget.grid_forget()
        except: pass

        try:
            self.titlelabel.grid_forget()
            self.titleentry.grid_forget()
            self.xlabellabel.grid_forget()
            self.xlabelentry.grid_forget()
            self.ylabellabel.grid_forget()
            self.ylabelentry.grid_forget()
        except: pass

        try:
            self.fontsizelabel.grid_forget()
            self.fontsizewidget.grid_forget()
        except: pass

        try:
            self.xlabelfontlabel.grid_forget()
            self.xlabelfontentry.grid_forget()
            self.xlabelfontunit.grid_forget()
            self.ylabelfontlabel.grid_forget()
            self.ylabelfontentry.grid_forget()
            self.ylabelfontunit.grid_forget()
            self.xaxisfontlabel.grid_forget()
            self.xaxisfontentry.grid_forget()
            self.xaxisfontunit.grid_forget()
            self.yaxisfontlabel.grid_forget()
            self.yaxisfontentry.grid_forget()
            self.yaxisfontunit.grid_forget()
            self.titlefontlabel.grid_forget()
            self.titlefontentry.grid_forget()
            self.titlefontunit.grid_forget()
            self.legendfontlabel.grid_forget()
            self.legendfontentry.grid_forget()
            self.legendfontunit.grid_forget()
        except: pass

        try:
            self.xlabellabel.grid_forget()
            self.xlabelentry.grid_forget()
            self.ylabellabel.grid_forget()
            self.ylabelentry.grid_forget()
            self.titlelabel.grid_forget()
            self.titleentry.grid_forget()
        except: pass

        try:
            self.MCavglabel.grid_forget()
            self.MCavgwidget.grid_forget()
            self.MCsdlabel.grid_forget()
            self.MCsdwidget.grid_forget()
            self.MCmaxminlabel.grid_forget()
            self.MCmaxminwidget.grid_forget()
        except: pass

        row = 7

        if self.sketch.get() != 'Hide sketch':

            self.typelabel.grid(    row = row, column = 2, columnspan = 2, sticky = 'WE', pady = 1)
            self.typewidget.grid(   row = row, column = 4, columnspan = 1, sticky = 'WE', pady = 1)

            row = row + 1

            self.linelabel.grid(    row = row, column = 2, columnspan = 2, sticky = 'WE', pady = 1)
            self.linewidget.grid(   row = row, column = 4, columnspan = 1, sticky = 'WE', pady = 1)

            row = row + 1

        self.plottextlabel.grid(    row = row,  column = 1, sticky = 'E', pady = 1)
        self.plottextwidget.grid(   row = row,  column = 2, sticky = 'WE', columnspan = 3, pady = 1)

        row = row + 1

        if self.textflag.get() == self.textflags[1]:

            self.titlelabel.grid( row = row, column = 1, columnspan = 1,sticky = 'E', pady = 2)
            self.titleentry.grid( row = row, column = 2, sticky = 'W', columnspan = 4, pady = 2)

            row = row + 1

            self.xlabellabel.grid( row = row, column = 1, columnspan = 1,sticky = 'E', pady = 2)
            self.xlabelentry.grid( row = row, column = 2, sticky = 'W', columnspan = 4, pady = 2)

            row = row + 1

            self.ylabellabel.grid( row = row, column = 1, columnspan = 1,sticky = 'E', pady = 2)
            self.ylabelentry.grid( row = row, column = 2, sticky = 'W', columnspan = 4, pady = 2)

            row = row + 1

        else:
            if self.type == 'Spatial profile':
                if self.variable == 'Material fraction':
                    self.title.set('Solid material fraction Profiles')
                    self.xlabel.set('Volumetric fraction')
                elif self.variable == 'Solid concentration':
                    self.title.set('Solid Concentration Profiles')
                    self.xlabel.set('Solid concentration ('+self.concunit[:-1]+'kg)')
                else:
                    if self.variable == 'Concentration':
                        self.title.set('Porewater Concentration Profiles')
                        self.xlabel.set('Porewater concentration ('+self.concunit+')')
                    if self.variable == 'Flux':
                        self.title.set('Flux Profiles')
                        self.xlabel.set('Flux (' + self.concunit[:-1] + self.lengthunit + u'\u00B2' + '/' + self.timeunit+')')
                    if self.variable == 'Total concentration':
                        self.title.set('Porespace Concentration Profiles')
                        self.xlabel.set('Porespace Concentration ('+self.concunit+')')
                self.ylabel.set('Depth ('+self.lengthunit+')')
            else:
                if self.variable == 'Material fraction':
                        self.title.set('Solid material fraction Time Profiles')
                        self.ylabel.set('Volumetric fraction ('+self.concunit+')')

                elif self.variable == 'Solid concentration':
                        self.title.set('Solid Concentration Time Profiles')
                        self.ylabel.set('Solid concentration ('+self.concunit[:-1]+'kg)')

                elif self.variable == 'Water concentration':
                        self.title.set('Overlying Water Concentration Time Profiles')
                        self.ylabel.set('Overlying water concentration ('+self.concunit+')')

                elif self.variable == 'Mass':
                        self.title.set('Accumulated Mass Time Profiles')
                        self.ylabel.set('Mass(' + self.concunit[:-1] + 'cm'+u'\u00B2'+'/yr)')
                else:
                    if self.variable == 'Concentration':
                        self.title.set('Porewater Concentration Time Profiles')
                        self.ylabel.set('Porewater concentration ('+self.concunit+')')
                    if self.variable == 'Flux':
                        self.title.set('Flux Time Profiles')
                        self.ylabel.set('Flux (' + self.concunit[:-1]+ self.lengthunit + u'\u00B2' + '/' + self.timeunit+')')
                    if self.variable == 'Total concentration':
                        self.title.set('Porespace Concentration Time Profiles')
                        self.ylabel.set('Porespace Concentration ('+self.concunit+')')
                self.xlabel.set('Time ('+self.timeunit+')')

        self.fontsizelabel.grid(    row = row, column = 1, sticky = 'E', pady = 1)
        self.fontsizewidget.grid(   row = row, column = 2, sticky = 'WE', columnspan = 3, pady = 1)

        row = row + 1

        if self.fontsize.get() == self.fontsizes[1]:

            self.titlefontlabel.grid( row = row, column = 2, columnspan = 1,sticky = 'E', pady = 1)
            self.titlefontentry.grid( row = row, column = 4, sticky = 'W', pady = 1)
            self.titlefontunit.grid(  row = row, column = 5, sticky = 'W', pady = 1)
            row = row + 1

            self.legendfontlabel.grid( row = row, column = 2, columnspan = 1,sticky = 'E', pady = 1)
            self.legendfontentry.grid( row = row, column = 4, sticky = 'W', pady = 1)
            self.legendfontunit.grid(  row = row, column = 5, sticky = 'W', pady = 1)
            row = row + 1

            self.xlabelfontlabel.grid( row = row, column = 2, columnspan = 1, sticky = 'E', pady = 1)
            self.xlabelfontentry.grid( row = row, column = 4, sticky = 'W', pady = 1)
            self.xlabelfontunit.grid(  row = row, column = 5, sticky = 'W', pady = 1)

            row = row + 1

            self.ylabelfontlabel.grid( row = row, column = 2, columnspan = 1,sticky = 'E', pady = 1)
            self.ylabelfontentry.grid( row = row, column = 4, sticky = 'W', pady = 1)
            self.ylabelfontunit.grid(  row = row, column = 5, sticky = 'W', pady = 1)

            row = row + 1

            self.xaxisfontlabel.grid(  row = row, column = 2, columnspan = 1, sticky = 'E', pady = 1)
            self.xaxisfontentry.grid(  row = row, column = 4, sticky = 'W', pady = 1)
            self.xaxisfontunit.grid(   row = row, column = 5, sticky = 'W', pady = 1)

            row = row + 1
            self.yaxisfontlabel.grid(  row = row, column = 2, columnspan = 1,sticky = 'E', pady = 1)
            self.yaxisfontentry.grid(  row = row, column = 4, sticky = 'W', pady = 1)
            self.yaxisfontunit.grid(   row = row, column = 5, sticky = 'W', pady = 1)

            row = row + 1

        else:
            self.titlefontsize.set(12)
            self.xaxisfontsize.set(8)
            self.yaxisfontsize.set(8)
            self.xlabelfontsize.set(10)
            self.ylabelfontsize.set(10)
            self.legendfontsize.set(9)

        if self.MCplots <> None:

            self.MCflaglabel.grid(    row = row,  column = 1, sticky = 'E', pady = 1)
            self.MCflagwidget.grid(   row = row,  column = 2, sticky = 'WE', columnspan = 3, pady = 1)

            row = row + 1

            if self.MCflag.get() == self.MCflags[1]:

                self.MCavglabel.grid(    row = row, column = 2, columnspan = 2, sticky = 'WE', pady = 1)
                self.MCavgwidget.grid(   row = row, column = 4, columnspan = 1, sticky = 'WE', pady = 1)

                row = row + 1

                self.MCsdlabel.grid(    row = row, column = 2, columnspan = 2, sticky = 'WE', pady = 1)
                self.MCsdwidget.grid(   row = row, column = 4, columnspan = 1, sticky = 'WE', pady = 1)

                row = row + 1

                self.MCmaxminlabel.grid(    row = row, column = 2, columnspan = 2, sticky = 'WE', pady = 1)
                self.MCmaxminwidget.grid(   row = row, column = 4, columnspan = 1, sticky = 'WE', pady = 1)

                row = row + 1

            else:
                self.MCavg.set(1.0)
                self.MCsd.set(0.2)
                self.MCmaxmin.set(0.1)

        row = row + 1

        self.focusbutton = self.okbutton

        self.master.geometry()
        self.master.center()

    def OK(self, event = None):
        """Finish and move on.  Checks that the number chemicals are less than the
        total number of chemicals in database."""

        error_flag = 0

        if self.xaxismin.get() > self.xaxismax.get():  error_flag = 1
        if self.yaxismin.get() > self.yaxismax.get():  error_flag = 1

        if self.master.window.top is not None: self.master.open_toplevel()
        elif error_flag == 1: self.warning()
        else: self.master.tk.quit()

    def warning(self):

        tkmb.showerror(title = self.version, message = 'The minimum axis must be smaller than the maximum axis value')
        self.focusbutton = None
        self.master.tk.lift()

    def Cancel(self):

        self.cancelflag = 1

        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()


class ResultExporter:

    def __init__(self, master, system, output):
    
        self.master    = master
        self.version   = system.version
        self.fonttype  = system.fonttype
        self.tframe    = Frame(master.tframe)
        self.frame     = Frame(master.frame)
        self.bframe    = Frame(master.bframe)
        self.top       = None
        self.flag      = 0
        
        self.system    = system
        self.output    = output
        self.filename  = self.master.master.window.graphpath

        self.cancelflag =  0

        self.types      = ['CSV file', 'Doc Report']
        self.chemicals  = [chemical.name for chemical in system.chemicals]
        self.MWs        = [chemical.MW   for chemical in system.chemicals]
        self.components = [component.name for component in system.components]

        if system.biomix == 1:  self.list = self.chemicals + self.components
        else:                   self.list = self.chemicals

        self.type      = StringVar(value = self.types[0])
        self.name      = StringVar(value = 'output')
        self.filepath  = StringVar(value = self.master.master.window.graphpath)
        self.chemical  = StringVar(value = self.chemicals[0])
        self.MW        = self.MWs[0]

    def make_widgets(self):

        self.instructions   = Label(self.frame, text = ' Please provides the following output information                    ')

        self.startcolumn    = Label(self.frame, text = ' ', width = 2)
        self.paracolumn     = Label(self.frame, text = ' ', width = 10)
        self.valuecolumn    = Label(self.frame, text = ' ', width = 20)
        self.browsecolumn   = Label(self.frame, text = ' ', width = 10)
        self.endcolumn      = Label(self.frame, text = ' ', width = 2)

        self.typelabel      = Label(self.frame, text = 'File type:')
        self.namelabel      = Label(self.frame, text = 'File name:')
        self.dirlabel       = Label(self.frame, text = 'File directory:')

        self.typewidget     = Label(self.frame, textvariable = self.type, width = 20, justify = 'center')
        self.namewidget     = Entry(self.frame, textvariable = self.name, width = 20, justify = 'center')
        self.csvext         = Label(self.frame, text = '.csv')
        self.docext         = Label(self.frame, text = '.doc')
        
        self.dirwidget      = Entry(self.frame, textvariable = self.filepath, width = 20)
        self.dirbutton      = Button(self.frame, text = 'Browse', command = self.browse_dir, width = 8)

        self.okbutton      = Button(self.frame, text = 'OK', width = 15, command = self.OK)
        self.cancelbutton  = Button(self.frame, text = 'Cancel', width = 15, command = self.Cancel)
        self.blank1        = Label(self.frame, text = ' ')
        self.blank2        = Label(self.frame, text = ' ')
        
        self.instructions.grid( row = 0, column = 0, columnspan = 6, padx = 8, sticky = 'W')
        
        self.startcolumn.grid(  row = 1, column = 0)
        self.paracolumn.grid(   row = 1, column = 1)
        self.valuecolumn.grid(  row = 1, column = 2)
        self.browsecolumn.grid( row = 1, column = 3)
        self.endcolumn.grid(    row = 1, column = 4)

        self.typelabel.grid(    row = 2, column = 1, sticky = 'E',  padx = 4, pady = 1)
        self.typewidget.grid(   row = 2, column = 2, sticky = 'WE', padx = 4, pady = 1)

        self.namelabel.grid(    row = 3, column = 1, sticky = 'E',  padx = 4, pady = 1)
        self.namewidget.grid(   row = 3, column = 2,                padx = 4, pady = 1)
        
        self.dirlabel.grid(     row = 4, column = 1, sticky = 'E',  padx = 4, pady = 1)
        self.dirwidget.grid(    row = 4, column = 2,                padx = 4, pady = 1)
        self.dirbutton.grid(    row = 4, column = 3, sticky = 'WE', padx = 4, pady = 1)

        self.click_type()
        
    def click_type(self, event = None):

        try:
            self.csvext.grid_forget()
            self.docext.grid_forget()
        except: pass

        row = 5
        
        if self.type.get() == self.types[0]:
            self.csvext.grid(    row = 3, column = 3, sticky = 'W', padx = 4, pady = 1)
            row = row + 1

        self.blank1.grid(row = row)
        row = row + 1
        self.okbutton.grid(row = row, columnspan = 11)
        row = row + 1
        self.cancelbutton.grid(row = row, columnspan = 11)
        row = row + 1
        self.blank2.grid(row = row)
        row = row + 1
        
        self.okbutton.bind('<Return>', self.OK)
        self.focusbutton = self.okbutton
        self.master.geometry()

    def browse_dir(self, event = None):

        if self.type.get() == self.types[0]:    filetype = [('CSV files',  '*.csv')]
        else:                                   filetype = [('Word files', '*.doc')]
        
        self.filename = tkfd.asksaveasfilename(initialdir=self.filepath.get(), initialfile = self.name.get(), title='Please define the directory to save the file', filetypes = filetype)

        if self.filename.count('.csv') > 0: self.filename = self.filename[0:-4]
        if self.filename.count('.doc') > 0: self.filename = self.filename[0:-4]

        if len(self.filename) > 0:
            i = 1
            while self.filename[-i] != '/':
                i = i + 1

            self.filepath.set(self.filename[0:-i+1])
            self.name.set(self.filename[-i+1:])

        self.frame.update()
        self.focusbutton = self.okbutton        
        self.master.geometry()
        
    def OK(self, event = None):
        """Finish and move on.  Checks that the number chemicals are less than the
        total number of chemicals in database."""

        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()
        
    def Cancel(self):
        
        self.cancelflag = 1
            
        if self.master.window.top is not None: self.master.open_toplevel()
        else: self.master.tk.quit()

def postprocess_data(system, output):
    """Shows the results of the simulation."""

    root = CapSimWindow(buttons = 2)
    root.make_window(Postprocess(root, system, output))
    root.tk.config(background = root.window.color)
    root.buttonframe.config(background = root.window.color)
    root.blank.config(background = root.window.color)
    root.mainloop()    
    main = root.main.get()
    root.destroy()
    
    return main

def postprocess_batch(type, filenames, systems, outputs):
    """Shows the results of the simulation."""

    top = CapSimWindow(buttons = 1)

    if type == 'Multiple Senarios':
        for system in systems:
            output = outputs[systems.index(system)]
            postprocess = Postprocess(top, system, output, batch_type = type)
            filename = postprocess.graphpath + '\\' + filenames[systems.index(system)]
            postprocess.make_csv_file(filename, 'Regular')
    else:
        output     = outputs[0].copy()

        output.O    = mean(array([outputs[i].O   for i in range(len(outputs))]), axis = 0)
        output.OL   = mean(array([outputs[i].OL  for i in range(len(outputs))]), axis = 0)
        output.C    = mean(array([outputs[i].C   for i in range(len(outputs))]), axis = 0)
        output.CL   = mean(array([outputs[i].CL  for i in range(len(outputs))]), axis = 0)
        output.Fi   = mean(array([outputs[i].Fi  for i in range(len(outputs))]), axis = 0)
        output.FiL  = mean(array([outputs[i].FiL for i in range(len(outputs))]), axis = 0)
        output.F    = mean(array([outputs[i].F   for i in range(len(outputs))]), axis = 0)
        output.q    = mean(array([outputs[i].q   for i in range(len(outputs))]), axis = 0)
        output.qL   = mean(array([outputs[i].qL  for i in range(len(outputs))]), axis = 0)
        output.W    = mean(array([outputs[i].W   for i in range(len(outputs))]), axis = 0)
        output.WL   = mean(array([outputs[i].WL  for i in range(len(outputs))]), axis = 0)
        output.qm   = mean(array([outputs[i].qm  for i in range(len(outputs))]), axis = 0)
        output.qmL  = mean(array([outputs[i].qmL for i in range(len(outputs))]), axis = 0)
        output.Cw   = mean(array([outputs[i].Cw  for i in range(len(outputs))]), axis = 0)
        output.FM   = mean(array([outputs[i].FM  for i in range(len(outputs))]), axis = 0)
        output.DM   = mean(array([outputs[i].DM  for i in range(len(outputs))]), axis = 0)
        output.M    = mean(array([outputs[i].M   for i in range(len(outputs))]), axis = 0)
        output.RM   = mean(array([outputs[i].RM  for i in range(len(outputs))]), axis = 0)

        output.O_std    = std(array([outputs[i].O   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.OL_std   = std(array([outputs[i].OL  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.C_std    = std(array([outputs[i].C   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.CL_std   = std(array([outputs[i].CL  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.Fi_std   = std(array([outputs[i].Fi  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.FiL_std  = std(array([outputs[i].FiL for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.F_std    = std(array([outputs[i].F   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.q_std    = std(array([outputs[i].q   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.qL_std   = std(array([outputs[i].qL  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.W_std    = std(array([outputs[i].W   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.WL_std   = std(array([outputs[i].WL  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.qm_std   = std(array([outputs[i].qm  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.qmL_std  = std(array([outputs[i].qmL for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.Cw_std   = std(array([outputs[i].Cw  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.FM_std   = std(array([outputs[i].FM  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.DM_std   = std(array([outputs[i].DM  for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.M_std    = std(array([outputs[i].M   for i in range(len(outputs))]), axis = 0, ddof = 1)
        output.RM_std   = std(array([outputs[i].RM  for i in range(len(outputs))]), axis = 0, ddof = 1)

        output.O_max    = nanmax(array([outputs[i].O   for i in range(len(outputs))]), axis = 0)
        output.OL_max   = nanmax(array([outputs[i].OL  for i in range(len(outputs))]), axis = 0)
        output.C_max    = nanmax(array([outputs[i].C   for i in range(len(outputs))]), axis = 0)
        output.CL_max   = nanmax(array([outputs[i].CL  for i in range(len(outputs))]), axis = 0)
        output.Fi_max   = nanmax(array([outputs[i].Fi  for i in range(len(outputs))]), axis = 0)
        output.FiL_max  = nanmax(array([outputs[i].FiL for i in range(len(outputs))]), axis = 0)
        output.F_max    = nanmax(array([outputs[i].F   for i in range(len(outputs))]), axis = 0)
        output.q_max    = nanmax(array([outputs[i].q   for i in range(len(outputs))]), axis = 0)
        output.qL_max   = nanmax(array([outputs[i].qL  for i in range(len(outputs))]), axis = 0)
        output.W_max    = nanmax(array([outputs[i].W   for i in range(len(outputs))]), axis = 0)
        output.WL_max   = nanmax(array([outputs[i].WL  for i in range(len(outputs))]), axis = 0)
        output.qm_max   = nanmax(array([outputs[i].qm  for i in range(len(outputs))]), axis = 0)
        output.qmL_max  = nanmax(array([outputs[i].qmL for i in range(len(outputs))]), axis = 0)
        output.Cw_max   = nanmax(array([outputs[i].Cw  for i in range(len(outputs))]), axis = 0)
        output.FM_max   = nanmax(array([outputs[i].FM  for i in range(len(outputs))]), axis = 0)
        output.DM_max   = nanmax(array([outputs[i].DM  for i in range(len(outputs))]), axis = 0)
        output.M_max    = nanmax(array([outputs[i].M   for i in range(len(outputs))]), axis = 0)
        output.RM_max   = nanmax(array([outputs[i].RM  for i in range(len(outputs))]), axis = 0)

        output.O_min    = nanmin(array([outputs[i].O   for i in range(len(outputs))]), axis = 0)
        output.OL_min   = nanmin(array([outputs[i].OL  for i in range(len(outputs))]), axis = 0)
        output.C_min    = nanmin(array([outputs[i].C   for i in range(len(outputs))]), axis = 0)
        output.CL_min   = nanmin(array([outputs[i].CL  for i in range(len(outputs))]), axis = 0)
        output.Fi_min   = nanmin(array([outputs[i].Fi  for i in range(len(outputs))]), axis = 0)
        output.FiL_min  = nanmin(array([outputs[i].FiL for i in range(len(outputs))]), axis = 0)
        output.F_min    = nanmin(array([outputs[i].F   for i in range(len(outputs))]), axis = 0)
        output.q_min    = nanmin(array([outputs[i].q   for i in range(len(outputs))]), axis = 0)
        output.qL_min   = nanmin(array([outputs[i].qL  for i in range(len(outputs))]), axis = 0)
        output.W_min    = nanmin(array([outputs[i].W   for i in range(len(outputs))]), axis = 0)
        output.WL_min   = nanmin(array([outputs[i].WL  for i in range(len(outputs))]), axis = 0)
        output.qm_min   = nanmin(array([outputs[i].qm  for i in range(len(outputs))]), axis = 0)
        output.qmL_min  = nanmin(array([outputs[i].qmL for i in range(len(outputs))]), axis = 0)
        output.Cw_min   = nanmin(array([outputs[i].Cw  for i in range(len(outputs))]), axis = 0)
        output.FM_min   = nanmin(array([outputs[i].FM  for i in range(len(outputs))]), axis = 0)
        output.DM_min   = nanmin(array([outputs[i].DM  for i in range(len(outputs))]), axis = 0)
        output.M_min    = nanmin(array([outputs[i].M   for i in range(len(outputs))]), axis = 0)
        output.RM_min   = nanmin(array([outputs[i].RM  for i in range(len(outputs))]), axis = 0)

        postprocess = Postprocess(top, systems[0], output, batch_type = type)

        #for chemical in systems[0].chemicals:
        filename = postprocess.graphpath + '\\' + filenames[0]
        postprocess.make_csv_file(filename, 'Monte Carlo Simulation')

    top.destroy()

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

    if round(dep, 8) > round(1.5*delz, 8):
        ans = interp(zpoint+dep, [dep]+z[jo+1:], Cn[jo:])

    else:
        for j in range(len(z)-1):
            if round(dep, 8) >= round(z[j] + 0.5 * delz,8) and round(dep,8) < round(z[j+1] + 0.5 * delz, 8):
                jj = j + 1
        ans = interp(zpoint+dep, [dep]+z[jj+1:], Cn[jj:])

    return ans
