#!/usr/bin/env python
#
#This script is used to make the window for the output options from CapSim.

from Tkinter             import Frame, Label, Entry, Checkbutton, OptionMenu, StringVar, IntVar, DoubleVar
from capsim_object_types import CapSimWindow

import os, tkMessageBox as tkmb

class InputOptions:
    """Makes a window for the user to specify output options."""

    def __init__(self, master, system):
        """Constructor method."""

        self.master      = master
        self.tframe      = Frame(master.tframe)
        self.frame       = Frame(master.frame)
        self.bframe      = Frame(master.bframe)
        self.fonttype    = system.fonttype
        self.version     = system.version
        self.top         = None

        self.adv         = system.adv
        self.timeunit    = system.timeunit

        self.batchfileoptions          = ['None', 'Multiple Senarios', 'Monte Carlo Simulation']

        self.inputname          = StringVar(value = 'input')
        self.batchfileoption    = StringVar(value = self.batchfileoptions[0])
        self.batchfilename      = StringVar(value = 'input')

        try:
            self.inputname.set(         system.cpsmfilename)
            self.batchfileoption.set(   system.batchfileoption)
            self.batchfilename.set(     system.batchfilename)
        except: pass

    def make_widgets(self):
        """Make the widgets."""

        self.instructions     = Label(self.frame, text = 'Please specify the simulation options for the system:  ')
        self.blank1           = Label(self.frame, text = ' ' * 2, width = 2)
        self.blank2           = Label(self.frame, text = ' ' * 2, width = 10)
        self.blank3           = Label(self.frame, text = ' ' * 2, width = 20)
        self.blank4           = Label(self.frame, text = ' ' * 2, width = 10)
        self.blank5           = Label(self.frame, text = ' ' * 2, width = 2)
        self.blank6           = Label(self.frame, text = ' ' * 2, width = 2)
        self.blank7           = Label(self.frame, text = ' ' * 2, width = 2)

        self.inputfilelabel   = Label(self.frame, text = 'Input parameter filename:')
        self.inputfileentry   = Entry(self.frame, justify = 'center', width = 14, textvariable = self.inputname)
        self.inputsuffix      = Label(self.frame, text = '.cpsm')

        self.inputoptionlabel = Label(self.frame, text = 'Batch file option:')
        self.inputoptionmenu  = OptionMenu(self.frame, self.batchfileoption, *self.batchfileoptions, command = self.updatewidgets)
        self.inputoptionmenu.config(width = 20)

        self.saveaslabel   = Label(self.frame, text = 'Batch filename:')
        self.saveasentry   = Entry(self.frame, justify = 'center', width = 14, textvariable = self.batchfilename)
        self.saveassuffix  = Label(self.frame, text = '.txt')

        #show the widgets that don't change on the grid

        self.instructions.grid(     row = 0, padx = 8, columnspan = 4, sticky = 'W')

        self.blank1.grid(           row = 1, column = 0)
        self.blank2.grid(           row = 1, column = 1)
        self.blank3.grid(           row = 1, column = 2)
        self.blank4.grid(           row = 1, column = 3)
        self.blank5.grid(           row = 1, column = 4)

        self.inputfilelabel.grid(   row = 2, column = 1, sticky = 'E', padx = 2)
        self.inputfileentry.grid(   row = 2, column = 2, sticky = 'WE',padx = 1)
        self.inputsuffix.grid(      row = 2, column = 3, sticky = 'W', padx = 2)

        self.inputoptionlabel.grid( row = 3, column = 1, sticky = 'E',  padx = 2)
        self.inputoptionmenu.grid(  row = 3, column = 2, sticky = 'W',  padx = 1, columnspan = 2)

        self.updatewidgets()


    def updatewidgets(self, default = None):

        try:
            self.saveaslabel.grid_forget()
            self.saveasentry.grid_forget()
            self.saveassuffix.grid_forget()
        except: pass

        if self.batchfileoption.get() == self.batchfileoptions[1] or self.batchfileoption.get() == self.batchfileoptions[2]:

            self.saveaslabel.grid(   row = 4, column = 1, sticky = 'E', padx = 2)
            self.saveasentry.grid(   row = 4, column = 2, sticky = 'WE',padx = 1)
            self.saveassuffix.grid(  row = 4, column = 3, sticky = 'W', padx = 2)

        else:
            self.blank6.grid(           row = 4, column = 0)

        self.blank7.grid(           row = 5, column = 0)

        self.focusbutton = None
        self.master.geometry()
        self.master.center()



def get_inputoptions(system, step):
    """Get the output options."""

    root = CapSimWindow(master = None, buttons = 1)
    root.make_window(InputOptions(root, system))
    root.mainloop()

    if root.main.get() == 1:    system = None
    else:
        if root.step.get() == 1:
            system.get_inputoptions(root.window)
        else:
            system.cpsmfilename     = 'input'
            system.batchfileoption  = 'None'
            system.batchfilename    = 'input'


    root.destroy()

    return system, step + root.step.get()
