#! /usr/bin/env python
#
#This file is used to get the flow and system parameters for the simulation.


from Tkinter                import Frame, Label, Button, Entry, Text, DoubleVar, OptionMenu, StringVar
from capsim_functions       import get_superfont
from capsim_object_types    import CapSimWindow
import tkMessageBox as tkmb

class KblEstimator:
    """Used to make a window to compute the value of "kbl" for Rivers."""

    def __init__(self, master, version, fonttype, bltype, blcoefs):
        """Constructor method."""

        self.version     = version
        self.fonttype    = fonttype
        self.tframe      = Frame(master.tframe)
        self.frame       = Frame(master.frame)
        self.bframe      = Frame(master.bframe)
        self.superfont   = get_superfont(self.fonttype)

        self.bltypes     = ['River', 'Lake']
        self.bltype      = StringVar(value = bltype)
        self.top         = None

        self.cancelflag  = 0

        #Parameters for Mass transfer coefficient in river
        self.vx          = DoubleVar(value = blcoefs['vx'])
        self.n           = DoubleVar(value = blcoefs['n'])
        self.hriver      = DoubleVar(value = blcoefs['hriver'])
        self.rh          = DoubleVar(value = blcoefs['rh'])
        self.nu          = DoubleVar(value = blcoefs['nu'])

        #Parameters for Mass transfer coefficient in lake
        self.rhoair      = DoubleVar(value = blcoefs['rhoair'])
        self.rhowater    = DoubleVar(value = blcoefs['rhowater'])
        self.vwind       = DoubleVar(value = blcoefs['vwind'])
        self.hlake       = DoubleVar(value = blcoefs['hlake'])
        self.llake       = DoubleVar(value = blcoefs['llake'])

    def make_widgets(self):
        """Makes the widgets for the window."""

        self.bgcolor       = self.frame.cget('bg')
        self.text          = Label(self.frame,  text = 'Please provide the following information regarding the benthic boundary layer:')
        self.blank1        = Label(self.frame,  text = '')
        self.blank2        = Label(self.frame,  text = '')
        self.blank3        = Label(self.frame,  text = '')
        self.blank4        = Label(self.frame,  text = '')

        self.paracolumn    = Label(self.frame, text = ' ', width = 20)
        self.valuecolumn   = Label(self.frame, text = ' ', width = 10)
        self.unitcolumn    = Label(self.frame, text = ' ', width = 10)

        self.paralabel     = Label(self.frame, text = 'Parameters')
        self.valuelabel    = Label(self.frame, text = 'Values')
        self.unitlabel     = Label(self.frame, text = 'Units')

        self.bllabel       = Label(self.frame, text = 'Overlying waterbody type')
        self.blwidget      = OptionMenu(self.frame, self.bltype, *self.bltypes, command = self.updateparameters)

        self.vxlabel       = Label(self.frame, text = 'River velocity:')
        self.vxentry       = Entry(self.frame, textvariable = self.vx, width = 8, justify = 'center')
        self.vxunits       = Label(self.frame, text = 'm/s')

        self.nlabel        = Label(self.frame, text = "Manning's n:")
        self.nentry        = Entry(self.frame, textvariable = self.n, width = 8, justify = 'center')
        self.nunits        = Label(self.frame, text = '-')

        self.hlabel        = Label(self.frame, text = 'River depth:')
        self.hentry        = Entry(self.frame, justify = 'center', width = 8, textvariable = self.hriver)
        self.hunits        = Label(self.frame, text = 'm')

        self.rhlabel       = Label(self.frame, text = 'Hydraulic radius:')
        self.rhentry       = Entry(self.frame, textvariable = self.rh, width = 8, justify = 'center')
        self.rhunits       = Label(self.frame, text = 'm')

        self.nulabel       = Label(self.frame, text = 'Kinematic viscosity:')
        self.nuentry       = Entry(self.frame, textvariable = self.nu, width = 8, justify = 'center')
        self.nuunits       = Label(self.frame, text = u'm\u00B2/s')

        self.rhoairlabel   = Label(self.frame, text = 'Density of air:')
        self.rhoairentry   = Entry(self.frame, textvariable = self.rhoair, width = 8, justify = 'center')
        self.rhoairunits   = Label(self.frame, text = 'g/L')

        self.rhowaterlabel = Label(self.frame, text = 'Density of water:')
        self.rhowaterentry = Entry(self.frame, textvariable = self.rhowater, width = 8, justify = 'center')
        self.rhowaterunits = Label(self.frame, text = 'g/L')

        self.vwindlabel    = Label(self.frame, text = 'Wind speed:')
        self.vwindentry    = Entry(self.frame, textvariable = self.vwind, width = 8, justify = 'center')
        self.vwindunits    = Label(self.frame, text = 'm/s')

        self.hlakelabel    = Label(self.frame, text = 'Lake depth:')
        self.hlakeentry    = Entry(self.frame, textvariable = self.hlake, width = 8, justify = 'center')
        self.hlakeunits    = Label(self.frame, text = 'm')

        self.lakelabel     = Label(self.frame,  text = 'Lake fetch:')
        self.lakeentry     = Entry(self.frame,  textvariable = self.llake, width = 8, justify = 'center')
        self.lakeunits     = Label(self.frame,  text = 'm')

        self.calcbutton  = Button(self.frame, text = 'Calculate', width = 20, command = self.calculate)
        self.exitbutton  = Button(self.frame, text = 'Cancel',    width = 20, command = self.exit)
        self.calcbutton.bind('<Return>', self.calculate)
        self.exitbutton.bind('<Return>', self.exit)

        self.text.grid(         row = 0, columnspan = 3, padx = 8)
        self.blank1.grid(       row = 1)

        self.paracolumn .grid(  row = 2, column = 0, sticky = 'WE')
        self.valuecolumn.grid(  row = 2, column = 1, sticky = 'WE')
        self.unitcolumn.grid(   row = 2, column = 2, sticky = 'WE')

        self.paralabel.grid(    row = 3, column = 0, sticky = 'E')
        self.valuelabel.grid(   row = 3, column = 1, sticky = 'WE')
        self.unitlabel.grid(    row = 3, column = 2, sticky = 'W')

        self.bllabel.grid(      row = 4, column = 0, sticky = 'E')
        self.blwidget.grid(     row = 4, column = 1, sticky = 'WE')

        self.blank2.grid(row = 10)
        self.calcbutton.grid(row = 11, columnspan = 3, pady = 1)
        self.exitbutton.grid(row = 12, columnspan = 3, pady = 1)

        self.focusbutton = self.calcbutton
        self.updateparameters()

    def updateparameters(self, event = None):

        try:
            self.vxlabel.grid_forget()
            self.vxentry.grid_forget()
            self.vxunits.grid_forget()
            self.nlabel.grid_forget()
            self.nentry.grid_forget()
            self.nunits.grid_forget()
            self.hlabel.grid_forget()
            self.hentry.grid_forget()
            self.hunits.grid_forget()
            self.rhlabel.grid_forget()
            self.rhentry.grid_forget()
            self.rhunits.grid_forget()
            self.nulabel.grid_forget()
            self.nuentry.grid_forget()
            self.nuunits.grid_forget()
        except: pass

        try:
            self.rhoairlabel.grid_forget()
            self.rhoairentry.grid_forget()
            self.rhoairunits.grid_forget()
            self.rhowaterlabel.grid_forget()
            self.rhowaterentry.grid_forget()
            self.rhowaterunits.grid_forget()
            self.vwindlabel.grid_forget()
            self.vwindentry.grid_forget()
            self.vwindunits.grid_forget()
            self.hlakelabel.grid_forget()
            self.hlakeentry.grid_forget()
            self.hlakeunits.grid_forget()
            self.lakelabel.grid_forget()
            self.lakeentry.grid_forget()
            self.lakeunits.grid_forget()
        except: pass

        if self.bltype.get() == self.bltypes[0]:
            self.vxlabel.grid(      row = 5, column = 0, sticky = 'E', padx = 4)
            self.vxentry.grid(      row = 5, column = 1)
            self.vxunits.grid(      row = 5, column = 2, sticky = 'W')
            self.nlabel.grid(       row = 6, column = 0, sticky = 'E', padx = 4)
            self.nentry.grid(       row = 6, column = 1)
            self.nunits.grid(       row = 6, column = 2, sticky = 'W')
            self.hlabel.grid(       row = 7, column = 0,  sticky = 'E', padx = 4)
            self.hentry.grid(       row = 7, column = 1)
            self.hunits.grid(       row = 7, column = 2,  sticky = 'W')
            self.rhlabel.grid(      row = 8, column = 0, sticky = 'E', padx = 4)
            self.rhentry.grid(      row = 8, column = 1)
            self.rhunits.grid(      row = 8, column = 2, sticky = 'W')
            self.nulabel.grid(      row = 9, column = 0, sticky = 'E', padx = 4)
            self.nuentry.grid(      row = 9, column = 1)
            self.nuunits.grid(      row = 9, column = 2, sticky = 'W')
        else:
            self.rhoairlabel.grid(  row = 5, column = 0, sticky = 'E',padx = 4)
            self.rhoairentry.grid(  row = 5, column = 1)
            self.rhoairunits.grid(  row = 5, column = 2, sticky = 'W')
            self.rhowaterlabel.grid(row = 6, column = 0, sticky = 'E',padx = 4)
            self.rhowaterentry.grid(row = 6, column = 1)
            self.rhowaterunits.grid(row = 6, column = 2, sticky = 'W')
            self.vwindlabel.grid(   row = 7, column = 0, sticky = 'E',padx = 4)
            self.vwindentry.grid(   row = 7, column = 1)
            self.vwindunits.grid(   row = 7, column = 2, sticky = 'W')
            self.hlakelabel.grid(   row = 8, column = 0, sticky = 'E',padx = 4)
            self.hlakeentry.grid(   row = 8, column = 1)
            self.hlakeunits.grid(   row = 8, column = 2, sticky = 'W')
            self.lakelabel.grid(    row = 9, column = 0, sticky = 'E', padx = 4)
            self.lakeentry.grid(    row = 9, column = 1)
            self.lakeunits.grid(    row = 9, column = 2, sticky = 'W')

    def calculate(self, event = None): 
        """Finish and move on."""

        self.frame.quit()

    def exit(self, event = None): 
        """Used to close the window without computing kbl."""

        self.cancelflag = 1

        self.frame.quit()


class TauEstimator:
    """Used to make a window to compute the value of "kbl" for Rivers."""

    def __init__(self, master, version, fonttype, timeunit, lengthunit, concunit, chemicals, taucoefs, BCs, U):
        """Constructor method."""

        self.master      = master
        self.version     = version
        self.fonttype    = fonttype
        self.tframe      = Frame(master.tframe)
        self.frame       = Frame(master.frame)
        self.bframe      = Frame(master.bframe)
        self.superfont   = get_superfont(self.fonttype)

        self.top         = None

        self.cancelflag  = 0

        self.timeunit    = timeunit
        self.lengthunit  = lengthunit
        self.concunit    = concunit

        self.Q           = DoubleVar(value = taucoefs['Q'])
        self.V           = DoubleVar(value = taucoefs['V'])
        self.h           = DoubleVar(value = taucoefs['h'])
        self.Qevap       = DoubleVar(value = taucoefs['Qevap'])
        self.doc         = DoubleVar(value = taucoefs['DOC'])

        self.U           = U

        self.chemicals   = chemicals

        self.h_factor    = 100
        self.lengthunits = [u'\u03BCm', 'cm', 'm']

        if self.lengthunit == self.lengthunits[0]:
            self.h_factor    = self.h_factor * 10000
        elif self.lengthunit == self.lengthunits[2]:
            self.h_factor    = self.h_factor/100

    def make_widgets(self):
        """Makes the widgets for the window."""

        self.bgcolor       = self.frame.cget('bg')
        self.text          = Label(self.frame,  text = 'Please provide the information to estimate the retention time of the overlying water body:')
        self.blank1        = Label(self.frame,  text = '')
        self.blank2        = Label(self.frame,  text = '')
        self.blank3        = Label(self.frame,  text = '')
        self.blank4        = Label(self.frame,  text = '')

        self.leftcolumn    = Label(self.frame, text = ' ', width = 15)
        self.paracolumn    = Label(self.frame, text = ' ', width = 10)
        self.valuecolumn   = Label(self.frame, text = ' ', width = 15)
        self.unitcolumn    = Label(self.frame, text = ' ', width = 10)
        self.rightcolumn   = Label(self.frame, text = ' ', width = 10)

        self.paralabel     = Label(self.frame, text = 'Parameters')
        self.valuelabel    = Label(self.frame, text = 'Values')
        self.unitlabel     = Label(self.frame, text = 'Units')

        self.qlabel       = Label(self.frame, text = 'Inflow rate:')
        self.qentry       = Entry(self.frame, textvariable = self.Q, width = 10, justify = 'center')
        self.qunits       = Label(self.frame, text = 'm' + u'\u00B3/' + self.timeunit)

        self.Vlabel       = Label(self.frame, text = 'Water body volume:')
        self.Ventry       = Entry(self.frame, textvariable = self.V, width = 10, justify = 'center')
        self.Vunits       = Label(self.frame, text = 'm' + u'\u00B3')

        self.hlabel        = Label(self.frame, text = 'Water body depth:')
        self.hentry        = Entry(self.frame, justify = 'center', width = 10, textvariable = self.h)
        self.hunits        = Label(self.frame, text = 'm')

        self.qevaplabel    = Label(self.frame, text = 'Water evaporation rate:')
        self.qevapentry    = Entry(self.frame, textvariable = self.Qevap, width = 10, justify = 'center')
        self.qevapunits    = Label(self.frame, text = 'm' + u'\u00B3/' + self.timeunit)

        self.doclabel      = Label(self.frame, text = 'Water body DOC concentration:')
        self.docentry      = Entry(self.frame, textvariable = self.doc, width = 10, justify = 'center')
        self.docunits      = Label(self.frame, text = 'mg/L')

        self.calcbutton  = Button(self.frame, text = 'OK',        width = 20, command = self.calculate)
        self.exitbutton  = Button(self.frame, text = 'Cancel',    width = 20, command = self.exit)
        self.calcbutton.bind('<Return>', self.calculate)
        self.exitbutton.bind('<Return>', self.exit)

        self.text.grid(         row = 0, columnspan = 6, padx = 8)

        self.paracolumn .grid(  row = 1, column = 1)
        self.valuecolumn.grid(  row = 1, column = 2)
        self.unitcolumn.grid(   row = 1, column = 3)

        self.paralabel.grid(    row = 2, column = 1, sticky = 'E')
        self.valuelabel.grid(   row = 2, column = 2)
        self.unitlabel.grid(    row = 2, column = 3, sticky = 'W')

        self.blank1.grid(       row = 3)

        self.qlabel.grid(       row = 4, column = 1, sticky = 'E')
        self.qentry.grid(       row = 4, column = 2, pady = 4)
        self.qunits.grid(       row = 4, column = 3, sticky = 'W')

        self.Vlabel.grid(       row = 5, column = 1, sticky = 'E',  pady = 2)
        self.Ventry.grid(       row = 5, column = 2, pady = 4)
        self.Vunits.grid(       row = 5, column = 3, sticky = 'W',  pady = 2)

        self.hlabel.grid(       row = 6, column = 1, sticky = 'E',  pady = 2)
        self.hentry.grid(       row = 6, column = 2, pady = 4)
        self.hunits.grid(       row = 6, column = 3, sticky = 'W',  pady = 2)

        self.qevaplabel.grid(   row = 7, column = 1, sticky = 'E',  pady = 2)
        self.qevapentry.grid(   row = 7, column = 2, pady = 4)
        self.qevapunits.grid(   row = 7, column = 3, sticky = 'W',  pady = 2)

        self.doclabel.grid(     row = 8, column = 1, sticky = 'E',  pady = 2)
        self.docentry.grid(     row = 8, column = 2, pady = 4)
        self.docunits.grid(     row = 8, column = 3, sticky = 'W',  pady = 2)

        self.blank2.grid(row = 12)
        self.calcbutton.grid(row = 13, column = 0, columnspan = 6, pady = 1)
        self.exitbutton.grid(row = 14, column = 0, columnspan = 6, pady = 1)

        self.focusbutton = self.calcbutton


    def calculate(self, event = None):
        """Finish and move on."""

        if self.V.get() <= 0:
            tkmb.showerror(title = self.version, message = 'The water body volume can not be zero or negative value')
            self.master.tk.lift()
        elif self.h.get() <= 0:
            tkmb.showerror(title = self.version, message = 'The water body depth can not be zero or negative value')
            self.master.tk.lift()
        elif self.Q.get() - self.Qevap.get() + self.V.get()/self.h.get()*self.U/self.h_factor < 0:
            tkmb.showerror(title = self.version, message = 'The water body will deplete, please change the flowrate or water evaporation rates')
            self.master.tk.lift()
        else: self.frame.quit()

    def exit(self, event = None):
        """Used to close the window without computing kbl."""

        self.cancelflag = 1

        self.frame.quit()


class TauChemicalEditor:
    """Used to make a window to compute the value of "kbl" for Rivers."""

    def __init__(self, master, version, fonttype, timeunit, lengthunit, concunit,  chemicals, BCs):
        """Constructor method."""

        self.version     = version
        self.fonttype    = fonttype
        self.tframe      = Frame(master.tframe)
        self.frame       = Frame(master.frame)
        self.bframe      = Frame(master.bframe)
        self.superfont   = get_superfont(self.fonttype)

        self.cancelflag  = 0

        self.concunit    = concunit
        self.timeunit    = timeunit
        self.lengthunit  = lengthunit
        self.chemicals   = chemicals

        self.Cinits  = []
        self.kevaps  = []
        self.kdecays = []
        for chemical in self.chemicals:
            self.Cinits.append(DoubleVar(value = BCs[chemical.name].Cw.get()))
            self.kevaps.append(DoubleVar(value = BCs[chemical.name].kevap))
            self.kdecays.append(DoubleVar(value = BCs[chemical.name].kdecay))

    def make_widgets(self):
        """Makes the widgets for the window."""

        self.bgcolor       = self.frame.cget('bg')
        self.text          = Label(self.frame,  text = 'Provide the initial water concentration for the water column:')
        self.paralabel     = Label(self.frame, width = 6,  text = 'Parameter' )
        self.unitlabel     = Label(self.frame, width = 6,  text = 'Unit' )
        self.leftcolumn    = Label(self.frame, text = ' ', width = 2)
        self.rightcolumn   = Label(self.frame, text = ' ', width = 2)
        self.blank1        = Label(self.frame,  text = '')
        self.blank2        = Label(self.frame,  text = '')
        self.blank3        = Label(self.frame,  text = '')
        self.blank4        = Label(self.frame,  text = '')

        self.Cinitlabel     = Label(self.frame,  width = 20,  text = 'Initial water concentration')
        self.kevaplabel     = Label(self.frame,  width = 20,  text = 'Evaporation rate coefficient')
        self.kdecaylabel    = Label(self.frame,  width = 20,  text = 'Decay rate coefficient')
        self.Cinitunitlabel = Label(self.frame,  width = 6,  text = self.concunit)
        self.kevapunitlabel = Label(self.frame,  width = 6,  text = self.lengthunit + '/' + self.timeunit)
        self.kdecayunitlabel= Label(self.frame,  width = 6,  text = self.timeunit + u'\u207B\u00B9')

        self.text.grid(         row = 0, columnspan = 11, padx = 8, sticky = 'W')
        self.blank1.grid(       row = 1)

        self.leftcolumn.grid(   row = 2, column = 0, padx = 4, sticky = 'WE')
        self.paralabel.grid(    row = 2, column = 1, padx = 4, sticky = 'WE')
        self.unitlabel.grid(    row = 2, column = 2, padx = 4, sticky = 'WE')

        self.blank3.grid(       row = 3)

        self.Cinitlabel.grid(   row = 4, column = 1, padx = 1, sticky = 'WE', pady = 2)
        self.kevaplabel.grid(   row = 5, column = 1, padx = 1, sticky = 'WE', pady = 2)
        self.kdecaylabel.grid(  row = 6, column = 1, padx = 1, sticky = 'WE', pady = 2)

        self.Cinitunitlabel.grid(   row = 4, column = 2, padx = 1, sticky = 'WE', pady = 2)
        self.kevapunitlabel.grid(   row = 5, column = 2, padx = 1, sticky = 'WE', pady = 2)
        self.kdecayunitlabel.grid(  row = 6, column = 2, padx = 1, sticky = 'WE', pady = 2)

        self.chemicallabels = []
        self.chemicalentrys = []
        self.evapentrys     = []
        self.decayentrys    = []

        column = 3


        for n in range(len(self.chemicals)):
            self.chemicallabels.append(Label(self.frame, text =self.chemicals[n].name, justify = 'center'))
            self.chemicalentrys.append(Entry(self.frame, textvariable = self.Cinits[n], justify = 'center', width = 8))
            self.evapentrys.append(Entry(self.frame, textvariable = self.kevaps[n], justify = 'center', width = 8))
            self.decayentrys.append(Entry(self.frame, textvariable = self.kdecays[n], justify = 'center', width = 8))
            self.chemicallabels[-1].grid(row = 2, column = column, padx = 4, pady = 1, sticky ='WE')
            self.chemicalentrys[-1].grid(row = 4, column = column, padx = 4, pady = 1)
            self.evapentrys[-1].grid(row = 5, column = column, padx = 4, pady = 1)
            self.decayentrys[-1].grid(row = 6, column = column, padx = 4, pady = 1)
            column = column + 1

        self.rightcolumn.grid(   row = 2, column = column, padx = 4, sticky = 'WE')

        self.okbutton       = Button(self.frame, text = 'OK',     width = 20, command = self.OK)
        self.cancelbutton   = Button(self.frame, text = 'Cancel', width = 20, command = self.cancel)
        self.okbutton.bind('<Return>',      self.OK)
        self.cancelbutton.bind('<Return>',  self.cancel)

        self.blank2.grid(row = 7)
        self.okbutton.grid(     row = 8, columnspan = len(self.chemicals) + 5, pady = 1)
        self.cancelbutton.grid( row = 9, columnspan = len(self.chemicals) + 5, pady = 1)

        self.focusbutton = self.okbutton

    def OK(self, event = None):
        """Finish and move on."""

        self.frame.quit()

    def cancel(self, event = None):
        """Used to close the window without computing kbl."""

        self.cancelflag = 1

        self.frame.quit()