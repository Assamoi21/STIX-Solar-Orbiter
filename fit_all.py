from tkinter import *
from tkinter import messagebox
from tkinter.filedialog import askopenfilename
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting, FittableModel, Parameter
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import custom_model
from pandas.plotting import register_matplotlib_converters
import numpy as np
import new_window
import background
import rebin_flux as rebin
from matplotlib.ticker import LogLocator, NullFormatter, AutoMinorLocator
import tkinter as tk
from matplotlib.widgets import SpanSelector

register_matplotlib_converters()


class Fitting:
    """
    Class to perform a spectrum fitting
    """

    # fname = 'solo_L1A_stix-sci-spectrogram' \
    #         '-2102140001_20210214T014006-20210214T015515_008648_V01.fits'
    # rname = 'stx_srm_2021feb14_0140_0155.fits'

    fname = 'solo_L1_stix-sci-xray' \
            '-spec_20230319T175504-20230320T000014_V02_2303197888-65462.fits'
    rname = 'stx_srm_2303197888.fits'

    # create a new window called 'SPEX Fit Options'
    def __init__(self, root):
        """Creates a new window, providing widgets to perform fitting analysis"""
        self.sender = None

        self.top2 = Toplevel()
        self.top2.title('SPEX Fit Options')  # title of the window
        self.top2.geometry("1000x600")  # size of the new window
        #Label(self.top2,
          #    text="Fit Options",  # place the text at the top of the window
          #    fg="red",  # in red
          #    font="Helvetica 12 bold italic").pack()  # with specific text font

        self.root = root
        self.hdul = None            # Opened file
        self.hdul2 = None            # Opened file
        self.name = None            # Name of the .fits file imported
        self.name2  = None #           # Name of the .fits file imported (response matrix)
        
        self.counts = None          # Matrix contaning the counts per band in function of time time
        self.counts_err = None      # Matrix contaning the error of the counts per band in function of time
        self.times = None           # Index of times for x axis
        self.time_del = None        # Time delay for the data
        self.energies = None         # Energy values for y axis
        self.e_low_det = None 
        self.e_high_det = None
        
        self.area = 6               # Area of the surface of detection of the telescope in cm²; used for the flux
        
        self.e_low_true = None
        self.e_high_true = None
        self.matrix = None

        self.background_start = None
        self.background_end = None

        self.data = None            # Converts self.counts in chosen unit (rate, counts or flux) and adds time index
        self.bkg = None             # Calculated background noise
        self.data_bkg = None        # Data without background noise

        self.bkg_start_index = []  # Liste des indices de début pour chaque canal
        self.bkg_end_index = []    # Liste des indices de fin pour chaque canal
        self.var_sep_times = IntVar(value=0)



        self.energy_min_var = tk.DoubleVar()
        self.energy_max_var = tk.DoubleVar()

        self.energy_min_var = tk.DoubleVar(value=0)
        self.energy_min2 = tk.OptionMenu(self.top2, self.energy_min_var, 0)
        #self.energy_min2.pack()

        self.energy_max_var = tk.DoubleVar(value=0)
        self.energy_max2 = tk.OptionMenu(self.top2, self.energy_max_var, 0)
        #self.energy_max2.pack()



        self.sepBkVar = IntVar()

        self.lbl1 = Label(self.top2, text="Choose Fit Function Model:", fg='blue',
                          font=("Helvetica", 11, "bold"))  # name the listbox
        self.lbl1.place(relx=0.07, rely=0.07)  # set the position on window

        self.lbl2 = Label(self.top2, text="Information:", fg='blue',
                          font=("Helvetica", 11, "bold"))  # name the scrollbar
        self.lbl2.place(relx=0.44, rely=0.07)  # set the position

        self.lbl3 = Label(self.top2, text="Choose the files and energy range:", fg='blue',
                          font=("Helvetica", 11, "bold"))  # name the scrollbar
        self.lbl3.place(relx=0.65, rely=0.07)  # set the position

        # self.lblFunc = Label(self.top2, text="Set function components: ")  # name the scrollbar
        # self.lblFunc.place(relx=0.73, rely=0.20)  # set the position

        # Spectrum: file name
        self.label_filename = Label(self.top2, text="Spectrum: ")
        self.label_filename.place(relx=0.65, rely=0.2, anchor=W)
        self.text_filename = Entry(self.top2, width=30)
        self.text_filename.place(relx=0.72, rely=0.2, anchor=W)
        if Fitting.fname:
            self.text_filename.insert(0, Fitting.fname)
            self.open_file(Fitting.fname)
        else:
            self.text_filename.insert(0, "No file chosen")
        
        self.btn_browse = Button(self.top2, text='Browse ->', command=self.open_file)
        self.btn_browse.place(relx=0.92, rely=0.2, anchor=W)

        # Response matrix: file name
        self.label_filename2 = Label(self.top2, text="Response: ")
        self.label_filename2.place(relx=0.65, rely=0.25, anchor=W)
        self.text_filename2 = Entry(self.top2, width=30)
        self.text_filename2.place(relx=0.72, rely=0.25, anchor=W)
        if Fitting.rname:
            self.text_filename2.insert(0, Fitting.rname)
            self.open_srm_file(Fitting.rname)
        else:
            self.text_filename2.insert(0, "No file chosen")
        
        self.btn_browse2 = Button(self.top2, text='Browse ->', command=self.open_srm_file)
        self.btn_browse2.place(relx=0.92, rely=0.25, anchor=W)

        self.fit_model = str()

        # self.X_Label = Label(self.top2, text="Energy range to fit: ")  # name "Energy range(s) to fit"
        # self.X_Label.place(relx=0.75, rely=0.32)  # locate

        # FIXME: Empty Window, surely to be removed later on
        def Set_Function():  # new window Set_Function definition
            """Creates a new window for "Set spec_data axis" part"""
            newwin = Toplevel(root)
            newwin.title('Function values')  # title of the window
            newwin.geometry("600x400")  # size of the new window
            display = Label(newwin, text="Choose function values: ", fg='blue', font=("Helvetica", 11, "bold"))
            display.place(relx=0.04, rely=0.07)

        self.lblFunc = Label(self.top2, text="Set function components: ")  # name the scrollbar
        self.lblFunc.place(relx=0.75, rely=0.30)

        self.Value_Button = Button(self.top2, text="Function value(s)",
                                   command=Set_Function)  # place a "Function value" button
        self.Value_Button.place(relx=0.75, rely=0.35, relheight=0.05, relwidth=0.13)

        # Energies range(s) to fit

        fname = Fitting.fname
        rname = Fitting.rname
        if fname is None or rname is None:  # if file not choosen, print
            print('Please, choose input file')

        else:
            # counts_file = fits.open(fname)
            # srm_file = fits.open(rname)

            # Initialisation des listes pour le fond
            n_channels = len(self.e_low_det)

            self.bkg_start_index = [0] * n_channels
            self.bkg_end_index = [len(self.times)] * n_channels


            usable_channels = np.arange(min(self.matrix.shape[1], len(self.e_low_det)))

            e_low_det = self.e_low_det[usable_channels]
            e_high_det = self.e_high_det[usable_channels]

            print("Usable channels: ", usable_channels)
            print("Energies law: ", e_low_det)
            print("Energies high: ", e_high_det)

            self.text_min_energy = Label(self.top2, text="Min energy")
            self.text_min_energy.place(relx=0.75, rely=0.45, anchor=N)
            self.text_max_energy = Label(self.top2, text="Max energy")
            self.text_max_energy.place(relx=0.85, rely=0.45, anchor=N)

            e_low_values = sorted(set(e_low_det)) 
            e_high_values = sorted(set(e_high_det)) 

            e_high_values = [e for e in e_high_values if e != float('inf') and e != float('-inf')]

            e_low_values_int = [int(e) for e in e_low_values if e != 0]
            e_high_values_int = [int(e) for e in e_high_values]

            self.energy_min_var = IntVar()
            self.energy_max_var = IntVar()

            self.energy_min_var.set(min(e_low_values_int))
            self.energy_max_var.set(max(e_high_values_int))

            # Créer les OptionMenu pour l'énergie minimale et maximale
            self.energy_min2 = OptionMenu(self.top2, self.energy_min_var, *e_low_values_int)
            self.energy_max2 = OptionMenu(self.top2, self.energy_max_var, *e_high_values_int)

            self.energy_min2.place(relx=0.75, rely=0.50, anchor=N)
            self.energy_max2.place(relx=0.85, rely=0.50, anchor=N)
      

        # ============== Main window description ==============

        self.lbox = Listbox(self.top2, selectmode=EXTENDED, highlightcolor='red', bd=4, selectbackground='grey')
        """ 
        On the left side of the 'SPEX Fit Options' window: place a list of text alternatives (listbox).
        The user can choose(highlight) one of the options.
        Options(functions):
        1) One Dimensional Power Law;
        2) 1-D Broken Power Law;
        3) Single Power Law Times an Exponetial
        """
        self.lbox.place(relx=0.05, rely=0.15, relheight=0.45, relwidth=0.25)

        self.scroll = Scrollbar(self.top2, command=self.lbox.yview)
        self.scroll.place(relx=0.3, rely=0.15, relheight=0.45, relwidth=0.02)
        self.lbox.config(yscrollcommand=self.scroll.set)

        # New frame at the bottom. Locate there 'Plot Units' and 'Do Fit' widgets
        self.frameFit = LabelFrame(self.top2, relief=RAISED,
                                   borderwidth=10)  # determine the border of the frame and size
        self.frameFit.place(relx=0.05, rely=0.63, relheight=0.25, relwidth=0.85)  # the frame position

        self.PlotUnits5 = Label(self.frameFit, text="Plot Units: ", fg='blue',
                                font=("Helvetica", 11, "bold"))  # lay out new text file
        self.PlotUnits5.place(relx=0.04, rely=0.4)

        # Add button for Units: Rate, Counts, Flux
        # Allows user to make a choice between three parameters
        self.Component_choicesFit = ('Rate', 'Counts', 'Flux')
        self.var = StringVar(self.frameFit)
        self.var.set(self.Component_choicesFit[0])
        self.selection = OptionMenu(self.frameFit, self.var, *self.Component_choicesFit)
        self.selection.place(relx=0.15, rely=0.38, relheight=0.23, relwidth=0.15)

        self.show_params_var = IntVar(value=1)  # Par défaut cochée
        self.show_params_check = Checkbutton(
            self.frameFit,
            text="Display parameters",
            variable=self.show_params_var
        )
        self.show_params_check.place(relx=0.35, rely=0.7)

        self.grid_var = IntVar(value=0) 
        self.grid_check = Checkbutton(
            self.frameFit,
            text="Show grid",
            variable=self.grid_var
        )
        self.grid_check.place(relx=0.55, rely=0.7)

        self.show_db_var = IntVar(value=0) 
        self.show_db_check = Checkbutton(
            self.frameFit,
            text="Data-Background",
            variable=self.show_db_var,
            command=self.on_background_check
        )

        self.show_db_check.place(relx=0.35, rely=0.5)

        self.show_photon_var = IntVar(value=0) 
        def on_photon_toggle():
            if self.show_photon_var.get():
                self.ask_photon_axes_scale()

        self.show_photon_check = Checkbutton(
            self.frameFit,
            text="Photon",
            variable=self.show_photon_var,
            command=on_photon_toggle  
        )

        self.show_photon_check.place(relx=0.55, rely=0.5)


        self.DoFit5_Button = Button(self.frameFit, text="Do Fit",
                                    command=self._selective_fit)  # place a "Do Fit" button
        self.DoFit5_Button.place(relx=0.70, rely=0.38, relheight=0.23, relwidth=0.15)  # locate

        self.refreshButton5 = Button(self.top2, text="Refresh")  # add Refresh button at the buttom
        # resets original view
        self.refreshButton5.place(relx=0.4, rely=0.94)

        """Scrollbar with information related to each function"""
        self.closeButton5 = Button(self.top2, text="Close", command=self.destroy5)  # add Close button
        # Close "Fit Options" window
        self.closeButton5.place(relx=0.5, rely=0.94)
        self.models = ['PowerLaw1D', 'BrokenPowerLaw1D','Single Power Law Times an Exponential', 'V_TH', 'V_TH + PowerLaw']  # function names
        for p in self.models:
            """On the right: place an 'entry text' Scrollbar widget (scrollbar) When user highlight the function, 
            displays the text information about function description and input parameters"""
            self.lbox.insert(END, p)
        self.lbox.bind("<<ListboxSelect>>", self.onSelect)
        self.list = {'PowerLaw1D': {'One dimensional power law model', '\n\n',
                                    'amplitude – model amplitude at the reference energy', '\n',
                                    'energy_data – reference energy', '\n', 'alpha – power law index'},
                     # if user choose PowerLaw1D, display
                     'BrokenPowerLaw1D': {'One dimensional power law model with a break', '\n\n',
                                          'amplitude - model amplitude at the break energy', '\n',
                                          'alpha 1 – power law index for energy_data<x_break', '\n',
                                          'alpha 2 – power law index for energy_data>x_break'},
                     # if user choose BrokenPowerLaw1D, display
                     'Gaussian': {'Single Gaussian function(high quality), width in sigma', '\n',
                                  'does not go through DRM', '\n',
                                  'This function returns the sum of Gaussian and ', '\n', '2nd order Polynomial',
                                  'amplitude - integrated intensity, mean - centroid', '\n', 'stddev - sigma'},
                     # if user choose Gaussian, display
                     'Polynomial': {'Polynomial function with offset in energy_data', '\n',
                                    'c0 - 0th order coefficient', '\n', 'c1 - 1st order coefficient', '\n',
                                    'c2 - 2nd order coefficient', '\n',
                                    'c3 - 3rd order coefficient', '\n', 'c4 - 4th order coefficient', '\n',
                                    'c5 - energy_data offset, such that function value at energy_data = c5 is C0 '},  # Polynomial
                     'Exponential': {'Exponential function', '\n', 't0 - Normalization', '\n',
                                     't1 - Pseudo temperature'},  # Exponential
                     'Single Power Law Times an Exponential': {'Multiplication of Single Power Law and Exponential',
                                                              '\n',
                                                              'p0 - normalization at epivot for power-law', '\n',
                                                              'p1 - negative power - law index', '\n',
                                                              'p2 - epivot (kEv) for power - law', '\n',
                                                              'e1 - normalization for exponential', '\n',
                                                              'e2 - pseudo temperature for exponential'},
                     # Single Power Law Times an Exponential
                     'Logistic Regression': {'Returns a sigmoid function', '\n'},  # Logistic Regression
                     'Lorentz': {'One dimensional Lorentzian model', '\n\n',
                                 'Amplitude correponds to peak value', '\n',
                                 'x_0 is the peak position (default value is 0)'},  # Lorentz Model
                     'Moffat': {'able to accurately reconstruct point spread functions', '\n',
                                'Moffat distribution'},  # Moffat model
                     'Voigt Profile': {'model computes the sum of Voigt function with a 2nd order polynomial', '\n',
                                       'amplitude centered at x_0 with the specified Lorentzian and Gaussian widths'},
                     # Voigt
                     'V_TH': {'Thermal Bremsstrahlung Model', '\n',
                                    'T - Temperature (keV)', '\n',
                                    'EM - Emission Measure (cm^-3)'},
                     'V_TH + PowerLaw': {'Addition of V_TH and Single Power Law', '\n',
                                    'T - Temperature (keV)', '\n',
                                    'EM - Emission Measure (cm^-3)', '\n',
                                    'Amplitude - Model amplitude at the reference energy', '\n',
                                    'Alpha - Power law index'
                                    }
                     }
                    
        self.list_selection = Listbox(self.top2, highlightcolor='red', bd=4)
        self.list_selection.place(relx=0.33, rely=0.15, relheight=0.45, relwidth=0.30)

        # print("counts: ", self.counts.shape)
        # print("counts_err: ", self.counts_err.shape)
        # print("times: ", self.times.shape)
        # print("time_del: ", self.time_del.shape)
        # print("e_low_det: ", self.e_low_det.shape)
        # print("e_high_det: ", self.e_high_det.shape)
        # print("matrix: ", self.matrix.shape)
        # print("e_low_true: ", self.e_low_true.shape)
        # print("e_high_true: ", self.e_high_true.shape)

    def open_file(self, file=None):
        """Reads the input data using Astropy library. It can be any extension. RHESSI .fits files are analysed. \n
        Parameters: \n
            file: if a file has already been opened previously (i.e. in background), automatically re-reads it instead
            of asking user to choose it again."""
        if file:
            self.name = file
        else:
            self.name = askopenfilename(initialdir=".",
                                        filetypes=(("FITS files", "*.fits"), ("All Files", "*.*")),
                                        title="Please Select Spectrum or Image File")
        self.text_filename.delete(0, 'end')
        Fitting.fname = self.name

        if self.name:  # If file has been chosen by user
            with fits.open(self.name) as hdul:
                self.hdul = hdul
                self.text_filename.insert(0, self.name)  # Displays the input file name in Entry box
            # Loading data
            data1 = self.load_data(self.name)
            self.times = data1['time']
            self.counts = data1['counts']
            self.counts_err = data1['counts_err']
            self.e_high_det = data1['e_high']
            self.e_low_det = data1['e_low']
            self.time_del = data1['timedel']
            self.update_energy_range()
        else:
            self.text_filename.insert(0, "No file chosen")  

    def open_srm_file(self, file=None):
        """Reads the input data using Astropy library. It can be any extension. RHESSI .fits files are analysed. \n
        Parameters: \n
            file: if a file has already been opened previously (i.e. in background), automatically re-reads it instead
            of asking user to choose it again."""
        if file:
            self.name2 = file
        else:
            self.name2 = askopenfilename(initialdir=".",
                                        filetypes=(("FITS files", "*.fits"), ("All Files", "*.*")),
                                        title="Please Select Spectrum or Image File")
        self.text_filename2.delete(0, 'end')
        Fitting.rname = self.name2

        if self.name2:  # If file has been chosen by user
            with fits.open(self.name2) as hdul:
                self.hdul2 = hdul
                self.text_filename2.insert(0, self.name2)  # Displays the input file name in Entry box
            # Loading data
            data = self.load_srm_data(self.name2)
            self.e_low_true = data['ENERG_LO']
            self.e_high_true = data['ENERG_HI']
            self.matrix = data['MATRIX']
            self.update_energy_range()
        else:
            self.text_filename2.insert(0, "No file chosen")

    def update_energy_range(self):
        if self.e_low_det is None or self.e_high_det is None or self.matrix is None:
            return

        usable_channels = np.arange(min(self.matrix.shape[1], len(self.e_low_det)))

        e_low_det = self.e_low_det[usable_channels]
        e_high_det = self.e_high_det[usable_channels]

        e_low_values = sorted(set(e_low_det))
        e_high_values = sorted(set(e_high_det))
        e_high_values = [e for e in e_high_values if e != float('inf') and e != float('-inf')]
        e_low_values_int = [int(e) for e in e_low_values if e != 0]
        e_high_values_int = [int(e) for e in e_high_values]

        self.energy_min_var.set(min(e_low_values_int))
        self.energy_max_var.set(max(e_high_values_int))

        menu = self.energy_min2["menu"]
        menu.delete(0, "end")
        for val in e_low_values_int:
            menu.add_command(label=val, command=lambda v=val: self.energy_min_var.set(v))

        menu = self.energy_max2["menu"]
        menu.delete(0, "end")
        for val in e_high_values_int:
            menu.add_command(label=val, command=lambda v=val: self.energy_max_var.set(v))

    @staticmethod
    def load_data(file):
        """Reads the Data and Header contents from input file. Loads the input file choosen in 'Select Input' section.
        Returns respectively a table containing datas, energies, dates and channels.\n
        Parameters: \n
            file: contains the data in a fits file."""
        hdulist = fits.open(file)   # Reads the data
        hdulist.info()              # Displays the content of the read file

        result = {}

        for hdu in hdulist:
            if not hasattr(hdu, 'columns'):
                continue
            colnames = hdu.columns.names

            # Time & Timedel
            if 'time' in colnames and 'timedel' in colnames:
                result['time'] = hdu.data['time']
                result['timedel'] = hdu.data['timedel']

            # Counts
            if 'counts' in colnames:
                result['counts'] = hdu.data['counts']
            if 'counts_comp_err' in colnames or 'counts_err' in colnames:
                err_col = 'counts_comp_err' if 'counts_comp_err' in colnames else 'counts_err'
                result['counts_err'] = hdu.data[err_col]

            # Triggers (optionnel)
            if 'triggers' in colnames:
                result['triggers'] = hdu.data['triggers']

            # Energy bins
            if 'e_low' in colnames and 'e_high' in colnames:
                result['e_low'] = hdu.data['e_low']
                result['e_high'] = hdu.data['e_high']

            # Version info
            if 'obt_start' in colnames and 'obt_end' in colnames:
                result['obt_start'] = hdu.data['obt_start']
                result['obt_end'] = hdu.data['obt_end']

        # Vérifications de base
        required_keys = ['counts', 'counts_err', 'e_low', 'e_high', 'time', 'timedel']
        for key in required_keys:
            if key not in result:
                print(f"⚠️  Attention : {key} non trouvé dans le FITS.")
        return result

        # return hdulist[2].data, hdulist[3].data, hdulist[0].header, hdulist[3].header
    
    @staticmethod
    def load_srm_data(file):
        """Reads the Data and Header contents from input file. Loads the input file choosen in 'Select Input' section.
        Returns respectively a table containing datas, energies, dates and channels.\n
        Parameters: \n
            file: contains the data in a fits file."""
        hdulist = fits.open(file)   # Reads the data
        hdulist.info()              # Displays the content of the read file

        result = {}

        for hdu in hdulist:
            if not hasattr(hdu, 'columns'):
                continue
            colnames = hdu.columns.names

            # Matrix
            if 'MATRIX' in colnames:
                result['MATRIX'] = hdu.data['MATRIX']

            # Energy bins
            if 'ENERG_LO' in colnames and 'ENERG_HI' in colnames:
                result['ENERG_LO'] = hdu.data['ENERG_LO']
                result['ENERG_HI'] = hdu.data['ENERG_HI']


        # Vérifications de base
        required_keys = ['MATRIX', 'ENERG_LO', 'ENERG_HI']
        for key in required_keys:
            if key not in result:
                print(f"⚠️  Attention : {key} non trouvé dans le FITS.")
        return result
        # return hdulist[1].data

    @staticmethod
    def editEnergy(p1):
        """Call new class to edit spec_data axis"""
        new_window.Set_Energy(p1)

    def onSelect(self, event):
        """Affiche les infos de la fonction sélectionnée dans la Listbox de droite."""
        try:
            # Récupère l’index de la sélection
            selected_index = self.lbox.curselection()[0]
            selected_name = self.lbox.get(selected_index)
            self.fit_model = selected_name

            # Réactive temporairement la Listbox info
            self.list_selection.config(state='normal')
            self.list_selection.delete(0, END)

            # Récupère et insère les infos correspondantes
            if selected_name in self.list:
                for line in sorted(list(self.list[selected_name])):
                    self.list_selection.insert(END, line)
            else:
                self.list_selection.insert(END, "Aucune information disponible.")

            # Désactive la Listbox info (pour la rendre non cliquable)
            self.list_selection.config(state='disabled')

        except Exception as e:
            print("Erreur dans onSelect :", e)


    def update_file_list(self, file_list):
        """Updating the frame (in information:) and adding new function description, related to the user choice"""
        self.list_selection.delete(0, END)
        for i in file_list:
            self.list_selection.insert(END, i)

    def findfiles(self, val):
        """Finding the information related to the function name"""
        self.sender = val.widget

    def destroy5(self):
        """Closing 'SPEX Fit Options' window"""
        self.top2.destroy()

    # function to calculate the flux
    def integrate_flux(e1, e2, model_func, n_points=10):
        energies = np.linspace(e1, e2, n_points)
        fluxes = model_func(energies)
        return np.trapezoid(fluxes, energies) / (e2 - e1)
    
    class ForwardFoldedPowerLaw(FittableModel):
        n_inputs = 1
        n_outputs = 1

        amplitude = Parameter(default=1e-2)
        alpha = Parameter(default=2.0)
        x_0 = 100.0  # énergie pivot en keV, fixe ici

        def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
            super().__init__(**kwargs)
            self.e_low_true = e_low_true
            self.e_high_true = e_high_true
            self.matrix = matrix
            self.exposure = exposure

        def evaluate(self, x, amplitude, alpha):
            model_func = lambda E: amplitude * (E / self.x_0) ** (-alpha)
            true_fluxes = np.array([
                Fitting.integrate_flux(e1, e2, model_func)
                for e1, e2 in zip(self.e_low_true, self.e_high_true)
            ])
            folded = np.dot(true_fluxes, self.matrix) / self.exposure
            return folded

    # === Forward Folded Broken Power Law ===
    class ForwardFoldedBrokenPowerLaw(FittableModel):
        n_inputs = 1
        n_outputs = 1

        amplitude = Parameter(default=1e-2)
        E_break = Parameter(default=10.0)
        alpha_1 = Parameter(default=2.0)
        alpha_2 = Parameter(default=3.0)

        def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
            super().__init__(**kwargs)
            self.e_low_true = e_low_true
            self.e_high_true = e_high_true
            self.matrix = matrix
            self.exposure = exposure

        def evaluate(self, x, amplitude, E_break, alpha_1, alpha_2):
            def model_func(E):
                return amplitude * np.where(E < E_break, (E / E_break)**(-alpha_1), (E / E_break)**(-alpha_2))

            true_fluxes = np.array([
                Fitting.integrate_flux(e1, e2, model_func)
                for e1, e2 in zip(self.e_low_true, self.e_high_true)
            ])
            folded = np.dot(true_fluxes, self.matrix) / self.exposure
            return folded

    # === Forward Folded VTH ===
    class ForwardFoldedVTH(FittableModel):
        n_inputs = 1
        n_outputs = 1

        # T = Parameter(default=10.0)      # Température en keV
        # EM = Parameter(default=1e49)     # Emission Measure en cm^-3

        EM = Parameter(default=1e48, bounds=(1e44, 1e52))
        T = Parameter(default=1.0, bounds=(0.1, 50.0))

        def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
            super().__init__(**kwargs)
            self.e_low_true = e_low_true
            self.e_high_true = e_high_true
            self.matrix = matrix
            self.exposure = exposure

        def evaluate(self, x, EM, T):
            # Constantes
            gff = 1.2
            A = 1.07e-42 * gff

            # Sécurité : éviter division par zéro ou valeurs négatives
            safe_T = max(1e-3, T)

            def thermal_model(E):
                return (A * EM) / (E * np.sqrt(safe_T)) * np.exp(-E / safe_T)

            true_fluxes = np.array([
                Fitting.integrate_flux(e1, e2, thermal_model)
                for e1, e2 in zip(self.e_low_true, self.e_high_true)
            ])

            folded = np.dot(true_fluxes, self.matrix) / self.exposure
            return folded

    # === Forward Folded Exponential Power Law ===
    class ForwardFoldedExpPowerLaw(FittableModel):
        n_inputs = 1
        n_outputs = 1

        p0 = Parameter(default=1.0)
        p1 = Parameter(default=-2.0)
        p2 = Parameter(default=20.0)
        e3 = Parameter(default=1.0)
        e4 = Parameter(default=10.0)

        def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
            super().__init__(**kwargs)
            self.e_low_true = e_low_true
            self.e_high_true = e_high_true
            self.matrix = matrix
            self.exposure = exposure

        def evaluate(self, x, p0, p1, p2, e3, e4):
            def model_func(E):
                safe_E = np.where(E <= 0, 1e-6, E)
                return (p0 * (safe_E / p2)**p1) * np.exp(e3 - safe_E / e4)

            true_fluxes = np.array([
                Fitting.integrate_flux(e1, e2, model_func)
                for e1, e2 in zip(self.e_low_true, self.e_high_true)
            ])
            folded = np.dot(true_fluxes, self.matrix) / self.exposure
            return folded
    
    # === Forward Folded VTH + Power Law ===
    class ForwardFoldedVTHPlusPowerLaw(FittableModel):
        n_inputs = 1
        n_outputs = 1

        # Paramètres VTH
        EM = Parameter(default=1e48, bounds=(1e44, 1e52))
        T = Parameter(default=1.0, bounds=(0.1, 50.0))

        # Paramètres Power Law
        amplitude = Parameter(default=1e-2)
        alpha = Parameter(default=2.0)
        x_0 = 100.0  # énergie pivot en keV

        def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
            super().__init__(**kwargs)
            self.e_low_true = e_low_true
            self.e_high_true = e_high_true
            self.matrix = matrix
            self.exposure = exposure

        def evaluate(self, x, EM, T, amplitude, alpha):
            # Constante de Gaunt
            gff = 1.2
            A = 1.07e-42 * gff
            safe_T = max(1e-3, T)

            def model_total(E):
                # Thermal component
                thermal = (A * EM) / (E * np.sqrt(safe_T)) * np.exp(-E / safe_T)
                # Power-law component
                power = amplitude * (E / self.x_0) ** (-alpha)
                return thermal + power

            # Intégration du flux photonique dans chaque bin SRM
            true_fluxes = np.array([
                Fitting.integrate_flux(e1, e2, model_total)
                for e1, e2 in zip(self.e_low_true, self.e_high_true)
            ])

            # Forward-folding via SRM
            folded = np.dot(true_fluxes, self.matrix) / self.exposure
            return folded

    def ask_photon_axes_scale(self):
        """Ouvre une popup centrée pour choisir les échelles X et Y du graphe photonique."""
        def confirm():
            self.photon_xscale = x_choice.get()
            self.photon_yscale = y_choice.get()
            popup.destroy()

        popup = Toplevel(self.top2)
        popup.title("Photon Plot Axes")

        # Taille désirée
        window_width = 400
        window_height = 200

        # Calculer la position centrée
        screen_width = popup.winfo_screenwidth()
        screen_height = popup.winfo_screenheight()
        pos_x = int((screen_width / 2) - (window_width / 2))
        pos_y = int((screen_height / 2) - (window_height / 2))

        popup.geometry(f"{window_width}x{window_height}+{pos_x}+{pos_y}")
        popup.resizable(False, False)

        # Interface utilisateur
        Label(popup, text="Choose axes scale for photon model:", font=("Helvetica", 11, "bold")).pack(pady=10)

        Label(popup, text="X axis scale:").pack()
        x_choice = StringVar(popup)
        x_choice.set("log")
        OptionMenu(popup, x_choice, "linear", "log").pack()

        Label(popup, text="Y axis scale:").pack()
        y_choice = StringVar(popup)
        y_choice.set("log")
        OptionMenu(popup, y_choice, "linear", "log").pack()

        Button(popup, text="Confirm", command=confirm, bg="#4CAF50", fg="white").pack(pady=10)

    def on_background_check(self):
        if self.show_db_var.get():  # Si check activé
            # self.select_background_range()
            messagebox.showwarning("Data-Background", "Not yet ready.")
            self.show_db_var.set(0)  # Désactiver le check
            return


    # def select_background_range(self):
    #     """Affiche une fenêtre pour sélectionner la plage d’énergie de fond"""

    #     # Données nécessaires
    #     background_counts = np.mean(self.counts, axis=0)
    #     e_low = self.e_low_det
    #     e_high = self.e_high_det
    #     edges = np.append(e_low, e_high[-1])
    #     exposure = np.mean(self.time_del)
    #     dE_det = np.diff(edges)

    #     fig, ax = plt.subplots(figsize=(10, 5))  # taille agrandie
    #     plt.subplots_adjust(bottom=0.2)

    #     # Counts
    #     mean_counts = background_counts

    #     # Rate
    #     rate = mean_counts / exposure

    #     # Flux (photons / s / cm² / keV)
    #     flux = rate / (self.area * dE_det)

    #     # Unit selection
    #     unit = self.var.get()
        
    #     if unit == 'Rate':
    #         self.data = rate
    #     elif unit == 'Counts':
    #         self.data = mean_counts
    #     elif unit == 'Flux':
    #         self.data = flux
    #     else:
    #         raise ValueError("Choose unit = 'rate', 'counts' ou 'flux'")

    #     ax.step(edges[:-1], self.data, where='mid', color='red', label='Data')

    #     def onselect(xmin, xmax):
    #         self.background_start = xmin
    #         self.background_end = xmax
    #         print(f"✔️ Background interval selected: [{xmin:.2f}, {xmax:.2f}] keV")
    #         plt.close(fig)

    #     span = SpanSelector(
    #         ax, onselect, 'horizontal', useblit=True,
    #         props=dict(alpha=0.5, facecolor='gray')
    #     )

    #     ax.set_title("Select Background Energy Range and Close Window")
    #     ax.set_xlabel("Energy (keV)")
    #     ax.set_ylabel("Counts")
    #     ax.grid(True)
    #     ax.legend()
    #     plt.show()
     
    # def get_data_bkg(self):
    #     """Soustrait le fond au signal, canal par canal et temps par temps."""
    #     self.data_bkg = np.zeros_like(self.counts)

    #     for band in range(len(self.e_low_det)):
    #         for time in range(len(self.times)):
    #             self.data_bkg[time, band] = self.counts[time, band] - self.bkg[time, band]
    #             if self.data_bkg[time, band] < 0:
    #                 self.data_bkg[time, band] = 0

    #     print("✅ Données - Fond calculées :", self.data_bkg.shape)

    # def get_bkg(self):
    #     """Calcule le fond pour chaque canal énergétique selon un intervalle temporel spécifié."""
    #     self.bkg = np.zeros_like(self.counts)

    #     for band in range(len(self.e_low_det)):
    #         # Par défaut : même méthode et intervalle pour tous les canaux
    #         if self.var_sep_times.get() == 1 and band != 0:
    #             self.bkg_start_index[band] = self.bkg_start_index[0]
    #             self.bkg_end_index[band] = self.bkg_end_index[0]

    #         start = self.bkg_start_index[band]
    #         end = self.bkg_end_index[band]

    #         mean_bkg = np.mean(self.counts[start:end, band])
    #         self.bkg[:, band] = mean_bkg

    #         # Mettre à zéro les valeurs négatives
    #         self.bkg[:, band][self.bkg[:, band] < 0] = 0

    #     print("✅ Background calculé pour tous les canaux :", self.bkg.shape)


    # def select_background_range(self):
    #     """Sélectionne la plage d’énergie pour le fond via graphique"""

    #     # data used
    #     background_counts = np.mean(self.counts, axis=0)
    #     e_low = self.e_low_det
    #     e_high = self.e_high_det
    #     edges = np.append(e_low, e_high[-1])
    #     exposure = np.mean(self.time_del)
    #     dE_det = np.diff(edges)

    #     fig, ax = plt.subplots(figsize=(10, 5))
    #     plt.subplots_adjust(bottom=0.2)

    #     # calculating data according to the selected unit
    #     mean_counts = background_counts
    #     rate = mean_counts / exposure
    #     flux = rate / (self.area * dE_det)

    #     unit = self.var.get()
    #     if unit == 'Rate':
    #         self.data = rate
    #         y_label = "Rate"
    #     elif unit == 'Counts':
    #         self.data = mean_counts
    #         y_label = "Counts"
    #     elif unit == 'Flux':
    #         self.data = flux
    #         y_label = "Flux"
    #     else:
    #         raise ValueError("Choose unit = 'rate', 'counts' ou 'flux'")

    #     ax.step(edges[:-1], self.data, where='mid', color='red', label='Data')
    #     ax.set_title("Select Background Energy Range and Close Window")
    #     ax.set_xlabel("Energy (keV)")
    #     ax.set_ylabel(y_label)
    #     ax.grid(True)
    #     ax.legend()

    #     def onselect(xmin, xmax):
    #         # Convertir la sélection en indices
    #         start_idx = np.searchsorted(edges[:-1], xmin, side='left')
    #         end_idx = np.searchsorted(edges[1:], xmax, side='right')

    #         self.background_channel_start = start_idx
    #         self.background_channel_end = end_idx

    #         print(f"✔️ Background interval selected: {xmin:.2f}–{xmax:.2f} keV")
    #         print(f"↪️ Channels: {start_idx} to {end_idx}")
    #         plt.close(fig)

    #     span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
    #                         props=dict(alpha=0.5, facecolor='gray'))

    #     plt.show()

    # def get_bkg(self):
    #     """calculate the background for each energy channel according to a specified time interval."""
    #     if not hasattr(self, 'background_channel_start') or not hasattr(self, 'background_channel_end'):
    #         print("❌ No background range selected. Use select_background_range() first.")
    #         return

    #     start = self.background_channel_start
    #     end = self.background_channel_end

    #     # mean background for all times
    #     # self.bkg = np.median(self.counts[:, start:end], axis=1).reshape(-1, 1)  # (time, 1)
    #     # self.bkg = np.repeat(self.bkg, self.counts.shape[1], axis=1)  # Broadcast (time, channel)

    #     # Moyenne sur les canaux sélectionnés pour chaque instant de temps
    #     background_per_time = np.median(self.counts[:, start:end], axis=1)  # shape (N_times,)

    #     # Étendre à tous les canaux (broadcasting)
    #     self.bkg = np.tile(background_per_time[:, np.newaxis], (1, self.counts.shape[1]))  # shape (N_times, N_channels)

    #     print(f"✅ Background calculated from channels {start} to {end} for all times.")

    # def get_data_bkg(self):
    #     """subtract the background from the signal, channel by channel and time by time."""
    #     self.data_bkg = self.counts - self.bkg
    #     self.data_bkg[self.data_bkg < 0] = 0  # change negative values to 0
    #     print("✅ Data - Background calculated.")
    #     print("Data - Background shape: ", self.data_bkg.shape)
    #     print(f"Background shape: {self.bkg.shape}")
    #     print(f"Data original: {self.counts.shape}")
    #     print("_______________________________")
    #     print("Data - Background shape: ", self.data_bkg)
    #     print(f"Background shape: {self.bkg}")
    #     print(f"Data original: {self.counts}")

    def select_background_range(self):
        """Permet de sélectionner une plage de temps pour estimer le fond (background)."""

        # Moyenne des counts sur tous les canaux pour chaque temps
        mean_counts_per_time = np.mean(self.counts, axis=1)
        times = self.times
        exposure = np.mean(self.time_del)

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(times, mean_counts_per_time/exposure, color='red', label='Mean Counts vs Time')

        ax.set_title("Select Time Interval for Background")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Mean Counts")
        ax.grid(True)
        ax.legend()

        def onselect(xmin, xmax):
            self.background_time_start = xmin
            self.background_time_end = xmax

            print(f"✔️ Background time interval selected: {xmin:.2f} to {xmax:.2f} seconds")
            plt.close(fig)

        span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                            props=dict(alpha=0.5, facecolor='gray'))

        plt.show()


    def get_bkg(self):
        """Calcule le fond pour chaque canal à partir d’un intervalle de temps sélectionné."""

        if not hasattr(self, 'background_time_start') or not hasattr(self, 'background_time_end'):
            print("❌ Aucun intervalle de temps sélectionné pour le fond.")
            return

        # Trouver les indices de temps correspondants
        time_start_idx = np.searchsorted(self.times, self.background_time_start)
        time_end_idx = np.searchsorted(self.times, self.background_time_end)

        # Sélection de l’intervalle
        selected_counts = self.counts[time_start_idx:time_end_idx + 1, :]  # shape (N_t_bg, N_channels)

        # Moyenne temporelle du fond
        self.bkg = np.mean(selected_counts, axis=0)  # shape (N_channels,)
        self.bkg = np.tile(self.bkg, (self.counts.shape[0], 1))  # shape (N_times, N_channels)

        print(f"✅ Background computed from time index {time_start_idx} to {time_end_idx} ({self.background_time_start:.2f}–{self.background_time_end:.2f} s)")

    def get_data_bkg(self):
        """Soustrait le fond au signal sur tous les temps et canaux."""
        self.data_bkg = self.counts - self.bkg
        self.data_bkg[self.data_bkg < 0] = 0
        print("✅ Data - Background calculated.")



    def _selective_fit(self):
        """Selection depending on Plot Units and Function Model
          Predefine Input Data in energy_data and spec_data
          We equate three components to rate_data, counts_data, flux_data. The value of energy_data is the same for all cases
          energy_data - independent variable, nominally energy in keV
          spec_data - Plot Unit"""
        
        selection = self.lbox.curselection()
        if not selection:
            messagebox.showwarning("No Model Selected", "Please select a fit model before clicking 'Do Fit'.")
            return
        
        # load chosen file in Select Input section
        fname = Fitting.fname
        rname = Fitting.rname
        if fname is None and rname is None:  # if file not choosen, print
            print('Please, choose input file')

        else:
            # --- Read Counts data ---
            # if self.show_db_var.get() and self.background_start is not None and self.background_end is not None:
            #     # Crée un masque pour le background
            #     bkg_mask = (self.e_low_det >= self.background_start) & (self.e_high_det <= self.background_end)
            #     background = np.mean(self.counts[:, bkg_mask], axis=0)

            #     # Attention à la forme de background vs counts
            #     background = np.mean(background)  # moyenne plate si nécessaire
            #     new_data = self.counts - background
            #     new_data[new_data < 0] = 0

            #     used_data = new_data
            #     absolute_name = "Data - Background"
            # else:
            #     used_data = self.counts
            #     absolute_name = "Data"

            if self.show_db_var.get():
                self.get_bkg()
                self.get_data_bkg()
                used_data = self.data_bkg
                absolute_name = "Data - Background"
            else:
                used_data = self.counts
                absolute_name = "Data"




            counts_all = np.mean(used_data, axis=0)
            counts_err_all = np.mean(self.counts_err, axis=0)
            exposure = np.mean(self.time_del)
            e_low_det_all = self.e_low_det
            e_high_det_all = self.e_high_det

            # Read SRM data
            e_low_true = self.e_low_true
            e_high_true = self.e_high_true
            matrix = self.matrix

            # --- use mask channel to avoid shape problem ---

            usable_channels = np.arange(min(matrix.shape[1], len(e_low_det_all)))

            counts = counts_all[usable_channels]
            counts_err = counts_err_all[usable_channels]
            e_low_det = e_low_det_all[usable_channels]
            e_high_det = e_high_det_all[usable_channels]

           
            # remove NaN values and negative values from counts and counts_err
            valid = (counts_err > 0) & np.isfinite(counts_err) & np.isfinite(counts)

            counts = counts[valid]
            counts_err = counts_err[valid]
            x_fake = np.zeros_like(counts)  # X fake for plotting must be same shape as counts

            matrix = matrix[:, valid]       # matrix is 2D array, so we need to remove the same channels from it
            e_low_det = e_low_det[valid]
            e_high_det = e_high_det[valid]

            edges_det = np.append(e_low_det, e_high_det[-1])
            dE_det = np.diff(edges_det)

            Edges_photon = np.append(e_low_true, e_high_true[-1])

            # --- Fitting avec LevMarLSQFitter ---

            # === fitting range ===
            fit_Emin = self.energy_min_var.get()  # keV
            fit_Emax = self.energy_max_var.get()  # keV
            # fit_Emin = 10.0  # keV
            # fit_Emax = 20.0  # keV

            # mask for fitting range
            fit_mask = (edges_det[:-1] >= fit_Emin) & (edges_det[1:] <= fit_Emax)

            x_fit = x_fake[fit_mask]
            counts_fit = counts[fit_mask]
            counts_err_fit = counts_err[fit_mask]
            matrix_fit = matrix[:, fit_mask]

            # Counts
            mean_counts = counts
            mean_counts_err = counts_err

            # Rate
            rate = mean_counts / exposure
            rate_err = mean_counts_err / exposure

            # Flux (photons / s / cm² / keV)
            flux = rate / (self.area * dE_det)
            flux_err = rate_err / (self.area * dE_det)

            # Unit selection
            unit = self.var.get()
            
            if unit == 'Rate':
                y_data = rate
                y_err = rate_err
                y_label = "Rate [counts / (s keV)]"
            elif unit == 'Counts':
                y_data = mean_counts
                y_err = mean_counts_err
                y_label = "Counts (Counts)'"
            elif unit == 'Flux':
                y_data = flux
                y_err = flux_err
                y_label = "Flux (Counts/s/cm²/keV)"
            else:
                raise ValueError("Choose unit = 'rate', 'counts' ou 'flux'")

            plt.figure()
            plt.step(edges_det[:-1], y_data, where='mid', label=f'{absolute_name} ({unit})', color='red')
            # plt.axvspan(self.e_low_det[self.background_channel_start], self.e_high_det[self.background_channel_end - 1], 
            #             color='gray', alpha=0.3, label="Background Interval")

            if self.lbox.curselection()[0] == 0:
                self.fit_model = 'Power Law'
                # Create model
                model_fit = Fitting.ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix_fit, exposure)

                # Fitting
                fitter = LevMarLSQFitter()
                fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                                    weights=1.0 / (counts_err_fit / exposure))

                # Parameters
                amplitude = fitted_model.amplitude.value
                alpha = fitted_model.alpha.value
                x_0 = fitted_model.x_0

                # Modèle complet pour affichage sur tout le domaine
                model_display = Fitting.ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix, exposure)
                model_display.amplitude = fitted_model.amplitude
                model_display.alpha = fitted_model.alpha

                # Calcul du modèle simulé complet
                rate_modeled_full = model_display(x_fake)

                if unit == 'Rate':
                    model_y = (rate_modeled_full / dE_det)
                elif unit == 'Counts':
                    model_y = (rate_modeled_full * exposure)
                elif unit == 'Flux':
                    model_y = (rate_modeled_full / (self.area * dE_det))
                else:
                    raise ValueError("Unit most be choose")
                
                plt.step(edges_det[:-1], model_y, where='mid',
                    label='Fitted Model', color='blue')
                
                if self.show_params_var.get():
                    # show model parameters on the plot
                    plt.text(0.05, 0.4,
                        f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
                        transform=plt.gca().transAxes,
                        fontsize=10,
                        verticalalignment='top',
                        bbox=dict(facecolor='white', alpha=0.7))
                
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Channel Energy (keV)")
                plt.ylabel(y_label)
                plt.title(f"Fitting on [{fit_Emin}, {fit_Emax}] keV")
                if self.grid_var.get():
                    plt.grid(True, which="both", ls="--", alpha=0.5)
                else:
                    plt.grid(False)
                plt.legend()
                plt.tight_layout()
                    
                if self.show_photon_var.get():
                    # --- Photon ---
                    model_func = lambda E: amplitude * (E / x_0)**(-alpha)
                    flux_photons = np.array([
                        Fitting.integrate_flux(e1, e2, model_func)
                        for e1, e2 in zip(e_low_true, e_high_true)
                    ])

                    plt.figure()
                    plt.step(Edges_photon[:-1], flux_photons, where='mid',
                            label='Photon model', color='green')
                    xscale = getattr(self, "photon_xscale", "log")
                    yscale = getattr(self, "photon_yscale", "log")
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    plt.xlabel("Energy (keV)")
                    plt.ylabel("Photon flux [photons / (s cm² keV)]")
                    plt.title("Photon Flux Model")
                    if self.grid_var.get():
                        plt.grid(True, which="both", ls="--", alpha=0.5)
                    plt.legend()
                    if self.show_params_var.get():
                        plt.text(0.05, 0.4,
                            f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7))
                    plt.tight_layout()
                
            elif self.lbox.curselection()[0] == 1:
                self.fit_model = 'Broken Power Law'

                model_fit = Fitting.ForwardFoldedBrokenPowerLaw(e_low_true, e_high_true, matrix_fit, exposure)
                model_fit.alpha_1.bounds = (0.1, 10)
                model_fit.alpha_2.bounds = (0.1, 10)
                model_fit.E_break.bounds = (1, 100)
                model_fit.amplitude.bounds = (1e-5, 1e3)


                # Fitting
                print("x_fit:", x_fit)
                print("counts_fit / exposure:", counts_fit / exposure)
                print("weights:", 1.0 / (counts_err_fit / exposure))

    
                fitter = LevMarLSQFitter()
                fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                                    weights=1.0 / (counts_err_fit / exposure))

                # Parameters
                amplitude = fitted_model.amplitude.value
                E_break = fitted_model.E_break.value
                alpha_1 = fitted_model.alpha_1.value  
                alpha_2 = fitted_model.alpha_2.value

                # Modèle complet pour affichage sur tout le domaine
                model_display = Fitting.ForwardFoldedBrokenPowerLaw(e_low_true, e_high_true, matrix, exposure)
                model_display.amplitude = fitted_model.amplitude
                model_display.alpha_1 = fitted_model.alpha_1
                model_display.alpha_2 = fitted_model.alpha_2
                model_display.E_break = fitted_model.E_break

                # Calcul du modèle simulé complet
                rate_modeled_full = model_display(x_fake)

                if unit == 'Rate':
                    model_y = (rate_modeled_full / dE_det)
                elif unit == 'Counts':
                    model_y = (rate_modeled_full * exposure)
                elif unit == 'Flux':
                    model_y = (rate_modeled_full / (self.area * dE_det))
                else:
                    raise ValueError("Unit most be choose")
                
                plt.step(edges_det[:-1], model_y, where='mid',
                    label='Fitted Model', color='blue')
                
                if self.show_params_var.get():
                    # show model parameters on the plot
                    plt.text(0.05, 0.4,
                        f"Broken Power Law:\n amplitude = {amplitude:.2e}\n E_break = {E_break:.2f} \n Alpha_1 = {alpha_1:.2e}\n Alpha_2 = {alpha_2:.2f} \n",
                        transform=plt.gca().transAxes,
                        fontsize=10,
                        verticalalignment='top',
                        bbox=dict(facecolor='white', alpha=0.7))
                
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Channel Energy (keV)")
                plt.ylabel(y_label)
                plt.title(f"Fitting on [{fit_Emin}, {fit_Emax}] keV")
                if self.grid_var.get():
                    plt.grid(True, which="both", ls="--", alpha=0.5)
                else:
                    plt.grid(False)
                plt.legend()
                plt.tight_layout()
                    
                if self.show_photon_var.get():
                    # --- Photon ---
                    model_func = lambda E: amplitude * np.where(E < E_break, (E / E_break)**(-alpha_1), (E / E_break)**(-alpha_2))
                    flux_photons = np.array([
                        Fitting.integrate_flux(e1, e2, model_func)
                        for e1, e2 in zip(e_low_true, e_high_true)
                    ])

                    plt.figure()
                    plt.step(Edges_photon[:-1], flux_photons, where='mid',
                            label='Photon model', color='green')
                    xscale = getattr(self, "photon_xscale", "log")
                    yscale = getattr(self, "photon_yscale", "log")
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    plt.xlabel("Energy (keV)")
                    plt.ylabel("Photon flux [photons / (s cm² keV)]")
                    plt.title("Photon Flux Model")
                    if self.grid_var.get():
                        plt.grid(True, which="both", ls="--", alpha=0.5)
                    plt.legend()
                    if self.show_params_var.get():
                        plt.text(0.05, 0.4,
                            f"Broken Power Law:\n amplitude = {amplitude:.2e}\n E_break = {E_break:.2f} \n Alpha_1 = {alpha_1:.2e}\n Alpha_2 = {alpha_2:.2f} \n",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7))
                    plt.tight_layout()

            elif self.lbox.curselection()[0] == 2:
                self.fit_model = 'Exponential Power Law'

                model_fit = Fitting.ForwardFoldedExpPowerLaw(e_low_true, e_high_true, matrix_fit, exposure)
                # Fixe des bornes pour éviter l’explosion
                model_fit.p0.bounds = (1e-3, 1e5)
                model_fit.p1.bounds = (-5, 5)
                model_fit.p2.bounds = (1e-2, 100)
                model_fit.e3.bounds = (-10, 10)
                model_fit.e4.bounds = (0.1, 100)

                # Fitting
                fitter = LevMarLSQFitter()
                fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                                    weights=1.0 / (counts_err_fit / exposure))

                # Parameters
                p0 = fitted_model.p0.value
                p1 = fitted_model.p1.value
                p2 = fitted_model.p2.value
                e3 = fitted_model.e3.value
                e4 = fitted_model.e4.value

                # Modèle complet pour affichage sur tout le domaine
                model_display = Fitting.ForwardFoldedExpPowerLaw(e_low_true, e_high_true, matrix, exposure)
                model_display.p0 = fitted_model.p0
                model_display.p1 = fitted_model.p1 
                model_display.p2 = fitted_model.p2
                model_display.e3 = fitted_model.e3
                model_display.e4 = fitted_model.e4

                # Calcul du modèle simulé complet
                rate_modeled_full = model_display(x_fake)

                if unit == 'Rate':
                    model_y = (rate_modeled_full / dE_det)
                elif unit == 'Counts':
                    model_y = (rate_modeled_full * exposure)
                elif unit == 'Flux':
                    model_y = (rate_modeled_full / (self.area * dE_det))
                else:
                    raise ValueError("Unit most be choose")
                
                plt.step(edges_det[:-1], model_y, where='mid',
                    label='Fitted Model', color='blue')
                
                if self.show_params_var.get():
                    # show model parameters on the plot
                    plt.text(0.05, 0.4,
                        f"Exponential Power Law:\n p0 = {p0:.2e}\n p1 = {p1:.2f} \n p2 = {p2:.2f} \n e3 = {e3:.2f} \n e4 = {e4:.2f} \n",
                        transform=plt.gca().transAxes,
                        fontsize=10,
                        verticalalignment='top',
                        bbox=dict(facecolor='white', alpha=0.7))
                
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Channel Energy (keV)")
                plt.ylabel(y_label)
                plt.title(f"Fitting on [{fit_Emin}, {fit_Emax}] keV")
                if self.grid_var.get():
                    plt.grid(True, which="both", ls="--", alpha=0.5)
                else:
                    plt.grid(False)
                plt.legend()
                plt.tight_layout()
                    
                if self.show_photon_var.get():
                    # --- Photon ---
                    model_func = lambda E: (p0 * (E / p2)**p1) * np.exp(e3 - E / e4)
                    flux_photons = np.array([
                        Fitting.integrate_flux(e1, e2, model_func)
                        for e1, e2 in zip(e_low_true, e_high_true)
                    ])

                    plt.figure()
                    plt.step(Edges_photon[:-1], flux_photons, where='mid',
                            label='Photon model', color='green')
                    xscale = getattr(self, "photon_xscale", "log")
                    yscale = getattr(self, "photon_yscale", "log")
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    plt.xlabel("Energy (keV)")
                    plt.ylabel("Photon flux [photons / (s cm² keV)]")
                    plt.title("Photon Flux Model")
                    if self.grid_var.get():
                        plt.grid(True, which="both", ls="--", alpha=0.5)
                    plt.legend()
                    if self.show_params_var.get():
                        plt.text(0.05, 0.4,
                            f"Exponential Power Law:\n p0 = {p0:.2e}\n p1 = {p1:.2f} \n p2 = {p2:.2f} \n e3 = {e3:.2f} \n e4 = {e4:.2f} \n",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7))
                    plt.tight_layout()

            elif self.lbox.curselection()[0] == 3:
                self.fit_model = 'VTH'

                model_fit = Fitting.ForwardFoldedVTH(e_low_true, e_high_true, matrix_fit, exposure)

                # Fitting
                fitter = LevMarLSQFitter()
                fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                                    weights=1.0 / (counts_err_fit / exposure))

                T = fitted_model.T.value
                EM = fitted_model.EM.value

                model_display = Fitting.ForwardFoldedVTH(e_low_true, e_high_true, matrix, exposure)
                model_display.T = fitted_model.T
                model_display.EM = fitted_model.EM

                rate_modeled_full = model_display(x_fake)

                if unit == 'Rate':
                    model_y = rate_modeled_full / dE_det
                elif unit == 'Counts':
                    model_y = rate_modeled_full * exposure
                elif unit == 'Flux':
                    model_y = rate_modeled_full / (self.area * dE_det)

                plt.step(edges_det[:-1], model_y, where='mid', label='Fitted VTH Model', color='blue')

                if self.show_params_var.get():
                    plt.text(0.05, 0.4,
                            f"V_TH Model:\n T = {T:.2f} keV\n EM = {EM:.2e} cm⁻³",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7))
                    
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Channel Energy (keV)")
                plt.ylabel(y_label)
                plt.title(f"Fitting on [{fit_Emin}, {fit_Emax}] keV")
                if self.grid_var.get():
                    plt.grid(True, which="both", ls="--", alpha=0.5)
                else:
                    plt.grid(False)
                plt.legend()
                plt.tight_layout()

                # === AFFICHAGE PHOTONIQUE ===
                if self.show_photon_var.get():
                    # Constantes physiques
                    gff = 1.2
                    A = 1.07e-42 * gff
                    k_B_keV = 8.617333262e-8  # erg/K in keV

                    T_keV = (T * k_B_keV) / 1.60218e-9  # conversion K -> keV

                    model_func = lambda E: (A * EM) / (E * np.sqrt(T)) * np.exp(-E / T_keV)
                    flux_photons = np.array([
                        Fitting.integrate_flux(e1, e2, model_func)
                        for e1, e2 in zip(e_low_true, e_high_true)
                    ])

                    plt.figure()
                    plt.step(Edges_photon[:-1], flux_photons, where='mid',
                            label='Photon model', color='green')
                    xscale = getattr(self, "photon_xscale", "log")
                    yscale = getattr(self, "photon_yscale", "log")
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    plt.xlabel("Energy (keV)")
                    plt.ylabel("Photon flux [photons / (s cm² keV)]")
                    plt.title("Photon Flux Model")
                    if self.grid_var.get():
                        plt.grid(True, which="both", ls="--", alpha=0.5)
                    plt.legend()
                    if self.show_params_var.get():
                        plt.text(0.05, 0.4,
                                f"V_TH Model:\n T = {T:.2e} K\n EM = {EM:.2e} cm⁻³",
                                transform=plt.gca().transAxes,
                                fontsize=10,
                                verticalalignment='top',
                                bbox=dict(facecolor='white', alpha=0.7))
                    plt.tight_layout()

            elif self.lbox.curselection()[0] == 4:
                self.fit_model = 'V_TH + Power Law'

                
                model_fit = Fitting.ForwardFoldedVTHPlusPowerLaw(e_low_true, e_high_true, matrix_fit, exposure)

                fitter = LevMarLSQFitter()
                fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                                    weights=1.0 / (counts_err_fit / exposure))

                # Paramètres du modèle
                EM = fitted_model.EM.value
                T = fitted_model.T.value
                amplitude = fitted_model.amplitude.value
                alpha = fitted_model.alpha.value

                # Création du modèle à afficher sur tout le domaine
                model_display = Fitting.ForwardFoldedVTHPlusPowerLaw(e_low_true, e_high_true, matrix, exposure)
                model_display.EM = fitted_model.EM
                model_display.T = fitted_model.T
                model_display.amplitude = fitted_model.amplitude
                model_display.alpha = fitted_model.alpha

                rate_modeled_full = model_display(x_fake)

                if unit == 'Rate':
                    model_y = rate_modeled_full / dE_det
                elif unit == 'Counts':
                    model_y = rate_modeled_full * exposure
                elif unit == 'Flux':
                    model_y = rate_modeled_full / (self.area * dE_det)

                plt.step(edges_det[:-1], model_y, where='mid', label='Fitted VTH Model', color='blue')

                if self.show_params_var.get():
                    if self.show_params_var.get():
                        plt.text(
                            0.05, 0.4,
                            f"V_TH + Power Law:\n"
                            f"T  = {T:.2e} keV\n"
                            f"EM = {EM:.2e} cm⁻³\n"
                            f"amplitude = {amplitude:.2e}\n"
                            f"alpha     = {alpha:.2f}",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7)
                        )
                    
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("Channel Energy (keV)")
                plt.ylabel(y_label)
                plt.title(f"Fitting on [{fit_Emin}, {fit_Emax}] keV")
                if self.grid_var.get():
                    plt.grid(True, which="both", ls="--", alpha=0.5)
                else:
                    plt.grid(False)
                plt.legend()
                plt.tight_layout()

                # === AFFICHAGE PHOTONIQUE ===
                if self.show_photon_var.get():
                    model_func = lambda E: (
                        (1.07e-42 * 1.2 * EM) / (E * np.sqrt(max(1e-3, T))) * np.exp(-E / T) +
                        amplitude * (E / 100.0) ** (-alpha)
                    )

                    flux_photons = np.array([
                        Fitting.integrate_flux(e1, e2, model_func)
                        for e1, e2 in zip(e_low_true, e_high_true)
                    ])


                    plt.figure()
                    plt.step(Edges_photon[:-1], flux_photons, where='mid',
                            label='Photon model', color='green')
                    xscale = getattr(self, "photon_xscale", "log")
                    yscale = getattr(self, "photon_yscale", "log")
                    plt.xscale(xscale)
                    plt.yscale(yscale)
                    plt.xlabel("Energy (keV)")
                    plt.ylabel("Photon flux [photons / (s cm² keV)]")
                    plt.title("Photon Flux Model")
                    if self.grid_var.get():
                        plt.grid(True, which="both", ls="--", alpha=0.5)
                    plt.legend()
                    if self.show_params_var.get():
                        plt.text(
                            0.05, 0.4,
                            f"V_TH + Power Law:\n"
                            f"T  = {T:.2e} keV\n"
                            f"EM = {EM:.2e} cm⁻³\n"
                            f"amplitude = {amplitude:.2e}\n"
                            f"alpha     = {alpha:.2f}",
                            transform=plt.gca().transAxes,
                            fontsize=10,
                            verticalalignment='top',
                            bbox=dict(facecolor='white', alpha=0.7)
                        )

                    plt.tight_layout()

            plt.show()