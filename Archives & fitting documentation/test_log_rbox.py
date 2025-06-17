from tkinter import *


class LogRbox(Radiobutton):
    def __init__(self, master=None, **kwargs):
        """This class uses the 'Radiobox' funtion included in the tkinter library. It allows display a pre-built line
         with two radioboxes to choose between linear and log scale."""
        # Variables to save text and buttons
        self.text = StringVar()
        self.rb1 = None
        self.rb2 = None

        # Recording variable to store the checked radiobox
        self.rec = StringVar()
        self.rec.set("linear")

        # Variable to return, containing 'linear' or 'log'
        self.scale = StringVar()

        # Initialisation of Radiobutton class

    def log_rbox(self, window, relx, rely, text_label, *args):
        """Displays the radioboxes allowing the user to choose between linear and logarithmic axis.
        Parameters:
            window: current window on which the radio boxes need to be displayed;
            relx: relative horizontal position of the North-West of the whole set on the window;
            rely: relative vertical position of the North-West of the whole set on the window;
            text_label: string variable containing the text displayed at the left of radioboxes.
        Returns scale parameter 'linear' or 'log' as a StringVar."""
        self.text = Label(window, text=text_label, fg='blue')
        self.text.place(relx=relx, rely=rely)

        self.rb1 = Radiobutton(window, text="Linear", variable=self.rec, value="linear",
                               command=self.radio)
        self.rb1.place(relx=relx + 0.2, rely=rely)
        self.rb2 = Radiobutton(window, text="Logarithmic", variable=self.rec, value="log",
                               command=self.radio)
        self.rb2.place(relx=relx + 0.35, rely=rely)

        return self.scale

    def radio(self):
        """Commands radiobox to save the scale in log_axis function.
        Saves in scalex a StringVar() 'linear' or 'log'."""
        self.scale = self.rec.get()
