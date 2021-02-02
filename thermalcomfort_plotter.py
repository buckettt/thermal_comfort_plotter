"""
===============
Embedding Comfort models in Tk Interface
===============

"""

import tkinter
from tkinter import Text, END, LEFT, Checkbutton
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.figure import Figure
import numpy as np
from mf_modules import comfortmodels
from mf_modules.ui import MFButton, MFLabelFrame, MFLabel, MFLabelBlack, MFHeader, MFOptionMenu, MFTk

VARS = ['DBT', 'MRT', 'RH', 'VEL', 'MET', 'CLO']
METRICS = ['PMV', 'PPD', 'SET']

class App_Window(MFTk):
    '''
    Tkinter app modules for thermal comfort plotter.

    '''

    def __init__(self,parent):
        MFTk.__init__(self,parent)
        self.title("Comfort Plotter")
        self.parent = parent
        self.initialise()

    def initcase(self):
        self.xrange = np.linspace(20, 30, 20)
        self.yrange = np.linspace(20, 30, 20)

        self.xaxis = tkinter.StringVar(self)
        self.yaxis = tkinter.StringVar(self)
        self.metric = tkinter.StringVar(self)
        

        self.xaxis.set("DBT")
        self.yaxis.set("MRT")
        self.metric.set("PMV")

        #Defaults###
        self.case={'DBT': 24,
              'MRT': 24,
              'RH': 90,
              'VEL': 0.1,
              'MET': 1.1,
              'CLO': 0.5,}

    def initialise(self):
        self.initcase()
        
        ### Banner ###
        MFHeader(self, text="Comfort Plotter").pack(fill=tkinter.BOTH)
        MFLabelBlack(self, text="An app to plot Comfort metrics for variations in different parameters. \n Change the values and hit plot to refresh the view.").pack(fill=tkinter.BOTH)

        ### Figure###
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        ##matplotlib toolbar###
        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()

        ###Buttons and controls###
        groupplot = MFLabelFrame(self, text="Plotting", padx=5, pady=5)
        groupplot.pack()

        button = MFButton(master=groupplot, text="Quit", command=self._quit)
        button2 = MFButton(master=groupplot, text="Plot", command=self.redraw)

        button.pack(side=tkinter.LEFT)
        button2.pack(side=tkinter.LEFT)
        
        popupMenu = MFOptionMenu(groupplot, self.metric, *METRICS)
        popupMenu.pack(side=tkinter.LEFT)

                ### Default/Other Parameters###
        groupdef = MFLabelFrame(self, text="Defaults", padx=5, pady=5)
        groupdef.pack()

        MFLabel(groupdef, text="Dry Bulb Temp:").pack(side=LEFT)
        self.dbtText = Text(groupdef, height=1, width=4)
        self.dbtText.pack(side=tkinter.LEFT)
        self.dbtText.insert(END, "24")
        MFLabel(groupdef, text="Radiant Temp:").pack(side=LEFT)
        self.mrtText = Text(groupdef, height=1, width=4)
        self.mrtText.pack(side=tkinter.LEFT)
        self.mrtText.insert(END, "24")
        MFLabel(groupdef, text="RH:").pack(side=LEFT)
        self.rhText = Text(groupdef, height=1, width=4)
        self.rhText.pack(side=tkinter.LEFT)
        self.rhText.insert(END, "50")
        MFLabel(groupdef, text="Air Speed:").pack(side=LEFT)
        self.velText = Text(groupdef, height=1, width=4)
        self.velText.pack(side=tkinter.LEFT)
        self.velText.insert(END, "0.2")
        MFLabel(groupdef, text="Met:").pack(side=LEFT)
        self.metText = Text(groupdef, height=1, width=4)
        self.metText.pack(side=tkinter.LEFT)
        self.metText.insert(END, "1.1")
        MFLabel(groupdef, text="Clothing:").pack(side=LEFT)
        self.cloText = Text(groupdef, height=1, width=4)
        self.cloText.pack(side=tkinter.LEFT)
        self.cloText.insert(END, "0.5")

        ### x-axis controlss###
        groupx = MFLabelFrame(self, text="x-axis", padx=5, pady=5)
        groupx.pack(side=tkinter.LEFT)

        popupMenu = MFOptionMenu(groupx, self.xaxis, *VARS)
        popupMenu.pack(side=tkinter.LEFT)

        MFLabel(groupx, text="Min:").pack(side=LEFT)
        self.xminText = Text(groupx, height=1, width=4)
        self.xminText.pack(side=tkinter.LEFT)
        self.xminText.insert(END, "20")

        MFLabel(groupx, text="Max:").pack(side=LEFT)
        self.xmaxText = Text(groupx, height=1, width=4)
        self.xmaxText.pack(side=tkinter.LEFT)
        self.xmaxText.insert(END, "30")

        ### y-axis controlss###
        groupy = MFLabelFrame(self, text="y-axis", padx=5, pady=5)
        groupy.pack(side=tkinter.LEFT)

        popupMenu = MFOptionMenu(groupy, self.yaxis, *VARS)
        popupMenu.pack(side=tkinter.LEFT)

        MFLabel(groupy, text="Min:").pack(side=LEFT)
        self.yminText = Text(groupy, height=1, width=4)
        self.yminText.pack(side=tkinter.LEFT)
        self.yminText.insert(END, "20")

        MFLabel(groupy, text="Max:").pack(side=LEFT)
        self.ymaxText = Text(groupy, height=1, width=4)
        self.ymaxText.pack(side=tkinter.LEFT)
        self.ymaxText.insert(END, "30")

        self.update()
        self.redraw()
        


    def _quit(self):
        self.quit()     # stops mainloop
        self.destroy()  # this is necessary on Windows to prevent

    def redraw(self):
        
        metric = self.metric.get()
        
        for i, m in enumerate(METRICS):
            if metric==m:
                metricid = i

        self.case={
              'DBT': float(self.dbtText.get("1.0", END)),
              'MRT': float(self.mrtText.get("1.0", END)),
              'RH': float(self.rhText.get("1.0", END)),
              'VEL': float(self.velText.get("1.0", END)),
              'MET': float(self.metText.get("1.0", END)),
              'CLO': float(self.cloText.get("1.0", END)),
              }

        print(self.case)

        self.xrange = np.linspace(float(self.xminText.get("1.0", END)), float(self.xmaxText.get("1.0", END)), 20)
        self.yrange = np.linspace(float(self.yminText.get("1.0", END)), float(self.ymaxText.get("1.0", END)), 20)

        vmin=[-3,0,10]
        vmax=[3,100,44]

        levels = [[-3,-2,-1,0,1,2,3],
                  [0,10,20,30,40,50,60,70,80,90,100],
                  [10,14.5,17.5,22.2,25.6,30,34.5,37.5,44]]


        self.ax.clear()
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        C = self.plot_data(self.case, metricid)
        PLOT = self.ax.pcolor(self.xrange, self.yrange, C,  cmap='RdBu_r', vmin=vmin[metricid], vmax=vmax[metricid])
        CS = self.ax.contour(self.xrange, self.yrange, C, cmap='gnuplot', levels=levels[metricid])

        fmt = {}
        strs = [['Cold (-3)', 'Cool (-2)', 'Slightly Cool (-1)', 'Neutral (0)', 'Slightly Warm (1)', 'Warm (2)', 'Hot (3)'],
                ['0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'],
                ['Cold','Cool','Slightly Cool','Comfortable','Slightly Warm','Warm','Hot','Very Hot']]
        
        
        for l, s in zip(CS.levels, strs[metricid]):
            fmt[l] = s

        self.ax.clabel(CS, inline=1, fontsize=10, fmt=fmt)
        self.ax.set_xlabel(self.xaxis.get())
        self.ax.set_ylabel(self.yaxis.get())

        self.fig.colorbar(PLOT)

        self.canvas.draw()



        print("redraw")


    def plot_data(self, case, metricid):
        C=np.zeros((len(self.xrange),len(self.yrange)))

        xaxis = self.xaxis.get()
        yaxis = self.yaxis.get()

        print(xaxis, yaxis)
        for i, x in enumerate(self.xrange):
            case[xaxis] = self.xrange[i]
            print(i,"/20")
            for j, y in enumerate(self.yrange):
                #print(case)
                case[yaxis] = self.yrange[j]
                #comfortmodels.pmvElevatedAirspeed(25.4,20,0.1,90,1.1,0.5,0.0)
                if metricid == 0:
                    C[j,i] = comfortmodels.pmvElevatedAirspeed(case['DBT'], case['MRT'], case['VEL'], case['RH'], case['MET'], case['CLO'], 0.0)['pmv']
                elif metricid == 1:
                    C[j,i] = comfortmodels.pmvElevatedAirspeed(case['DBT'], case['MRT'], case['VEL'], case['RH'], case['MET'], case['CLO'], 0.0)['ppd']
                elif metricid == 2:
                    C[j,i] = comfortmodels.pierceSET(case['DBT'], case['MRT'], case['VEL'], case['RH'], case['MET'], case['CLO'], 0.0)
                
        return C


if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()
