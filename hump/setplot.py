
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
import os
from numpy import zeros
import matplotlib.pyplot as plt


add_qref = True   # add plot of reference solution?
if not os.path.isdir('_output_qref'):
    print "*** Did not find directory _output_qref"
    print "*** See README for how to create reference solution"
    add_qref = False
    qrefdir = None
else:
    qrefdir = os.path.abspath('_output_qref')


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    plotdata.clearfigures()  # clear any old figures,axes,items data

    

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-8,16] #'auto'
    plotaxes.ylimits =  [0.8, 1.5] #'auto'
    plotaxes.title = 'Density'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    # Reference solution:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.outdir = qrefdir
    plotitem.show = add_qref
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    
    # Plot outline of interface
    def aa(current_data):
      from pylab import linspace,plot,annotate,text
      #gcs = 2.0/200.0
      #x = [-1.3,-1.3,1.3,1.3]
      #y = [-1000000,10000000,1000000,-1000000]
      x = [0.0, 0.0]
      y = [-1000000, 1000000]
      #y[:] = [xx - gcs for xx in y]
      plot(x,y,'k',linewidth=2.0)
      #plot(-8.0, 180000, 'vk', markersize=10) 
      #plot(-2.0, 180000, 'vk', markersize=10) 
      #plot(0.0, 180000, 'vk', markersize=10) 
      #plot(2.0, 180000, 'vk', markersize=10)
      text(2.0,285000,'Water',fontweight='bold',fontsize=20)
      #text(-0.8,285000,'PS',fontweight='bold',fontsize=20)
      text(-8.0,285000,'Air',fontweight='bold',fontsize=20)
      annotate('Gauges: \n1', xy=(-8.0, 250000), xytext=(-8, 270000),
      arrowprops=dict(facecolor='black', shrink=0.05),
      )
      annotate('2', xy=(-2.0, 250000), xytext=(-2.0, 270000),
      arrowprops=dict(facecolor='black', shrink=0.05),
      )
      #annotate('3', xy=(0.0, 250000), xytext=(0, 270000),
            #arrowprops=dict(facecolor='black', shrink=0.05),
            #)
      annotate('4', xy=(2.0, 250000), xytext=(2, 270000),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

    plotaxes.afteraxes = aa
    

    # Figure for Momentum
    plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-8,16] #'auto'
    plotaxes.ylimits = [-100,150] #'auto'
    plotaxes.title = 'Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    # Reference solution:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.outdir = qrefdir
    plotitem.show = add_qref
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    
    # Plot interface
    plotaxes.afteraxes = aa
    

    # Figure for Energy
    plotfigure = plotdata.new_plotfigure(name='Energy', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-8,16] #'auto'
    plotaxes.ylimits = [200000,400000] #'auto'
    plotaxes.title = 'Energy'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 2
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    # Reference solution:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.outdir = qrefdir
    plotitem.show = add_qref
    plotitem.plot_var = 2
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    
    # Plot interface
    plotaxes.afteraxes = aa
    
    # Figure for Pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-8.5,16]#[-8.5,16] #'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.ylimits = [90000,300000] #'auto'
    #plotaxes.ylimits = [101305,101385] #'auto'
    plotaxes.title = 'Pressure'
    
    plotitem = plotaxes.new_plotitem(plot_type='1d')

    def Pressure(current_data):
        q = current_data.q   # solution when this function called
        aux = current_data.aux
        gamma = aux[0,:]
        gamma1 = aux[0,:] - 1.0
        pinf = aux[1,:]
        omega = aux[2,:]
        rho = q[0,:]           # density
        mom = q[1,:]           # momentum
        ene = q[2,:]           # energy
        P = gamma1*(ene - 0.5*mom*mom/rho)/(1.0 - omega*rho)
        P = P - gamma*pinf 
        return P

    plotitem.plot_var = Pressure  # defined above
    plotitem.plotstyle = '-o'
    plotitem.color = 'r'
    
    # Plot interface
    plotaxes.afteraxes = aa
    
    
     # Figure for Pressure Gauge
    #plotfigure = plotdata.new_plotfigure(name='PressureGauge', figno=4)

    ## Set up for axes in this figure:
    #plotaxes = plotfigure.new_plotaxes()
    #plotaxes.xlimits = 'auto'
    #plotaxes.ylimits = [0,4] #'auto'
    #plotaxes.title = 'Pressure Gauge'
    
    #plotitem = plotaxes.new_plotitem(plot_type='1d')

    #Open clawpack simulation output for Pressure Gauge data
    f = open("./_output/a-pgauge.dat",'r')
    data = [line.split() for line in f]
    f.close()
    time = zeros(len(data))
    PressureG = zeros(len(data))
    #Assign the data to the correponsding arrays
    for i in range(len(data)):
      time[i] = float(data[i][0])
      PressureG[i] = float(data[i][1])
      
    print max(PressureG)
    plt.plot(time,PressureG)
    plt.title('Pressure at Gauge 4 (x=2.0)')
    plt.xlabel('Time in seconds')
    plt.ylabel('Pressure in Pa (N/m^2)')
    plt.ylim(80000,300000)
    plt.xlim(0.0,0.053)
    plt.savefig('./_plots/gauge/Pressure.png')
    
    #plotitem.plot_var = PressureG  # defined above
    #plotitem.plotstyle = '-o'
    #plotitem.color = 'r'

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
