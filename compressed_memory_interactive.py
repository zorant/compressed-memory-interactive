''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve sliders.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders
in your browser.
'''
import numpy as np
from scipy.misc import factorial
import pylab as pl

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput
from bokeh.plotting import figure

def compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_tau, dtau):

    N = buff_len+2*k
    
    #Create power-law growing Taustarlist and corresponding s
    alpha = (Taustar_max/Taustar_min)**(1./buff_len)-1
    pow_vec = np.arange(-k,buff_len +k) #-1
    Taustarlist = Taustar_min * (1+alpha)**pow_vec
    
    s = k/Taustarlist
    
    #Create DerivMatrix
    DerivMatrix = np.zeros((N,N))
    for i in range(1,N-1):
        DerivMatrix[i, i-1] = -(s[i+1]-s[i])/(s[i]-s[i-1])/(s[i+1] - s[i-1])
        DerivMatrix[i, i] = ((s[i+1]-s[i])/(s[i]- s[i-1])-(s[i]-s[i-1])/(s[i+1]-s[i]))/(s[i+1] - s[i-1])
        DerivMatrix[i, i+1] = (s[i]-s[i-1])/(s[i+1]-s[i])/(s[i+1] - s[i-1])
      
    #Create Time vector and input signal 
    time_vec = np.arange(0,length_tau,dtau)
    f_tau = np.concatenate((np.zeros(3850), np.ones(100), np.zeros(4000), np.ones(100), np.zeros(1950)))
    t = np.zeros((N,int(length_tau/dtau)))
    T = np.zeros((N,int(length_tau/dtau)))
    
    #Run the main loop 
    tau_index = 0
    for tau in pl.frange(dtau,length_tau-dtau,dtau):
        tau_index = tau_index+1
        t[:,tau_index] = t[:,tau_index-1]+((-s.T*t[:,tau_index-1]+f_tau[tau_index])*dtau)
        t_diff = np.dot(np.linalg.matrix_power(DerivMatrix, k), t[:,tau_index])
        L1 = (-1)**k*s**(k+1)
        L2 = t_diff/factorial(k)
        T[:,tau_index] = L1.T*L2.T
        
    return t, T, time_vec
    
    
#Initialize parameters
buff_len = 5
k = 4
Taustar_min = 1
Taustar_max = 10 
length_tau = 10
dtau = 0.001

t, T, time_vec = compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_tau, dtau)
    
# Set up plot
plot = figure(plot_height=400, plot_width=400, title="small big t test",
              tools="crosshair,pan,reset,save,wheel_zoom",
               x_range=[0, length_tau - dtau])

for i in range(k,k+buff_len-1,1):
   # Set up data
   x = time_vec
   y = T[i, :]
   source = ColumnDataSource(data=dict(x=x, y=y))
   plot.line('x', 'y', source=source, line_width=2, line_alpha=0.6)
   #plot.multi_line('x', 'y', source=source, line_width=2, line_alpha=0.6)

#x = time_vec
#y = T[20, :]

#source = ColumnDataSource(data=dict(x=x, y=y))
#plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
    
# Set up widgets
text = TextInput(title="title", value='big T')
k_slider = Slider(title="k", value=k, start=1, end=10, step=1)
Taustar_min_slider = Slider(title="tstr min", value=Taustar_min, start=1, end=4, step=1)
Taustar_max_slider = Slider(title="tstr max", value=Taustar_max, start=6, end=15, step=1)
buff_len_slider = Slider(title="number of tstr nodes", value=buff_len , start=3, end=100, step=1)

# Set up callbacks
def update_title(attrname, old, new):
    plot.title.text = text.value

text.on_change('value', update_title)

def update_data(attrname, old, new):

    # Get the current slider values
    k = k_slider.value
    Taustar_min = Taustar_min_slider.value
    Taustar_max = Taustar_max_slider.value
    buff_len = buff_len_slider.value

    t, T, time_vec = compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_tau, dtau)
    
    for i in range(k,k+buff_len-1,1):
        # Set up data
        x = time_vec
        y = T[i, :]
        source.data = dict(x=x, y=y)
        #source = ColumnDataSource(data=dict(x=x, y=y))
        plot.line('x', 'y', source=source, line_width=2, line_alpha=0.6)
    # Set up data
    #x = time_vec
    #y = T[20, :]
    #source.data = dict(x=x, y=y)

for w in [k_slider, Taustar_min_slider, Taustar_max_slider, buff_len_slider]:
    w.on_change('value', update_data)


# Set up layouts and add to document
inputs = widgetbox(text, k_slider, Taustar_min_slider, Taustar_max_slider, buff_len_slider)

curdoc().add_root(row(inputs, plot, width=800))
curdoc().title = "Sliders"





