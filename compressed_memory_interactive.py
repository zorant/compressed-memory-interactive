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
import numpy.matlib
from scipy.misc import factorial
import pylab as pl

from bokeh.io import curdoc
from bokeh.layouts import row, widgetbox, layout
from bokeh.models.widgets import Slider, TextInput
from bokeh.core.properties import field
from bokeh.models import (
    ColumnDataSource, HoverTool, SingleIntervalTicker, Slider, Button, Label,
    CategoricalColorMapper,
)
from bokeh.palettes import Spectral6
from bokeh.plotting import figure

def compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_time, dtime):

    N = buff_len+2*k
    Nt = int(length_time/dtime)
    
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
    time_vec = np.arange(0,length_time,dtime)
    f = np.concatenate((np.zeros(2850), np.ones(100), np.zeros(4500), np.ones(100), np.zeros(2450)))
    t = np.zeros((N,Nt))
    T = np.zeros((N,Nt))
    
    #Run the main loop 
    time_index = 0
    for time in pl.frange(dtime,length_time-dtime,dtime):
        time_index = time_index+1
        t[:,time_index] = t[:,time_index-1]+((-s.T*t[:,time_index-1]+f[time_index])*dtime)
        t_diff = np.dot(np.linalg.matrix_power(DerivMatrix, k), t[:,time_index])
        L1 = (-1)**k*s**(k+1)
        L2 = t_diff/factorial(k)
        T[:,time_index] = L1.T*L2.T
        
    return f, t, T, time_vec, Taustarlist
    
def compute_source(time_limit, T, t, f):
    source_dict = dict(ff=np.matlib.repmat(f[0:int(time_limit/dtime):100],buff_len, 1).tolist(), tt=t[k:-k,0:int(time_limit/dtime):100].tolist(), TT=T[k:-k,0:int(time_limit/dtime):100].tolist(), xx=np.matlib.repmat(np.arange(0,time_limit,dtime*100), buff_len, 1).tolist())
    #import pdb; pdb.set_trace()
    return source_dict
    
def compute_source_vec(time_limit, T, Taustarlist):
    Taustarlist = -Taustarlist
    source_dict = dict(TTaustarlist=Taustarlist[k:-k].tolist(), TT=T[k:-k:,int(time_limit/dtime)].tolist())
    #import pdb; pdb.set_trace()
    return source_dict
    
#Initialize parameters
buff_len = 25
k = 8
Taustar_min = 1
Taustar_max = 10 
length_time = 10
dtime = 0.001

f, t, T, time_vec, Taustarlist = compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_time, dtime)
    
# Set up plot
plotLI = figure(plot_height=150, plot_width=500, title="Leaky integrators",
              tools="", logo = None, toolbar_location = None,
               x_range=[0, length_time - dtime], y_range=[-0.01, 0.1])
               
plotTC = figure(plot_height=150, plot_width=500, title="Time cells",
              tools="", logo = None, toolbar_location = None,
               x_range=[0, length_time - dtime], y_range=[-0.01, 0.2])
plotf = figure(plot_height=150, plot_width=500, title="Input signal",
              tools="", logo = None, toolbar_location = None,
               x_range=[0, length_time - dtime], y_range=[-0.01, 1.1])
plotr = figure(plot_height=150, plot_width=500, title="Memory representation",
              tools="", logo = None, toolbar_location = None,
               x_range=[-Taustar_max, 0], y_range=[-0.01, 0.2])

source = ColumnDataSource(compute_source(0, T, t, f))
plotLI.multi_line('xx', 'tt', source=source)
plotTC.multi_line('xx', 'TT', source=source)
plotf.multi_line('xx', 'ff', source=source) #this can be improved since only one dimension is needed, right now f is constructed from repmat
#import pdb; pdb.set_trace()
source_vec = ColumnDataSource(compute_source_vec(0, T, Taustarlist))
plotr.line('TTaustarlist', 'TT', source=source_vec) 

label = Label(x=1.1, y=18, text=str(0), text_font_size='70pt', text_color='#eeeeee')
plotLI.add_layout(label)
    
# Set up widgets
#text = TextInput(title="title", value='big T')
k_slider = Slider(title="k", value=k, start=1, end=10, step=1)
Taustar_min_slider = Slider(title="tstr min", value=Taustar_min, start=1, end=4, step=1)
Taustar_max_slider = Slider(title="tstr max", value=Taustar_max, start=6, end=15, step=1)
buff_len_slider = Slider(title="number of tstr nodes", value=buff_len , start=3, end=100, step=1)
time_limit_slider = Slider(title="Time", start=0., end=10.-dtime, value=0, step=1)
button_update = Button(label='Update', width=60)


# Set up callbacks
#def update_title(attrname, old, new):
 #   plot.title.text = text.value


#text.on_change('value', update_title)


def update_data():
    # Get the current slider values
    k = k_slider.value
    Taustar_min = Taustar_min_slider.value
    Taustar_max = Taustar_max_slider.value
    buff_len = buff_len_slider.value
    f, t, T, time_vec, Taustarlist = compute_t_T(k, Taustar_max, Taustar_min, buff_len, length_time, dtime)
    source.data = compute_source(length_time, T, t, f)
    source_vec.data = compute_source_vec(length_time, T, Taustarlist)
    return t, T, time_vec

def animate_update():
    time_limit = time_limit_slider.value + 1
    if time_limit >= length_time-dtime:
        time_limit = time_limit -1
        curdoc().remove_periodic_callback(animate_update)
        button_play.label = '► Play'
    label.text = str(time_limit)
    source.data = compute_source(time_limit, T, t, f)
    source_vec.data = compute_source_vec(time_limit, T, Taustarlist)
    print(time_limit)
    time_limit_slider.value = time_limit 


def time_limit_slider_update(attrname, old, new):
    time_limit = time_limit_slider.value
    label.text = str(time_limit)
    source.data = compute_source(time_limit, T, t, f)
    source_vec.data = compute_source_vec(length_time, T, Taustarlist)
        
        
time_limit_slider.on_change('value', time_limit_slider_update)


def animate():
    if button_play.label == '► Play':
        button_play.label = '❚❚ Pause'
        curdoc().add_periodic_callback(animate_update, 200)
    else:
        button_play.label = '► Play'
        curdoc().remove_periodic_callback(animate_update)


button_play = Button(label='► Play', width=60)
button_play.on_click(animate)

button_update.on_click(update_data)


#for w in [k_slider, Taustar_min_slider, Taustar_max_slider, buff_len_slider]:
#    w.on_change('value', update_data)

# Set up layouts and add to document
inputs = widgetbox(k_slider, Taustar_min_slider, Taustar_max_slider, buff_len_slider, button_update)

layout = layout([
    [plotf], 
    [plotLI], 
    [plotTC],
    [plotr],
    [time_limit_slider, button_play],
], sizing_mode='fixed')

curdoc().add_root(row(inputs, layout, width=800))
#curdoc().add_root(row(layout, width=800))
curdoc().title = "Compressed memory"
