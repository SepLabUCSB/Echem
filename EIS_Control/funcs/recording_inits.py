import matplotlib.pyplot as plt
import os
import time
from datetime import date, datetime

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def init_recording(Rec, new_time_plot, n_plots, save):
    # Rec: Recorder object
    
    Rec.update_config_file()
    
    current_range = Rec.current_range.get('1.0', 'end')
    current_range = float(current_range)
    
    init_plots(Rec, new_time_plot, n_plots)
        
    vdiv1, vdiv2, voffset1, voffset2, sara, frame_time = read_scope_params(Rec)
    
    _, freqs = Rec.get_waveform()
            
    return locals()




    



def init_plots(Rec, new_time_plot, n_plots):
    # Initialize plots
    
    
    Rec.ax.set_xscale('linear')
    Rec.ax.clear()
    Rec.ax.set_xscale('linear')
    Rec.ax2.clear()
    
    # line1: Z data
    # line2: phase data
    # line3: Z fit
    # line4: phase fit
    Rec.line1, = Rec.ax.plot([],[], 'o', color=colors[0])
    Rec.line2, = Rec.ax2.plot([],[], 'x', color=colors[1])
    Rec.line3, = Rec.ax.plot([],[], '-', color=colors[0])
    Rec.line4, = Rec.ax2.plot([],[], '-', color=colors[1])
    
    # clear time plot
    if new_time_plot:
        Rec.init_time_plot(n_plots)
#        Rec.line5, = Rec.timeax.plot([],[], 'o', color='k')
                  
    
    Rec.ax.set_xscale('log')
    Rec.ax.set_xlabel('Frequency/ Hz')
    Rec.fig.tight_layout()
    Rec.canvas.draw_idle()



def read_scope_params(Rec):
    # Connect to scope
    inst = Rec.rm.open_resource(Rec.scope.get())
         
    
    # Set and record some scope parameters
    inst.write('TRMD AUTO')
    inst.write('MSIZ 70K')
    inst.write('TDIV 100MS')
    inst.write('TRMD STOP')
    
    vdiv1       = float(inst.query('C1:VDIV?')[8:-2])
    voffset1    = float(inst.query('C1:OFST?')[8:-2])
    vdiv2       = float(inst.query('C2:VDIV?')[8:-2])
    voffset2    = float(inst.query('C2:OFST?')[8:-2])
    sara        = float(inst.query('SARA?')[5:-5])
    tdiv        = float(inst.query('TDIV?')[5:-2])
    frame_time  = 14*tdiv
    
    return vdiv1, vdiv2, voffset1, voffset2, sara, frame_time



def init_save(Rec, save_path):
    # Create metadata file
    meta_file = os.path.join(save_path, '0000_Metadata.txt')
    
    with open(meta_file, 'w') as f:
        f.write('Waveform Vpp (mV): '+ str(Rec.waveform_vpp.get('1.0', 'end')))
        f.write('Waveform: '+ str(Rec.waveform.get())) 
        s_t = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
        f.write('\nStart time: %s'%s_t)
        
    f.close()
    
    # Start time list file
    time_file = os.path.join(save_path, '0000_time_list.txt')
    
    # Start fits file
    fits_file = os.path.join(save_path, '0000_fits.txt')
    
    # Start mean (DC) current file
    DC_file = os.path.join(save_path, '0000_DC_currents.txt')
    
    d = {'time_file': time_file,
         'meta_file': meta_file,
         'fits_file': fits_file,
         'DC_file'  : DC_file,
         'save_path': save_path}
    
    return d
    
    
