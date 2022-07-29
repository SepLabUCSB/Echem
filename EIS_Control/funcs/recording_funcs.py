import time
from array import array
import numpy as np
import pandas as pd


        

def record_frame(Rec, inst, frame_time, recording_params,
                 start_time, frame, save, process_last=True):
    
    frame_start_time = time.time()
    
    # Starts recording
    inst.write('TRMD AUTO')
    
    if process_last:
        process_frame(Rec, frame - 1, save)
        
    while time.time() - frame_start_time < frame_time:
        Rec.root.after(1)
    
    inst.write('TRMD STOP')
            
    volts1, volts2 = read_data(inst, recording_params)
    
    ft = transform_data(volts1, volts2, recording_params, 
                        frame_start_time)
    
    # Calculate Z
    V = ft.CH1data
    I = -ft.CH2data * recording_params['current_range']
    
    
    if Rec.potentiostat.get() != 'Autolab':
        I = -I
    
    Z       = V/I
    phase   = np.angle(Z, deg=True)
    freqs   = recording_params['freqs']
    
    df = pd.DataFrame({'freqs': ft.freqs,
                       'Z': Z,
                       'phase': phase})
    df = df[df['freqs'].isin(freqs)]
    
    # Apply reference correction
    if Rec.ref_corr_var.get():
        Z_corr, phase_corr = Rec.get_correction_values()
        
        if Z_corr == 0:
            return
        
        df['Z']     = df['Z'] / Z_corr
        df['phase'] = df['phase'] - phase_corr
    
    # Save to FTData object
    ft.freqs    = df['freqs'].to_numpy()
    ft.Z        = df['Z'].to_numpy()
    ft.phase    = df['phase'].to_numpy()
    ft.waveform = Rec.waveform.get()
    
    
    return ft



def read_data(inst, recording_params):
    
    vdiv1    = recording_params['vdiv1']
    vdiv2    = recording_params['vdiv2']
    voffset1 = recording_params['voffset1']
    voffset2 = recording_params['voffset2']
    
    # Get CH 1 data
    inst.write('C1:WF? DAT2')
    trace1 = inst.read_raw()
    wave1 = trace1[22:-2]
    adc1 = np.array(array('b', wave1))
    
    # Get CH 2 data
    inst.write('C2:WF? DAT2')
    trace2 = inst.read_raw()
    wave2 = trace2[22:-2]
    adc2 = np.array(array('b', wave2))
    
    # Convert to voltages
    volts1 = adc1*(vdiv1/25) - voffset1 
    volts2 = adc2*(vdiv2/25) - voffset2           
    
    return volts1, volts2



def transform_data(volts1, volts2, recording_params,
                   frame_start_time, sample_time=1):
    
    sara            = recording_params['sara']
    current_range   = recording_params['current_range']
    
    
    times = np.array([(1/sara*i) for i in range(len(volts1))])
    
    # Only Fourier transform first sample_time s
    end = None
    if sample_time:
        end = np.where(times == times[0] + sample_time)[0][0]
        
    freqs = sara*np.fft.rfftfreq(len(volts1[:end]))[1:]
    ft1   =      np.fft.rfft(volts1[:end])[1:]
    ft2   =      np.fft.rfft(volts2[:end])[1:]
    
    mean_I = np.mean(volts2 * current_range)
    Vpp = max(volts1) - min(volts1)
    
    ft = FourierTransformData(time  = frame_start_time,
                            freqs   = freqs,
                            CH1data = ft1,
                            CH2data = ft2,
                            Vpp     = Vpp,
                            mean_I  = mean_I)
    
    return ft
    


def process_frame(Rec, frame, save, save_path,
                  DC_file, time_file, update_time_plot):
    
    params = None
    
    if Rec.fit.get():
        fit_frame(Rec, frame)
        if Rec.circuit.get() == 'Randles_adsorption':
            Rct = Rec.ft[frame].params["R2"]
            Cad = Rec.ft[frame].params["Q2"]
            ket = 1/(2*Rct*Cad)
            print(f'Rct: {Rct}, Cad: {Cad}, ket: {ket}')
    
        params = Rec.ft[frame].params
        
        
    # Plot this result to figure canvas
    d = Rec.ft[frame]
    Z = np.abs(d.Z)
    phase = d.phase
    
    
    plot_Z      = Rec.plot_Z.get()
    plot_phase  = Rec.plot_phase.get()
    
    # Determine which plot to make
    if plot_Z:
        if not plot_phase:
            Rec.line1.set_xdata(d.freqs)          
            Rec.line1.set_ydata(Z)
            if hasattr(Rec.ft[frame], 'fits'):
                Rec.line3.set_xdata(d.freqs)
                Rec.line3.set_ydata(np.abs(Rec.ft[frame].fits))
            Rec.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
            Rec.ax.set_ylabel('|Z|/ $\Omega$')
            Rec.ax2.set_yticks([])
    
    if plot_phase:
        if not plot_Z:
            Rec.line2.set_xdata(d.freqs)
            Rec.line2.set_ydata(phase)
            if hasattr(Rec.ft[frame], 'fits'):
                Rec.line4.set_xdata(d.freqs)
                Rec.line4.set_ydata(np.angle(Rec.ft[frame].fits))
            Rec.ax2.set_ylim(min(phase)-10, max(phase)+10)
            Rec.ax2.set_ylabel('Phase/ $\degree$')
            Rec.ax2.set_yticks([])
        
    if plot_Z and plot_phase:
        Rec.line1.set_xdata(d.freqs)
        Rec.line2.set_xdata(d.freqs)
        Rec.line1.set_ydata(Z)
        Rec.line2.set_ydata(phase)
        if hasattr(Rec.ft[frame], 'fits'):
            Rec.line3.set_xdata(d.freqs)
            Rec.line4.set_xdata(d.freqs)
            Rec.line3.set_ydata(np.abs(Rec.ft[frame].fits))
            Rec.line4.set_ydata(np.angle(Rec.ft[frame].fits, deg=True))
        Rec.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
        Rec.ax2.set_ylim(min(phase)-10, max(phase)+10)
        Rec.ax.set_ylabel('|Z|/ $\Omega$')
        Rec.ax2.set_ylabel('Phase/ $\degree$')
        
        
    # Draw the plot
    if (plot_Z or plot_phase):
        Rec.fig.tight_layout()
        Rec.ax.set_xticks([1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
        Rec.ax.set_xlim(0.7*min(d.freqs), 1.5*max(d.freqs))
        Rec.fig.canvas.draw()
        Rec.fig.canvas.flush_events()
    
    if update_time_plot:
        Rec.update_time_plot(d.time, d.freqs, d.Z, d.phase, params,
                             ax=None)
    
    if save:
        # Add frame time to time list
        with open(time_file, 'a') as f:
            f.write(str(d.time) + '\n')
            f.close()
        with open(DC_file, 'a') as f:
            f.write(str(d.mean_I) + '\n')
            f.close()
        # Save frame as tab separated .txt
        Rec.save_frame(frame, d.freqs, np.real(d.Z),
                    np.imag(d.Z), save_path)
        
    
    
    return d.time, d.freqs, Z, phase, params



def fit_frame(Rec, frame, average=1, n_iter=25,
              starting_guess=None, **kwargs):
    
    if average < 0:
        print('Average must be greater than 0 frames!')
        return
    
    
    if frame == 0:
        bounds, starting_guess, params = Rec.initialize_circuit()
                
    elif frame > 0:
        starting_guess = Rec.ft[frame-1].params
        bounds = {param:[val/2, val*2] for param, val in
                  Rec.ft[frame-1].params.items()}
    
    
    if average == 1:
        Z     = Rec.ft[frame].Z
        freqs = Rec.ft[frame].freqs
        
    elif average > 1:
        x = frame - average
        y = frame + 1
        Z     = np.mean(Rec.ft[x:y].Z, axis=0)
        freqs = np.mean(Rec.ft[x:y].freqs, axis=0)
    
    
    # Perform fit
    print('Real time fitting not currently supported!')
    # DataFile = EIS_fit.DataFile(file='', circuit=Rec.circuit.get(), 
    #                     Z=Z, freqs=freqs, bounds=bounds)
    
    # DataFile.ga_fit(n_iter = n_iter, starting_guess = starting_guess, **kwargs)
    # DataFile.LEVM_fit(timeout = 0.4) # Needs short timeout
    #                                  # to not interfere with data
    #                                  # collection
    
    # # Save fit parameters, if fit was successful
    # try:
    #     self.ft[frame].params = DataFile.params # R, C parameters
    #     self.ft[frame].fits   = DataFile.fits   # Fitted Z vs freq
        
    #     if save:
    #         with open(fits_file, 'a') as f:
    #             if frame == 0:
    #                 f.write('time,')
    #                 for key, _ in self.ft[frame].params.items():
    #                     f.write(key + ',')
    #                 f.write('\n')
                
    #             f.write(str(self.ft[frame].time) + ',')
    #             for key, val in self.ft[frame].params.items():
    #                 f.write(str(val) + ',')
    #             f.write('\n')
    #             f.close()
    
    # except:
    #     pass        
    
    

class FourierTransformData:
    
    def __init__(self, time, freqs, CH1data, CH2data, Z = None,
                 phase = None, waveform = None, Vpp = None,
                 mean_I = None):
        self.time = time
        self.freqs = freqs
        self.CH1data = CH1data
        self.CH2data = CH2data
        
        self.Z = Z
        self.phase = phase
        self.waveform = waveform
        self.Vpp = Vpp
        self.mean_I = mean_I
        
        self.params = None
        
        
        