import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from io import StringIO
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')


# filedir = r'C:\Users\bwroe\Desktop\Analysis folder'



default_params = {
    'out_dir': r'C:\Users\bwroe\Desktop\Analysis folder\Output',
    'csv_dir': r'C:\Users\bwroe\Desktop\git\Brian_scripts\HEKA\csv',
    
    'number_of_freqs': 31,
    'highest_freq' : '1k',
    'lowest_freq' : '1',
    'filter_setting': '30k',
    'correct_filter': True,
    'filter_correct_mode': False,
    
    'sample_time' : 10,
    'AUTOLAB': 1,
    'autolab_current_range': 1e-6,    
    'plot_Nyquist' : 1,
    'plot_Bode' : 1,
    
    'save_Excel' : 0,
    'save_ASCII' : 0,
        }



class FT_EIS:
    
    def __init__(self, file, params):
        
        self.file = file
        
        self.out_dir = params['out_dir']
        self.csv_dir = params['csv_dir']
        
        self.number_of_freqs = params['number_of_freqs']
        self.highest_freq = params['highest_freq']
        self.lowest_freq = params['lowest_freq']
        self.filter_setting = params['filter_setting']
        self.correct_filter = params['correct_filter']
        self.filter_correct_mode = params['filter_correct_mode']
        
        self.sample_time = params['sample_time']
        self.AUTOLAB = params['AUTOLAB']
        self.autolab_current_range = params['autolab_current_range']
        self.plot_Nyquist = params['plot_Nyquist']
        self.plot_Bode = params['plot_Bode']
        
        self.save_Excel = params['save_Excel']
        self.save_ASCII = params['save_ASCII']
        
        self.Nyquist = None
        self.Bode = None
        
        
        self.run_file()




    def run_file(self):
        self.get_metadata()
        self.get_data()
        self.FourierTransform()
        
        if self.plot_Nyquist:
            self.make_Nyquist()
        
        if self.plot_Bode:
            self.make_Bode()
        
        if self.save_Excel:
            if not self.filter_correct_mode:
                self.WriteToExcel()
            
        if self.save_ASCII:
            if not self.filter_correct_mode:
                self.WriteToASCII()
                
        if self.filter_correct_mode:
            self.filter_correct()
        

    ##########################################################
    #####                                                #####
    #####           DATA EXTRACTION FUNCTIONS            #####
    #####                                                #####
    ##########################################################
    
    def get_metadata(self):
        # Record acquisition date, metadata line, user comment from PATCHMASTER,
        # and the starting line for the actual dataset
        comment = str()
        startline = 5
        with open(self.file) as f:
            for lineno, line in enumerate(f):
                
                if lineno == 1:
                    meta = line
                    self.date = line.split()[3].replace('/', '-')
                
                if lineno > 1:
                    if line.startswith('"Index"'):
                        startline = lineno
                        break
                    else:
                        comment = comment + ' ' + line
                
                if lineno > startline:
                    break
        
        self.metadata = pd.Series({'meta':meta,
                          'comment': comment}) 
        self.startline = startline
        
    
    
    
    def get_data(self):
        
        d = dict()
        
        # Use StringIO for faster file reading, and handling multiple
        # PATCHMASTER Series in a single file
        
        s = StringIO()
        
        def isfloat(x):
            try:
                float(x)
                return True
            except:
                return False
        
        with open(self.file) as f:
            for line in f:
                if isfloat(line.split(',')[0]):
                    #skip rows which don't start with the index number (1, 2, ...)
                    s.write(line)
    
        
        s.seek(0) #return to top of StringIO
        
        
        df = pd.read_csv(s, names=(
            "index", "time", "current", "time2", "voltage"), 
            engine='c', dtype=float)
        df = df.drop(['time2', 'index'], axis=1)
        
        
        if self.AUTOLAB == 1:
            # Correct for *-1 factor
            # Convert AD voltage -> current
            df['current'] = -df['current']*self.autolab_current_range
        
        
        
        # Sample data by sample_time
        samplefreq = np.round(1/(df.iloc[1,0] - df.iloc[0,0]))
        points_per_cycle = self.sample_time * samplefreq
        cycles = np.arange(1, int(len(df)/points_per_cycle) + 1)
        
        
        for i in cycles:
            time_lower = (i-1)*self.sample_time
            time_upper = i*self.sample_time
            d[i] = df.iloc[int(samplefreq*time_lower) : 
                            int(samplefreq*time_upper), :] 
        
        print('Sampling Frequency: %d kHz, %s s, %d cycles.' %(
            samplefreq/1000, int(len(df)/samplefreq), len(cycles)))    
    
        self.d = d
    
    
    def P2R(self, radii, angles):
        # Polar to rectangular coordinates
        return radii * np.exp(1j*angles*np.pi/180)
    
    
    def FT(self, x):
        # Helper function to fourier transform and discard 0 freq.
        return np.fft.rfft(x)[1:]
    
    
    def FourierTransform(self):
        '''
        Main function to fourier transform dict of time-domain data
        
        Parameters
        ----------
        d : dict
            d[1], d[2], ... = time domain dataframes corresponding
            to spectra 1, 2, ...
    
        Returns
        -------
        ft : dict
            Dict of fourier transformed spectra.
    
        '''
        
        # Get frequencies   
        try:
            file = self.csv_dir + r'\f_%s_%s_%sfreqs.csv' %(self.highest_freq, 
                                                            self.lowest_freq,
                                                            self.number_of_freqs)
                

            f = pd.read_csv(file)
        
        except FileNotFoundError:
            print('\nERROR: File not found:')
            print(file)
            print('''Check that highest_freq, lowest_freq, and number_of_freqs
                  are correct''')
            import sys
            sys.exit()
        
        
        freq_array = f['frequency'].to_numpy()
          
        
        samplefreq = np.round(1/(self.d[1].iloc[1,0] - self.d[1].iloc[0,0]))
        freqs = np.fft.rfftfreq(self.d[1]['voltage'].size)[1:]
        freqs = np.round(freqs*samplefreq, 3)
        
        
        
        # Fourier transforms
        ft = dict()
        
        for i in self.d:
            ft[i] = pd.DataFrame(
                {'f':freqs,
                    'V':self.FT(self.d[i]['voltage']),
                    'I':self.FT(self.d[i]['current'])}
                )
            
            # Keep only applied frequencies
            ft[i] = ft[i][ft[i]['f'].isin(freq_array)]
            
            ft[i]['Z'] = np.absolute(ft[i]['V']/ft[i]['I'])
            ft[i]['phase'] = np.angle(ft[i]['V']/ft[i]['I'], deg=True)
            
            # Correct for filter transfer function
            if self.correct_filter == True and self.filter_correct_mode == False and self.AUTOLAB == False:
                
                corr_file = self.csv_dir + r'\c_%sBessel_%s.xlsx' %(
                    self.filter_setting, self.highest_freq)
                
                corr_df = pd.read_excel(corr_file)
                Z_corr = corr_df['Z_factor'].to_numpy()
                phase_corr = corr_df['phase_factor'].to_numpy()
                
                ft[i]['Z'] = ft[i]['Z']/Z_corr
                ft[i]['phase'] = ft[i]['phase'] - phase_corr
            
            # Re-calculate complex impedance
            ft[i]['Z'] = self.P2R(ft[i]['Z'], ft[i]['phase'])
            ft[i]['re'] = np.real(ft[i]['Z'])
            ft[i]['im'] = np.imag(ft[i]['Z'])
            
        self.ft = ft
    
    
    def get_units(self, ft):
        max_Z = max(abs(ft[1]['Z']))
        
        if max_Z <= 1000:
            return ('', 1)
        
        if max_Z > 1e3 and max_Z <= 1e6:
            return('k', 1e3)
        
        if max_Z > 1e6:
            return('M', 1e6)
    
    
    
    def make_Nyquist(self, add_comment=True):
        fig, ax = plt.subplots()
        prefix, factor = self.get_units(self.ft)
    
        
        if len(self.ft) == 1:
            for i in self.ft:
                ax.plot(self.ft[i]['re']/factor, -self.ft[i]['im']/factor, 'o-')
        else: 
            colors = plt.cm.plasma(np.linspace(0.2,0.8, len(self.ft)))
            
            for i in self.ft:
                # First cycle is often noisy, don't plot it if possible
                if i > 1:
                    ax.plot(self.ft[i]['re']/factor, -self.ft[i]['im']/factor, 
                            'o-', color=colors[i-1])
        
        unit = prefix + '$\Omega$'
        ax.set_xlabel("Z'/ " + unit)
        ax.set_ylabel("Z''/ " + unit)
        
        xmax = ax.get_xlim()[1]
        ymax = ax.get_ylim()[1]
        
        if xmax > ymax:
            ax.set_xlim(0, )
            ax.set_ylim(ax.get_xlim())
            
        else:
            ax.set_ylim(0, )
            ax.set_xlim(ax.get_ylim())
        
        if add_comment == True:
            plt.text(0.1, -0.2, self.metadata['comment'], transform=plt.gcf().transFigure)
        plt.show()
        
        self.Nyquist = fig
    
    
        
    def make_Bode(self, add_comment=True):
        fig, ax1 = plt.subplots()
        prefix, factor = self.get_units(self.ft)
        ax1.set_xscale('log')
        ax2 = ax1.twinx()
        if len(self.ft) == 1:
            for i in self.ft:
                ax1.plot(self.ft[i]['f'], np.absolute(self.ft[i]['Z'])/factor, '.-')
                ax2.plot(self.ft[i]['f'], self.ft[i]['phase'], 'x')
        else:
            colors = plt.cm.plasma(np.linspace(0.2,0.8, len(self.ft)))
            for i in self.ft:
                # First cycle is often noisy, don't plot it if possible
                if i > 1:
                    ax1.plot(self.ft[i]['f'], np.absolute(self.ft[i]['Z'])/factor, 
                             '-', color=colors[i-1])
                    ax2.plot(self.ft[i]['f'], self.ft[i]['phase'], 'x',
                             color=colors[i-1])
        unit = prefix + '$\Omega$'
        ax1.set_xlabel("Frequency/ Hz")
        ax1.set_ylabel("|Z|/ " + unit)
        ax2.set_ylabel('Phase/ $\degree$')
        if add_comment == True:
            plt.text(0.1, -0.2, self.metadata['comment'], transform=plt.gcf().transFigure)
        plt.show()
        
        self.Bode = fig
    
    
    
    
    
    ##########################################################
    #####                                                #####
    #####           FILE SAVING FUNCTIONS                #####
    #####                                                #####
    ##########################################################
     
    def createFolder(self, directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
     
        
     
    def WriteToExcel(self) : 
        #saves DataFrame as Excel sheet
        
        file     = self.file
        ft       = self.ft
        out_dir  = self.out_dir
        date     = self.date
        metadata = self.metadata
        
        fileOut = date + "_" + file[:-4] +'.xlsx'
        writer = pd.ExcelWriter(out_dir+'\\'+fileOut, engine = 'xlsxwriter')
        
        
        for num in ft:
            if num == 1:
                ft[num].to_excel(writer, sheet_name = 'Impedance', index=False, 
                                   startrow=0, startcol=0, columns = ['f', 're', 'im', 'phase'], 
                                   header = ['Frequency/ Hz', 'realz' + str(num), 'imagz' + str(num), 'phase' + str(num)]) 
            else:
                ft[num].to_excel(writer, sheet_name = 'Impedance', index=False, 
                                   startrow=0, startcol=(3*num - 2), columns = ['re', 'im', 'phase'], 
                                   header = ['realz' + str(num), 'imagz' + str(num), 'phase' + str(num)])
         
        metadata.to_excel(writer, sheet_name = 'Metadata', index = False, header = False, startcol = 0)
        writer.save()
        print('Saved as ' + fileOut)
    
    
    
    def WriteToASCII(self):  
        #saves DataFrames(s) as ASCII for MEISP
        
        file     = self.file.split('\\')[-1]
        ft       = self.ft
        out_dir  = self.out_dir
        date     = self.date
        metadata = self.metadata
        
        folder_name = date + "_" + file[:-4]
        folder = os.path.join(out_dir, folder_name)
        print(folder)
        self.createFolder(folder)
        
        
        for num in ft:
            if num > 0 and num < 10:
                FileOutName = '000%ss.txt' %str(num)
            elif num >= 10 and num < 100:
                FileOutName = '00%ss.txt' %str(num)
            elif num >= 100 and num < 1000:
                FileOutName = '0%ss.txt' %str(num) 
            elif num >= 1000 and num < 10000:
                FileOutName = '%ss.txt' %str(num)
            elif num >= 10000:
                print('Need to write too many files!')
            ft[num].to_csv(folder +'\\'+FileOutName, columns = ['f', 're', 'im'],
                             header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], sep = '\t', index = False,
                             encoding='ascii')
    
        meta_file = folder +'\\'+ '0000_METADATA.txt'
        metadata.to_csv(meta_file, sep = ' ', index=False)
        
        
        if self.Nyquist:
            self.Nyquist.savefig(folder +'\\0000_Nyquist', dpi=300)
            
        if self.Bode:
            self.Bode.savefig(folder + '\\0000_Bode', dpi=300)
        
        print('Saved as ' + folder_name)
        
        
        
    
    
    
    def filter_correct(self):
        
        ftdf = self.ft
        max_freq = self.highest_freq
        filter_setting = self.filter_setting
        date = self.date
        csvdir = self.csv_dir
        
        correctdf = pd.DataFrame({'frequency': ftdf[1]['f'],
                                        'realz': ftdf[1]['re'],
                                        'imagz': ftdf[1]['im'],
                                        'phase_factor': ftdf[1]['phase']})
        correctdf['Z'] = np.sqrt(correctdf['realz']**2 + correctdf['imagz']**2)
        correctdf['Z_factor'] = correctdf['Z']/correctdf.iloc[0, 4]

        # Save the filter correction file
        main_folder = r'C:\Users\BRoehrich\Desktop\HEKA python'
        fileOut = 'c_%sBessel_%s.xlsx'%(filter_setting, max_freq)
        writer = pd.ExcelWriter(csvdir+'\\'+fileOut, engine = 'xlsxwriter')
        correctdf.to_excel(writer, index=False)
        writer.save()
        print('Saved filter correction as %s.'%fileOut)    
        
        # Save filter correction file to log
        fileOut_hist = date +'_%sfilter_%s.xlsx' %(filter_setting, max_freq)
        fileOutdir_hist = r'C:\Users\BRoehrich\Desktop\HEKA python\Filter corrections'
        writer2 = pd.ExcelWriter(fileOutdir_hist+'\\'+fileOut_hist, engine = 'xlsxwriter')
        correctdf.to_excel(writer2, index=False)
        writer2.save()










# if __name__ == '__main__':
#     l = []
#     n = 1
#     for file in os.listdir(filedir):
#         if file.endswith('.asc'):
#             print('File %s: ' % n, file)
            
#             f = os.path.join(filedir,file)
            
#             df = FT_EIS(f, default_params)
        

