


import sys
import numpy as np
import h5py

from mpi4py import MPI


from proc.algorithms.Converter import Converter


from para.MasterWorker import MasterWorker




class Coinc(MasterWorker):
         

    def __init__(self,config):

        super(Coinc, self).__init__(map_func=self.process_data,config=config)
        
        

 
        self.config = config
        self.init_params(config)




        print('Starting worker: {0}.'.format(self.mpi_rank))
        sys.stdout.flush()

        return

   
    def init_params(self,config):
    
        
        
        
        self.input_folder = config['Input']['folder']
        self.input_runs = config['Input']['runs'].split(',') 
        
        self.atts = config['Input']['atts'].split(',') 
        
        
        self.dn1 = config['Input']['dn1']
        self.dn2 = config['Input']['dn2']
        
        self.dn_plseng = config['Input']['dn_plseng']
        self.bunch_offset = int(config['Input']['bunch_offset'])
        
                
        self.output_folder = config['Output']['folder']
        self.output_name1 = config['Output']['name1']

        
        self.signals = config['Signals']['names'].split(',')
        
        self.sigdict = {}
        self.sigdict['maxnumpls'] = int(config['Signals']['maxnumpls'])
        self.sigdict['maxtime'] = int(config['Signals']['maxtime'])
        self.sigdict['sample_interval'] = float(config['Signals']['sample_interval'])
        
        self.paramCFD = {}
        self.paramHit = {}
        
            
        for sig in self.signals:
            self.paramCFD[sig] = {}
            for key in config[sig]:
                if key == 'polarity' or key == 'id':
                    self.paramCFD[sig][key] = config[sig][key]
                else:
                    self.paramCFD[sig][key] = float(config[sig][key])
         
         
        for key in config['Reconstruction']:

            self.paramHit[key] = float(config['Reconstruction'][key])       
            
        
        
 


    def process_data(self,files,files_plseng,atts_plseng,mpi_rank):
    
        svfn = self.output_folder+self.output_name1+'_wk'+str(mpi_rank)+'.h5'
        conv = Converter(len(files),files,files_plseng,atts_plseng,self.sigdict,self.dn1,self.dn2,self.dn_plseng,self.bunch_offset,self.paramCFD,self.paramHit,svfn)
        conv.convert()
        conv.savedata()
    

        return 
        
