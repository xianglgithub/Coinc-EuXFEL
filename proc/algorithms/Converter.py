

from proc.algorithms.HitFinder import HitFinder
from proc.algorithms.AcqirisPeakFinder import AcqirisPeakFinder
from proc.algorithms.AcqirisPeakFinderTg import AcqirisPeakFinderTg
import numpy as np
import h5py
import itertools


class Converter:
    def __init__(self, numInputHf,seqlst,seqlst_plseng,atts_plseng,sigdict,dn1,dn2,dn_plseng,bunch_offset,paramCFD,paramHit,hfsvn):
        
        
        #self.triggerT = triggerT
        
        tgnum = sigdict['maxnumpls']
        
        self.numPls = -1
        
        #self.plsSpacing = np.mean(self.triggerT[1:]-self.triggerT[:-1])
        
        self.wt = np.arange(0,sigdict['maxtime'],sigdict['sample_interval'])
        
       # self.plsID = np.zeros(len(self.wt),np.int)
        
       # for i, trigger in enumerate(self.triggerT):
       #     if i < self.numPls-1:
       #         self.plsID[(self.wt>=trigger) & (self.wt<self.triggerT[i+1])] = i
       #     else:
       #         self.plsID[self.wt>=trigger] = i
        
        self.dn1 = dn1
        self.dn2 = dn2
        self.dn_plseng = dn_plseng
        self.bunch_offset = bunch_offset
        
        self.cnmap = {}
        #old:self.cnmap = {'mcp': '4_B', 'u1':'2_D','u2':'1_B','v1':'3_D','v2':'2_B','w1':'3_B','w2':'1_D'} 
       # self.cnmap = {'mcp': '4_B', 'u1':'1_B','u2':'2_D','v1':'3_D','v2':'2_B','w1':'3_B','w2':'1_D'}  
        
        self.channels = ['mcp','u1','u2','v1','v2','w1','w2']
        
        for ch in self.channels:
            self.cnmap[ch] = paramCFD[ch]['id']
            
        
        self.seqlst = seqlst
        self.seqlst_plseng = seqlst_plseng
        self.atts_plseng = atts_plseng
        
        self.hf = HitFinder(paramHit)
        
        self.pf = {}
        
        self.hfsv = h5py.File(hfsvn,'w')
       # self.numpks = {}
       # self.lens = {}
        
        self.pf_tg = AcqirisPeakFinderTg(paramCFD['tg'])
        
        self.dset = {}
        
       
        self.dset['totalNumPls'] = self.hfsv.create_dataset('totalNumPls',(1,),dtype=np.int) 
        self.dset['totalNumPls'][0] = 0
        ##
        self.dset['plseng'] = self.hfsv.create_dataset('plseng',(numInputHf*500*tgnum,),dtype=np.int)
        ##
        for k in self.cnmap.keys():
            #(self, params,numPls,plsID,triggerT):
            self.pf[k] = AcqirisPeakFinder(paramCFD[k])

          #  self.numpks[k] = 0
           # self.lens[k] = 0
            self.dset['peaks_'+k] = self.hfsv.create_dataset('peaks_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,))     
            
            self.dset['numPeaksPerPls_'+k] = self.hfsv.create_dataset('numPeaksPerPls_'+k, 
                                                                      (numInputHf*500*tgnum,),
                                                                      dtype=np.int)  
            
            
            
        for k in ['u','v','w']:
            self.dset['sum_'+k] = self.hfsv.create_dataset('sum_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,)) 
            
            self.dset['sub_'+k] = self.hfsv.create_dataset('sub_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,)) 
                            
            
                                                                                                          
        
        for k in ['uv','uw','vw']:
            self.dset['sum_'+k] = self.hfsv.create_dataset('sum_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,)) 
            
            self.dset['sub_'+k] = self.hfsv.create_dataset('sub_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,))             
            
            self.dset['hits_t_'+k] = self.hfsv.create_dataset('hits_t_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,))  
            self.dset['hits_x_'+k] = self.hfsv.create_dataset('hits_x_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,))  
            self.dset['hits_y_'+k] = self.hfsv.create_dataset('hits_y_'+k, (0,), maxshape=(None,),
                            dtype=np.float, chunks=(10**4,))              
            
            self.dset['numHitsPerPls_'+k] = self.hfsv.create_dataset('numHitsPerPls_'+k, 
                                                                     (numInputHf*500*tgnum,),
                                                                     dtype=np.int)                  

            
    def setTg(self,num_tg,tg):
        self.numPls = num_tg
        self.triggerT = tg
        for k in self.cnmap.keys():
            self.pf[k].setTg(num_tg,tg)
            
        
        
        
    
    def convertOneF(self,seq,seq_plseng,att_plseng):
        df = h5py.File(seq,'r')
        df_plseng = h5py.File(seq_plseng,'r')
         
        self.peaks = {}
        self.hits = {}
        
        trainNum = df[self.dn1+self.cnmap['u1']+self.dn2].shape[0]
        
        for i in range(trainNum):
            #t0 = time.time()
            wf_tg = df[self.dn1+'4_D'+self.dn2][i]
            tg = self.pf_tg.cfd(wf_tg,self.wt)
            num_tg = len(tg)
            if num_tg == 0:
                continue
            if num_tg != self.numPls:
                self.setTg(num_tg,tg)
                
            ##
            self.dset['plseng'][self.dset['totalNumPls'][0]:(self.dset['totalNumPls'][0]+self.numPls)] = (df_plseng[self.dn_plseng][i][self.bunch_offset:(self.bunch_offset+self.numPls)])*att_plseng
            ##
                            
            for k in self.cnmap.keys():
                
                
                wf = df[self.dn1+self.cnmap[k]+self.dn2][i]
                self.peaks[k] = self.pf[k].cfd(wf,self.wt)
                #pks = self.pf[k].cfd(wf,self.wt)
                
                for j, triger in enumerate(self.triggerT):
                
                    npk = len(self.peaks[k][j])
                  #  self.numpks[k] += npk
                
                
                    self.dset['numPeaksPerPls_'+k][self.dset['totalNumPls'][0]+j] = npk
                    if npk > 0:
                        self.dset['peaks_'+k].resize(self.dset['peaks_'+k].shape[0]+npk,axis=0)
                        self.dset['peaks_'+k][-npk:] = self.peaks[k][j]
                        
            
            for j, triger in enumerate(self.triggerT):
            
                self.hf.FindHits(self.peaks['mcp'][j],self.peaks['u1'][j],self.peaks['u2'][j],
                                 self.peaks['v1'][j],self.peaks['v2'][j],self.peaks['w1'][j],
                                 self.peaks['w2'][j])
                
                for k in ['u','v','w']:
                    ns = len(self.hf.sumf[k])
                    if ns>0:
                        self.dset['sum_'+k].resize(self.dset['sum_'+k].shape[0]+ns,axis=0)
                        self.dset['sum_'+k][-ns:] = self.hf.sumf[k]
                        
                        self.dset['sub_'+k].resize(self.dset['sub_'+k].shape[0]+ns,axis=0)
                        self.dset['sub_'+k][-ns:] = self.hf.subf[k]                        
                    
                
                for k in ['uv','uw','vw']:
                    
                    
                    nh = len(self.hf.Xf[k])
                    
                    self.dset['numHitsPerPls_'+k][self.dset['totalNumPls'][0]+j] = nh
                    if nh>0:
                    #print(k,':',nh)
                    
                   # print('1x:',self.dset['hits_x_'+k].shape,self.hf.Xf[k].shape)
                   # print('1y:',self.dset['hits_y_'+k].shape,self.hf.Yf[k].shape)
                    
                        self.dset['hits_x_'+k].resize(self.dset['hits_x_'+k].shape[0]+nh,axis=0)
                        self.dset['hits_x_'+k][-nh:] = self.hf.Xf[k]
                    

                        
                        self.dset['hits_y_'+k].resize(self.dset['hits_y_'+k].shape[0]+nh,axis=0)
                        self.dset['hits_y_'+k][-nh:] = self.hf.Yf[k]
                    
                   # print('2x:',self.dset['hits_x_'+k].shape,self.hf.Xf[k].shape)
                   # print('2y:',self.dset['hits_y_'+k].shape,self.hf.Yf[k].shape)                        

                        self.dset['hits_t_'+k].resize(self.dset['hits_t_'+k].shape[0]+nh,axis=0)
                        self.dset['hits_t_'+k][-nh:] = self.hf.Tf[k]    
                        
          
            self.dset['totalNumPls'][0] = self.dset['totalNumPls'][0] + self.numPls  
                                        

                                        
                                
            #t1 = time.time()

        df.close()
                        
                        
                
    def convert(self):
        for i, seq in enumerate(self.seqlst):
            self.convertOneF(seq,self.seqlst_plseng[i],self.atts_plseng[i])
            print(seq+' processed!********************')
            
                
    def savedata(self):
        
        self.hfsv.close()
    
    
