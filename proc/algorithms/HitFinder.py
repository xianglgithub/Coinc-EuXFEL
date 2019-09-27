

import numpy as np


class HitFinder:

    def __init__(self, params):
                                #parameter initializations
        self.uRunTime = params['runtime_u']
        
        self.vRunTime = params['runtime_v']
        
        self.wRunTime = params['runtime_w']
        
        self.uTSumAvg = params['tsum_avg_u']
        self.uTSumLow = self.uTSumAvg - params['tsum_hw_u']
        self.uTSumHigh = self.uTSumAvg + params['tsum_hw_u']
        
        self.vTSumAvg = params['tsum_avg_v'] 
        self.vTSumLow = self.vTSumAvg - params['tsum_hw_v']
        self.vTSumHigh = self.vTSumAvg + params['tsum_hw_v']
          
        self.wTSumAvg = params['tsum_avg_w']
        self.wTSumLow = self.wTSumAvg - params['tsum_hw_w']
        self.wTSumHigh = self.wTSumAvg + params['tsum_hw_w']
        
        self.f_u = params['f_u']
        self.f_v = params['f_v']
        self.f_w = params['f_w']
        
         
        self.sqrt3 = np.sqrt(3.)
        
    
        
        
    def FindHits(self, McpSig, u1Sig, u2Sig, v1Sig, v2Sig, w1Sig, w2Sig):
        
        
        
        t1u = (-self.uRunTime+2*McpSig+self.uTSumAvg)/2
        t2u = (self.uRunTime+2*McpSig+self.uTSumAvg)/2
            
        t1v = (-self.vRunTime+2*McpSig+self.vTSumAvg)/2
        t2v = (self.vRunTime+2*McpSig+self.vTSumAvg)/2
            
        t1w = (-self.wRunTime+2*McpSig+self.wTSumAvg)/2
        t2w = (self.wRunTime+2*McpSig+self.wTSumAvg)/2   
        
        
        self.Xf = {}
        self.Yf = {}
        self.Tf = {}
        
        self.subf = {}
        self.sumf = {}
        
        for k in ['u','v','w']:
            self.subf[k] = np.array([])
            self.sumf[k] = np.array([])        
        
        for k in ['uv','uw','vw']:
            self.Xf[k] = np.array([])
            self.Yf[k] = np.array([])
            self.Tf[k] = np.array([])

    
               
        for i_McpT, McpT in enumerate(McpSig):
           
            u1 = u1Sig[(u1Sig>t1u[i_McpT]) & (u1Sig<t2u[i_McpT])]
            u2 = u2Sig[(u2Sig>t1u[i_McpT]) & (u2Sig<t2u[i_McpT])]
            v1 = v1Sig[(v1Sig>t1v[i_McpT]) & (v1Sig<t2v[i_McpT])]
            v2 = v2Sig[(v2Sig>t1v[i_McpT]) & (v2Sig<t2v[i_McpT])]      
            w1 = w1Sig[(w1Sig>t1w[i_McpT]) & (w1Sig<t2w[i_McpT])]
            w2 = w2Sig[(w2Sig>t1w[i_McpT]) & (w2Sig<t2w[i_McpT])]                          
            
            u1u2_sum = u1[:,np.newaxis] + u2[np.newaxis,:] - 2*McpT 
            v1v2_sum = v1[:,np.newaxis] + v2[np.newaxis,:] - 2*McpT 
            w1w2_sum = w1[:,np.newaxis] + w2[np.newaxis,:] - 2*McpT 
            
            
            #time sum conditions
            
          #  u1_ind, u2_ind = np.where(True | (u1u2_sum>self.uTSumLow) | (u1u2_sum<self.uTSumHigh))
          #  v1_ind, v2_ind = np.where(True | (v1v2_sum>self.vTSumLow) | (v1v2_sum<self.vTSumHigh))
          #  w1_ind, w2_ind = np.where(True | (w1w2_sum>self.wTSumLow) | (w1w2_sum<self.wTSumHigh))
            
            u1_ind, u2_ind = np.where((u1u2_sum>self.uTSumLow) & (u1u2_sum<self.uTSumHigh))
            v1_ind, v2_ind = np.where((v1v2_sum>self.vTSumLow) & (v1v2_sum<self.vTSumHigh))
            w1_ind, w2_ind = np.where((w1w2_sum>self.wTSumLow) & (w1w2_sum<self.wTSumHigh))
                    
            
            
            sub_u = u1[u1_ind]-u2[u2_ind]
            sub_v = v1[v1_ind]-v2[v2_ind]
            sub_w = w1[w1_ind]-w2[w2_ind]
            
            sub_uf = sub_u*self.f_u/2
            sub_vf = sub_v*self.f_v/2
            sub_wf = sub_w*self.f_w/2
            
            sum_u = u1[u1_ind]+u2[u2_ind] - 2*McpT 
            sum_v = v1[v1_ind]+v2[v2_ind] - 2*McpT 
            sum_w = w1[w1_ind]+w2[w2_ind] - 2*McpT 
            

            self.subf['u'] = np.concatenate([self.subf['u'],sub_u],axis=0)
            self.subf['v'] = np.concatenate([self.subf['v'],sub_v],axis=0)
            self.subf['w'] = np.concatenate([self.subf['w'],sub_w],axis=0)
            
            self.sumf['u'] = np.concatenate([self.sumf['u'],sum_u],axis=0)
            self.sumf['v'] = np.concatenate([self.sumf['v'],sum_v],axis=0)
            self.sumf['w'] = np.concatenate([self.sumf['w'],sum_w],axis=0)            
            
            Xuv = sub_uf[:,np.newaxis] + 0*sub_vf[np.newaxis,:]
            Yuv = (sub_uf[:,np.newaxis] - 2*sub_vf[np.newaxis,:])/self.sqrt3      
            
            Xuw = sub_uf[:,np.newaxis] + 0*sub_wf[np.newaxis,:]
            Yuw = (2*sub_wf[np.newaxis,:] - sub_uf[:,np.newaxis])/self.sqrt3
            
            Xvw = sub_vf[:,np.newaxis] + sub_wf[np.newaxis,:]
            Yvw = (sub_wf[np.newaxis,:] - sub_vf[:,np.newaxis])/self.sqrt3           
            
            self.Xf['uv'] = np.concatenate([self.Xf['uv'],np.ravel(Xuv)],axis=0)
            self.Yf['uv'] = np.concatenate([self.Yf['uv'],np.ravel(Yuv)],axis=0)
            self.Tf['uv'] = np.concatenate([self.Tf['uv'],McpT*np.ones(len(np.ravel(Xuv)))],axis=0)
            
            self.Xf['uw'] = np.concatenate([self.Xf['uw'],np.ravel(Xuw)],axis=0)
            self.Yf['uw'] = np.concatenate([self.Yf['uw'],np.ravel(Yuw)],axis=0)   
            self.Tf['uw'] = np.concatenate([self.Tf['uw'],McpT*np.ones(len(np.ravel(Xuw)))],axis=0)
            
            self.Xf['vw'] = np.concatenate([self.Xf['vw'],np.ravel(Xvw)],axis=0)
            self.Yf['vw'] = np.concatenate([self.Yf['vw'],np.ravel(Yvw)],axis=0)   
            self.Tf['vw'] = np.concatenate([self.Tf['vw'],McpT*np.ones(len(np.ravel(Xvw)))],axis=0)
