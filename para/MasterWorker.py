
# Adapted from OnDA's parallelization script to use the Maxwell cluster.


import sys
import mpi4py.MPI
#import mpi4py
import math
import time
import datetime

import glob

import numpy as np




class MasterWorker(object):

    NOMORE = 998
    DIETAG = 999
    DEADTAG = 1000

    def __init__(self, map_func,config):

        debug = False


        self._buffer = None


        self.mpi_rank = mpi4py.MPI.COMM_WORLD.Get_rank()
        self.mpi_size = mpi4py.MPI.COMM_WORLD.Get_size()
        if self.mpi_rank == 0:
            self.role = 'master'
        else:
            self.role = 'worker'



        self.map = map_func

        
        self.files = []
        self.files_plseng = []
        
        self.atts_plseng = []
        
 


        
        if self.role == 'worker':
            self.num_lost_events_evt = 0        
            self.num_lost_events_time = 0
            self.num_lost_events_data = 0            
            
        if self.role == 'master':
            self.num_lost_events_evt = 0                    
            self.num_lost_events_time = 0
            self.num_lost_events_data = 0 
            
            self.num_reduced_events = 0
            self.num_nomore = 0

        return

    def shutdown(self, msg='Reason not provided.'):

        print ('Shutting down: {0}'.format(msg))

        if self.role == 'worker':
            self._buffer = mpi4py.MPI.COMM_WORLD.send(dest=0, tag=self.DEADTAG)
            mpi4py.MPI.Finalize()
            sys.exit(0)

        if self.role == 'master':
          #  self.save_func()
            try:
                for nod_num in range(1, self.mpi_size()):
                    mpi4py.MPI.COMM_WORLD.isend(0, dest=nod_num,
                                                tag=self.DIETAG)
                num_shutdown_confirm = 0
                while True:
                    if mpi4py.MPI.COMM_WORLD.Iprobe(source=mpi4py.MPI.ANY_SOURCE, tag=0):
                        self._buffer = mpi4py.MPI.COMM_WORLD.recv(source=mpi4py.MPI.ANY_SOURCE, tag=0)
                    if mpi4py.MPI.COMM_WORLD.Iprobe(source=mpi4py.MPI.ANY_SOURCE, tag=self.DEADTAG):
                        num_shutdown_confirm += 1
                    if num_shutdown_confirm == self.mpi_size() - 1:
                        break
              #  self.save_func(self.num_lost_events_time,self.num_lost_events_data)                  
                mpi4py.MPI.Finalize()
            except Exception:
              #  self.save_func(self.num_lost_events_time,self.num_lost_events_data)              
                mpi4py.MPI.COMM_WORLD.Abort(0)
            sys.exit(0)
        return

    def start(self, verbose=False):

        if self.role == 'worker':
            seqlst_DA1 = []
            seqlst_DA2 = []
            att = []
        
            for irun, run in enumerate(self.input_runs):
                #fd = self.path+'r'+run+'/'
                #fd = '/gpfs/exfel/exp/SQS/201921/p002412/raw/r0077/'
                fd = self.input_folder+run+'/'
                tDA1 = sorted(glob.glob(fd+'*DA01*'))
                #tDA2 = sorted(glob.glob(fd+'*DA02*'))
                tDA2 = [da1[:-14]+'DA02'+da1[-10:] for da1 in tDA1]
                seqlst_DA1 = seqlst_DA1 + tDA1
                seqlst_DA2 = seqlst_DA2 + tDA2
                att = att + [float(self.atts[irun])]*len(tDA1)   
                
 
                
            
            
            req = None


            for i in range(len(seqlst_DA2)):
               # if i==57:
               #     break
                if i % (self.mpi_size-1) == (self.mpi_rank-1):
                    self.files.append(seqlst_DA2[i])
                    self.files_plseng.append(seqlst_DA1[i])
                    self.atts_plseng.append(att[i])
                    
                    
                    
                
                
                # Check if a shutdown message is coming from the server
            if mpi4py.MPI.COMM_WORLD.Iprobe(source=0, tag=self.DIETAG):
                self.shutdown('Shutting down RANK: {0}.'.format(self.mpi_rank))
                


        #    result = 0

                # send the mapped event data to the master process
         #   if req:
          #      req.Wait()  # be sure we're not still sending something
          #  req = mpi4py.MPI.COMM_WORLD.isend(result, dest=0, tag=0)

            self.map(self.files,self.files_plseng,self.atts_plseng,self.mpi_rank)


            # When all events have been processed, send the master a
            # dictionary with an 'end' flag and die
            end_dict = {'end': True}
            if req:
                req.Wait()  # be sure we're not still sending something
            mpi4py.MPI.COMM_WORLD.isend((end_dict, self.mpi_rank), dest=0, tag=0)
            mpi4py.MPI.Finalize()
            sys.exit(0)

        # The following is executed on the master
        elif self.role == 'master':

            if verbose:
                print ('Starting master.')

            
            
            # Loops continuously waiting for processed data from workers
            while True:

                try:

                    buffer_data = mpi4py.MPI.COMM_WORLD.recv(
                        source=mpi4py.MPI.ANY_SOURCE,
                        tag=0)
                        
                    if 'end' in buffer_data[0].keys():
       
                                              
                        print ('Finalizing {0}'.format(buffer_data[1]))
                        self.num_nomore += 1
                        if self.num_nomore == self.mpi_size - 1:

                            print('All workers finished jobs.')
                            
                                                                  
                            print('Shutting down.')
                            self.end_processing()

                            mpi4py.MPI.Finalize()
                            sys.exit(0)
                        continue

                   # self.reduce(buffer_data)


                except KeyboardInterrupt as e:
                    print ('Recieved keyboard sigterm...')
                    print (str(e))
                    print ('shutting down MPI.')
                    
                    self.shutdown()
                    print ('---> execution finished.')
                    sys.exit(0)

        return

    def end_processing(self):
        print('Processing finished.')

        pass
