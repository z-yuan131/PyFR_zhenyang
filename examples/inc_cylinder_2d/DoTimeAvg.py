#!/usr/bin/env python
# coding: utf-8

# In[6]:


from pyfr.readers.native import NativeReader
from pyfr.inifile import Inifile

import numpy as np
import h5py
from collections import defaultdict


# In[13]:


class time_average(object):
    def __init__(self):
        # Series time
        start  = 0  #60
        end    = 75   #90
        self.dt = dt = 5     #5

        tt = np.arange(start, end, dt)
        self.time = list()
        for i in range(len(tt)):
            self.time.append("{:.2f}".format(tt[i]))
        
        
    def load_and_avg(self):
        
        self.soln_mean = defaultdict()
        temp = defaultdict()
        for t in self.time:
        
            name_prefix = f'inc_cylinder_2d-{t}.pyfrs'
            
            soln = NativeReader(name_prefix)
            
            for k in soln:
                if k.split('_')[0] == 'soln':
                    try:
                        self.soln_mean[k] += soln[k]
                    except KeyError:
                        self.soln_mean[k] = soln[k]
                        self.mesh_uuid = soln['mesh_uuid']
                        self.cfg = Inifile(soln['config'])
                        self.stats = Inifile(soln['stats'])
                        
        self.write_out()
        
            
    def write_out(self):
        
        out_putname = f'avg_soln.pyfrs'
        
        f = h5py.File(out_putname,'w')
        f['config'] = self.cfg.tostr()
        f['mesh_uuid'] = self.mesh_uuid
        for key in self.soln_mean:
            f[key] = self.soln_mean[key]/len(self.time)
            print(self.soln_mean[key].shape)
        f['stats'] = self.stats.tostr()
        f.close()
                            
            
time_average().load_and_avg()


# In[ ]:




