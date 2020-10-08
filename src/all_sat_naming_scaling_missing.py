#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 09:58:03 2018

@author: vincentch
"""
import numpy as np

class AllSatNamingScalingMissing():
    """ Gathers all static parameters configurations for each source dataset.
    TODO: collect and verify all sats values.
            now: ok for msg_alb values.
    """
    def __init__(self):
#        self.is_h5zip_needed={
#                        'vgt_alb_cgls':False,
#                        'vgt_alb_c3s':False,
#                        'eps_alb':True,
#                        'msg_alb_bb':True,
#                        'msg_alb_bb':True,
#                        'modis6_alb':False,
#                        'probav_alb':False
#                      }
        self.preproc_fct={
                        'vgt_alb_cgls':'open_mfdataset_preproc_for_CGLS_VGT',
                        'vgt_alb_c3s':'',
                        'eps_alb':'open_mfdataset_preproc_for_LSAF_msg_eps',
                        'msg_alb_bb':'open_mfdataset_preproc_for_LSAF_msg_eps',
                        'msg_alb_sp':'open_mfdataset_preproc_for_LSAF_msg_eps',
                        'msg_alb_sp_kernels':'open_mfdataset_preproc_for_LSAF_msg_eps',
                        'msg_alb_sp_cov':'open_mfdataset_preproc_for_LSAF_msg_eps',
                        'modis6_alb':'open_mfdataset_preproc_for_modis6',
                        'probav_alb':''
                      }
        self.in_default_scaling = {
                        'vgt_alb_cgls':{'bsa_shwv':1/10000,'bsa_vis':1/10000,'bsa_nir':1/10000,'wsa_shwv':1/10000,'wsa_vis':1/10000,'wsa_nir':1/10000,'bsa_shwv_err':1/10000,'bsa_vis_err':1/10000,'bsa_nir_err':1/10000,'wsa_shwv_err':1/10000,'wsa_vis_err':1/10000,'wsa_nir_err':1/10000,'brdf_alb_qual_flag':1,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1},
                        'vgt_alb_c3s':{'bsa_vis':1,'bsa_shwv':1,'bsa_nir':1,'wsa_vis':1,'wsa_shwv':1,'wsa_nir':1,'bsa_vis_err':1,'bsa_shwv_err':1,'bsa_nir_err':1,'wsa_vis_err':1,'wsa_shwv_err':1,'wsa_nir_err':1,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1},
                        'eps_alb':{'bsa_vis':1/10000,'bsa_shwv':1/10000,'bsa_nir':1/10000,'wsa_vis':1/10000,'wsa_shwv':1/10000,'wsa_nir':1/10000,'bsa_vis_err':1/10000,'bsa_shwv_err':1/10000,'bsa_nir_err':1/10000,'wsa_vis_err':1/10000,'wsa_shwv_err':1/10000,'wsa_nir_err':1/10000,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1},
                        'msg_alb_bb':{'bsa_vis':1/10000,'bsa_shwv':1/10000,'bsa_nir':1/10000,'wsa_vis':1/10000,'wsa_shwv':1/10000,'wsa_nir':1/10000,'bsa_vis_err':1/10000,'bsa_shwv_err':1/10000,'bsa_nir_err':1/10000,'wsa_vis_err':1/10000,'wsa_shwv_err':1/10000,'wsa_nir_err':1/10000,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1}, # ok verified
                        'msg_alb_sp':{'bsa_sp':1/10000,'wsa_sp':1/10000,'bsa_sp_err':1/10000,'wsa_sp_err':1/10000,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1}, # ok verified
                        'msg_alb_sp_kernels':{'brdf_alb_qual_flag':1,'age_info':1,'nmod':1, 'k0':1/10000,'k1':1/10000,'k2':1/10000}, # ok verified
                        'msg_alb_sp_cov':{'c00':1/1000,'c01':1/1000,'c02':1/1000,'c11':1/1000,'c12':1/1000,'c22':1/1000},
                        'modis6_alb':{'bsa_vis':1,'bsa_shwv':1,'bsa_nir':1,'wsa_vis':1,'wsa_shwv':1,'wsa_nir':1,'bsa_vis_err':1,'bsa_shwv_err':1,'bsa_nir_err':1,'wsa_vis_err':1,'wsa_shwv_err':1,'wsa_nir_err':1,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1},
                        'probav_alb':{'bsa_vis':1,'bsa_shwv':1,'bsa_nir':1,'wsa_vis':1,'wsa_shwv':1,'wsa_nir':1,'bsa_vis_err':1,'bsa_shwv_err':1,'bsa_nir_err':1,'wsa_vis_err':1,'wsa_shwv_err':1,'wsa_nir_err':1,'brdf_alb_qual_flag':1,'age_info':1,'nmod':1}
                      }
        self.in_default_offset = {
                        'vgt_alb_cgls':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0},
                        'vgt_alb_c3s':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0},
                        'eps_alb':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0},
                        'msg_alb_bb':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0}, # ok verified
                        'msg_alb_sp':{'bsa_sp':0,'wsa_sp':0,'bsa_sp_err':0,'wsa_sp_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0}, # ok verified
                        'msg_alb_sp_kernels':{'brdf_alb_qual_flag':0,'age_info':0,'nmod':0, 'k0':0,'k1':0,'k2':0}, # ok verified
                        'msg_alb_sp_cov':{'c00':0,'c01':0,'c02':0,'c11':0,'c12':0,'c22':0},
                        'modis6_alb':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0},
                        'probav_alb':{'bsa_vis':0,'bsa_shwv':0,'bsa_nir':0,'wsa_vis':0,'wsa_shwv':0,'wsa_nir':0,'bsa_vis_err':0,'bsa_shwv_err':0,'bsa_nir_err':0,'wsa_vis_err':0,'wsa_shwv_err':0,'wsa_nir_err':0,'brdf_alb_qual_flag':0,'age_info':0,'nmod':0}
                      }
        self.in_default_fillvalue = {
                        'vgt_alb_cgls':{'bsa_vis':65535,'bsa_shwv':65535,'bsa_nir':65535,'wsa_vis':65535,'wsa_shwv':65535,'wsa_nir':65535,'bsa_vis_err':65535,'bsa_shwv_err':65535,'bsa_nir_err':65535,'wsa_vis_err':65535,'wsa_shwv_err':65535,'wsa_nir_err':65535,'brdf_alb_qual_flag':None,'age_info':None,'nmod':None},
                        'vgt_alb_c3s':{'bsa_vis':65535,'bsa_shwv':65535,'bsa_nir':65535,'wsa_vis':65535,'wsa_shwv':65535,'wsa_nir':65535,'bsa_vis_err':65535,'bsa_shwv_err':65535,'bsa_nir_err':65535,'wsa_vis_err':65535,'wsa_shwv_err':65535,'wsa_nir_err':65535,'brdf_alb_qual_flag':None,'age_info':None,'nmod':None},
                        'eps_alb':{'bsa_vis':-1,'bsa_shwv':-1,'bsa_nir':-1,'wsa_vis':-1,'wsa_shwv':-1,'wsa_nir':-1,'bsa_vis_err':-1,'bsa_shwv_err':-1,'bsa_nir_err':-1,'wsa_vis_err':-1,'wsa_shwv_err':-1,'wsa_nir_err':-1,'brdf_alb_qual_flag':-1,'age_info':-1,'nmod':-1},
                        'msg_alb_bb':{'bsa_vis':-1,'bsa_shwv':-1,'bsa_nir':-1,'wsa_vis':-1,'wsa_shwv':-1,'wsa_nir':-1,'bsa_vis_err':-1,'bsa_shwv_err':-1,'bsa_nir_err':-1,'wsa_vis_err':-1,'wsa_shwv_err':-1,'wsa_nir_err':-1,'brdf_alb_qual_flag':-1,'age_info':-1,'nmod':-1},  # ok verified
                        'msg_alb_sp':{'bsa_sp':-1,'wsa_sp':-1,'bsa_sp_err':-1,'wsa_sp_err':-1,'brdf_alb_qual_flag':-1,'age_info':-1,'nmod':-1}, # ok verified
                        'msg_alb_sp_kernels':{'brdf_alb_qual_flag':-1,'age_info':-1,'k0':-32768,'k1':-32768,'k2':-32768}, # ok verified
                        'msg_alb_sp_cov':{'c00':-32768,'c01':-32768,'c02':-32768,'c11':-32768,'c12':-32768,'c22':-32768},
                        'modis6_alb':{'bsa_vis':32767,'bsa_shwv':32767,'bsa_nir':32767,'wsa_vis':32767,'wsa_shwv':32767,'wsa_nir':32767,'bsa_vis_err':None,'bsa_shwv_err':None,'bsa_nir_err':None,'wsa_vis_err':None,'wsa_shwv_err':None,'wsa_nir_err':None,'brdf_alb_qual_flag':255,'age_info':None,'nmod':None},
                        'probav_alb':{}
                      }
#        self.in_default_not_nan_flag_values = {
#                        'vgt_alb_cgls':list(np.where(np.arange(256) & 0b11111100101 == 0b00000000000)[0]),
#                        'vgt_alb_c3s':[],
#                        'eps_alb':list(np.where(np.arange(256) & 0b10001111 == 0b00001001)[0]),  # (cf LSASAF ETAL PUM, pp18-19, https://landsaf.ipma.pt/GetDocument.do?id=644)
#                        'msg_alb':list(np.where(np.arange(256) & 0b10001111 == 0b00000101)[0]) ,# (cf LSASAF MDAL PUM p39, https://landsaf.ipma.pt/GetDocument.do?id=615)
#                        'modis6_alb':[0,1],
#                        'probav_alb':[]
#                      }
        self.in_default_valid_range_before_scaling_and_offset = {   ### TO ADAPT TO REALISTIC DEFAULT VALUES
                        'vgt_alb_cgls':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]},# a verifier
                        'vgt_alb_c3s':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]} ,# a verifier
                        'eps_alb':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]} ,
                        'msg_alb_bb':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]} ,  
                        'msg_alb_sp':{'bsa_sp': [0, 10000],'wsa_sp': [0, 10000],'bsa_sp_err': [0, 10000],'wsa_sp_err': [0, 10000],'brdf_alb_qual_flag': [None,None],'age_info': [None,None],'nmod': [None,None]},
                        'msg_alb_sp_kernels':{'brdf_alb_qual_flag': [None,None],'age_info': [None,None],'nmod': [None,None],'k0': [None,None],'k1': [None,None],'k2': [None,None]},
                        'msg_alb_sp_cov':{'c00': [None,None],'c01': [None,None],'c02': [None,None],'c11': [None,None],'c12': [None,None],'c22': [None,None]},
                        'modis6_alb':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]} ,# a verifier
                        'probav_alb':{'bsa_vis':[0,10000],'bsa_shwv':[0,10000],'bsa_nir':[0,10000],'wsa_vis':[0,10000],'wsa_shwv':[0,10000],'wsa_nir':[0,10000],'bsa_vis_err':[0,10000],'bsa_shwv_err':[0,10000],'bsa_nir_err':[0,10000],'wsa_vis_err':[0,10000],'wsa_shwv_err':[0,10000],'wsa_nir_err':[0,10000],'brdf_alb_qual_flag':[None,None],'age_info':[None,None]} # a verifier
                      }
        self.in_default_snow_flag_values = {
                        'vgt_alb_cgls':list(np.where(np.arange(256) & 0b11111100111 == 0b00000000010)[0]),
                        'vgt_alb_c3s':list(np.where(np.arange(256) & 0b11111100111 == 0b00000000010)[0]),
                        'eps_alb':list(np.where(np.arange(256) & 0b10101111 == 0b00101001)[0]),  # (cf LSASAF ETAL PUM, pp18-19, https://landsaf.ipma.pt/GetDocument.do?id=644)
                        'msg_alb_bb':list(np.where(np.arange(256) & 0b10101111 == 0b00100101)[0]),# (cf LSASAF MDAL PUM p39, https://landsaf.ipma.pt/GetDocument.do?id=615)
                        'msg_alb_sp':list(np.where(np.arange(256) & 0b10101111 == 0b00100101)[0]),# (cf LSASAF MDAL PUM p39, https://landsaf.ipma.pt/GetDocument.do?id=615)
                        'msg_alb_sp_kernels':list(np.where(np.arange(256) & 0b10101111 == 0b00100101)[0]),# (cf LSASAF MDAL PUM p39, https://landsaf.ipma.pt/GetDocument.do?id=615)
                        'msg_alb_sp_cov':list(np.where(np.arange(256) & 0b10101111 == 0b00100101)[0]),# (cf LSASAF MDAL PUM p39, https://landsaf.ipma.pt/GetDocument.do?id=615)
                        'modis6_alb':[],
                        'probav_alb':[]
                      }
#        self.in_default_clim_flag_values = {
#                        'vgt_alb_cgls':[],
#                        'vgt_alb_c3s':[],
#                        'eps_alb':[],
#                        'msg_alb':[],
#                        'modis6_alb':[],
#                        'probav_alb':[]
#                      }
        # here shall gather only the names of the variables that are indeed in the input files
        self.input_varnames = {
                        'vgt_alb_cgls':{'AL-DH-BB':'bsa_shwv','AL-DH-VI':'bsa_vis','AL-DH-NI':'bsa_nir','AL-BH-BB':'wsa_shwv','AL-BH-VI':'wsa_vis','AL-BH-NI':'wsa_nir','AL-DH-BB-ERR':'bsa_shwv_err','AL-DH-VI-ERR':'bsa_vis_err','AL-DH-NI-ERR':'bsa_nir_err','AL-BH-BB-ERR':'wsa_shwv_err','AL-BH-VI-ERR':'wsa_vis_err','AL-BH-NI-ERR':'wsa_nir_err','AL-BH-QFLAG':'brdf_alb_qual_flag','AL-DH-QFLAG':'brdf_alb_qual_flag','NMOD':'nmod'},
                        'vgt_alb_c3s':{},
                        'eps_alb':{'AL-VI-DH':'bsa_vis','AL-BB-DH':'bsa_shwv','AL-NI-DH':'bsa_nir','AL-VI-BH':'wsa_vis','AL-BB-BH':'wsa_shwv','AL-NI-BH':'wsa_nir','AL-VI-DH-ERR':'bsa_vis_err','AL-BB-DH-ERR':'bsa_shwv_err','AL-NI-DH-ERR':'bsa_nir_err','AL-VI-BH-ERR':'wsa_vis_err','AL-BB-BH-ERR':'wsa_shwv_err','AL-NI-BH-ERR':'wsa_nir_err','Q-Flag':'brdf_alb_qual_flag','Z_Age':'age_info'},
                        'msg_alb_bb':{'AL-VI-DH':'bsa_vis','AL-BB-DH':'bsa_shwv','AL-NI-DH':'bsa_nir','AL-BB-BH':'wsa_shwv','AL-VI-DH-ERR':'bsa_vis_err','AL-BB-DH-ERR':'bsa_shwv_err','AL-NI-DH-ERR':'bsa_nir_err','AL-BB-BH-ERR':'wsa_shwv_err','Q-Flag':'brdf_alb_qual_flag','Z_Age':'age_info'},
                        'msg_alb_sp':{'AL-SP-DH':'bsa_sp','AL-SP-BH':'wsa_sp','AL-SP-DH-ERR':'bsa_sp_err','AL-SP-BH-ERR':'wsa_sp_err','Q-Flag':'brdf_alb_qual_flag','Z_Age':'age_info'},
                        'msg_alb_sp_kernels':{'Q-Flag':'brdf_alb_qual_flag','Z_Age':'age_info','K0':'k0','K1':'k1','K2':'k2'},
                        'msg_alb_sp_cov':{'C00': 'c00','C01': 'c01','C02': 'c02','C11': 'c11','C12': 'c12','C22': 'c22'},
                        'modis6_alb':{},
                        'probav_alb':{}
                      }
        # used to force dtype & minise the memory load when reading the data & storing it in dataset
        self.input_variables_dtypes = {
                        #'vgt_alb_cgls':{'AL-DH-BB':np.float32,'AL-DH-VI':np.float32,'AL-DH-NI':np.float32,'AL-BH-BB':np.float32,'AL-BH-VI':np.float32,'AL-BH-NI':np.float32,'AL-DH-BB-ERR':np.float32,'AL-DH-VI-ERR':np.float32,'AL-DH-NI-ERR':np.float32,'AL-BH-BB-ERR':np.float32,'AL-BH-VI-ERR':np.float32,'AL-BH-NI-ERR':np.float32,'AL-BH-QFLAG':np.int16,'AL-DH-QFLAG':np.int16,'NMOD':np.int8},
                        #'vgt_alb_c3s':{},
                        'eps_alb':{'AL-VI-DH':np.int16,'AL-BB-DH':np.int16,'AL-NI-DH':np.int16,'AL-VI-BH':np.int16,'AL-BB-BH':np.int16,'AL-NI-BH':np.int16,'AL-VI-DH-ERR':np.int16,'AL-BB-DH-ERR':np.int16,'AL-NI-DH-ERR':np.int16,'AL-VI-BH-ERR':np.int16,'AL-BB-BH-ERR':np.int16,'AL-NI-BH-ERR':np.int16,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_bb':{'AL-VI-DH':np.int16,'AL-BB-DH':np.int16,'AL-NI-DH':np.int16,'AL-BB-BH':np.int16,'AL-VI-DH-ERR':np.int16,'AL-BB-DH-ERR':np.int16,'AL-NI-DH-ERR':np.int16,'AL-BB-BH-ERR':np.int16,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_sp':{'AL-SP-DH':np.int16,'AL-SP-BH':np.int16,'AL-SP-DH-ERR':np.int16,'AL-SP-BH-ERR':np.int16,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_sp_kernels':{'Q-Flag':np.int16,'Z_Age':np.int8,'K0':np.int16,'K1':np.int16,'K2':np.int16},
                        'msg_alb_sp_cov':{'C00':np.int16,'C01':np.int16,'C02':np.int16,'C11':np.int16,'C12':np.int16,'C22':np.int16},
                        #'modis6_alb':{},
                        #'probav_alb':{}
                      }
        # used to force dtype & minise the memory load when storing the scaled data in dataset
        # could probably be enough to use float16 for albedo data precision
        self.scaled_variables_dtypes = {
                        'vgt_alb_cgls':{'AL-DH-BB':np.float32,'AL-DH-VI':np.float32,'AL-DH-NI':np.float32,'AL-BH-BB':np.float32,'AL-BH-VI':np.float32,'AL-BH-NI':np.float32,'AL-DH-BB-ERR':np.float32,'AL-DH-VI-ERR':np.float32,'AL-DH-NI-ERR':np.float32,'AL-BH-BB-ERR':np.float32,'AL-BH-VI-ERR':np.float32,'AL-BH-NI-ERR':np.float32,'AL-BH-QFLAG':np.int16,'AL-DH-QFLAG':np.int16,'NMOD':np.int8},
                        'vgt_alb_c3s':{},
                        'eps_alb':{'AL-VI-DH':np.float32,'AL-BB-DH':np.float32,'AL-NI-DH':np.float32,'AL-VI-BH':np.float32,'AL-BB-BH':np.float32,'AL-NI-BH':np.float32,'AL-VI-DH-ERR':np.float32,'AL-BB-DH-ERR':np.float32,'AL-NI-DH-ERR':np.float32,'AL-VI-BH-ERR':np.float32,'AL-BB-BH-ERR':np.float32,'AL-NI-BH-ERR':np.float32,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_bb':{'AL-VI-DH':np.float32,'AL-BB-DH':np.float32,'AL-NI-DH':np.float32,'AL-BB-BH':np.float32,'AL-VI-DH-ERR':np.float32,'AL-BB-DH-ERR':np.float32,'AL-NI-DH-ERR':np.float32,'AL-BB-BH-ERR':np.float32,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_sp':{'AL-SP-DH':np.float32,'AL-SP-BH':np.float32,'AL-SP-DH-ERR':np.float32,'AL-SP-BH-ERR':np.float32,'Q-Flag':np.int16,'Z_Age':np.int8},
                        'msg_alb_sp_kernels':{'Q-Flag':np.int16,'Z_Age':np.int8,'K0':np.float32,'K1':np.float32,'K2':np.float32},
                        'msg_alb_sp_cov':{'C00':np.float32,'C01':np.float32,'C02':np.float32,'C11':np.float32,'C12':np.float32,'C22':np.float32},
                        'modis6_alb':{},
                        'probav_alb':{}
                      }