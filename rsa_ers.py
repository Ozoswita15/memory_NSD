from ccxt_config import *
from nilearn.image import get_data
import warnings
# import rsa_utils

rois = ['vtc']
encoded = 'FE'
stim = 'CC'
memory_target_ses = 2

def ccx_ers (sub,rois,encoded,stim,memory_target_ses):
    
    con = encoded + '_' + stim
    subj = bids_meta(sub)
    
    for roi in rois:
        
        try:
            mask = get_data(f'{subj.masks}/{roi}.nii.gz') #load in the mask
        except:
            mask = get_data(f'{group_masks}/{roi}.nii.gz')
    
        if encoded == 'FC':
            # loading in what we need 
                        # loading in what we need 
            
            beta_mem_ses2= get_data(f'{subj.betas}/ses-2_task-memory_deriv_bmap.nii.gz')
            try:
                beta_mem_ses3= get_data(f'{subj.betas}/ses-3_task-memory_deriv_bmap.nii.gz')
            except:
                beta_mem_ses3=None
            beta_fear_ses1= get_data(f'{subj.betas}/ses-1_task-fear_deriv_bmap.nii.gz')
            beta_fearext_ses1= get_data(f'{subj.betas}/ses-1_task-fearextinction_deriv_bmap.nii.gz')
                
            events_mem_ses2 = pd.read_csv(path(subj.timing,f'ses-2/func/{subj.fsub}_ses-2_task-memory_events.tsv'),sep='\t')
            try:
                events_mem_ses3= pd.read_csv(path(subj.timing,f'ses-3/func/{subj.fsub}_ses-3_task-memory_events.tsv'),sep='\t')
            except:
                events_mem_ses3=None
            events_fear_ses1= pd.read_csv(path(subj.timing,f'ses-1/func/{subj.fsub}_ses-1_task-fear_events.tsv'),sep='\t')   
            events_fearext_ses1= pd.read_csv(path(subj.timing,f'ses-1/func/{subj.fsub}_ses-1_task-fearextinction_events.tsv'),sep='\t')

            weights_mem_ses2= get_data(f'{subj.rsa}/ses-2_task-memory/{con}_all.nii.gz')
            try:
                weights_mem_ses3= get_data(f'{subj.rsa}/ses-3_task-memory/{con}_all.nii.gz')
            except:
                weights_mem_ses3 = None
            weights_fear_ses1= get_data(f'{subj.rsa}/ses-1_task-fear/{con}_early.nii.gz')
            weights_fearext_ses1= get_data(f'{subj.rsa}/ses-1_task-fearextinction/{con}_all.nii.gz')
            
            # select specific stim in events
            events_mem_ses2 = events_mem_ses2[events_mem_ses2.trial_type == con].reset_index(drop=True)
            try:
                events_mem_ses3= events_mem_ses3[events_mem_ses3.trial_type == con].reset_index(drop=True)
            except:
                events_mem_ses3= None
            events_fear_ses1 = events_fear_ses1[events_fear_ses1.trial_type == con].reset_index(drop=True)
            events_fearext_ses1= events_fearext_ses1[events_fearext_ses1.trial_type == con].reset_index(drop=True)
            
            try:
                events_mem = pd.concat([events_mem_ses2,events_mem_ses3]).reset_index(drop=True)
            except:
                events_mem = events_mem_ses2
            events_fear = pd.concat([events_fear_ses1,events_fearext_ses1]).reset_index(drop=True)
            
            #select specific betas
            data_mem_ses2 = beta_mem_ses2[:,:,:,events_mem_ses2[events_mem_ses2.trial_type == con].index]
            try:
                data_mem_ses3 = beta_mem_ses3[:,:,:,events_mem_ses3[events_mem_ses3.trial_type == con].index]
            except:
                data_mem_ses3 = None
            
            data_fear_ses1 = beta_fear_ses1[:,:,:,events_fear_ses1[events_fear_ses1.trial_type == con].index]
            data_fearext_ses1 = beta_fearext_ses1[:,:,:,events_fearext_ses1[events_fearext_ses1.trial_type == con].index]
            
            # weighting our selected betas
            data_mem_ses2  *= weights_mem_ses2[:, :, :, np.newaxis]
            try:
                data_mem_ses3  *= weights_mem_ses3[:, :, :, np.newaxis]
            except:
                data_mem_ses3 = None  
            data_fear_ses1  *= weights_fear_ses1[:, :, :, np.newaxis]
            data_fearext_ses1  *= weights_fearext_ses1[:, :, :, np.newaxis]

            # get memory into one beta series
            
            events_mem = events_mem.sort_values('cs')
            events_mem_index = events_mem.index
            try:
                data_mem = np.concatenate((data_mem_ses2 ,data_mem_ses3), axis=3)
            except:
                data_mem = data_mem_ses2
            data_mem = data_mem[:, :, :, events_mem_index]
            
            # get fear into one beta series
            events_fear = events_fear.sort_values('cs')
            #events_fear_index = events_fear.index
            
            data_fear = np.concatenate((data_fear_ses1 ,data_fearext_ses1), axis=3)
            #data_ext = data_ext[:, :, :, events_ext_index]

            # select memory session
            
            data_mem = data_mem[:,:,:,events_mem[events_mem.ses == memory_target_ses].index]
            masked_events_mem = events_mem[events_mem.ses == memory_target_ses]
            data_fear = data_fear[:,:,:,events_fear[events_fear['cs'].isin(masked_events_mem['cs'])].index]
        
            # mask
            
            x=apply_mask(mask,data_mem)
            y=apply_mask(mask,data_fear)


        elif encoded == 'FE':
            
            # loading in what we need 
            
            beta_mem_ses2= get_data(f'{subj.betas}/ses-2_task-memory_deriv_bmap.nii.gz')
            try:
                beta_mem_ses3= get_data(f'{subj.betas}/ses-3_task-memory_deriv_bmap.nii.gz')
            except:
                beta_mem_ses3=None
            beta_fearext_ses1= get_data(f'{subj.betas}/ses-1_task-fearextinction_deriv_bmap.nii.gz')
            beta_ext_ses1= get_data(f'{subj.betas}/ses-1_task-extinction_deriv_bmap.nii.gz')
                
            events_mem_ses2 = pd.read_csv(path(subj.timing,f'ses-2/func/{subj.fsub}_ses-2_task-memory_events.tsv'),sep='\t')
            try:
                events_mem_ses3= pd.read_csv(path(subj.timing,f'ses-3/func/{subj.fsub}_ses-3_task-memory_events.tsv'),sep='\t')
            except:
                events_mem_ses3=None
            events_fearext_ses1= pd.read_csv(path(subj.timing,f'ses-1/func/{subj.fsub}_ses-1_task-fearextinction_events.tsv'),sep='\t')
            events_ext_ses1= pd.read_csv(path(subj.timing,f'ses-1/func/{subj.fsub}_ses-1_task-extinction_events.tsv'),sep='\t')   

            weights_mem_ses2= get_data(f'{subj.rsa}/ses-2_task-memory/{con}_all.nii.gz')
            try:
                weights_mem_ses3= get_data(f'{subj.rsa}/ses-3_task-memory/{con}_all.nii.gz')
            except:
                weights_mem_ses3 = None
            weights_fearext_ses1= get_data(f'{subj.rsa}/ses-1_task-fearextinction/{con}_all.nii.gz')
            weights_ext_ses1= get_data(f'{subj.rsa}/ses-1_task-extinction/{con}_late.nii.gz')
            
            # select specific stim in events
            events_mem_ses2 = events_mem_ses2[events_mem_ses2.trial_type == con].reset_index(drop=True)
            try:
                events_mem_ses3= events_mem_ses3[events_mem_ses3.trial_type == con].reset_index(drop=True)
            except:
                events_mem_ses3= None

            events_fearext_ses1= events_fearext_ses1[events_fearext_ses1.trial_type == con].reset_index(drop=True)
            events_ext_ses1 = events_ext_ses1[events_ext_ses1.trial_type == con].reset_index(drop=True)
            
            try:
                events_mem = pd.concat([events_mem_ses2,events_mem_ses3]).reset_index(drop=True)
            except:
                events_mem = events_mem_ses2
            events_ext = pd.concat([events_ext_ses1,events_fearext_ses1]).reset_index(drop=True)
            
            #select specific betas
            data_mem_ses2 = beta_mem_ses2[:,:,:,events_mem_ses2[events_mem_ses2.trial_type == con].index]
            try:
                data_mem_ses3 = beta_mem_ses3[:,:,:,events_mem_ses3[events_mem_ses3.trial_type == con].index]
            except:
                data_mem_ses3 = None

            data_fearext_ses1 = beta_fearext_ses1[:,:,:,events_fearext_ses1[events_fearext_ses1.trial_type == con].index]
            data_ext_ses1 = beta_ext_ses1[:,:,:,events_ext_ses1[events_ext_ses1.trial_type == con].index]
            
            # weighting our selected betas
            data_mem_ses2  *= weights_mem_ses2[:, :, :, np.newaxis]
            try:
                data_mem_ses3  *= weights_mem_ses3[:, :, :, np.newaxis]
            except:
                data_mem_ses3 = None  
            data_fearext_ses1  *= weights_fearext_ses1[:, :, :, np.newaxis]
            data_ext_ses1  *= weights_ext_ses1[:, :, :, np.newaxis]

            # get memory into one beta series
            
            events_mem = events_mem.sort_values('cs')
            events_mem_index = events_mem.index
            try:
                data_mem = np.concatenate((data_mem_ses2 ,data_mem_ses3), axis=3)
            except:
                data_mem = data_mem_ses2
            data_mem = data_mem[:, :, :, events_mem_index]
            
            # get ext into one beta series
            events_ext = events_ext.sort_values('cs')
            #events_ext_index = events_ext.index
            
            data_ext = np.concatenate((data_ext_ses1 ,data_fearext_ses1), axis=3)
            #data_ext = data_ext[:, :, :, events_ext_index]

            # select memory session
            
            data_mem = data_mem[:,:,:,events_mem[events_mem.ses == memory_target_ses].index]
            masked_events_mem = events_mem[events_mem.ses == memory_target_ses]
            data_ext = data_ext[:,:,:,events_ext[events_ext['cs'].isin(masked_events_mem['cs'])].index]
        
            # mask
            
            x=apply_mask(mask,data_mem)
            y=apply_mask(mask,data_ext)
            
        # rsa

        mat = np.arctanh(np.corrcoef(x,y)) #run the correlations
        mat[np.eye(mat.shape[0]).astype(bool)] = np.nan #correct the eye/diagonal to be nan
        mat = mat[x.shape[0]:,:x.shape[0]] #get the lower left square of between comparisons
        mat = np.diagonal(mat) # get diagonal to get ERS
        df = pd.DataFrame(columns = ['sub',"roi",'memory_ses','encoded','trial_type','rsa','type'])
        df.loc[df.shape[0]] = [None,None,None,None,None,None,None]
        df['sub'] = sub
        df['roi'] = roi
        df['memory_ses'] = memory_target_ses
        df['encoded'] = encoded
        df['trial_type'] = stim
        df['rsa'] = mat.mean()
        df['type'] = 'ers'

        df.to_csv(f'{DATA}/rsa_out/{subj.fsub}_{roi}_encoding-{encoded}-{stim}_retrieval-ses-{memory_target_ses}_ers.csv',index=False)
        print(f'done with {subj}!')