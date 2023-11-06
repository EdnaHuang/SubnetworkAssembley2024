#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import subprocess
from datetime import date
import sys


# In[2]:

tracker_path = '/work2/08623/abbajpai/stampede2/luproject/lulab_data/Tracker'
df = pd.read_excel(os.path.join(tracker_path,"06022022_Job_Run_tracker.xlsx"))


# In[15]:


def create_job(job_file, t, job_directory, matrix_type, matlab_file):
    
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH -J Lulab\n")
        fh.writelines("#SBATCH -p skx-normal\n")
        fh.writelines("#SBATCH -o output_%j.txt\n")
        fh.writelines("#SBATCH -e error_%j.err\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=$USER@iu.edu\n")
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --ntasks-per-node=1\n")
        fh.writelines("#SBATCH --cpus-per-task=48\n")
        fh.writelines(f'#SBATCH --time={t}\n')
        # fh.writelines("#SBATCH --mem=100GB\n")
        fh.writelines("#SBATCH -A TG-BIO220061\n")

        fh.writelines("module load matlab\n")
        fh.writelines(f'matlab -batch "{matlab_file}({matrix_type}, {job_directory})"')
        
    fh.close()    
    #os.system("sbatch %s" %job_file)
    sbatch_command = "sbatch %s" %job_file
    sbatch_response = subprocess.getoutput(sbatch_command) 
    job_id = sbatch_response.split(' ')[-1].strip()
    
    return job_id


# In[16]:


# Job Submission

path = '/work2/08623/abbajpai/stampede2/luproject/community_detection/comm_detect'
job_counter = 0
n = int(sys.argv[1]) # edit to decide number of jobs to launch next 
for i in range(0,4737):
    
    if pd.isna(df.iloc[i,10]) == True: # Job has not been created for this row
        
        job_directory = "'" + df.iloc[i,2] + "'"       
        
           
        # Time based on matrix size (check hyperparameter study)
        mat_size = df.iloc[i,7]
        matrix = df.iloc[i,9]

        # Select model type
        if df.iloc[i,8] == 'Asymmetric':
            matlab_file = "Modol_script_asymmetry"          
        else:
            matlab_file = "Modol_script_uniform" 
        
        # Select time assignment 
        if mat_size >= 3000:
            t = '16:00:00'
        elif mat_size in range(2000,3000):                
            t = '6:30:00'
        elif mat_size in range(1200,2000):    
            t = '3:00:00'
        else:
            t = '1:30:00'   	

	# Select matrix type
        matrix_type = "'" + df.iloc[i,9] + "'"
        
        # create job
        job_file = os.path.join(path,"run_batch_" + str(i) + ".sh")
        job_id = create_job(job_file, t, job_directory, matrix_type, matlab_file)
        df.iloc[i,10] = job_id   
        print(job_id)
        job_counter +=1
            
    if job_counter == n: 
        print(f'{n} jobs submitted')
        break


# In[24]:


df.to_excel(str(date.today()) + "_Job_Run_tracker.xlsx", index = False) # a daily copy

# In[ ]:


# df.to_excel(os.path.join(tracker_path, "06022022_Job_Run_tracker.xlsx"), index = False) # update to main tracker

# get all job details


print("Initiating Tracker Update")


subprocess.getoutput('sacct -p --starttime=2022-06-02 --format=jobid,exit,maxrss,Elapsed,Start,Submit>Job_details.csv') 

job_details = pd.read_csv('Job_details.csv')
print("....")

'''Cleaning the Job details file'''

'Cleaning Job Details file'
job_details[['JobID', 'ExitCode','MaxRSS','Elapsed','Start','Submit']] = job_details["JobID|ExitCode|MaxRSS|Elapsed|Start|Submit|"].str.split('|', 5, expand=True)
job_details.replace('', np.nan, inplace=True)
job_details['Elapsed'] = job_details['Elapsed'].apply(lambda x: x[:8])
job_details['Submit'] = job_details['Submit'].apply(lambda x: x[:19])
job_details[['JobID', 'Type']] = job_details["JobID"].str.split('.', 1, expand=True)


'Additional processing on Job Details '
#Start and Submit date is on different lines'

# batch row
job_details_batch = job_details[job_details.Type =='batch']
job_details_batch = job_details_batch.drop(columns = ['Start', "Submit"])
job_details_batch['JobID'] = job_details_batch['JobID'].map(int)

# non-batch row
job_details_non_batch = job_details[job_details.Type !='batch']
job_details_non_batch = job_details_non_batch.drop(columns = ['Elapsed', "ExitCode", "MaxRSS", 'Type','JobID|ExitCode|MaxRSS|Elapsed|Start|Submit|'])
job_details_non_batch['JobID'] = job_details_non_batch['JobID'].map(int)

# merge
job_details = job_details_batch.merge(job_details_non_batch, on ='JobID', how='inner')

# calculate runtime
job_details = job_details.drop(columns = ["JobID|ExitCode|MaxRSS|Elapsed|Start|Submit|", "Type"])
job_details.reset_index(drop=True, inplace =True)
job_details[['Hrs', 'Mins', 'Seconds']] = job_details['Elapsed'].str.split(':', 3, expand=True)

# Merge with tracker
df.loc[df['job_Id'].notnull(), 'job_Id'] = df.loc[df['job_Id'].notnull(), 'job_Id'].map(int)
df = df.merge(job_details, how='left', left_on ='job_Id', right_on='JobID')

if len(df)>4879:
    print(f'Merged record count: {len(df)}')	
    print('error in merge')


'''Updating the final tracker'''
for index, row in df.iterrows():
    # if job id is assigned on tracker and it is executed
    if (pd.isna(df.iloc[index,10]) == False) & (pd.isna(df.iloc[index,22])==False): 
        
        if df.iloc[index,20]=='0:0': 
            
            # Runtime
            df.iloc[index,11] = df.iloc[index,22]
            if pd.isna(df.iloc[index,21]) == False:
            # Memory conversion to GB from Kilobytes(K)/ Megabytes(M)
                if df.iloc[index,21][-1] =='K':
                    df.iloc[index,12] = np.round(float(df.iloc[index,21][:-1])/(1024*1024),2) 
                else:
                    df.iloc[index,12] = np.round(float(df.iloc[index,21][:-1])/(1024),2) # Megabytes
            else:
                df.iloc[index,12] = 0	
                   
            # Elapsed Time (mins)
            df.iloc[index,13] = int(df.iloc[index,25])*60 + int(df.iloc[index,26]) 
            df.iloc[index,17] = df.iloc[index,23]
            df.iloc[index,18] = df.iloc[index,24]
            
        else:
            print(f'Failed Job ID: {df.iloc[index,10]}')  
                
df = df.drop(columns = ['JobID','ExitCode','MaxRSS','Elapsed','Hrs','Mins','Seconds','Start','Submit'])            
            
print("....")

df.to_excel(os.path.join(tracker_path, "06022022_Job_Run_tracker.xlsx"), index = False) # update to main tracker
df[['job_Id','batch_Id', 'runtime','memory','runtime_mins']].to_csv("Job_Run_quick_view.csv", sep= "|", index = False)

print("Tracker Update Complete")
