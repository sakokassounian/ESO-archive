#Imported external functions
from tools.table_extraction import eso_import, decode_pandas, search_simbad_names, \
get_simbad_names, query_gaia_id, download_raw, download_process, cross_match

#required packages for this code
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from astropy.io import fits
from astropy.table import Table
import numpy as np
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime

#%%#---------------------------- USER SECTION ---------------------------------------------------

"""
Check instrument documentation first for processed and raw data guidelines.
'eso_import' searches in VO_identifier_pro, waveband_pro, waveband_raw
to do further searches using the 'group_by_dict' function to generate new
columns which can be used for filtering. 

Be carefull when cleaning and downloading for more than one instrument. 
If both instruments don't have similar parameters many rows might be 
cleaned out and excluded. 

---> Before downloading visualize the following variables:
    - top
    - accepted_pro
    - accepted_raw_unprocessed
"""



#################################
### DATA GATHERING PARAMETERS ###
#################################

#Select search type
search_type = 'by_name' #'by_coordinates' or 'by_name'

#Target name to be later used
TARGET_SESAME_NAME = "M 67"

#specific target coordinate
target_ra = 132.841544
target_dec = 11.777780


# square window search size in degrees 
RESOLVE_WINDOW =  1/60  

#select instrument name which matches the database
INSTRUMENT_NAMES = ['UVES']

#choose all or as many as you want: options are -> ['ACQUISITION','CALIB','SCIENCE']
#Note that if you choose only 'CALIB','ACQUISITION', or both, the script will
#will not retreive the processed data and cross match it with the raw files  
PRODUCT_TYPE_RAW = ['SCIENCE']


#observation dates
"""
If you pass 'None' or 'False' to the observation dates, all of the dates will 
be returned for the data and they will not be filtered. Use the same format
as the dates below
"""

#obs_night_beg = '1800-01-01T00:00:00.0000'
#obs_night_end = '3000-12-31T12:59:59.9999'

obs_night_beg = None
obs_night_end = None

#cross match search radius
rad = 5 #arcsec

#save_file_anotation
note = str(datetime.now()) 
#when the code is ran the summaries are saved by the running datetime
#the user can also change it

##################################
### DATA FILTERING  PARAMETERS ###
##################################

#generate timestamps based on dates passed
obs_night_beg = pd.to_datetime(obs_night_beg,format='%Y-%m-%dT%H:%M:%S.%f')
obs_night_end = pd.to_datetime(obs_night_end,format='%Y-%m-%dT%H:%M:%S.%f')


#the dictionary that will create columns based on these list of strings 
PRODUCT_TYPE_PROCESS = {'reduced':['SFLX','SRED'],
                        'stacked':['SOBF','SOBR']}

#the dictionary that will create columns based on these list of strings 
WAVEBAND={'red':['FLAMES','DIC1R','DIC2R','_RED'],
          'blue': ['DIC1B','DIC2B','BLUE','_BLU']}

activate_post_cleaning = True

#selection criteria based on the above created columns
process_product_selection = ['reduced','stacked']
wave_selection = ['red']


#Items for extracting data from fits files
fits_items = ['ESO INS GRAT1 WLEN',
              'ESO INS GRAT2 WLEN', 
              'ESO DET WIN1 BINX' ,
              'ESO DET WIN1 BINY' ,
              'ESO INS SLIT3 WID' ,
              'ESO INS SLIT1 WID']
mos_table_name = 'OzPoz_table'
mos_type =  'U                   ' #they have a lot of spaces 

##################################
### DATA SAVING  PARAMETERS ###
##################################

#creating saving directory
save_directory_global = "/home/sako/Desktop/"+ TARGET_SESAME_NAME
save_directory_raw = save_directory_global +"/Raw/"
save_directory_process = save_directory_global +"/Processed/"

if os.path.exists(save_directory_global): 
    print(f'[INFO]:{save_directory_global} exists')
else: 
    print(f'[INFO]:creating {save_directory_global}')
    os.mkdir(save_directory_global)
    os.mkdir(save_directory_process)
    os.mkdir(save_directory_raw)
    
 

#%% ------------------------- IMPLEMENTATION ---------------------------------

#####################
### Import data  ####
#####################
if search_type == 'by_name': coordinates = SkyCoord.from_name(TARGET_SESAME_NAME)
if search_type =='by_coordinates': coordinates = SkyCoord(target_ra*u.deg, target_dec*u.deg, frame='icrs')

complete_query=eso_import(coordinates,
               RESOLVE_WINDOW,
               instrument_name = INSTRUMENT_NAMES,
               product_type_raw = PRODUCT_TYPE_RAW,
               waveband = WAVEBAND,
               product_type_process = PRODUCT_TYPE_PROCESS)

if len(complete_query) == 0: 
    print("[INFO]: Exit code")
    sys.exit()

if (len(PRODUCT_TYPE_RAW)>0) and ('SCIENCE' not in PRODUCT_TYPE_RAW) :
    complete_query = complete_query[complete_query.product_type_raw.notna()]
    



        
#create ra & dec coupling
unique_stars = pd.DataFrame({'RA':list(complete_query.RA_raw)+ 
                 list(complete_query.RA_pro),'DEC':list(complete_query.DEC_raw)+  
                 list(complete_query.DEC_pro)})
    
unique_stars.drop_duplicates(subset=['RA','DEC'],inplace=True)
unique_stars['cluster']= len(unique_stars)*[TARGET_SESAME_NAME]


plt.figure(figsize=(10,10))
plt.plot(complete_query.RA_raw, complete_query.DEC_raw,'r*', markersize=20,label='raw')
plt.plot(complete_query.RA_pro, complete_query.DEC_pro,'b.',label='processed')
plt.legend()
plt.xlabel('RA',fontsize=15)
plt.ylabel('DEC',fontsize=15)
plt.title(f"N = {len(unique_stars)} unique products")
plt.show()



#%%

#########################
### CLEAN UP for UVES ###
#########################

#cleaning the retrieved data 
if activate_post_cleaning:
         

    """ As the PK values are generated based on the observation time for each raw
        file, as a result we can recreate the values using the string pk values
        and convert them to datetime. ADQL query doesn't support the SQL 'convert'
        funtion, as a result we use pandas to apply the filtering"""

    #generate timestamp using the PK values
    obs_datetime=[]
    for i in complete_query.pk:
        temp = i.split(".")
        obs_datetime.append(temp[1]+'.'+temp[2])    
    complete_query["observed_datetime"] =pd.to_datetime(obs_datetime,format='%Y-%m-%dT%H:%M:%S.%f')    
        
    
    if obs_night_beg or obs_night_end:
        top = complete_query[ (complete_query.observed_datetime >= obs_night_beg) & (complete_query.observed_datetime<=obs_night_end) ]   
    else: top = complete_query.copy()    


    #clean processed data
    accepted_pro = top[(top.VO_identifier_pro.notna())]
    accepted_pro = accepted_pro[accepted_pro.product_type_pro.isin(process_product_selection)] #select processed product types
    accepted_pro = accepted_pro[accepted_pro.waveband_pro.isin(wave_selection)] #select waveband arm
    accepted_pro['processing'] = len(accepted_pro)*['process']
    accepted_pro.index=range(len(accepted_pro))
    
    #clean raw data
    accepted_raw = top[top.VO_identifier_raw.notna() ]
    accepted_raw = accepted_raw[accepted_raw.waveband_raw.isin(wave_selection)] #select waveband arm
    accepted_raw['processing'] = len(accepted_raw)*['raw']



    #creating sets
    pro = set(accepted_pro.pk)
    raw =  set(accepted_raw.pk)
    raw_unprocess = raw - pro.intersection(raw) #  A - A and B
    

    #select the raw files which are not processed
    accepted_raw_unprocess = top[top.pk.isin(list(raw_unprocess)) ]
    accepted_raw_unprocess['processing'] = len(accepted_raw_unprocess)*['raw']
    accepted_raw_unprocess.index = range(len(accepted_raw_unprocess))
    
    
    
 
#%%

#########################
### DOWNLOAD RAW DATA ###
#########################    

accepted_raw_unprocess=accepted_raw_unprocess.iloc[:4]

print("[INFO]: Downloading Raw data")
counter=1
status_raw=[]
for i in range(len(accepted_raw_unprocess)):
    url_fits = accepted_raw_unprocess.access_url_raw.iloc[i]
    try:
        name=download_raw(url_fits,save_directory_raw)
        print(f"[INFO]: Downloading ->{counter}/{len(accepted_raw_unprocess)} -> SUCCESS" )
        status_raw.append(name)

    except KeyboardInterrupt:
        sys.exit()
        
    except: 
        print(f"[INFO]: Downloading ->{counter}/{len(accepted_raw_unprocess)} -> FAILED: No access"  ) 
        status_raw.append("Not Downloaded")
        
    counter = counter + 1 
    
    

accepted_raw_unprocess['local_file_name']= status_raw


#RAW DATA
info={ k: []  for k in fits_items}
for i  in accepted_raw_unprocess.local_file_name:
    if i!="Not Downloaded":
        file = fits.open(i)
        header=file[0].header
        for j in fits_items:
            try:
                temp = info[j]
                temp.append(header[j])
                info[j] = temp
            
            except:
                temp = info[j]
                temp.append(np.nan)
                info[j] = temp

            
        
for i in fits_items:          
    accepted_raw_unprocess[i] = info[i]  



#%% --- save raw data --- 

accepted_raw_unprocess.to_csv(save_directory_global+'/raw_summary_'+note+'.csv')


#%%

###############################
### DOWNLOAD PROCESSED DATA ###
###############################    
print("[INFO]: Downloading Processed data")
status_pro=[]
counter=1
for i in range(len(accepted_pro)):
    id_fits = accepted_pro.access_url_pro.iloc[i].split("=")[1]
    try:
        name=download_process(id_fits,save_directory_process)
        print(f"[INFO]: Downloading ->{counter}/{len(accepted_pro)} -> SUCCESS" ) 
        status_pro.append(name)
    except: 
        print(f"[INFO]: Downloading ->{counter}/{len(accepted_pro)} -> FAILED: No access"  ) 
        status_pro.append("Not Downloaded")
        
    counter = counter + 1    

print("\n")

accepted_pro['local_file_name']= status_pro
if accepted_pro.shape[0] ==0 :  print("NO processed data")


#Processed data
info={ k: []  for k in fits_items}
for i  in accepted_pro.local_file_name:
    if i!="Not Downloaded":
        file = fits.open(i)
        header=file[0].header
        for j in fits_items:
            try:
                temp = info[j]
                temp.append(header[j])
                info[j] = temp
            
            except:
                temp = info[j]
                temp.append(np.nan)
                info[j] = temp

            
for i in fits_items:          
    accepted_pro[i] = info[i]  

#%% --- save processed data --

accepted_pro.to_csv(save_directory_global+'/process_summary_'+note+'.csv')


#%%
#####################################
### MOS fits file data extraction ###
#####################################
mos = accepted_raw_unprocess[accepted_raw_unprocess.technology_raw=="MOS"]
results=pd.DataFrame()

if len(mos)>0:
    for index,name in enumerate(mos.local_file_name):
        
        try: 
            fits_file = fits.open(name)
            t = Table.read(fits_file,hdu=mos_table_name,format='fits') 
            t = decode_pandas(t.to_pandas())
            t ['local_file_name'] = len(t)*[name]   
            t ['pk'] = len(t)*[mos.pk.iloc[index]]    
            results=results.append(t)
            print(f'{name} -> extraction successful')
        except:
            print(f'{name} -> FAILED')
                    
    
    accepted_mos = results[results.TYPE == mos_type]
    accepted_mos.to_csv(save_directory_global+'/MOS_summary_'+note+'.csv')



else: 
    accepted_mos = pd.DataFrame(data=None,columns=['pk','RA','DEC','local_file_name'])
    print("\n")
    print("[INFO]: No MOS files retrieved")




#%%

##########################################
### Cross match with external catalogs ###
##########################################


#cleaning data to be sent for cross match
mos_clean = accepted_mos[['pk','RA','DEC','local_file_name']]
mos_clean['processing'] = len(mos_clean)*['raw_MOS']
pro_clean = accepted_pro[['pk','RA_pro','DEC_pro','processing','local_file_name']].rename({'RA_pro':'RA','DEC_pro':'DEC'},axis=1)
raw_clean = accepted_raw_unprocess[['pk','RA_raw','DEC_raw','processing','local_file_name']]\
    [accepted_raw_unprocess.technology_raw!="MOS"].rename({'RA_raw':'RA','DEC_raw':'DEC'},axis=1)


data = mos_clean.append(pro_clean)
data = data.append(raw_clean) 

#%%
if len(data) == len(data.drop_duplicates())  : print('[INFO]: No duplicates')
else: print('[INFO]: Cautions! ->  Duplicates exits')


#get galactic coordinates
coords = SkyCoord(data.RA*u.deg, data.DEC*u.deg, frame='icrs')
coords = coords.galactic
data['l_eso'] = coords.l.to_value()
data['b_eso'] = coords.b.to_value()


#Get Simbad catalog
#TODO: Some of the processed stars seem to be off target -> increased proximity"""
print('[INFO]: Cross match with Simbad')
df_sim=cross_match(data,'simbad','RA','DEC',rad).rename({'ra':'RA_simbad','dec':'DEC_simbad','angDist':'angDist_simbad'},axis=1) 

#Exclude non-star values
df_sim[df_sim.other_types.str.contains("*",regex=False)]


#gather simbad names
df_sim_names = df_sim[['main_id']].drop_duplicates().dropna()

names_all = get_simbad_names(df_sim_names.main_id)
names_HD = search_simbad_names(names_all,'HD')
names_HR = search_simbad_names(names_all,'HR')
names_2MASS = search_simbad_names(names_all,'2MASS')

df_sim_names['HD_names']=[ names_HD[i]  for i in df_sim_names.main_id  ]
df_sim_names['HR_names']=[names_HR[i]  for i in df_sim_names.main_id  ]
df_sim_names['2MASS_names']=[names_2MASS[i]  for i in df_sim_names.main_id  ]

df_sim = pd.merge(df_sim,df_sim_names,how='outer',on=['main_id'])
df_sim = pd.merge(df_sim,data,how='outer',on=['pk','l_eso','b_eso','RA','DEC','local_file_name','processing'])


#gather UBVRI
Simbad.add_votable_fields('flux(U)','flux(I)')
df_UBVRI = Simbad.query_objects(list(df_sim.main_id.dropna().unique()))
df_UBVRI = decode_pandas(df_UBVRI.to_pandas())[["MAIN_ID","FLUX_U","FLUX_I"]].rename({"MAIN_ID":"main_id",'FLUX_U':'U','FLUX_I':'I'},axis=1)
df_sim = pd.merge(df_sim,df_UBVRI,how='outer',on=['main_id'])


df_sim.to_csv(save_directory_global +f"/simbad_{rad}_arc.csv")


#%%----------2MASS -------------
print('[INFO]: Cross match with 2MASS')
df_2mass=cross_match(data,'vizier:II/246/out','RA','DEC',rad)       
df_2mass = pd.merge(df_2mass,data,how='outer',on=['pk','l_eso','b_eso','RA','DEC','local_file_name','processing']) 

df_2mass.to_csv(save_directory_global +f"/2mass_{rad}_arc.csv")

          


#%%----------GAIA edr3 -------------
print('[INFO]: Cross match with Gaia EDR3')
df_gaia_pre=cross_match(data,'vizier:I/350/gaiaedr3','RA','DEC', rad)       
df_gaia_pre = pd.merge(df_gaia_pre,data,how='outer',on=['pk','l_eso','b_eso','RA','DEC','local_file_name','processing']) 
df_gaia_pre = df_gaia_pre[['angDist','pk','l_eso','b_eso','RA','DEC','local_file_name','processing','source_id']]
df_gaia_complete=query_gaia_id(list(df_gaia_pre.source_id.dropna().astype(int).unique()))

df_gaia = pd.merge(df_gaia_complete,df_gaia_pre,how='outer',on=["source_id"])


df_gaia.to_csv(save_directory_global +f"/gaia_{rad}_arc.csv")




#%%
"""
maxer=0
maxer_id=None
for index,value in df_sim.groupby(["pk","RA","DEC"]):
    if len(value)>maxer : 
        maxer = len(value)
        maxer_id = index
        print(value[["angDist_simbad","pk","RA","DEC","main_type"]].to_string())
        jack=input("press Enter")
"""          


#%% ___ save_result in csvs ___

#fix mistakes

"""
#For raw
raw=glob.glob(save_directory_raw+'*')
file_names=[]
for pk in accepted_raw_unprocess.pk:
    for i in raw:
        if pk in i : file_names.append(i)

accepted_raw_unprocess.local_file_name=file_names + (len(accepted_raw_unprocess)-len(file_names))*["Not Downloaded"]

#For process
process=glob.glob(save_directory_process+'*')
file_names=[]
for pk in accepted_pro.access_url_pro:
    for i in process:
        if pk.split("=")[-1] in i : file_names.append(i)

accepted_pro.local_file_name=file_names + (len(accepted_pro)-len(file_names))*["Not Downloaded"]

"""

