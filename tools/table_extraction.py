#!/usr/bin/env python
# Version: 2019-02-13


from pyvo.dal import tap
from astropy.io import fits
from astropy.table import Table
from astropy import units as u

import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import numpy as np
from math import sqrt, ceil


from pathlib import Path
import unlzw3
import os
import warnings
import urllib
from astroquery.xmatch import XMatch

#To call all available data
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
Simbad.ROW_LIMIT= -1 
Vizier.ROW_LIMIT= -1 
Gaia.ROW_LIMIT= -1 

#___________________ Suplimentary functions____________________________________


def unzip(zip_directory):
        
    uncompressed_data = unlzw3.unlzw(Path(zip_directory))
    bytes_ = uncompressed_data
    f = open(zip_directory+".fits", "wb")
    f.write(bytes_)
    f.close()
    
    
def pd_to_table(df):
    t = Table()
    for i in df.columns:
        t[i] = df[i]
    return t


def unpack_sql(arr):
    start="("
    for i in range(len(arr)):
        start=start+ f"'{arr[i]}',"
    end = start[:-1] + ")"
    return end

def unpack_sql_num(arr):
    start="("
    for i in range(len(arr)):
        start=start+ f"{arr[i]},"
    end = start[:-1] + ")"
    return end
    

        
def standardize(arr):
    return (arr - np.nanmean(arr))/np.nanstd(arr)    


def normalize(arr):
    return (arr - np.nanmin(arr))/(np.nanmax(arr) - np.nanmin(arr))    


def decode_pandas(dataframe):
    
    """
    Alot of data retrieved from astropy comes in 'Table' form and the strings
    are in a byte structure; e.g.: x = b'hello world' 
    
    To visualize the data in pandas and be able to apply string operations,
    this function decodes the strings into a UTF-8 format
    """
    
    df = dataframe.copy()
    for i in df.columns:
        if type(df[i].iloc[0]) == bytes : df[i] = df[i].str.decode('UTF-8')
    return df 



def get_datettime(time_string):

    try: 
        temp = time_string.split('-')
        a=temp[0][-4:]
        b=temp[1]
        c=temp[2][:15]
        return a +'-'+ b + '-' + c 
    except: return "" 
    
    
    
def group_by_dict(arr,dic):
    """
    Parameters
    ----------
    arr : list or array
        The array which contains a set of elements in which we will aim to group
        based on a specific mapping based by the "dic" dictionary

        
    dic : dict
        The dictionary which contains the mapping criteria
        
    Returns
    -------
    a list which includes the groupping results based on the passed parameters


    Example
    -------
    arr = ["Sarkis_Kassounian","James_Gordon","Jhon_Smith","Michael_Jordan","William_Jordan"]
    dic = {"family 1": ["Smith","Gordon"],
           "family 2":["Jordan"]}

    >>> search_dict(arr,dic)
    >>> [nan, 'family 1', 'family 1', 'family 2', 'family 2']

    """
    
    arr_ret=[]
    for arr_item in arr:
        temp=0
        for dic_key in dic.keys():
            for value in dic[dic_key]:
                if (not pd.isna(arr_item)) and (value in arr_item):    
                    temp=dic_key
                    arr_ret.append(dic_key)
                if temp != 0 : break
            
        if temp==0: arr_ret.append(np.nan)
    
    return arr_ret    
    


def download_process(dp_id,saving_directory):    
    download_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?%s&eso_download=file" % dp_id    
    now = str(datetime.now()).split('.')
    new_file = saving_directory+ now[0]+'-'+now[1]+"_"+ dp_id +'.fits'   
    urllib.request.urlretrieve(download_url,new_file)    
    return new_file

def download_raw(url_raw,saving_directory):    
    now = str(datetime.now()).split('.')
    dp_id = url_raw.split("/")[-1]
    new_file = saving_directory+ now[0]+'-'+now[1]+"_"+ dp_id
    urllib.request.urlretrieve(url_raw,new_file)    
    unzip(new_file)
    os.remove(new_file)
    return new_file+'.fits'


def cross_match(pandas_df,catalog_name,ra_str,dec_str,proximity):
    tb = pd_to_table(pandas_df)
    table = XMatch.query(cat1=tb,
                      cat2=catalog_name,
                      max_distance= proximity * u.arcsec, colRA1=ra_str,
                      colDec1=dec_str)
    
    warnings.simplefilter("ignore")
    return_df=decode_pandas(table.to_pandas())     

    return return_df 


def get_simbad_names(simbad_main_ids):
    #get all possible values of names in simbad database
    name_dic={}
    for name in simbad_main_ids:
        query = Simbad.query_objectids(name)
        name_list = list(query['ID'])
        name_dic[name] = name_list
        
    return name_dic    


def search_simbad_names(simbad_names_dic,prefix):
    temp={}
    for name in simbad_names_dic.keys():
        for name_list in simbad_names_dic[name]:
            if prefix in name_list.split(" ") : temp[name] = name_list
            
        try:
            temp[name]
        except:
            temp[name] = np.nan
    
    return temp




#__________________________ Main functions____________________________________


def get_NN_coord(coo_1,coo_2,K):
    """

    Parameters
    ----------
    coo_1 : list or an array which contains a group of tuples
        The tuples correspond to the RA and DEC of each star and will be
        cross matched with that of the coo_2 values.
    coo_2 : list or an array which contains a group of tuples
        This contains the list of tuples which we will cross match
    K : Integer
        K is the number of nearest neighbours that will be created.

    Returns
    -------
    df : Pandas dataframe
        It will contain the results of the nearest neighbours:
            - Euclidean distance (d)
            - nearest RA 
            - nearest DEC
            - idx is the index in coo_2 nearest to an item in coo_1

    """
    df=pd.DataFrame(data=None)
    
    for i in range(len(coo_1)):        

        distances=[]
        for j in range(len(coo_2)):
            d = np.sqrt((coo_1[i][0] - coo_2[j][0])**2 + (coo_1[i][1] - coo_2[j][1])**2)
            distances.append(d)
        
        top=np.argsort(distances)

        items={f'd{nn}':distances[top[nn]] for nn in range(K)}
        items.update({f'RA{nn}': coo_2[top[nn]][0] for nn in range(K) })
        items.update({f'DEC{nn}': coo_2[top[nn]][1] for nn in range(K) })
        
        items.update({f'idx{nn}':top[nn] for nn in range(K)})    
        
        df=df.append(pd.DataFrame(data=items,index=[coo_1[i]]))
        
    return df   

        
    
def extract_eso_data(position,window_size,table,instrument = None,data_category=None,max_rows=30000000):
    """
    Description
    -----------
    This function sends queries in the eso archive base using ADQL, to retrieve
    tables about raw or processed data. 
    
    
    Parameters
    ----------
    position : Skycoord object
        'Skycoordinate' object which contains the right accension and declination.

    window_size : Float
        window size of the search square. Standard ADQL queries has a parameter which
        generates a circle. ESO archive doesn't support that circle generation. 
        As a result the query is done as a square.
        
    table : str
        The name of the table found in the ESo-archive database which will be
        queried.

    instrument : list of str, optional
        a list of strings which includes the names of the instruments which we 
        are intrested to query for in the database. 
        Note: Make sure how the names are in database. There exists some
        abbreviations for instrument names which should be used.
        The default is None. As a result all of the available values are
        retrieved. 

    data_category : list of str, optional
        Some instruments have specific product or data types. Dealing with it
        has the same logic as the 'instrument' variable above.        
        The default is None.  As a result all of the available values are
        retrieved. 

    max_rows : int, optional
        Number of maximum retrieved rows. The default is 30 million, since 
        that is the largest number of existing values. Test it by running the
        follwing query which counts how many rows you have:
            
                SELECT COUNT(*)
                from required_table
            
        
    Returns
    -------
    res_df : pandas dataframe
        The function will return a pandas dataframe or an empty set of the 
        retrived table after decoding from the byte strings to UTF-8, using the
        decode_pandas() function. 
    """

    ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
    tapobs = tap.TAPService(ESO_TAP_OBS)
     

    ra_min = position.ra.to_value() - window_size
    ra_max = position.ra.to_value() + window_size 
    dec_min = position.dec.to_value() - window_size
    dec_max = position.dec.to_value() + window_size

    #select ra and dec
    if table == 'dbo.ssa' : 
        ra_col = 's_ra'
        dec_col = 's_dec'
        instrument_col = 'COLLECTION'
        
        
    #select ra and dec
    if table== 'ivoa.ObsCore': 
        ra_col = 's_ra'
        dec_col = 's_dec'
        instrument_col = 'instrument_name'
        
    
    
    if table=='dbo.raw' : 
        ra_col = 'ra'
        dec_col = 'dec'
        instrument_col = 'instrument'
        data_category_col = 'dp_cat'
    
    
    
    #top = "TOP %d" % (3)
    query=f"""SELECT *
    from {table}
    where {ra_col}  between {ra_min} and {ra_max}
    and {dec_col} between {dec_min} and {dec_max}
    """
     
    #if we add dp_cat = 'SCIENCE' we only get science fits files 
    
        
    if instrument:    query=query + f" and  {instrument_col}   in {unpack_sql(instrument)}"                    
    if data_category: query=query + f" and {data_category_col} in {unpack_sql(data_category)}"    
    if table == 'ivoa.ObsCore' and bool(instrument): query=query + f" and obs_collection in {unpack_sql(instrument)}" 
            


    #Apply query
    try: 
        res = tapobs.search(query=query,maxrec=max_rows)
        res_table = res.to_table()
    except: 
        print("[INFO]: No results")        
        res_df = []
        res_table = [] 
    
    
    if len(res_table) != 0: 
        warnings.simplefilter("ignore")
        res_df=decode_pandas(res_table.to_pandas())
    
    else: 
        print("[INFO]: No results")
        res_df=[]    
    
    
    return res_df





def eso_import(target_pos,search_window,
               instrument_name=[],product_type_raw=[],
               waveband={},product_type_process={}):
    """
    Disclaimer:
    -----------    
    This function is optimzed for specific types of instruments.
    It is very important to note that the eso-archive raw and processed database
    doesn't conform to a special database management system.
    
    Therefore if the results are not satisfactory, it is advised to 
    first download the entire instrument database without any conditions and columns
    renaming and redit based on your personal requirement,
    visualize the database without any filters, and then change the function
    scripts based on the values of columns if deemed neccesary. 
    
    
    Parameters
    ----------
    target_coord : string (obligatory)
         'Skycoordinate' object which contains the right accension and declination.

    search_radius : float
        The value of the search window size in degrees.

    instrument_name : list or array (obligatory)
        The list of instruments which the user should pass so that the query would
        retrieve them. If nothing is passed the query retrieves all available 
        instrument data. It is passed as a list since it is sent to the query
        to optimize the query time. 
        Ex: instrument_name = ["UVES","GIRAFFE"]

    waveband : dict (optional)
        The user passes a mapping dictionary of strings which represent the arm
        colors of the instruments. Even wavelenghth values which later will appear in 
        the final dataframe and be used for cleaning the pandas dataframe. 

        Ex: waveband = {'RED':['FLAMES','DIC1R','DIC2R','RED'],
                        'BLU': ['DIC1B','DIC2B','BLUE']}
        
        This does a search in the 'obs_creator_did' columns of the processed data,
        hence it is advised to just use it for arm colors or any other existing
        exact numerical value matching. 

        Note: To get the information of instruments and their indications
        visit http://www.eso.org/sci/observing/phase3/data_streams.html ,
        https://www.eso.org/observing/dfo/quality/UVES/pipeline/settings.html 
        (for UVES), and other associated documents that can describe the data
        and the insturments mentioned in the identifiers and names of columns. 
        
    product_type_process : dict (optional)
        The user should pass a mapping dictionary of strings which represent
        the set of product 
        
        s to be used. The idea here is that these products 
        are selected from a set of strings and not cleaned out from the query.
        This gives you the ability to download the required data and manually 
        clean it later on. You can use this column for any other additional
        quantification that you want regarding search for data type. 
        This does a search in the 'obs_creator_did' columns of the processed data. 

        Ex:  product_type_process = {'reduced': ['SRED', 'SFLX']}
        
        Note: To get the information of instruments and their indications
        visit http://www.eso.org/sci/observing/phase3/data_streams.html

    product_type_raw : list or array (optional)
        The list of instruments which the user should pass so that the query would
        retrieve them. If nothing is passed the query retrieves all available 
        instrument raw data. It is passed as a list since it is sent to the query
        to optimize the query time. 
        
        Ex:  product_type_raw = ["SCIENCE","CALIB"]
        
        
    Returns
    -------
    A pandas dataframe including the results of the query after being editted
    and combined between raw and processed eso archive data. 

    """    
    
            

    
    ################
    ### RAW DATA ### 
    ################     
    #get raw data and clean it
    print(f'[INFO]: Target = \n {target_pos}')
    print(f'[INFO]: Radius = {search_window}')
    
    df_raw = extract_eso_data(target_pos,search_window,'dbo.raw',instrument_name,product_type_raw)
    
    #In case nothing is returned and if raw is empty 
    if len(df_raw)==0: 
        print("[INFO]: No DATA in ESO")
        return []
    
    #generate homogenous raw identifier (primary key)
    pk_raw = [i.split('.')[1:] for i in list(df_raw.dp_id)]
    pk_raw = [ i[0]+'.'+i[1] for i in pk_raw ]  
    df_raw['pk']= df_raw.instrument +"."+pk_raw


    #Generate arm(waveband) identifier for easy cleaning for UVES
    df_raw['waveband_raw']= group_by_dict(list(df_raw.origfile),waveband)

    #select columns
    df_raw = df_raw[['instrument','pk','tpl_start','datalink_url','access_url','dp_id','ra', \
                         'dec','dp_tech',"dp_cat","dp_type","obs_mode",'waveband_raw','origfile']].rename(columns={'ra': 'RA_raw',\
                         'dec':'DEC_raw','access_url':'access_url_raw',\
                        'dp_tech':'technology_raw',\
                        'datalink_url':"VO_identifier_raw",'dp_cat':'product_type_raw',"dp_type":'data_type_raw','instrument':'instrument_raw',"origfile":"name_raw_file"})


                                                                                                                            
    ######################
    ### processed DATA ### 
    ######################  
        
    df_pro = extract_eso_data(target_pos,search_window,'ivoa.ObsCore',instrument_name)
    
    
    
    #create columns of additional info    
    if df_pro.shape[0]>0:
        df_pro['waveband_pro'] = group_by_dict(list(df_pro.obs_creator_did),waveband) 
        df_pro['product_type_pro'] = group_by_dict(list(df_pro.obs_creator_did),product_type_process)
        
        
        pks=[i.split("?")[1] for i in df_pro.obs_creator_did]
        pks_clean = [get_datettime(i) for i in pks]   
        
        pk_pro=[]
        for index,value in enumerate(pks_clean):
            if value =="": pk_pro.append(df_pro.instrument_name.iloc[index]+'.'+ "".join(df_pro.obs_creator_did.iloc[index].split("?")[1]))    
            else:  pk_pro.append(df_pro.instrument_name.iloc[index] +'.' +value) 
        df_pro['pk']= pk_pro

    else:         
        df_pro['waveband_pro']= []
        df_pro['product_type_pro']=[]
        df_pro['pk']= []
        
        
        
        
        
    print(f'[INFO]: RAW retrieved = {len(df_raw)}')
    print(f'[INFO]: PROCESSED retrieved = {len(df_pro)}')
        
    

    
    
    df_pro = df_pro[['pk','obs_creator_did','access_url','s_ra','s_dec', \
                              'snr','obstech','waveband_pro','product_type_pro','instrument_name']].rename(columns={'s_ra': 'RA_pro',\
                             's_dec':'DEC_pro','access_url':'access_url_pro',\
                            'obstech':'technology_pro',\
                            'obs_creator_did':"VO_identifier_pro",'instrument_name':"instrument_pro"})
                                                                

    
                                                                                                                                        
    #merge based on dp_id to know which has a complete set of data
    merged=pd.merge(df_raw,df_pro,on='pk',how='outer')
    merged=merged.sort_values(by='pk')[['tpl_start','pk','snr','instrument_pro','RA_pro','DEC_pro','instrument_raw','RA_raw','DEC_raw',\
            'VO_identifier_raw','VO_identifier_pro','access_url_pro','access_url_raw', \
            'technology_raw','technology_pro',"product_type_raw","product_type_pro",'waveband_pro','waveband_raw',"data_type_raw","obs_mode","name_raw_file"]]
    
    merged.index=range(len(merged))
    
    #merged["observed_datetime"] =pd.to_datetime([i.split(".")[-1] for i in merged.pk],format='%Y-%m-%dT%H:%M:%S.%f')
     
      
    
    
    return merged




#%%

def extract_hst_data(position,window_size,table,max_rows=500000):
    url = "http://vao.stsci.edu/HSCTAP/tapservice.aspx"
    tapobs = tap.TAPService(url)
     

    ra_min = position.ra.to_value() - window_size
    ra_max = position.ra.to_value() + window_size 
    dec_min = position.dec.to_value() - window_size
    dec_max = position.dec.to_value() + window_size
    
    ra_col = 'ra'
    dec_col = 'dec'
    
    
#top = "TOP %d" % (3)
    query=f"""SELECT *
    from {table}
    where {ra_col}  between {ra_min} and {ra_max}
    and {dec_col} between {dec_min} and {dec_max}
    """


    #Apply query
    try: 
        res = tapobs.search(query=query,maxrec=max_rows)
        res_table = res.to_table()
    except: 
        print("[INFO]: No results")        
        res_df = []
        res_table = [] 
    
    
    if len(res_table) != 0: 
        warnings.simplefilter("ignore")
        res_df=decode_pandas(res_table.to_pandas())
    
    else: 
        print("[INFO]: No results")
        res_df=[]    
    
    
    return res_df



#%%

def visualize_fits_spectra_eso(file_directory):
    """
    This functions visualizes the fits files of the processed spectra
    downloaded from the ESO archive. The output is a subplots structure 
    of the data found in the fits file. The x-axis is taken as the wavelength
    which has the title "WAVE". 

    Parameters
    ----------
    file_directory : str
        Directory of the fits file that you want to visualize
    """
    
    df_spectra=pd.DataFrame(data=None)
    fits_file = fits.open(file_directory)
    

    #extract header for use
    #main_header = fits_file[0].header
    spectra_header = fits_file[1].header
    
    #extract spectra 
    spectra = fits_file[1].data[0]
    
    #Build dataframe
    for i in range(len(spectra)):
        df_spectra[spectra_header[f"TTYPE{i+1}"]] = list(spectra[i])
    
    df_spectra.index=df_spectra['WAVE']
    df_spectra.drop(['WAVE'],axis=1,inplace=True)

    dim = ceil(sqrt(len(spectra)))
    df_spectra.plot(subplots=True,layout=(dim,dim),figsize=(20,20),sharex=False)
    plt.show()




def query_gaia_id(ids, table='gaiaedr3.gaia_source'):
    """
    Queries the original gaia database and retrieves all possible data
    
    """
        
    query=f"""SELECT *
              from {table}
              where source_id in {unpack_sql_num(ids)}
           """
                      
    job = Gaia.launch_job_async(query)
    
    warnings.simplefilter("ignore")
    df= decode_pandas(job.get_results().to_pandas())

    return df



def demo_query(table):
    ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
    tapobs = tap.TAPService(ESO_TAP_OBS)
     

    query=f"""SELECT TOP 100 *
    from {table}
    """
  

    #Apply query
    try: 
        res = tapobs.search(query=query)
        res_table = res.to_table()
    except: 
        print("[INFO]: No results")        
        res_df = []
        res_table = [] 
    
    
    if len(res_table) != 0: 
        warnings.simplefilter("ignore")
        res_df=decode_pandas(res_table.to_pandas())
    
    else: 
        print("[INFO]: No results")
        res_df=[]    
	

    return res_df		



#%%



"""

def remove_outlier_2delta(DF,exceptions):
    
    df=DF.copy()
    
    for col in df.columns:
        if col in exceptions: continue
        std1 = df[col].quantile(0.05)
        std2 = df[col].quantile(0.95)

        df = df[(df[col] > std1) & (df[col] < std2)]

    return df    
    

#%%


"""



#%% ___________ Cross match _____________
"""  
#%%astropy cross-matching 
#coo_raw = SkyCoord(raw['ra'], raw['dec'])
#coo_pro = SkyCoord(process['s_ra'], process['s_dec'])
#idx, d2d, d3d = coo_pro.match_to_catalog_sky(coo_raw,nthneighbor=1)
"""

#%%______________K-NN selection______________

"""
K_NN = 5

#K-NN cross matching using the entire database 
coo_pro = list(set([*zip(process['s_ra'], process['s_dec'])])) # remove duplicates of processed data
coo_raw = list(set([*zip(raw['ra'], raw['dec'])])) # remove duplicates of raw data

idx_df=get_NN(coo_pro,coo_raw,K_NN) 
idx_df.where(idx_df[[f'd{i}' for i in range(K_NN)]]<2/3600,inplace=True)


df_excluded = raw_df.copy()

coo_exc=[]
for i in range(len(idx_df)):
    for nn in range(K_NN):
        if not np.isnan(idx_df[f'd{nn}'].iloc[i]):
            coo_exc.append((idx_df[f'RA{nn}'].iloc[i],idx_df[f'DEC{nn}'].iloc[i]))
            
idx_exc=[]            
for ra,dec in coo_exc:            
    temp = raw_df[(raw_df.ra==ra) & (raw_df.dec==dec)]
    idx_exc=idx_exc+list(temp.index)
    
df_excluded.drop(list(set(idx_exc)),inplace=True) 

    
#%%
plt.figure(figsize=(10,10))
#plt.plot(raw_all['ra'],raw_all['dec'],'go',markersize=15,label=f'raw -> all: {len(raw_all)}')
plt.plot(df_excluded['ra'],df_excluded['dec'],'r*',markersize=20,label=f'raw_excluded: {len(df_excluded)}')
plt.plot(process_df['s_ra'],process_df['s_dec'],'b.',markersize=5,label=f'processed: {len(process)}')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.legend()
plt.show()
"""






"""
###########################################################################
##### This section was in the visualization for extract the data ##########
###########################################################################
    if saving_directory:
        import urllib
        print("Downloading the non-proprietary files (%d files)" % (len(res.to_table())))
        track_files=[]        
        n_files = 1
        for row in res_table:
            # Change from previous version: the access_url points to the datalink and not directly to the file download.
            # Given a dp_id of a public file, the link to download it is constructed as follows:
            download_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?%s&eso_download=file" % row["dp_id"].decode()
        
            now = str(datetime.now()).split('.')
            new_file = saving_directory+ now[0]+'-'+now[1]+"_"+ row["dp_id"].decode()
        
        
            urllib.request.urlretrieve(download_url,new_file)
            unzip(new_file)
            
        
            # informing user:
            print(f"[INFO]: Downloading and unzipping {n_files}/{len(res_table)} : {now[0]}-{now[1]}_{row['dp_id'].decode()}")        
            n_files += 1
        res_table['local_file_tracker'] = track_files   
        
        
    
    if visualize_fits:
        zip_files=glob.glob(saving_directory+'*UVES*')
        
        
        unzip(zip_files[0])
        
        
        fits_file = fits.open(zip_files[0]+'.fits')
        im = fits_file[0].data
        
        
        #gamma level change to improve visualization
        plt.imshow(im**0.1,cmap='nipy_spectral')
        plt.show()
        os.remove(zip_files[0]+'.fits')
        
"""        




"""
################################
### GAIA EDR3 data gathering ###                                       
################################


gaia_df = extract_gaia_data(target_pos,search_radius,'gaiaedr3.gaia_source')
    
gaia_select = gaia_df[['ra','dec','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error']]



#%%
#generate nearest neighbours
crossed = get_NN([*zip(to_cross.RA_raw,to_cross.DEC_raw)],[*zip(gaia_select.ra,gaia_select.dec)],1)

gaia_crossed = gaia_select.loc[crossed.idx0]
gaia_crossed['name_id']=list(to_cross.name_id)


descion=pd.merge(final,gaia_crossed,on='name_id',how='inner')

descion.to_csv(f"/home/sarkis/Desktop/Thesis files/git_scripts/thesis_codes/production_codes/fixtures/VO school/{target_name}.csv")
  
#%%___________ Plot catalog coordinates___________

plt.figure(figsize=(10,10))
plt.quiver(descion.ra, descion.dec, descion.pmra, descion.pmdec,scale=200,alpha=1,label='gaia proper motion')
plt.scatter(descion.ra,descion.dec,c=descion.parallax,label=f'{target_name} after ESO/GAIA cuts')
plt.colorbar().set_label(label='GAIA Parallax',size=15)
plt.legend()
plt.xlabel('RA',fontsize=15)
plt.ylabel('DEC',fontsize=15)
plt.tight_layout()
plt.savefig(f"/home/sarkis/Desktop/Thesis files/git_scripts/thesis_codes/production_codes/fixtures/VO school/{target_name}.jpg",format='jpg',dpi=300)
plt.show()

"""










        
