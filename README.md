# ESO-archive
When I was an Astronomy PhD student, one of my duties was to prepare an ETL pipeline for extracting data from Astronomical large surveys. 

This script queries and extracts raw and/or processed spectroscopic metadata and the associated FITS files. It creates new sets of metadata for pre-processing and downloads the required files from http://archive.eso.org/cms.html. 

The examples given in http://archive.eso.org/programmatic/ are not practical but they are the inspiration to these codes. 

The ESO archive is not organized in an efficient way to match processed and raw data and ERDs do not exit. As an astronomer you can always use TOPCAT to do the same query of metadata and understand the basics of schema of the database but you can't download the FITS files.  

These scripts will do that for you and download the raw files which don't have any processed files.

Follow the steps in the jupyter notebook and later you can automatically run the script from the normal python script.  
