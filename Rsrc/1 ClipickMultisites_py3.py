#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# enable debugging
import cgitb
cgitb.enable()
import csv
from io import StringIO # import StringIO if using Python 2
import urllib
import urllib.request # import urllib2 if using Python 2
import urllib.parse
import time

# Needed daily DATA:
# Column 0: Year
# Column 1: Month
# Column 2: Day
# Column 3: Average temperature (degrees C)
# Column 3: Radiation (MJ / sq meter)
# Column 3: Precipitation (mm)

#CONTROL VARIABLES OF THIS ALGORTIHM:
  
with open('Sites_Lon_Lat.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        print(row[1],row[2],row[3])

         
with open('Sites_Lon_Lat.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    next(readCSV)
    for row in readCSV:
        Longitude= row[2]
        Latitude = row[3]
        StartYear 	        = 1990 # Min=1951, Max=2100
        StartMonth 	        = 1
        StartDay 	        = 1
        EndYear             = 2015
        EndMonth 	        = 12
        EndDay 		        = 31
        IPCCAssessmentReport = 4 # either 4 or 5
        Dataset              = 'METO-HC_HadRM3Q0_A1B_HadCM3Q0_DM_25km' # if IPCCAssessmentRepo    
        outFileName   = row[1]+'.csv' # tip you can build the name of the file to be according to the dates extracted
        outFileHandle = open(outFileName, 'w')
        start_time = time.time() # this is facultative, just to calculate timming of retrieval
        pars = {}
        pars['lat']        = Latitude
        pars['lon']       = Longitude
        pars['fmt']       = 'csv' # either csv, htmltable
        pars['tspan']   = 'd'# d=daily; m =monthly
        pars['sd']        = StartDay #'01'
        pars['sm']       = StartMonth #'01'
        pars['sy']         = StartYear
        pars['ed']        = EndDay
        pars['em']       = EndMonth
        pars['ey']        = EndYear
        pars['dts']       = Dataset# Beware of dates for extraction
        pars['ar']         = IPCCAssessmentReport # either 4 or 5
        pars['mod']     = "hisafe" # either yieldsafe or hisafe
        url                    = 'http://www.isa.ulisboa.pt/proj/clipick/climaterequest_fast.php'
        
        url_pars           = urllib.parse.urlencode(pars)
        full_url           = url + '?' + url_pars
        print ("Request made to " + full_url)
        response           = urllib.request.urlopen(full_url) #urllib2.urlopen(full_url) if using Python 2
        the_page           = response.read().decode('utf-8') # just response.read() if using Python 2
 
        f                  = StringIO(the_page)
        reader             = csv.reader(f, delimiter=',')
        
        # CEATE AN ARRAY FROM THE REQUESTED CSV OUTPUT
        result=[]
        for row in reader:
            result.append(row)
        print ("Output is being written in " + outFileName)
        for i in result:
            outFileHandle.write(",".join(i) + "\n")
        outFileHandle.flush()
        outFileHandle.close()
        end_time = time.time()
        print ("Processed in " + str(round(end_time - start_time,0)) + " seconds")
        print ("Output stored in " + outFileName)
        print ("done")
        
# WRITE IT DOWN IN THE OUTPUT FILE
'''
the daily output comes as 
yieldsafe : Day, Month, Year, tas, rss, pr
hisafe    : Day, Month, Year, tasmax,tasmin,hursmax,hursmin,rss,pr,wss
in AR5 datasets, there are no min and max relative humidity (at the time of this deliverable).
therefore, at the moment, hisafe format for AR5 are as follows:
hisafe: : Day, Month, Year, tasmax,tasmin,hurs,evspsbl,rsds,pr,sfcWind

Currently Valid variables are (nomenclature as ENSEMBLES, and CORDEX project):
    "pr"     : precipitation
    "tas"    : Average Temperature
    "tasmin" : Minimum Temperature
    "tasmax" : Maximum temperature
    "rss"    : Radiation
    "evspsbl": Evaporation
    "hurs"   : Relative Humidity
    "hursmax": Maximum Relative humidity
    "hursmin": Minimum Relative humidity
    "wss"    : Wind Speed (ensembles)
    "sfcWind": Wind Speed (cordex)
'''

