import cgitb
cgitb.enable()
import csv
import StringIO
import urllib2
import urllib
import time

# test longitude 8.385980, test latitude 48.255170
def getWeatherData(Longitude, Latitude, StartYear = 1951, StartMonth = 1, StartDay = 1, EndYear = 2100, EndMonth = 12, EndDay = 31):
  IPCCAssessmentReport = 4 # either 4 or 5
  Dataset             = 'METO-HC_HadRM3Q0_A1B_HadCM3Q0_DM_25km' # if IPCCAssessmentReport =4 use METO-HC_HadRM3Q0_A1B_HadCM3Q0_DM_25km. If IPCCAssessmentReport =5 use either knmihistorical, knmievaluation, knmircp45, knmircp85

  # checking that the arguments are valid 
  if StartYear < 1951 or StartYear > 2100:
    raise ValueError("StartYear must be between 1951 and 2100.")
  
  if StartMonth < 1 or StartMonth > 12:
    raise ValueError("Invalid StartMonth value.")
    
  if StartDay < 1 or StartDay > 31:
    raise ValueError("Invalid StartDay value.")
    
  if EndYear < 1951 or EndYear > 2100:
    raise ValueError("EndYear must be between 1951 and 2100.")
    
  if EndMonth < 1 or EndMonth > 12:
    raise ValueError("Invalid EndMonth value.")
    
  if EndDay < 1 or EndDay > 31:
    raise ValueError("Invalid EndDay value.")

  # convert argument values to int  
  StartYear = int(StartYear)
  StartMonth = int(StartMonth)
  StartDay = int(StartDay)
  EndYear = int(EndYear)
  EndMonth = int(EndMonth)
  EndDay = int(EndDay)
    
#  for testing    
#  Longitude = 8
#  Latitude = 48
#  StartYear = 1951 
#  StartMonth = 1 
#  StartDay = 1 
#  EndYear = 1951 
#  EndMonth = 1 
#  EndDay = 31
 
  #the file to export the output:
#  outFileName   = 'Clipick/outputs/ClipickExportedData.csv' # tip you can build the name of the file to be according to the dates extracted
#  outFileHandle = open(outFileName, 'w')
 
  start_time = time.time() # this is facultative, just to calculate timing of retrieval
 
  # Build the HTTP REQUEST
  pars = {}
  pars['lat']        = Latitude
  pars['lon']        = Longitude
  pars['fmt']        = 'csv' # either csv, htmltable
  pars['tspan']      = 'd'# d=daily; m =monthly
  pars['sd']         = StartDay #'01'
  pars['sm']         = StartMonth #'01'
  pars['sy']         = StartYear
  pars['ed']         = EndDay
  pars['em']         = EndMonth
  pars['ey']         = EndYear
  pars['dts']        = Dataset# Beware of dates for extraction
  pars['ar']         = IPCCAssessmentReport # either 4 or 5
  pars['mod']        = "hisafe" # either yieldsafe or hisafe
  url                = 'http://www.isa.ulisboa.pt/proj/clipick/climaterequest_fast.php'
  url_pars           = urllib.urlencode(pars)
  full_url           = url + '?' + url_pars
  print "Request made to " + full_url
  response           = urllib2.urlopen(full_url)
  the_page           = response.read()
 
  f                  = StringIO.StringIO(the_page)
  reader             = csv.reader(f, delimiter=',')
   
  # CEATE AN ARRAY FROM THE REQUESTED CSV OUTPUT
  result=[]
  for row in reader:
      result.append(row)
      
      # WRITE THE RESULTS IN THE FILE AND CLOSE IT
#  print "Output is being written in " + outFileName
#  for i in result:
#      outFileHandle.write(",".join(i) + "\n")
  #outFileHandle.flush()
#  outFileHandle.close()
 
  #Facultative...
  end_time = time.time()
  print "Processed in " + str(round(end_time - start_time,0)) + " seconds"
#  print "Output stored in " + outFileName
  print "done"
  return result