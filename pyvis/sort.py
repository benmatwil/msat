import numpy as np

def sortdate(files):
  files = np.array(files)
  dates = files.copy()
  
  for i in xrange(len(files)):
    date = files[i][14:22]
    dates[i] = date[:4]+'-'+date[4:6]+'-'+date[6:]

  dates = np.array(dates, dtype=np.datetime64)

  isort = np.argsort(dates)

  return files[isort], dates[isort]

# def sortdate(files)
#   files = np.array(files)
#   years = np.zeros(len(files), dtype=np.int)
#   months = np.copy(years)
#   days = np.copy(years)

#   for i in xrange(len(files)):
#     date = files[i][9:17]
#     years[i] = int(date[0:4])
#     months[i] = int(date[4:6])
#     days[i] = int(date[6:])

#   dates = np.core.records.fromarrays([years, months, days], names='years, months, days')

#   sorted = np.argsort(dates, order=('years', 'months','days'))

#   return files[sorted]
  

#   a = np.recarray(10, dtype=[('x',np.int32),('y',np.float64,3)]
