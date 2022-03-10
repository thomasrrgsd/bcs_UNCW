"""

  - BCS Analysis Configuration file

  - Written for UNIX

  - All file paths are absolute, but you could use relative paths if
    you are comfortable with file structures.

  - The order of list are always 0T, 0M, 0B, 3T, 3M, 3B, ... 18T, 18M, 18B.

  - If you do not have the files for a setting, then leave it as an empty
    list.

"""

import sys

class configuration:

  """

  Add and edit settings here. Each setting is accessed in the following format:
  
       config = configuration() # Create configuration object
       config.savepath          # Access each field by '.'
       config.tof_spectra[4]    # List are accessed like arrays
       
  """

  # File path where root files will be saved to.
  savepath = './'

  # List of all gas in trees.
  tree_in  = ['./nme_888.root',
              './nme_999.root',
              './nme_111.root',
              './nme_222.root']

  # List of all gas out trees.
  tree_out = ['./nme_333.root',
              './nme_444.root']

  # Time of Flight (TOF) spectra paths.
  # If you do not have any spectra, leave list empty. tof_spectra = []
  tof_spectra = []

  # Path where your cutpath root file is saved.
  cutpath  = './'

  # Names of cut objects within the cutpath root file.
  cut_names = []


  """

  Do not touch the functions below.

  """

  def __init__(self): # Empty constructor.
    pass

  def __str__(self):  # Print dunder method. Print(configuration)
    return "Nothing yet"   
