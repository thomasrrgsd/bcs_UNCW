#################################
#				#
# PyROOT BCS Script             #
# 			       	#
#				#
################################# 

import logging
import ROOT
import sys
import select
import time

sys.path.append(".")
from config import configuration

# Create a basic logger and a file/stdout logger.
log1 = logging.getLogger('console_log')
log1.setLevel(logging.INFO)

log2 = logging.getLogger('file_log')
log2.setLevel(logging.INFO)

format = "%(asctime)s: %(message)s"
formatter = logging.Formatter(format, datefmt="%H:%M:%S")

fh = logging.FileHandler(filename="temp.log")
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)

log1.addHandler(ch)
log2.addHandler(fh)
log2.addHandler(ch)

config = configuration()

#print(config.savepath)
#print(config.tree_in[0])
#print(config.tree_in[1])
#print(config.tree_in[3])
#print(config)

#log1.info("this is a log1 message")
#log2.info("this is a log2 message")

menu_options = ['  [0] Generate PH versus PSD histograms',
                '  [1] Generate PH versus PSD with cuts shown',
                '  [2] Generate TOF spectra',
                '  [3] Fit function to combined top, middle, and bottom TOF spectra',
                '  [4] Quit Script']

print(
"""
Welcome to the 25 MeV TUNL analysis of the Summer 2018 Argon run. The filled sample cell
runs are 3681, 3682, 3687, and 3688. The empty runs are 3609, 3691, 3692, 3693, and 3694.
Select an option below and then select the angle to be analyzed. Make sure you have set up
all the paths in the configuration file. IF YOU HAVE NOT PROPERLY SET UP configuration.py,
THEN QUIT THE SCRIPT AND DO SO NOW.
"""
)

for i in menu_options:
  print(i)

while (True):
  
  try:
    option = int(input("\nEnter your choice: "))
  except:
    print("Invalid input, please enter a number.")
    continue

  if (option == 1):
    print("Option 1 selected")
  elif (option == 2):
    print("Option 2 selected")
  elif (option == 3):
    print("Option 3 selected")
    sys.exit(0)
  else:
    print("Invalid option. Please try again.")


sys.exit(0) # Exit code 0 means the script finished as intended.

