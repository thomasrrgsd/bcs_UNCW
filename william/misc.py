#################################
#				#
# Miscellaneous functions       #
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

def menu_call (log1, log2):

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

  return choice, angle, root_file_name
