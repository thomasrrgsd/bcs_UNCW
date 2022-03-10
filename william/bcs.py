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

print(config.savepath)
print(config.tree_in[0])
print(config.tree_in[1])
print(config.tree_in[3])
print(config)

log1.info("this is a log1 message")
log2.info("this is a log2 message")

