#################################
#				#
# William T. Jarratt		#
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

format = "%(asctime)s: %(message)s"
file_handler = logging.FileHandler(filename="temp.log")
stdout_handler = logging.StreamHandler(sys.stdout)
handlers = [file_handler, stdout_handler]

logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
file_handler.setFormatter(logging.Formatter(format, datefmt="%H:%M:%S"))
logging.getLogger().addHandler(file_handler)
#logging.getLogger().addHandler(stdout_handler)

config = configuration()

print(config.savepath)
print(config.tree_in[0])
print(config.tree_in[1])
print(config.tree_in[3])
print(config)

logging.info("this is a message")

