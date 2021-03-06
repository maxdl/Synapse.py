# -*- coding: utf-8 -*-

import os.path
import sys

title = "Synapse"
author = "Max Larsson"
version = "1.2.1"
date = ("May", "31", "2018")
email = "max.larsson@liu.se"
homepage = "www.liu.se/medfak/forskning/larsson-max/software"
if hasattr(sys, 'frozen'):
    if '_MEIPASS2' in os.environ: 
        path = os.environ['_MEIPASS2']
    else:
        path = sys.argv[0]
else:
    path = __file__
app_path = os.path.dirname(path)
icon = os.path.join(app_path, "syn.ico")
