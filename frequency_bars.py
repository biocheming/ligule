""" Samuel D. Lotz - Lab of Dr. Alex Dickson - Michigan State
University 2015-2-16

This script takes multiple pandas.Series of all the same size and
makes a bar plot corresponding to the numbers in the
arrays. pandas.Series.name is used for labelling bars on legend.
Indices for each Series must be the same and is used for the bar
labels.

"""


import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd

# the frequencies for the bar chart
df = pd.concat(sys.argv[1:], axis=1)
