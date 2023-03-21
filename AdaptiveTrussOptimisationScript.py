# -*- coding: utf-8 -*-

import numpy as np
import matplotlib

blankLength = 0.520; # length of the blank in m
blankHeight = 0.250; # height of the blank in m

lengthDisc = 0.01; # discretisation for length in m
heightDisc = 0.01; # discretisation for height in m

numOfLengthNodes = blankLength/lengthDisc;
numOfHeightNodes = blankHeight/heightDisc;

l = np.zeros()
