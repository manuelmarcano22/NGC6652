#!/bin/sh



cd q1/
python vimos_bias.py
esorex vmbias --StackMethod=Average bias.sof

python vimos_calib.py
esorex vmifucalib calib.sof

python ifustandard.py
esorex vmifustandard ifustandard.sof

python ifuscience.py
esorex vmifuscience --CalibrateFlux=true ifuscience.sof

