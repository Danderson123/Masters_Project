#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 11:25:05 2020

@author: danielanderson
"""


import subprocess

f = str((open('D39V.fasta','r').read())) 
# Set up the echo command and direct the output to a pipe
output = subprocess.Popen(['prokka', f], stdout = subprocess.PIPE).communicate()[0]
