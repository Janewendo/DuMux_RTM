# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:37 2024

@author: janewendu
"""
# from dune.grid import reader, leafGridView
# 
# domain2d = (reader.dgf, "leo2d.dgf")
# 
# # specify dimgrid since it cannot be easily extracted from the dgf file
# ugView = leafGridView(domain2d, dimgrid=2)
# ugView.plot(figsize=(5,5))
import dune.grid
print(dir(dune.grid))

from dune.grid import structuredGrid, plotGrid  # Import the correct functions

# Specify the domain as a 2D grid from the DGF file
domain2d = structuredGrid(reader.dgf("leo2d.dgf"), dimgrid=2)

# Create a grid view and plot it
plotGrid(domain2d, figsize=(5, 5))