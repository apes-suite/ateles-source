# Utilizing VTK's Lagrange Cells


This branch aims at the utilization of Lagrange Cells in
VTK as described on:
https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/

This allows for polynomial representation of data up to a polynomial degree
of 10.

The sorting of the nodes can be found in:

https://gitlab.kitware.com/vtk/vtk/-/blob/master/Common/DataModel/vtkHigherOrderHexahedron.cxx#L601

The cell type for High-Order Hexahedra is 72.
