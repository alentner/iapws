# iawps::IF97
### A Water Properties Python Library
==================================

A Python package for providing water and steam properties based on 
the 2007 revised release of the Industrial Formulation 1997 from the 
International Association for the Properties of Water and Steam.

---------------------------------------------------------------------

_Insert descriptions of the following: overview of the library, basic usage examples, etc._

Not implemented equations list:
* Supplementary ... Metastable-Vapor 	: Region 2
* Backward Equation v(p, h)				: Region 3
* Backward Equation v(p, T)				: Region 3
* Basic Equation g(p, T)       			: Region 5
* Backward Equation T(p, s)				: Region 5
* Backward Equation T(p, h)				: Region 5

Not implemented unit testing list:
* Wrapper consistancy (partials, sat, ext ...)
* Bacward equations using (P, s)
* Backward Partials w.r.t. h in regions 1 and 2, 
                           and dhdp in region 4
* Region identifications