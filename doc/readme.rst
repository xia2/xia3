xia3
----

============
Introduction
============

Want to design xia3, doing so by means of mocking up a user manual to figure out what we actually want to have. Use case 0: processing some prion data from Diamond Light Source i03.

==========
Use Case 0
==========

Start with a data set consisting of 48 sweeps of 100 x 0.1 degree
images. Want to:

* integrate data in P1
* look at the unit cells, decide on one (P1) and rerun with this
  target
* discover symmetry, indexing ambiguity
* refine in correct symmetry, maybe reprocess with symmetry applied
* sort together
* assess isomorphism
* scale
* assess isomorphism

In the first cut, however then we want to also handle these ideas:

* add another data set
* remove a data set

API something like

xia3.process /path/to/images # will remember data sets in data000.json
xia3.collate_unit_cells
xia3.process triclinic_unit_cell=a,b,c,al,be,ga
xia3.cosym
xia3.apply_symmetry
xia3.reindex_and_refine
xia3.scale
xia3.cluster
xia3.exclude foo_01_####.cbf
xia3.exclude bar_02_####.cbf:101:200
xia3.scale

xia3.add /path/to/more # will update data000.json to data001.json
xia3.scale

