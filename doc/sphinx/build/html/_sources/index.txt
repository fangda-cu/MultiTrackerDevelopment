.. Giaro documentation master file, created by
   sphinx-quickstart on Thu Oct 28 10:54:10 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Figaro
=================================================

What is Figaro?
-------------------

Figaro is Weta's hair dynamics simulation.  It takes the air created in Barbershop and adds dynamic details to create hair that moves realistically.  

.. toctree::
   :maxdepth: 2
   
   creating_rods.rst
   attaching_collision_objects.rst
   animating_rods.rst	
   node_reference.rst
   
   

Needs
------

The latest need is wmFigaro-0043_M2010_64


Rods
---------

Figaro uses rods to solve the dynamics for the hair. To create the dynamics, you need to create rods.  The rods then drive the hair when it is animated. 

=============================================	============================
.. image:: images/example_rodsuncreated.png 	.. image:: images/example_rodscreated.png
*Initial furset*.				*Rods created. The rods are shown in green.* 
=============================================	============================

You can control the details of the rods, including:

* 	how many rods are created.
* 	the collision model to use.
* 	rod details such as size, stiffness (lateral and torsional), damping, and so on. 

Collisions
-----------

Figaro supports both *standard* and *volumetric* collisions. 

Standard collisions can be either:

* 	Edge/edge - this calculates collisions by checking the entire edge of each hair. This gives accurate collisions, but takes time.
* 	Point/triangle - this calculates collisons by sampling points along the hair. This is less accurate but much faster.

Volumetric collisions are quite different from standard collisions.  Volumetric collisions calculate the collision by dividing the worldspace into cubes, calculating density within those cubes, and using this density to figure out probable collisions. Volumetric collisions are not as accurate as standard collisions, but are fast to calculate and give quite naturalistic effects. Volumetric collisions are particularly good at interactions between groups of hair, for example when one set of hair is pulling away from another. 




Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

