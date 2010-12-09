.. index:: wmFigaro node, nodes:wmFigaro, nodes:wmFigRodNode, wmFigRodNode, gravity, threads, start time, FPS, frames per second, steps, collisions, self-collisions, volumetric collisions, flip, 



Node attributes
===========================


.. _attributes_wmFigaroNode:

wmFigaro node attributes
--------------------------

=======================================	===========
Attribute				Definition
=======================================	===========
Solver Type 				Sets the simulation type. At present, leave this as **Symmetric Implicit Euler**, as this is the most stable choice.
Enabled 				Turns the simulation on and off. 
Number of Threads			The number of threads to use in the simulation. Most of the simulation calculations are threaded, except self-collisions (if they are on).  
Start Time 				The frame that the simulation resets at. Everything gets rebuilt at this frame.
Frames Per Second 			The animation speed the scene is meant to run at. For example, if this is 24, the simulation assumes the scene is at 24fps. 
					If you this setting does not match the scene as a whole, the simulation can look very strange. 
Sub Steps 				The minimum number of steps to take per frame. The simulation is adaptively stepped when it needs to be. The lower the number of steps, the faster the simulation runs, so this sets the maximum speed of the simulation. 
					Be careful setting this to low values, as it can make the simulation too "bouncy". 
Gravity 				Sets the gravity. By default, this attribute is connected to all rod nodes so they pick up the world gravity from the Figaro node.
Plastic Deformations 			Ignore this at present.
**Solver Tolerances** section		These settings control the tolerances of the Newtonian solver used to calculate the rod movement. The solver works iteratively until it gets a result within these tolerances. Use higher values to get faster but less accurate results, lower values for more accuracy at the expense of computation time.
**Clumping** section 			Ignore this at present. 
**Collisions**	section		
Object Collisions Enabled		Turns the collision handling for meshes on and off. 
Self Collision Penalty Forces		Controls whether the self collision penalty is on. If this is on, all rods found to be too close to each other at the start of a time step are pushed apart slightly.
Full Self Collisions 			Controls whether to use the full iterative edge/edge collision detection and response. If this is off, Figaro does not calculate self-collisions. 
Full Self Collision Iterations		When **Full Self Collisions** is on, this sets the number of iterations of edge/edge self collisions to do. Higher values take longer to calculate, but if the value is too low, the self-collision may not be evaluated properly. 
Full Self Collision COR			When **Full Self Collisions** is on, this sets the coefficient of restitution to use in edge/edge collisions. At 1, the collision is very elastic (things bounce back with the entire collision force).  At 0, the objects stop dead when they collide (no bounce at all).  
**Volumetric Collisions**	
Volumetric Collisions			Turns on the volumetric collision handling. 
Grid DX  				The edge length of the grid cells used for calculating the volumetric collision.
Target Edge Density  			The density desired in each cell. Numbers smaller than the hair density will pull hairs in, larger will puff them out.
Volumetric Radius			The radius of influence of a hair on the grid. This should usually be at least the actual radius of the individual hair. 
Flip					Sets how grid velocities are interpolated to particles. At 0, you get a viscous effect. At 1, you get much more free results, but some collisions may be missed.
Slip 					Affects how much of the grid velocity is used. At 0, Figaro uses the grid velocities; at 1, Figaro ignores volumetric collisions entirely.
Separation Condition			Debug information, ignore.  
Display Grid etc 			Debug information, ignore. 
**Timings** section			Debug information, ignore. 
**Drawing** Section			Debug information, ignore. 
**XML Output** section 			Debug information, ignore. 
=======================================	===========



.. _attributes_wmFigRodNode:

wmFigRodNode attributes
-------------------------

======================================	================
Attribute				Definition
======================================	================
**Input Resampling** section		
Vertex Spacing	 			Set this to 0 to use the vertex spacing from the original model. Otherwise, Figaro resamples the input. This produces much more stable results than using input from the original models (which may have too many vertices, or vertices that are not well spaced out).
Minimum Rod Length			Any rods shorter than this will not be simulated.
**Rod Parameters** section
Cache Frame 				Turn this on to write out the rod positions to disk for later playback.
Read from Cache				Turn this on to read from the cached file rather than running the simulation. This needs you to have previously run the simulation with *Cache Frame* turned on. 
Lock first edge to input		Locks the segment of the rod in place. If this is off, the entire rod (all segments) can move; if it's on, the bottom segment of the rod stays fixed in position as the rest of the rod deforms.
Percentage of Barbershop strands 	The percentage of Barbershop hairs in the furset to simulate.
Gravity 				Sets the gravity. Use this if you want to use a value different from the global gravity for the rods in this node. 
Minor Radius 				The minor radius of the rods.
Major Radius 				The major radius of the rods.
Youngs Modulus 				Sets how stiff the rod is (how much it resists bending). xxx sensible values? 
Shear Modulus 				How much the rod resists twisting along its length. xxx sensible values? 
Internal Damping 			How much internal resistence the rod has to resist oscillating. If you find your hair is too floppy, try increasing the **Youngs Modulus** and this value. 
Density 				The density of the rod, defined in cm^3. 
					This is not defined as mass per vertex (as with Maya hair) because it means the behaviour changes depending on the number of vertices. The hair behaves consistently even if you change the number of vertices.
Mass Damping 				Amount of damping to apply (this is a standard sim system damping).
Simulation Set 				If you want to only simulate some of the rods, enter the indices of those rods.  
					*This may be broken at present.* 
**Drawing** section		
Draw 3DRod				Draws the rods in a 3d view, rather than as a series of CVs linked by lines. 
Draw Scale				When **Draw 3DRod** is on, this sets the width scale for the 3d representations of the rods. Increase this to make the rods thicker onscreen. 
Cache Path 				The location to read or write cache data from/to.
======================================	================


.. _attributes_collision_node:

Collision node attributes
-------------------------

===============================	============
Attribute			Definition
===============================	============
Friction 			The amount of friction to use during collisions. This makes the rods stick rather than slide on the mesh.
Coefficient of Restitution	The amount of bounce to use during collisions. This makes the rods stick or bounce as they hit the mesh.
Separation Strength 		This scales the amount of separation force (how much the two objects "don't want" to collide) in the initial stage of the collision. Higher values give more bounce out, but less stability. 
Thickness			The distance from the surface where things are considered in collision. 
Edge Collisions			Turn this on to calculate full edge/edge collisions rather than just point-triangle.
Draw Collision Data		Debug, ignore. 
=============================== ============

