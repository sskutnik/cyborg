CyBORG: A physics-based reactor module for Cyclus
==============

CyBORG (CYclus-Based ORiGen) is an `Origen-powered <https://scale.ornl.gov>`_ reactor analysis module for the `Cyclus nuclear fuel cycle simulator <https://fuelcycle.org>`_. It uses the `Origen 6.2 API <https://wawiesel.github.io/OrigenAPI-Demo/>`_ and reactor data libraries to perform physics-based reactor depletion calculations from within Cyclus. 

CyBORG uses Origen to dynamically calculate depleted fuel recipes using Origen, interpolating Origen reactor data libraries to the problem-specific conditions including initial enrichment and cycle burnup while also supporting user-specified interpolation parameters for future fuel and reactor types. CyBORG calls the Origen solver API directly to calculate depleted fuel compositions, which then become the output recipe for the CyBORG reactor archetype.

Known issues
~~~~~~~~~~~~

CyBORG is designed to only update the fuel recipe when it needs to; i.e., if the irradiation conditions (power, discharge burnup, initial fuel composition) are unchanged, CyBORG needs to only call Origen once to calculate the discharge fuel recipe. At present, the most computationally-expensive component of this is actually in the reactor data library disk I/O and interpolation; here again, these are one-time costs per instance. CyBORG caches recipes by the input conditions (i.e., initial enrichment, burnup, etc.), so that once any instance of the archetype in a simulation has calculated a discharge fuel recipe, it is available to all instances. 

Nonetheless, as a high-fidelity reactor simulation archetype, it is substantially more computationally-expensive for the first irradiation cycle compared to a strictly recipe-based model (such as the Cycamore reactor archetype). We are currently looking into ways to decrease the disk I/O cost for performing the first reactor data library interpolation.

Getting CyBORG
==============

To build and use CyBORG, you will need to download the source code. 

You can get CyBORG either by using git to
`clone the repository <https://github.com/sskutnik/cyborg.git>`_ or by
`downloading the zip file <https://github.com/sskutnik/cyborg/archive/develop.zip>`_.
Let's put this code in a ``cyborg`` directory and go into it.

**Getting cyborg via git:**

.. code-block:: bash

    $ git clone https://github.com/sskuntik/cyborg.git cyborg
    $ cd cyborg

**Getting cyborg via zip:**

.. code-block:: bash

    $ curl -L https://api.github.com/repos/sskutnik/cyborg/zipball > cyborg.zip
    $ unzip tutorial.zip
    $ mv sskutnik-cyborg-* tutorial
    $ cd cyborg


------------

Building and installing CyBORG
==============================

Installation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~

To install CyBORG, you'll need the following:

* SCALE v6.2.1 
* `Cyclus <https://github.com/cyclus/cyclus>`_ (with `Cycamore <https://cyclus.github.com/cyclus/cycamore>`_)

Note that the SCALE 6.2.1 release libraries are compiled using GNU gcc 4.8.3; 
thus, if compiling Cyclus and CyBORG with a similar version of gcc, you should
be able to build directly against the SCALE shared libraries. However, if using
another version (e.g., gcc 4.9+ or 5.0+), you will need to rebuild SCALE from 
source. Please consult the SCALE README for instructions on building SCALE.

CyBORG doesn't require any additional third-party libraries beyond what is 
required for Cyclus and SCALE; as long as you can build these two packages 
on your target system, you can build CyBORG.

Building CyBORG
~~~~~~~~~~~~~~~

CyBORG includes a convenient install script ``install.py``, which will try
to locate all of the dependencies CyBORG needs to build, including the 
appropriate SCALE shared object libraries used by ORIGEN.

When building, you will likely want to point CyBORG's installer to the default
location where it can find the tagged ORIGEN reactor data libraries. Do this
by using the ``orglib_root`` flag with the ``install.py`` script, i.e.:

.. code-block:: bash

   cyborg $ python install.py --orglib_root=/path/to/tagged/origen/libraries

"Tagged" Origen libraries
~~~~~~~~~~~~~~~~~~~~~~~~~

CyBORG uses the newest ``Obiwan`` package in SCALE 6.2 for reactor data library
interpolation, based on reading reactor data library descriptor information 
used for interpolation (e.g., enrichment, burnup, etc.) directly off of the 
library. As of SCALE 6.2.1, this information is not directly included on the 
libraries released with SCALE; thus, it must be added by the user. 

A simple script has been developed which uses information on the 
``arpdata.txt`` text-based database to automatically add tag information to
Origen reactor data libraries. 

**TODO:** 
   - Add information on the library tagging script
   - Add the auto-tagging script to the repository

------------

Using CyBORG for physics-based depletion analysis
=================================================

CyBORG supports the following **numeric** arguments to describe the fuel burnup

   - ``power_cap``: Reactor thermal power (in MW)
   - ``assem_size``: Mass of a single fuel assembly (in kg)
   - ``n_assem_core``: Total number of assemblies in the core
   - ``n_assem_batch``: Number of fuel assemblies per batch. Defaults to ``n_assem_core`` (i.e., single-batch core)
   - ``n_assem_fresh``: Minimum number of fresh fuel assemblies to keep in storage if possible (default: 0)
   - ``n_assem_spent``: Number of discharged fuel assemblies that can be stored in the reactor (default: 1000000000)
   - ``cycle_time``: Length of a full irradiation cycle (excluding refueling time), in Cyclus time steps
   - ``refuel_time``: Length of a refueling (down) time, in Cyclus time steps
   - ``reactor_lifetime``: Lifetime of the reactor in the simulation, in Cyclus time steps (default: 400)

In addition, CyBORG takes the following **string** type arguments:

   - ``fuel_type``: Reactor fuel type (UOX, MOX, or "other". Default: UOX). Used for determining fissile information to extract from the input recipe for interpolation (i.e., U-235 content for UOX, Pu-239 and Pu fraction for MOX, and nothing for "other").
   - ``assembly_type``: Origen reactor data library to use for assembly design. (default: "w17x17").
   - ``spent_fuel``: Name of the spent fuel commodity generated by this CyBORG reactor instance (default: "spent_fuel")
   - ``power_name'': Name of the power commodity the reactor produces (default: "power")
   - ``lib_path``: Path to the ORIGEN reactor data libraries. Defaults to the value set by the ``ORGLIB_ROOT`` flag when building CyBORG.

Finally, CyBORG can take the following **complex** argument types (i.e., nested XML data):

   - ``fuel_recipes`` One or more string values corresponding to fuel recipe names accepted by the reactor
   - ``fuel_incommods`` One or more string values corresponding to input commodity names that this reactor will bid for
   - ``fuel_prefs`` Real-valued fuel preferences - one value per recipe / incommodity. Defaults to 1.0 for all preferences if not specified.
   - ``core_power_frac`` List of double values (one per cycle) to indicate the core power fraction for each batch of assemblies; used to specify non-uniform burnups between cycles. Valid values are in (0,1). Unnormalized values are renormalized to sum to 1. Number of entries must equal the number of batches as determined by ``n_assem_core`` / ``n_assem_batch``
   - ``tags`` Tag/value pairs for interpolation parameters to be used for problem-dependent library interpolation (i.e., expressed as ``<tags><tag> tagName </tag> <value> tagValue </value> ... </tags>``). Tags must be present on the specified Origen library.
   
Usage examples
==============

Examples of how to use CyBORG can be found in the `inputs <https://github.com/sskutnik/cyborg/tree/develop/inputs>`_ directory. These illustrate how CyBORG can be configured for use within a Cyclus simulation to generated depleted fuel recipes which are passed back into the Cyclus simulation.

