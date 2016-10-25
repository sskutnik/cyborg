.. _hello_world:

Hello, Cyborg!
==============

CyBORG (CYclus-Based ORiGen) is an [Origen-powered](https://scale.ornl.gov) 
reactor analysis module for Cyclus. It uses the Origen 6.2 API and reactor 
data libraries to perform physics-based reactor depletion calculations from 
within Cyclus. 

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

#Building and installing CyBORG

##Installation requirements

To install CyBORG, you'll need the following:

* SCALE v.6.2.1 
* Cyclus (with Cycamore)

Note that the SCALE 6.2.1 release libraries are compiled using GNU gcc 4.8.3; 
thus, if compiling Cyclus and CyBORG with a similar version of gcc, you should
be able to build directly against the SCALE shared libraries. However, if using
another version (e.g., gcc 4.9+ or 5.0+), you will need to rebuild SCALE from 
source. Please consult the SCALE README for instructions on building SCALE.

CyBORG doesn't require any additional third-party libraries beyond what is 
required for Cyclus and SCALE; as long as you can build these two packages 
on your target system, you can build CyBORG.

##Building CyBORG

CyBORG includes a convenient install script ``install.py``, which will try
to locate all of the dependencies CyBORG needs to build, including the 
appropriate SCALE shared object libraries used by ORIGEN.

When building, you will likely want to point CyBORG's installer to the default
location where it can find the tagged ORIGEN reactor data libraries. Do this
by using the ``orglib_root`` flag with the ``install.py`` script, i.e.:

.. code-block:: bash

   cyborg $ python install.py --orglib_root=/path/to/tagged/origen/libraries

##"Tagged" Origen libraries

CyBORG uses the newest ``Obiwan`` package in SCALE 6.2 for reactor data library
interpolation, based on reading reactor data library descriptor information 
used for interpolation (e.g., enrichment, burnup, etc.) directly off of the 
library. As of SCALE 6.2.1, this information is not directly included on the 
libraries released with SCALE; thus, it must be added by the user. 

A simple script has been developed which uses information on the 
``arpdata.txt`` text-based database to automatically add tag information to
Origen reactor data libraries. 

**TODO:** 
- [ ] Add information on the library tagging script
- [ ] Add the auto-tagging script to the repository

------------

#Using CyBORG for physics-based depletion analysis


