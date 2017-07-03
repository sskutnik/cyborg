Example CyBORG inputs & validation cases
==============

The following directories include a series of validation cases for the depletion performance of the CyBORG module against standalone Origen calculations. These include the following cases:

- UOX, PWR (Westinghouse 17x17)
- UOX, BWR (GE 10x10-8)
- MOX, BWR (GE 10x10-8)


For each case, two Origen validation inputs are provided: one using the older ARP interpolation module and the newer Obiwan module. CyBORG relies on the same interpolation method as Obiwan; ergo this case is considered a more representative validation.

Details on the specifications of each case are as follows.

UOX (PWR)
~~~~~~~~~~~

The UOX-PWR case consists of a Westinghouse 17x17 PWR type assembly, with an initial U-235 enrichment of 4.2%. The assembly is irradiated for 3 12-month cycles at a specific power of 33.333 MW/MTU (W/g) for a total discharge burnup of 33,000 MWd/MTU.

For each batch, an inter-cycle decay time of 1 month is assumed.

UOX (BWR)
~~~~~~~~~~~

The UOX-BWR case consists of a GE 10x10 BWR assembly type with UOX fuel, with an initial average U-235 enrichment of 4.2% and an average moderator density of 0.55 g/cc. The assembly is irradited for 3 12-month cycles at a specific power of 33.333 MW/MTU, for a total discharge burnup of 33,000 MWd/MTU.

For each batch, an inter-cycle decay time of 1 month is assumed.

MOX (BWR)
~~~~~~~~~~~

The MOX-BWR case consists of a GE 10x10-8 BWR assembly type with MOX fuel, with an initial Pu-239 fraction of 62% and a plutonium enrichment of 4.5%, with an average moderator density of 0.55 g/cc. Like the prior cases, the assembly is irradiated for 3 12-month cycles at a specific power of 33.333 MW/MTU, for a total discharge burnup of 33,000 MWd/MTU.

For each batch, an inter-cycle decay time of 1 month is assumed.

The Cyclus representation of this case assumes a simple source of plutonium shipped directly to the reactor as MOX fuel (in order to simplify the validation). A more physically realistic instance would include, for example, irradiation of UOX fuel, followed by reprocessing, stream blending, and UOX fuel fabrication. However, the intent of this case is to provide a simple and clear validation of the Origen depletion case, and thus the input is simplified as much as possible.
