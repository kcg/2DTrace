2DTrace
=======

2DTrace is a two-dimensional ray tracing program/module written in Python.

Dependencies
------------

* numpy
* scipy (optimize, misc)
* pylab/matplotlib (for the visualizer)

Usage
-----

2DTrace provides a class structure to define ray sources, objects and detectors.
The ray is then traced from the source until it is lost (leaves the "world") or until it hits a detector.
So far there is no absorption implemented, however this could be easily done.
The modular concept allows the addition of different objects. So far, a simple object defined by two analythic, mathematical functions (and its derivatives) is implemented.

Warning
-------

The module is still in an alpha satus! Beware of bugs!
