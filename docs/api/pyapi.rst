Python API
##########

This page lists the public Python API.


High-Level Classes
==================

.. autoclass:: cigsegy.SegyNP
   :members:

.. autoclass:: cigsegy.SegyWriter
   :members:


Factory Functions
=================

.. autofunction:: cigsegy.textual_header

.. autofunction:: cigsegy.metaInfo

.. autofunction:: cigsegy.fromfile

.. autofunction:: cigsegy.collect

.. autofunction:: cigsegy.tofile

.. autofunction:: cigsegy.to_npy

.. autofunction:: cigsegy.create_by_sharing_header

.. autofunction:: cigsegy.get_trace_keys

.. autofunction:: cigsegy.create


Tools
=====

.. automodule:: cigsegy.tools
   :members:


Transform
=========

.. automodule:: cigsegy.transform
   :members:


Interp
======

.. automodule:: cigsegy.interp
   :members:


Plot
====

.. automodule:: cigsegy.plot
   :members:


Pysegy
======

``Pysegy`` is the low-level C++ binding.  Prefer ``SegyNP`` and
``SegyWriter`` unless you need direct access to the backend.

.. autoclass:: cigsegy.Pysegy
   :members:
   :undoc-members:
