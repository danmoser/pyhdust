pyhdust
========
**Analysis tools for multi-technique astronomical data and** *hdust* **models**.

|Tests| 

.. |Tests| image:: https://github.com/danmoser/pyhdust/actions/workflows/main.yaml/badge.svg
   :target: https://github.com/danmoser/pyhdust/actions/workflows/main.yaml

Full documentation at `pyhdust.readthedocs.io <http://pyhdust.readthedocs.io>`_.

**pyhdust** should be independent of plataform (Linux, Mac, Windows) and compatible with Python3 (3.6+).

If you make use of **pyhdust** in your work, please cite Section 2.1.4 of my thesis (Faes 2015: `arXiv <https://arxiv.org/abs/1512.06094>`_, `ADS <https://ui.adsabs.harvard.edu/abs/2015PhDT........60F>`_, `BibTeX <https://ui.adsabs.harvard.edu/abs/2015PhDT........60F/exportcitation>`_).


How to install/uninstall
--------------------------
I strongly suggest users to use the ``pip`` plataform:

.. code:: bash

    pip install pyhdust

**pyhdust** dependencies are constantly evolving. If you face an installation problem, check the `versioning page <https://pyhdust.readthedocs.io/versioning.html>`_.

If your are not the root of the system, add the flag ``--user`` to the command above. Dependencies (as filters efficiency and stellar evolution models) are included in the subdirectories.

.. warning::

    Never combine ``sudo`` with ``--user``! Otherwise you will face critical permission problems for your packages!

.. note:: 

    To use the **pyhdust** scripts, the binaries path of your pip installation directory must be in system ``PATH``. If you don't find them, adapt the following command to your ``$HOME/.bashrc``:

    .. code:: bash

        PATH=$PATH:~/.local/bin/


**pyhdust** requires numpy, six and astropy modules. Optionally, it makes use of: 

    - emcee
    - matplotlib
    - pandas
    - pIDLy
    - scipy
    - wget
    - xmltodict

To only **update** the package:

.. code:: bash

    pip install -U --no-deps pyhdust

``-U`` forces the upgrade and ``--no-deps`` do not reinstall the dependent packages. 

For other options, consult the `full documentation <http://pyhdust.readthedocs.io>`_.


How to use the modules
-------------------------
To make use of all routines, the suggestion is to import them as follows:

.. code:: python

    import pyhdust as hdt
    import pyhdust.beatlas as bat
    import pyhdust.fieldstars as fls
    import pyhdust.images as img
    import pyhdust.input as inp
    import pyhdust.interftools as intt
    import pyhdust.jdcal as jdcal
    import pyhdust.poltools as polt
    import pyhdust.phc as phc
    import pyhdust.rotstars as rot
    import pyhdust.singscat as sst
    import pyhdust.spectools as spt
    import pyhdust.stats as stt


License
-----------
The code is free, available under the terms of the `GNU GPL v3.0 license <https://github.com/danmoser/pyhdust/blob/master/LICENSE>`_.
