pyhdust
========
**Analysis tools for multi-technique astronomical data and** *hdust* **models**.

Full documentation at `pyhdust.readthedocs.io <http://pyhdust.readthedocs.io>`_.

About *hdust* code: Carciofi & Bjorkman (`2006 <http://adsabs.harvard.edu/abs/2006ApJ...639.1081C>`_, `2008 <http://adsabs.harvard.edu/abs/2008ApJ...684.1374C>`_).

**pyhdust** should be independent of plataform (Linux, Mac, Windows) and compatible with any version of Python (superior to 2.7).


How to install/uninstall
--------------------------
I strongly suggest users to use the ``pip`` plataform:

.. code:: bash

    pip install pyhdust

If your are not the root of the system, add the flag ``--user`` to the command above. Dependencies (as filters efficiency and stellar evolution models) are included in the subdirectories.

.. warning::

    Never combine ``sudo`` with ``--user``! Otherwise you will face critical permission problems for your packages!

.. note:: 

    To use the **pyhdust** scripts, the binaries path of your pip installation directory must be in system ``PATH``. If you don't find them, adapt the following command to your ``$HOME/.bashrc``:

    .. code:: bash

        PATH=$PATH:~/.local/bin/


**pyhdust** requires numpy and six modules. Optionally, it makes use of: 

    - matplotlib
    - pyfits
    - emcee
    - pIDLy
    - scipy
    - pyqt_fit

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
