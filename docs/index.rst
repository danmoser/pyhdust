.. HDUST Python Tools documentation master file

Welcome to **pyhdust** documentation!
=============================================================
**Analysis tools for multi-technique astronomical data and** *hdust* **models**. About the *hdust* code, see Carciofi & Bjorkman (`2006 <http://adsabs.harvard.edu/abs/2006ApJ...639.1081C>`_, `2008 <http://adsabs.harvard.edu/abs/2008ApJ...684.1374C>`_).

Pyhdust is currently at **version 1.5.11**.

**pyhdust** should be independent of plataform (Linux, Mac, Windows) and compatible with any version of Python (3.6+).

.. toctree::
    :titlesonly:
    :caption: Package contents:

    pyhdust
    bcd
    beatlas
    fieldstars
    hdrpil
    images
    input
    interftools
    jdcal
    oifits
    phc
        phc_list
    poltools
    releases
    rotstars
    singscat
    spectools
    stats
    tabulate
    triangle

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

    - astropy
    - emcee
    - matplotlib
    - pandas
    - pIDLy
    - scipy
    - wget
    - xmltodict

To only *update* the package:

.. code:: bash

    pip install -U --no-deps pyhdust

``-U`` forces the upgrade and ``--no-deps`` do not reinstall the dependent packages.

To uninstall it:

.. code:: bash

    pip uninstall pyhdust


Alternative procedures
^^^^^^^^^^^^^^^^^^^^^^^^^^^
See the :doc:`Alternative download/installation <alternative>` page.


How to use the modules
--------------------------
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


How to contribute to the module
----------------------------------
See the :doc:`How to contribute to pyhdust <contribute>` page.


Indices and tables
--------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* `Source code index`_

.. _`Source code index`: _modules/index.html


License
----------

The code is free, available under the terms of the `GNU GPL v3.0 license`_.

.. _GNU GPL v3.0 license:
	https://github.com/danmoser/pyhdust/blob/master/LICENSE

.. raw:: html

    <a href="http://github.com/danmoser/pyhdust/" target="_blank"><img
    style="position: fixed; top: 0; right: 0; border: 0;"
    src="https://github.blog/wp-content/uploads/2008/12/forkme_right_orange_ff7600.png"
    alt="Fork me on GitHub"></a>
