Alternative downloads and installations
=========================================

Alternative downloads
---------------------------
git
^^^^
The package is available under **git** platform (https://github.com/danmoser/pyhdust/)

.. _https://github.com/danmoser/pyhdust/: https://github.com/danmoser/pyhdust/

.. code-block:: bash

    git clone https://github.com/danmoser/pyhdust.git

The command above saves the package in *pyhdust* subfolder. 

To **update** the package with *git*:

.. code-block:: bash

    git pull


zipball
^^^^^^^^
If neither *pip* or *git* are not installed in your system, you can catch the current version as a zipball:

.. code-block:: bash

    curl -L https://github.com/danmoser/pyhdust/zipball/master > pyhdust.zip
    unzip pyhdust.zip


Alternative installations
---------------------------
In the downloaded pyhdust folder, type:

.. code-block:: bash

    python setup.py install

If your are not the root of the system, add the flag ``--user`` to the command
above.

**Alternatively**, to import the modules declare in your *~/.bashrc* file:

.. code-block:: bash
    
    export PYTHONPATH=$PYTHONPATH:~/pyhdust/

where *~/pyhdust/* is the package installation directory.

.. note::

    To uninstall pyhdust from the alternative installations, using *git* you need to remove its folder manually. Using the zipball installation, you can use the following commands:

    .. code-block:: bash

        python setup.py install --record files.txt
        cat files.txt | xargs rm -rf

End
----
Go back to the :doc:`Main Page <index>`.
