How to contribute to pyhdust
=========================================
References:
    - [PT-BR] http://pythonclub.com.br/como-fazer-fork-clone-push-pull-request-no-github.html
    - [PT-BR] http://rogerdudler.github.io/git-guide/index.pt_BR.html

1. Access `github page of pyhdust <https://github.com/danmoser/pyhdust>`_ and make a *fork*.

2. In *download*, copy the url and paste it in the following command in your terminal (be sure to do it in the right local path):

.. code:: bash

    git clone https://github.com/GITUSER/pyhdust.git


3. If you had already downloaded it, sync it with the original source to avoid conflicts:

.. code:: bash

    git remote add upstream https://github.com/danmoser/pyhdust.git  # you one need to do this in the first time
    git fetch upstream
    git checkout master
    git merge upstream/master

4. Make your changes/contribuitions.

5. To see them, type inside your local folder:

.. code:: bash

    git status

6. The are two ways for accepting the changes:

.. code:: bash

    git add pyhdust/interftools.py  # example file
    # or
    git add .
    # The later, inclues all the newly modified files 

7. Save your actions in the github repository:

.. code:: bash
    
    git commit -m "Message describing what you have done"

    # if it is your first time, you need to configure git (i.e., provide your user name and your email on github):
    git config --global user.email "you@example.com"
    git config --global user.name "Your Name"

8. Send the changes to your github account:

.. code:: bash
    
    git push

9. Do *pull request* in the `github page of pyhdust <https://github.com/danmoser/pyhdust>`_.


Example script
------------------
.. code:: bash

    #!/bin/bash

    # This is a script to submit pyhdust updates

    # synchronize with original source
    # git remote add upstream https://github.com/danmoser/pyhdust.git
    git fetch upstream
    git checkout master
    git merge upstream/master
    echo
    echo Perform the changes now
    echo

    # submit changes
    git add .
    git commit -m "commit message"
    git push
    echo
    echo If you are not the owner of this repository, please open a PULL REQUEST at github
    echo

End
----
Go back to the :doc:`Main Page <index>`.
