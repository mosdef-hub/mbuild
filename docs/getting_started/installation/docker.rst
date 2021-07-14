Using mBuild with Docker
========================

Docker and other containerization technologies allow entire applications
and their dependencies to be packaged and distributed as images. This
simplifies the installation process for the user and substantially
reduces platform dependence (e.g., different compiler versions, libraries,
etc). This section is a how-to guide for using mBuild with docker.

Prerequisites
-------------
A docker installation on your machine. This
`Docker installation documentation <https://docs.docker.com/get-docker/>`_ has instructions to get docker running on your machine.
If you are not familiar with docker, the Internet is full of good tutorials like these from
`Docker curriculum <https://docker-curriculum.com/>`_ and
`YouTube <https://www.youtube.com/watch?v=zJ6WbK9zFpI&feature=youtu.be>`_.

Jupyter Quick Start
-------------------
After you have a working docker installation, use the following command to
start a Jupyter notebook with mBuild and all the required dependencies:

.. code-block:: bash

    $ docker pull mosdef/mbuild:latest
    $ docker run -it --name mbuild -p 8888:8888 mosdef/mbuild:latest

If no command is provided to the container (as above), the container starts a
``jupyter-notebook`` at the (container) location ``/home/anaconda/data``.
To access the notebook, paste the notebook URL into a web browser on your
computer. When you are finished, you can use control-C to exit the notebook
as usual. The docker container will exit upon notebook shutdown.

.. warning::

    Containers by nature are ephemeral, so filesystem changes (e.g., adding
    a new notebook) only persists until the end of the container's lifecycle.
    If the container is removed, any changes or code additions will not persist.
    See the section below for persistent data.

.. note::

    The ``-it`` flags connect your keyboard to the terminal running in the
    container. You may run the prior command without those flags, but be
    aware that the container will not respond to any keyboard input. In
    that case, you would need to use the ``docker ps`` and ``docker kill``
    commands to shut down the container.


Persisting User Volumes
-----------------------
If you are using mBuild from a docker container and need access to data
on your local machine or you wish to save files generated in the container,
you can mount user volumes in the container. User volumes will provide a way
to persist filesystem changes made to a container regardless of the container
lifecycle. For example, you might want to create a directory called
``mbuild-notebooks`` on your local system, which will store all of your mBuild
notebooks/code. In order to make that accessible from within the container
(where the notebooks will be created/edited), use the following steps:

.. code-block:: bash

    $ mkdir mbuild-notebooks
    $ cd mbuild-notebooks/
    $ docker run -it --name mbuild --mount type=bind,source=$(pwd),target=/home/anaconda/data -p 8888:8888  mosdef/mbuild:latest

You can easily mount a different directory from your local machine by changing
``source=$(pwd)`` to ``source=/path/to/my/favorite/directory``.

.. note::

    The ``--mount`` flag mounts a volume into the docker container. Here we
    use a ``bind`` mount to bind the current directory on our local filesystem
    to the ``/home/anaconda/data`` location in the container. The files you see
    in the ``jupyter-notebook`` browser window are those that exist on your
    local machine.

.. warning::

    If you are using the container with jupyter notebooks you should use
    the ``/home/anaconda/data`` location as the mount point inside the container;
    this is the default notebook directory.

Running Python scripts in the container
---------------------------------------
Jupyter notebooks are a great way to explore new software and prototype
code. However, when it comes time for production science, it is often
better to work with python scripts. In order to execute a python script
(``example.py``) that exists in the current working directory of your
local machine, run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/mbuild:latest "python data/test.py"

Note that once again we are ``bind`` mounting the current working directory
to ``/home/anaconda/data``. The command we pass to the container is
``python data/test.py``. Note the prefix ``data/`` to the script; this is because
we enter the container in the home folder (``/home/anaconda``), but our script
is located under ``/home/anaconda/data``.

.. warning::

    Do not bind mount to ``target=/home/anaconda``. This will cause errors.


If you don't require a Jupyter notebook, but just want a Python interpreter,
you can run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/mbuild:latest python

If you don't need access to any local data, you can of course drop the
``--mount`` command:

.. code-block:: bash

    $ docker run mosdef/mbuild:latest python


Different mBuild versions
-------------------------
Instead of using ``latest``, you can use the image ``mosdef/mbuild:stable``
for most recent stable release of mBuild.

Cleaning Up
-----------
You can remove the *container* by using the following command.

.. code-block:: bash

    $ docker container rm mbuild

The *image* will still exist on your machine. See the tutorials at the
top of this page for more information.

.. warning::

    You will not be able to start a second container with the same name
    (e.g., mbuild), until the first container has been removed.

.. note::

    You do not need to name the container `mbuild` as shown in the above
    examples (``--name mbuild``). Docker will give each container a name
    automatically. To see all the containers on your machine, run
    ``docker ps -a``.
