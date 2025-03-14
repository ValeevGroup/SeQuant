Building the Documentation
================================================

SeQuant uses a combination of `Doxygen <https://www.doxygen.nl/>`_, `Sphinx <https://www.sphinx-doc.org/en/master/>`_ and `Breathe <https://breathe.readthedocs.io/en/latest/>`_ to build documentation. For building the documentation locally, configure CMake with `-DSEQUANT_BUILD_DOCS=ON`.

Dependencies
------------------------------------------------
Building the documentation requires several Python packages. You can install them using pip. The required packages are listed in the ``requirements.txt`` file in the ``doc/sphinx`` directory. You can install them using the following command:

.. code-block:: bash

    pip install -r doc/sphinx/requirements.txt

This will install the required packages, including Sphinx, Breathe, Exhale, and others.
Now you can build the documentation using CMake. The documentation will be generated in the ``doc/html`` directory. You can open the ``index.html`` file in your web browser to view the documentation.

.. code-block:: bash

    mkdir build
    cd build
    cmake -DSEQUANT_BUILD_DOCS=ON ..
    cmake --build . --target doc-sequant

See :doc:`../internal/docs_advanced` for more information on the documentation structure and how to add new documentation.
