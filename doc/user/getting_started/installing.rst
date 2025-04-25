Installation Guide
==================

Short Version
-------------
To build and install :code:`SeQuant`, run the following commands in the source directory:

.. code-block:: bash

    cmake -B build -S .
    cmake --build build
    cmake --build build --target install


Prerequisites
-------------

Mandatory
~~~~~~~~~

* `CMake <https://cmake.org/>`_ 3.27 or later
* `C++20 compatible compiler <https://en.cppreference.com/w/cpp/compiler_support#cpp20>`_
* `Boost <https://www.boost.org/>`_ 1.81 or later
* `Range-V3 <https://github.com/ericniebler/range-v3.git>`_ 0.12.0 or later
* `libperm <https://github.com/Krzmbrzl/libPerm>`_

.. note:: SeQuant can download and build Boost if configured with :code:`Boost_FETCH_IF_MISSING=ON`, but the use of Boost provided by the system package manager is recommended. The following non-header-only Boost libraries are required, hence Boost must be configured/built:

    * :code:`Boost.Regex`
    * :code:`Boost.Locale`


Optional
~~~~~~~~
* `TiledArray <https://github.com/ValeevGroup/tiledarray.git>`_ - for building coupled-cluster evaluation tests
* `Eigen3 <http://eigen.tuxfamily.org/>`_ - for building coupled-cluster integration tests

.. note:: If not found, SeQuant can download and build Range-V3, libperm and TiledArray.


Useful CMake Variables
----------------------

.. list-table::
   :widths: 20 10 70
   :header-rows: 1

   * - CMake Variable
     - Default
     - Description
   * - `CMAKE_CXX_COMPILER <https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html#variable:CMAKE_%3CLANG%3E_COMPILER>`_
     -
     - Specifies the C++ compiler to use.
   * - `CMAKE_PREFIX_PATH <https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html>`_
     -
     - This semicolon-separated list specifies search paths for dependencies (Boost, Range-V3, etc.).
   * - `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html>`_
     -
     - The installation path for SeQuant.
   * - SEQUANT_TESTS
     - `BUILD_TESTING <https://cmake.org/cmake/help/latest/variable/BUILD_TESTING.html>`_
     - Enables test targets, e.g. ``check-sequant``.
   * - SEQUANT_EVAL_TESTS
     - OFF
     - Enables SeQuant evaluation tests using ``TiledArray`` and ``BTAS``.
   * - SEQUANT_MIMALLOC
     - OFF
     - Use `mimalloc <https://github.com/microsoft/mimalloc>`_ for fast memory allocation.
   * - SEQUANT_BUILD_DOCS
     - OFF
     - Enables building of the documentation. See :ref:`documentation-guide` for detailed information.
   * - SEQUANT_PYTHON
     - OFF
     - Enables building of Python bindings.
   * - SEQUANT_USE_SYSTEM_BOOST_HASH
     - ON
     - Use system Boost for hashing? Set to OFF to make hashing independent of Boost, thus value-portable
   * - SEQUANT_IWYU
     - OFF
     - Whether to use the `include-what-you-use <https://github.com/include-what-you-use/include-what-you-use>`_ tool (if found)
   * - Boost_FETCH_IF_MISSING
     - OFF
     - If set to ON, SeQuant will download and build Boost if it is not found by ``find_package(Boost ...)``; this is not recommended.
   * - ENABLE_TBB
     - OFF
     - Enable TBB as an optional prerequisite for C++'s parallel STL


Configuring and Building
------------------------

To configure and build SeQuant, you can use various CMake variables to customize the build process. These variables can be set using the :code:`-D` flag when running the :code:`cmake` command. For example:

.. code-block:: bash

    cmake -B build -S . -D<VARIABLE_NAME>=<VALUE>

Now you can build SeQuant running the following command in the source directory:

.. code-block:: bash

    cmake --build build -S .
    cmake --build build --target check-sequant # for testing
    cmake --build build --target install
