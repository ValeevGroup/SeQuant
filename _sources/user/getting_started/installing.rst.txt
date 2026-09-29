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

* `CMake <https://cmake.org/>`_ 3.28 or later
* `C++20 compatible compiler <https://en.cppreference.com/w/cpp/compiler_support#cpp20>`_
* `Boost <https://www.boost.org/>`_ 1.81 or later (1.85 or later with MSVC)
* `Range-V3 <https://github.com/ericniebler/range-v3.git>`_ 0.12.0 or later
* `Eigen3 <http://eigen.tuxfamily.org/>`_ 3.0 or later
* `libperm <https://github.com/Krzmbrzl/libPerm>`_
* `polymorphic_variant <https://github.com/Krzmbrzl/polymorphic_variant>`_
* `Utfcpp <https://github.com/nemtrif/utfcpp>`_ 4.0 or later
* `DTL <https://github.com/cubicdaiya/dtl>`_ 1.12 or later - only used for unit tests
* `CLI11 <https://github.com/CLIUtils/CLI11>`_ 2.0 or later - only used for external interface
* `spdlog <https://github.com/gabime/spdlog>`_ - only used for external interface
* `nlohmann::json <https://github.com/nlohmann/json>`_ 3.0 or later - only used for external interface

.. note:: SeQuant can download and build Boost if configured with :code:`Boost_FETCH_IF_MISSING=ON`, but the use of Boost provided by the system package manager is recommended. The following non-header-only Boost libraries are required, hence Boost must be configured/built:

    * :code:`Boost.Regex`
    * :code:`Boost.Locale`


Optional
~~~~~~~~
* `TiledArray <https://github.com/ValeevGroup/tiledarray.git>`_ - for building coupled-cluster evaluation tests

.. note:: If not found, SeQuant can download and build all dependencies other than CMake and the C++ compiler, provided `git <https://git-scm.com/>`_
   is available on the system.


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
   * - SEQUANT_SKIP_LONG_TESTS
     - OFF
     - Skip long test cases within unit and integration tests. ``Debug`` configurations always skip them, regardless of
       this option.
   * - SEQUANT_BTAS
     - OFF
     - SeQuant will look for (or build) `BTAS tensor library <https://github.com/ValeevGroup/BTAS>`_ and enable its use as an evaluation backend.
   * - SEQUANT_TILEDARRAY
     - OFF
     - SeQuant will look for (or build) `TiledArray tensor framework <https://github.com/ValeevGroup/TiledArray>`_ and enable its use as an evaluation backend.
   * - SEQUANT_TAPP
     - OFF
     - SeQuant will look for (or build)  `TAPP <https://github.com/TAPPorg/reference-implementation>`_ and enable its use as an evaluation backend.
   * - SEQUANT_BENCHMARKS
     - ON
     - Enable SeQuant benchmarks.
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
     - Enable `TBB <https://en.wikipedia.org/wiki/Threading_Building_Blocks>`_ as an optional prerequisite for C++'s `parallel STL
       <https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t>`_
   * - SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
     - ON
     - If set to `OFF` the default context must be initialized and manipulated from single thread only (most users will want to do this).
   * - SEQUANT_ASSERT_BEHAVIOR
     - ``IGNORE`` in ``Release`` and ``MinSizeRel`` mode, ``ABORT`` otherwise
     - Controls how assertions within SeQuant's code are handled. Valid options are ``ABORT``, ``THROW`` and ``IGNORE``. The latter disables
       assertions, whereas the former keep them active and either abort the program or throw an exception on violation respectively.
       On the first configure of a build directory this also seeds the corresponding option of the dependencies SeQuant builds
       from source: ``TA_ASSERT_POLICY`` of TiledArray (which in turn seeds ``BTAS_ASSERT_POLICY`` of the BTAS it builds), or
       ``BTAS_ASSERT_POLICY`` of BTAS when SeQuant builds BTAS itself (``SEQUANT_TILEDARRAY=OFF``). As with any cached option,
       an explicit ``-DTA_ASSERT_POLICY=...``/``-DBTAS_ASSERT_POLICY=...`` wins, and a later change of ``SEQUANT_ASSERT_BEHAVIOR``
       does not re-seed them: set them explicitly, or use a fresh build directory.
   * - SEQUANT_LTO
     - context-dependent — see Description
     - Controls whether SeQuant will be built with `link-time optimizations (LTO) <https://en.wikipedia.org/wiki/Link-time_optimization>`_ in
       non-debug builds (`CMAKE_BUILD_TYPE` != `Debug`). Left unset, it defaults per target type: `ON` for static libraries if the compiler
       supports "fat" LTO objects, `OFF` otherwise; `ON` for other target types (shared libs, executables, etc.); `OFF` for all targets when
       SeQuant is consumed as a subproject. Set to `ON` or `OFF` explicitly to override the per-target-type default for all SeQuant targets;
       an explicit setting is honored regardless of whether SeQuant is the top-level project.


Configuring and Building
------------------------

To configure and build SeQuant, you can use various CMake variables to customize the build process. These variables can be set using the :code:`-D` flag when running the :code:`cmake` command. For example:

.. code-block:: bash

    cmake -B build -S . -D<VARIABLE_NAME>=<VALUE>

Now you can build SeQuant by running the following command in the source directory:

.. code-block:: bash

    cmake --build build
    cmake --build build --target check-sequant # for testing
    cmake --build build --target install


Windows
----------

SeQuant needs more stack space than the 1 MB that Windows reserves per thread by default (Linux and macOS typically
provide 8 MB). SeQuant's own executables are therefore linked with :code:`/STACK:8388608`, and executables using
SeQuant must be linked with an equivalent stack reserve too, e.g. via
:code:`target_link_options(<target> PRIVATE /STACK:8388608)`. This reserve is the default for every thread of the
program; threads created with an explicit stack size that call into SeQuant need a comparable size.
