.. _release-guide:

Release Guide
=============

This document describes the steps for making a new SeQuant release.

Version Numbering
-----------------

SeQuant uses `calendar versioning <https://calver.org/>`_: ``YY.N`` for a release and ``YY.N.P`` for a patch release,
with an optional prerelease identifier (e.g., ``alpha.1``, ``beta.2``, ``rc.1``). ``YY`` is the two-digit year, ``N``
counts the releases within that year starting from 0, and ``P`` counts the patch releases of ``YY.N`` starting from 1.
E.g., the third release of 2026 is ``26.2``, its first patch release is ``26.2.1``, and the first release of 2027 is
``27.0``. Releases up to ``2.2.0`` used semantic versioning; ``2.1.0`` and ``2.2.0``, both released in 2026, count as
the first two releases of that year, so the first calendar-versioned release is ``26.2``. The canonical version is
defined in the top-level ``CMakeLists.txt`` and propagated automatically to all other locations via CMake template
substitution.

Files to Update
---------------

Only **two files** require manual version updates:

1. **CMakeLists.txt** (lines ~37--40) — the source of truth:

   .. code-block:: cmake

      set(SEQUANT_YEAR_VERSION 26)   # two-digit year
      set(SEQUANT_YMAJOR_VERSION 2)  # zero-based release counter within the year
      set(SEQUANT_PATCH_VERSION 0)   # patch level; 0 is omitted from the version string
      set(SEQUANT_PRERELEASE_ID )

   Clear ``SEQUANT_PRERELEASE_ID`` for a final release, or set it to e.g. ``alpha.1`` for a pre-release.
   ``SEQUANT_PATCH_VERSION`` stays ``0`` for a ``YY.N`` release and is set to ``P`` for the ``YY.N.P`` patch release.

2. **CITATION.cff** (lines ~5--6):

   .. code-block:: yaml

      version: 26.2
      date-released: 2026-11-15

   Update both the version string and the release date.

Auto-populated Files (no manual edits needed)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following files use ``@SEQUANT_VERSION@`` (or similar) placeholders that CMake fills in at configure time:

- ``SeQuant/version.hpp.in`` — C++ version macros
- ``cmake/sequant-config.cmake.in`` — installed CMake package config
- ``doc/.doxygen/Doxyfile.in`` — Doxygen project number
- ``doc/.sphinx/conf.py.in`` — Sphinx docs version

Release Checklist
-----------------

1. **Create a release branch** (if not already on one), e.g. the 26.2 release and its 26.2.P patch releases would
   live on the v26.2.x branch:

   .. code-block:: bash

      git checkout -b vYY.N.x master

2. **Bump versions** in ``CMakeLists.txt`` and ``CITATION.cff`` as described above.

3. **Commit the version bump**:

   .. code-block:: bash

      git add CMakeLists.txt CITATION.cff
      git commit -m "bump version to YY.N[.P]"

4. **Verify the build** — ensure the project configures, builds, and tests pass:

   .. code-block:: bash

      cmake -B build -S .
      cmake --build build
      ctest --test-dir build

5. **Tag the release**:

   .. code-block:: bash

      git tag -a vYY.N[.P] -m "SeQuant YY.N[.P]"

6. **Push the branch and tag**:

   .. code-block:: bash

      git push origin vYY.N.x
      git push origin --tags

7. **Create a GitHub release** from the tag at
   `github.com/ValeevGroup/SeQuant/releases <https://github.com/ValeevGroup/SeQuant/releases>`_.
   Use the auto-generated release notes as a starting point.

8. **Post-release**: on ``master``, bump the version to the next development pre-release
   (``YY.N+1``, or ``YY+1.0`` if the next release will be the first of a new year, with ``SEQUANT_PRERELEASE_ID`` set
   to ``alpha``) and update ``CITATION.cff`` accordingly. A patch release ``YY.N.P`` is made on the ``vYY.N.x``
   branch and tagged ``vYY.N.P``; it does not change the version on ``master``.
