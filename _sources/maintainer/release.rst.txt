.. _release-guide:

Release Guide
=============

This document describes the steps for making a new SeQuant release.

Version Numbering
-----------------

SeQuant uses `semantic versioning <https://semver.org/>`_: ``MAJOR.MINOR.MICRO`` with an optional prerelease identifier
(e.g., ``alpha.1``, ``beta.2``, ``rc.1``). The canonical version is defined in the top-level ``CMakeLists.txt``
and propagated automatically to all other locations via CMake template substitution.

Files to Update
---------------

Only **two files** require manual version updates:

1. **CMakeLists.txt** (lines ~37--40) — the source of truth:

   .. code-block:: cmake

      set(SEQUANT_MAJOR_VERSION 2)
      set(SEQUANT_MINOR_VERSION 3)
      set(SEQUANT_MICRO_VERSION 0)
      set(SEQUANT_PRERELEASE_ID )

   Clear ``SEQUANT_PRERELEASE_ID`` for a final release, or set it to e.g. ``alpha.1`` for a pre-release.

2. **CITATION.cff** (lines ~5--6):

   .. code-block:: yaml

      version: 2.3.0
      date-released: 2026-03-15

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

1. **Create a release branch** (if not already on one), e.g. the 2.3.Z releases would live on the v2.3.x branch:

   .. code-block:: bash

      git checkout -b vX.Y.x master

2. **Bump versions** in ``CMakeLists.txt`` and ``CITATION.cff`` as described above.

3. **Commit the version bump**:

   .. code-block:: bash

      git add CMakeLists.txt CITATION.cff
      git commit -m "bump version to X.Y.Z"

4. **Verify the build** — ensure the project configures, builds, and tests pass:

   .. code-block:: bash

      cmake -B build -S .
      cmake --build build
      ctest --test-dir build

5. **Tag the release**:

   .. code-block:: bash

      git tag -a vX.Y.Z -m "SeQuant X.Y.Z"

6. **Push the branch and tag**:

   .. code-block:: bash

      git push origin vX.Y.Z
      git push origin --tags

7. **Create a GitHub release** from the tag at
   `github.com/ValeevGroup/SeQuant/releases <https://github.com/ValeevGroup/SeQuant/releases>`_.
   Use the auto-generated release notes as a starting point.

8. **Post-release**: on ``master``, bump the version to the next development pre-release
   (e.g., ``X.Y+1.0`` with ``SEQUANT_PRERELEASE_ID`` set to ``alpha``) and update ``CITATION.cff`` accordingly.
