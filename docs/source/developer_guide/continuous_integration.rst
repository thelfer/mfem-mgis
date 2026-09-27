.. _mfem_mgis_developer_guide_continuous_integration:

======================
Continuous integration
======================

.. contents::
    :depth: 2
    :local:

How it works
============

The continuous integration runs on GitHub Actions. It builds and tests
mfem-mgis in two complementary ways:

- ``cmake.yml`` installs the dependencies with Spack, then builds mfem-mgis
  with CMake in ``Release``, ``Debug`` and ``Coverage`` modes. It also builds
  and tests `mfem-mgis-examples <https://github.com/latug0/mfem-mgis-examples>`_
  and `mm-opera-hpc <https://github.com/rprat-pro/mm-opera-hpc>`_.
- ``spack.yml`` builds and tests mfem-mgis as a Spack package, as a user
  would install it.

Two kinds of caches keep the runs short:

- the dependencies are stored as Spack binaries in the GitHub Container
  Registry;
- the compilations are stored by ccache in the GitHub Actions cache.

Only the runs on ``master`` write to the Spack binary caches. Pull requests
only read them.

A run with empty caches builds all the dependencies, which takes about 1h30.
The next runs take about 15 minutes, because they find the dependencies and
most of the compilations in the caches.

Workflows
=========

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Workflow
     - Triggers
     - Role
   * - ``cmake.yml``
     - pushes and pull requests on ``master``, every night, manual
     - CMake builds, tests, examples, OperaHPC and coverage report
   * - ``spack.yml``
     - pushes and pull requests on ``master``, every night, releases, manual
     - Spack builds and tests of the checked out sources
   * - ``clean-build-cache.yml``
     - first day of every month, manual
     - empties the two Spack binary caches
   * - ``doxygen.yml`` and ``sphinx.yml``
     - pushes and pull requests on ``master``
     - build the documentation, published on pushes

``cmake.yml``
-------------

The ``run_tests`` job runs once per build type. It uses Spack 1.2.2 and the
``develop`` branch of the Spack packages. The compilers are restricted to
GCC 12.4. TFEL and MGIS are built from their ``master`` branches, without
Python. The ``Debug`` job also uses a debug build of MFEM.

The ``Coverage`` job uploads its report as the ``code-coverage-report``
artifact.

``spack.yml``
-------------

The ``build`` job covers two versions of Spack, with and without MPI:

- Spack 1.1.0 with the ``releases/v2025.11`` packages. They do not provide
  mfem-mgis, which comes from the
  `mfem-mgis Spack repository <https://github.com/rprat-pro/spack-repo-mfem-mgis>`_.
- Spack 1.2.2 with the ``develop`` packages, which provide mfem-mgis.

The job first creates a Spack environment, where ``spack develop`` points
mfem-mgis to the checked out sources. It installs the dependencies and, on
``master``, pushes them to the cache. ``spack install --test=root`` then builds
and tests mfem-mgis. Only the hash of mfem-mgis differs from a regular
installation, so the dependencies still come from the cache.

Spack binary caches
===================

Each workflow has its own cache, as their dependencies differ:

- ``ghcr.io/<owner>/mfem-mgis-buildcache`` for ``cmake.yml``;
- ``ghcr.io/<owner>/mfem-mgis-spack-buildcache`` for ``spack.yml``.

``<owner>`` is the owner of the repository. A fork thus has its own caches.

Reusing the binaries
--------------------

Spack only reuses a binary when its hash matches. The configuration must thus
stay the same from one run to the next. The ``setup-build-cache`` action sets
it:

- the ``x86_64_v3`` target, which all the runners support;
- padded installation paths, so that the binaries can be relocated;
- no reuse of the installed packages, so that the latest versions are used.

The ``master`` branches of TFEL and MGIS are pinned to their latest commit.
Otherwise, an older binary of ``tfel@master`` could be reused.

Any change of the configuration, of the Spack version or of the specs changes
the hashes. The first run after such a change rebuilds the dependencies.

Filling and cleaning
--------------------

1. Each job of a run on ``master`` pushes the dependencies it installed.
   mfem-mgis itself is never pushed.
2. Each job also records the hashes of all the packages it used.
3. Once all the jobs are done, ``update_build_cache`` removes the packages
   that no job used, but only if all the jobs succeeded. It then updates the
   index of the cache.
4. On the first day of every month, ``clean-build-cache.yml`` empties the
   caches. The next nightly runs fill them again.

Pull requests, manual runs on other branches and releases never write to
these caches.

Spack cannot prune a cache stored in a container registry.
``spack buildcache prune`` only supports local, S3 and Google Cloud Storage
mirrors. The ``clean-build-cache`` action thus deletes the package versions
with the GitHub API.

Visibility
----------

The first run on ``master`` creates the two caches as private packages. They
should be made public in their package settings. Public packages are free,
whereas private ones count against the storage quota of the owner.

Once public, the caches can also be used by the continuous integration of
other repositories, such as mfem-mgis-examples and mm-opera-hpc. These
repositories can only read them. Only the workflows of mfem-mgis can write to
them.

ccache
======

The compilations of mfem-mgis, of its tests, of the examples and of OperaHPC
go through ccache. In ``spack.yml``, only the compilation of mfem-mgis does.

- ``setup-ccache`` installs ccache. It restores the latest cache saved under
  the same key: the build type, or the Spack version and the MPI variant.
- ``save-ccache`` removes the entries that the job did not use. It then saves
  the cache under a new key.

GitHub removes the caches unused for 7 days and limits a repository to 10 GB.
A pull request reads the caches of ``master``. It saves its own caches, which
only this pull request can read.

The runners have different processors. The ccache of Ubuntu 24.04 ignores
what ``-march=native`` stands for, so it could reuse an object built for
another processor. ``cmake.yml`` thus builds mfem-mgis with
``-Denable-portable-build=ON``, which removes ``-march=native``. This caused no
measurable slowdown.

Safety
======

Nobody outside the project can fill or corrupt the caches.

- The steps writing to the Spack binary caches only run on ``master``. Only
  the users with write access to the repository can push to ``master``.
- A pull request from a fork runs with a read-only token. It cannot write to
  the Spack binary caches, even if it modifies the workflows. No workflow uses
  ``pull_request_target``, the event that would give it a write token.
- The size of the Spack binary caches is bounded. Each run on ``master``
  removes the unused packages, and the caches are emptied every month.
- A pull request can save ccache entries, but only this pull request can read
  them. They never reach ``master``.
- All the ccache entries share the 10 GB of the repository. Beyond this limit,
  GitHub removes the least recently used entries. A pull request can thus at
  worst slow down the next runs. It cannot make them wrong.

Common operations
=================

Empty the Spack binary caches
   Run the ``Clean the build caches`` workflow from the Actions tab.

Refresh the Spack binary caches
   Run ``Build with Cmake and Run Examples`` or ``Spack`` manually on
   ``master``.

Empty the ccache caches
   Delete them in the Caches page of the Actions tab, or run
   ``gh cache delete --all``.

Understand a slow job
   In the step installing the dependencies, ``fetching from build cache``
   marks a reused binary and ``no binary available`` a built one. The
   ``Save ccache`` step shows the ccache statistics of the job.

Rerun a failed job
   Reruns are safe. The artifacts of the previous attempt are overwritten.

Rules to keep
=============

- Never write to the Spack binary caches from a pull request.
- Never enable the ``Send write tokens to workflows from pull requests``
  setting of the repository.
- Never use ``-march=native`` in the compilations cached by ccache.
- Keep the Spack configuration shared by the workflows in
  ``setup-build-cache``.
- The nightly runs are scheduled at 01:17 UTC. GitHub delays the scheduled
  runs the most at the start of an hour. They may still start a few hours
  late.
