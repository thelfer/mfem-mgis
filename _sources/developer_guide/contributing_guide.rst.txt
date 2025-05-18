==================
Contributing guide
==================


Contributions to mfem-mgis are greatly appreciated. Please take a moment to review this document in order to make the contribution process easy and effective for everyone involved. Following these guidelines helps to communicate that you respect the time of the developers managing and developing this open source project. In return, they should reciprocate that respect in addressing your issue or assessing patches and features.

This part includes the rules written in the `CONTRIBUTING.md` file.

Using the issue tracker
-----------------------

The issue tracker is the preferred channel for bug reports, features requests and submitting pull requests, but please respect the following restrictions:

    Please do not use the issue tracker for personal support requests (contact directly the authors ``tfel-contact@cea.fr``

Bug reports
-----------

A bug is a demonstrable problem that is caused by the code in the repository. Good bug reports are extremely helpful - thank you!

Guidelines for bug reports:

  * Use the GitHub issue search: check if the issue has already been reported.

  * Check if the issue has been fixed: try to reproduce it using the latest master or development branch in the repository.

  * Isolate the problem: ideally create a [reduced test case].

A good bug report shouldn't leave others needing to chase you up for more information. Please try to be as detailed as possible in your report. What is your environment? What steps will reproduce the issue? What compiler(s) and OS experience the problem? What would you expect to be the outcome? All these details will help people to fix any potential bugs.

Example:

.. code-block:: text

   Short and descriptive example bug report title

   A summary of the issue, versions of MFEM, MGIS, TFEL/MFront used and the OS/compiler environment in which it occurs. If suitable, include the steps required to reproduce the bug.

   - This is the first step
   - This is the second step
   - Further steps, etc.

   Any other information you want to share that is relevant to the issue being reported. This might include the lines of code that you have identified as causing the bug, and potential solutions (and your opinions on their merits).



Feature requests
----------------


Feature requests are welcome. But take a moment to find out whether your idea fits with the scope and aims of the project. It's up to you to make a strong case to convince the project's developers of the merits of this feature. Please provide as much detail and context as possible.
Pull requests

Good pull requests - patches, improvements, new features - are a fantastic help. They should remain focused in scope and avoid containing unrelated commits.

Please ask first before embarking on any significant pull request (e.g. implementing features, refactoring code, porting to a different language), otherwise you risk spending a lot of time working on something that the project's developers might not want to merge into the project.

Please adhere to the coding conventions used throughout a project (indentation, accurate comments, etc.) and any other requirements (such as test coverage).

Adhering to the following this process is the best way to get your work included in the project. 

Please, before any merge request of a new feature, make sure that non-regression tests are performed with the command: 

.. code-block:: bash

   cd build && make check

.. note::

   New features must be accompanied by one or more ``non-regression tests`` in tests repository and a ``documentation`` in the /docs/source/developer_guide repository in sphinx format. 

.. warning:: 

  By submitting a patch, you agree to allow the project owners to license your work under the the terms of the LGPL License.
  
Before creating a pull request, please follow these steps to ensure that your developments are up to date with the master branch:

1. Fork the project, clone your fork, and configure the remotes:

.. code-block:: bash

    # Clone your fork of the repo into the current directory
    git clone https://github.com/<your-username>/mfem-mgis.git
    # Navigate to the newly cloned directory
    cd mfem-mgis
    # Assign the original repo to a remote called "upstream"
    git remote add upstream https://github.com/thelfer/mfem-mgis.git

2. If you cloned a while ago, get the latest changes from upstream:

.. code-block:: bash

   git checkout master
   git pull upstream master

3. Create a new topic branch (off the main project development branch) to contain your feature, change, or fix:

.. code-block:: bash

   git checkout -b <topic-branch-name>

4. Commit your changes in logical chunks. Please adhere to these git commit message guidelines or your code is unlikely be merged into the main project. Use Git's interactive rebase feature to tidy up your commits before making them public.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

.. code-block:: bash

   git pull [--rebase] upstream master

6. Check that non-regression tests are still running. 

.. code-block:: bash

   cd build && make check

7. Push your topic branch up to your fork:

.. code-block:: bash

   git push origin <topic-branch-name>

8. Open a Pull Request with a clear title and description.

.. note::

  You can accompany your pull request with a request to appear as a contributor to the MFEM-MGIS project. 
