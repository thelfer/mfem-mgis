.. _simulation_class_overview:

========================================
Overview of the simulation :code:`class`
========================================

.. contents::
    :depth: 3
    :local:

A simulation tries to compute the evolution of a physical system over
several temporal sequences. Temporal sequences are defined by the user.

An example of definition of temporal sequences would typically be:

.. code:: python

   temporal_sequences = [0, 2.5, 4]

which defines a first temporal sequence ``[0, 2.5]`` and a second one
``[2.5, 4]``. By defining such sequences, the user impose to |MFEM/MGIS| a
start and an end simulation time, but also impose that the intermediate
times (``2.5`` in the example above) will be among the simulation time
steps.

Temporal sequences are necessarily contiguous.
 
The management of time steps inside temporal sequences is handled by |MFEM/MGIS|.

The main role of a simulation is to manage the time steps in the
temporal sequences. To do so, it relies on the coupling scheme and:

- A time increment computer, whose role is to compute the time increment
  at the beginning of the next time step. The default time increment
  computer returns the difference between the end of the temporal
  sequence and the current time step.
- A time step validator, whose role is to validate a time step after the
  convergence of the coupling scheme. If the time step is invalidated,
  the time step validator shall provide a new time step. The default
  time step validator only relies on user-defined functions.
- A convergence failure handler, whose role is to propose a default time
  step in case of divergence of the coupling scheme. The default
  convergence handler just takes the current time step and divides it by
  two.

A simulation also:

- intializes the physical system if it not initialized at the beginning
  of the simulation.
- calls user-defined initialization tasks
- calls the post-processings defined by the physical system or some
  user-defined post-processing tasks.

Main parameters
===============

- :param:`Times`: list of times defining the temporal sequences.
- :param:`AllowSubStepping`: boolean stating if sub stepping in case of
  convergence failure is allowed. The default value of this parameter is
  `true`.
- :param:`IndependentTemporalSequences`: boolean stating if the temporal
  sequences are independent. Currently, this boolean only choose if the
  last time increment of the previous temporal sequence is taken into
  account to bound the first time increment of the current temporal
  sequence (bounding occurs if this boolean is false. The default value
  of this parameter is `true`.
- :param:`MinimalTimeIncrement`: minimal time increment. If the time
  increment decreases below this value, the simulation is stopped.
- :param:`MaximalTimeIncrement`: maximal time increment. The time step
  can't be greater than this value: the time is set to this value if its
  computation leads to a greater value.
- :param:`MaximumNumberOfFailuresPerTemporalSequence`: maximum number of
  failures or time step rejections allowed within each temporal sequence.
- :param:`TimeIncrementComputer`: strategy used to determine the next time step.
- :param:`TimeStepValidator`:  strategy used to determine the next time step.
- :param:`ConvergenceFailureHandler`: strategy used to determine how a
  divergence of the resolution shall be handled.
- :param:`LimitTimeIncrementIncrease`: boolean stating if the current
  estimate of the next time step can be greater than the previous time
  increment multiplied by the
  :param:`MaximalTimeIncrementRelativeIncrease` parameter.
- :param:`MaximalTimeIncrementRelativeIncrease`: coefficient used to
  determined the maximum ratio between the next time step and the
  previous one. The default value of this parameter is `1.1`, allowing a
  :math:`10\%` increase of the time step compared to the previous one.
- :param:`LimitTimeIncrementDecrease`: boolean stating if the current
  estimate of the next time step can be lower than the previous time
  increment multiplied by the
  :param`:`MaximalTimeIncrementRelativeDecrease` parameter.
- :param:`MaximalTimeIncrementRelativeDecrease`: coefficient used to
  determined the minimum ratio between the next time step and the
  previous one. The default value of this parameter is `0.2`, allowing
  the next time step to be 5 times smaller than the previous one.
- :param:`BalanceTimeIncrement`: boolean stating if time increments
  shall be balanced to avoid small time increments at the end of the
  temporal sequence.
- :param:`TimeIncrementBalancer`: strategy used to balance the time
  increments to avoid small time increments at the end of the temporal
  sequence.
- :param:`MaximumNumberOfTimeSteps`: maximum number of time steps
  allowed per call to the `run` method.

  See Section :ref:`simulation_ending` for details.
- :param:`NumberOfTimeStepsBetweenPostProcessings`: the number of time
  steps between two post-processing times marked as *explicitly
  requested by the user*.
  
  Note that the parameter :param:`NumberOfTimeStepsBetweenPostProcessings`
  has no effect if it is greater than the
  :param:`MaximumNumberOfTimeSteps` parameter.

  See Section :ref:`simulation_post_processings` for details.
- :param:`TimeBetweenPostProcessings`: time between two post-processing
  times marked as *explicitly requested by the user*.

  See Section :ref:`simulation_post_processings` for details.

- :param:`Monitors`: list of simulation monitors.

Time steps management
=====================

Determining the next time increment at the beginning of a time step
-------------------------------------------------------------------

Let :math:`t` be the time at the beginning of the current time step and
:math:`t_{e}` the end of the current temporal sequence and
:math:`\Delta\,t^{(p)}` the previous time step.

The next time increment is defined as follows:

- A first candidate :math:`\Delta\, t_{1}` is given by calling the
  method of the `getNextTimeIncrement` of the coupling scheme. This time
  increment is the minimum increment returned by the
  `getNextTimeIncrement` of all coupling items.
- A second candidate :math:`\Delta\, t_{2}` is given by calling the time
  increment computer.
- :math:`\Delta\, t_{3}` is the minimum of :math:`\Delta\, t_{1}` and
  :math:`\Delta\, t_{2}`, :math:`t_{e}-t` and :math:`\Delta\,
  t_{\textrm{max}}` if the maximal time increment :math:`\Delta\,
  t_{\textrm{max}}` has been defined (see the
  :param:`MaximalTimeIncrement`).
- if :math:`\Delta\, t_{3}` is greater :math:`\Delta\,t_{(p)}` and if
  the :param:`LimitTimeIncrementIncrease` parameter is `true`,
  :math:`\Delta\, t_{3}` is limited by
  :math:`\alpha_{\textrm{i}}\,\Delta\,t_{(p)}` where
  :math:`\alpha_{\textrm{i}}` is the value of the
  :param:`MaximalTimeIncrementRelativeIncrease` parameter.
- if :math:`\Delta\, t_{3}` is lower :math:`\Delta\,t_{(p)}` and if the
  :param:`LimitTimeIncrementDecrease` parameter is `true`,
  :math:`\Delta\, t_{3}` is limited by
  :math:`\alpha_{\textrm{d}}\,\Delta\,t_{(p)}` where
  :math:`\alpha_{\textrm{i}}` is the value of the
  :param:`MaximalTimeIncrementRelativeDecrease` parameter.
- if :math:`\Delta\, t_{3}` is lower than :math:`\Delta\,
  t_{\textrm{min}}`, the simulation is stopped.

The time increment :math:`\Delta\, t_{3}` can be balanced to avoid the
small time steps at the end of the time step (see the
:param:`BalanceTimeIncrement` and :param:`TimeIncrementBalancer`
parameters).

Determining the next time increment in case of convergence failure
------------------------------------------------------------------

In case of convergence failure, the convergence failure handler is
called to propose an estimate of a new time increment
:math:`\Delta\,t_{e}`.

The next time increment is computed as follows:

.. math::

   \Delta\,t  = \max(\Delta\,t_{e}, \alpha_{\textrm{i}}\, \Delta\,t_{c});

where :math:`\alpha_{\textrm{i}}` is the value of the
:param:`MaximalTimeIncrementRelativeIncrease` parameter and
:math:`\Delta\,t_{c}` the current time increment (corresponding to the
time step which lead to a failure).

The default convergence failure handler divides the current time step by
\(2\).

Determining the next time increment in case of invalidation of the time step
----------------------------------------------------------------------------

The time step validator must propose an estimate of a new time increment
:math:`\Delta\,t_{e}` when it invalidates the current time step.

The next time increment is computed as follows:

.. math::

   \Delta\,t  = \max(\Delta\,t_{e}, \alpha_{\textrm{i}}\, \Delta\,t_{c});

where :math:`\alpha_{\textrm{i}}` is the value of the
:param:`MaximalTimeIncrementRelativeIncrease` parameter and
:math:`\Delta\,t_{c}` the current time increment (corresponding to the
rejected time step).

The default time step validator just calls user defined validators. If
no user defined validator is defined, the time step is always validated.

.. _simulation_ending:

How simulation is stopped
-------------------------

Simulation may fail du to non-convergence of the coupling scheme and
excessive sub-steppings.

Simulation may stop on success for two reasons:

- the last temporal sequence has been completed,
- the maximum number of time steps has been reached.

In the second case, the :cxx:`run` method can be called again and again,
up to the point where the last temporal sequence has been completed.

To know, if last temporal sequence has been completed, one may check the
output of the :cxx:`getTimes` methods:

- if the returned array contains only one value (corresponding the the
  value of the time at the end of the last temporal sequence), then the
  last temporal sequence has been completed,
- if the returned array contains more than one value, one may safely call
  the :cxx:`run` method again.

.. _simulation_post_processings:

Conditions for post-processings to be executed
----------------------------------------------

Every time a post-processing is called, it recieves a boolean stating if
the post-processing time as been *explicitly requested by the user*.

Lightweight post-processings (:cxx:`Curves` for instance) may ignore
this boolean by default.

This boolean is mostly important for heavy post-processsings in terms of
memory, disk usage or computations, such as the :cxx:`VTKExport`
post-processing.

By convention, every post-processing shall expose a parameter named
boolean :param:`AllTimeSteps` which allows to select if the
post-processing must be executed at the end of each time step or only at
at times *explicitly requested by the user*.

The difference between lightweight and heavy post-processings is only
the default value of this parameter (:cxx:`true` for lightweight
post-processings, :cxx:`false` for heavy post-processings).

Post-processings times *explicitly requested by the user*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The end of each temporal sequence is always a post-processing times.

Two parameters allows to select additional post-processing times:

- :param:`NumberOfTimeStepsBetweenPostProcessings` the number of time
  steps between two post-processing times marked as *explicitly
  requested by the user*.
- :param:`TimeBetweenPostProcessings` between two post-processings
  marked as *explicitly requested by the user*.

The criteria (to mark a post-processing time as *explicitly requested by
the user*) associated with those parameters rely with some internal
counters which are reset at each call to the :cxx:`run` method.

Note also that those criteria are taken evaluated independently, i.e. if
the criterion associated with
:param:`NumberOfTimeStepsBetweenPostProcessings` is satisfied, it does
not reset the internal timer associated with the criterion associated
with the :param:`TimeBetweenPostProcessings`.


Note that:

- each end of a temporal sequence is always flaged as *explicitly
  requested by the user*, independently of the
  :param:`NumberOfTimeStepsBetweenPostProcessings` parameter,

Appendix
========

Description of the default time increment balancer
--------------------------------------------------

Let :math:`\Delta\, t_{c}` be the current estimate of the time
increment, :math:`t` be the time at the beginning of the current time
step and :math:`t_{e}` the end of the current temporal sequence.

Let :math:`q` be equal to :math:`\lfloor (t_{e}-t)/\Delta\, t_{c} \rfloor`
(where :math:`\lfloor x \rfloor` is the integral part of :math:`x`) and
:math:`r` be equal to :math:`(t_{e}-t)/\Delta\, t_{c}-q`. Note that `r`
is positive and lower than :math:`1`.

Assuming that the time increment remains constant, the number of time
steps will be:

1. :math:`q` if :math:`t_{e}-t` is exactly proportional to
   :math:`\Delta\, t_{c}`. This case can only be verified approximatly by
   verifying that :math:`r` is close enough to :math:`0`. The default
   time increment balancer compares `r` to a value
   :math:`r_{\textrm{min}}` whose default value is :math:`5\,10^{-2}` and
   which can be modified by a parameter named named
   :param:`MinimalRelativeRemainder` .
2. :math:`q+1` otherwise. However, this may lead to a final time step
   that can be much smaller than the other. To avoid this, :math:`r` is
   compared to a parameter named :param:`MaximalRelativeRemainder` whose
   default value is :math:`0.8`.

The default time increment balancer thus selects the time step as follows:

- If :math:`r<r_{\textrm{min}}`, the time increment is chosen as :math:`(t_{e}-t)/q`.
- If :math:`r>r_{\textrm{max}}`, the time increment is chosen as :math:`\Delta\, t_{c}`.
- If :math:`r_{\textrm{min}} \leq r \leq r_{\textrm{max}}`, the time
  increment is chosen as :math:`(t_{e}-t)/(q+1)`.

