.. _mfem_mgis_developer_guide_parameter:

=======================
Parameter Handling
=======================

Rationale
---------

The ``Parameter``, ``Parameters``, and ``ParametersValidator`` classes provide a 
type-safe and structured way to handle configuration parameters in ``mfem-mgis``.

- ``Parameter`` is a variant type that can hold various scalar and complex types 
  (booleans, integers, reals, strings, vectors of parameters, dictionaries, and functions).
- ``Parameters`` is a map associating names to ``Parameter`` values.
- ``ParametersValidator`` is a helper class to declare and validate parameters in a 
  structured manner, ensuring type safety and consistency.

The use of ``ParametersValidator`` is the **recommended** way to validate parameters,
though its adoption is still limited in the codebase. The older ``checkParameters`` 
functions are deprecated in favor of ``ParametersValidator``.

Basic Usage
-----------

Declaring Parameters
~~~~~~~~~~~~~~~~~~~~

Parameters are typically passed as dictionaries (``Parameters``):

.. code-block:: c++

   auto params = Parameters{
       {"IterativeMethod", "Newton"},
       {"Theta", 0.5},
       {"MaxIterations", 100}
   };

A ``Parameter`` can hold various types:

.. code-block:: c++

   Parameter p1 = 42;                  // integer
   Parameter p2 = 3.14;                // real
   Parameter p3 = "value";             // string
   Parameter p4 = true;                // boolean
   Parameter p5 = Parameters{};        // nested dictionary
   Parameter p6 = std::vector<Parameter>{1, 2, 3};  // list

Accessing Parameters
~~~~~~~~~~~~~~~~~~~~

Use the ``get`` functions to access parameter values. Two modes are available:

- **Non-throwing mode**: Uses a ``Context`` to report errors.
- **Throwing mode**: Uses the ``throwing`` attribute for constructors or functions 
  where no context is available.

.. code-block:: c++

   // Non-throwing mode (recommended)
   auto ctx = Context{};
   const auto value = get<int>(ctx, params, "MaxIterations");
   if (isInvalid(value)) {
       // Handle error
       return ctx.registerErrorMessage("invalid parameter");
   }

   // Throwing mode (for constructors or functions with attributes::Throwing)
   const auto value = get<int>(throwing, params, "MaxIterations");

The ``contains`` function checks if a parameter exists:

.. code-block:: c++

   if (contains(params, "Theta")) {
       // Parameter exists
   }

Validating Parameters
~~~~~~~~~~~~~~~~~~~~~

The ``ParametersValidator`` class provides a fluent interface to declare and validate 
parameters:

.. code-block:: c++

   auto validator = ParametersValidator{}
       .add<std::string>("Material", "Name of the material")
       .add<int>("MaxIterations", {.required = true})
       .add<real>("Theta", "Time integration parameter", {.required = true})
       .addIncompatibleParametersList({"OptionA", "OptionB"});

   // Validate parameters
   auto ctx = Context{};
   if (!validator.validate(ctx, params)) {
       // Handle validation error
       return ctx.registerErrorMessage("invalid parameters");
   }

The ``add`` method can:

- Declare allowed keys with optional descriptions.
- Specify if a parameter is required.
- Restrict the parameter to specific types using template arguments.
- Add custom validators.

The ``addIncompatibleParametersList`` method declares mutually exclusive parameters.

Predefined Validators
~~~~~~~~~~~~~~~~~~~~~

The ``ParametersValidator`` class provides predefined validators:

- ``addStrictlyPositiveIntegerCheck``: Ensures a parameter is a strictly positive integer.

Custom validators can be added using the ``add`` method with a validator function:

.. code-block:: c++

   auto validator = ParametersValidator{}
       .add("CustomParameter", [](Context& ctx, const Parameter& p) noexcept -> bool {
           if (!is<int>(p)) {
               return ctx.registerErrorMessage("parameter must be an integer");
           }
           const auto value = get<int>(throwing, p);
           if (value < 0 || value > 100) {
               return ctx.registerErrorMessage("parameter must be between 0 and 100");
           }
           return true;
       });

Error Handling Rules
--------------------

1. **Non-throwing mode**: Functions that accept a ``Context`` must **never** throw. 
   Errors are reported by returning ``false`` or an invalid result (e.g., ``std::nullopt``, 
   ``InvalidResult``). The ``Context`` accumulates error messages.

2. **Throwing mode**: Functions marked with ``attributes::Throwing`` may throw exceptions. 
   This mode is reserved for constructors or functions where no ``Context`` is available.

3. **Error propagation**: Use the ``|`` operator to propagate errors from a ``Context``:

   .. code-block:: c++

       auto or_raise = ctx.getThrowingFailureHandler();
       validator.validate(ctx, params) | or_raise;

4. **Error messages**: Provide clear and descriptive error messages. Include the parameter 
   name and the expected type or constraints.

Best Practices
--------------

1. **Use ``ParametersValidator``**: Prefer ``ParametersValidator`` over the deprecated 
   ``checkParameters`` functions for new code.

2. **Validate early**: Validate parameters as early as possible, ideally at the beginning 
   of functions or constructors.

3. **Document parameters**: Provide descriptions for parameters to improve error messages 
   and documentation.

4. **Type safety**: Use template arguments to restrict parameter types when possible.

5. **Required parameters**: Mark required parameters explicitly using the ``required`` 
   option in ``AddArguments``.

6. **Incompatible parameters**: Use ``addIncompatibleParametersList`` to enforce mutual 
   exclusivity between parameters.

Examples
--------

Validating a set of solver parameters:

.. code-block:: c++

   auto validator = ParametersValidator{}
       .add<std::string>("Solver", "Name of the solver", {.required = true})
       .add<int>("MaxIterations", "Maximum number of iterations")
       .add<real>("Tolerance", "Convergence tolerance", {.required = true})
       .addIncompatibleParametersList({"Verbose", "Silent"});

   auto ctx = Context{};
   if (!validator.validate(ctx, params)) {
       return ctx.registerErrorMessage("invalid solver parameters");
   }

Validating a strictly positive integer:

.. code-block:: c++

   auto validator = ParametersValidator{}
       .addStrictlyPositiveIntegerCheck("NumberOfSteps", {.required = true});

   auto ctx = Context{};
   if (!validator.validate(ctx, params)) {
       return false;
   }

Using custom validators:

.. code-block:: c++

   auto validator = ParametersValidator{}
       .add("CustomValue", [](Context& ctx, const Parameter& p) noexcept -> bool {
           if (!is<real>(p)) {
               return ParametersValidator::reportUnmatchedTypeError(ctx, "CustomValue");
           }
           const auto value = get<real>(throwing, p);
           if (value <= 0 || value > 1) {
               return ctx.registerErrorMessage("CustomValue must be in (0, 1]");
           }
           return true;
       }, "A custom value between 0 and 1");
