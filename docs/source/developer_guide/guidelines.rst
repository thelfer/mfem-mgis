.. _mfem_mgis_developer_guide_guidelines:

============
Guidelines
============

This document outlines the coding guidelines for developing in ``mfem-mgis``.
These rules ensure consistency, maintainability, and robustness of the codebase.

General Rules
-------------

1. **C++ Standard**: The project uses **C++20**. All code must conform to this standard.

2. **Code Style**: Follow the project's ``clang-format`` configuration (see ``.clang-format``).
   Run ``clang-format -i $(find . -name "*xx")`` to format your code.

3. **Header Guards**: All header files must use the standard header guard pattern:

   .. code-block:: c++

      #ifndef LIB_MFEMMGIS_<PATH>_HXX
      #define LIB_MFEMMGIS_<PATH>_HXX
      
      // code
      
      #endif /* LIB_MFEMMGIS_<PATH>_HXX */

4. **Includes**: Use forward declarations where possible to reduce compile-time dependencies.
   Include headers in the following order:

   - Standard library headers
   - Third-party library headers (e.g., MFEM, MGIS)
   - Project headers

5. **Namespaces**: All code must be in the ``mfem_mgis`` namespace or a nested namespace.

6. **Doxygen Comments**: All public classes, methods, and functions must have Doxygen comments.
   Avoid redundant ``\brief`` sections if they duplicate the description.

   .. code-block:: c++

      /*!
       * \brief Brief description.
       * 
       * Detailed description if necessary.
       * 
       * \param[in] arg: Description of the argument.
       * \return Description of the return value.
       */

7. **Error Messages**: Error messages must be clear, descriptive, and include relevant context
   (e.g., parameter names, expected types, or constraints).

Error Handling
--------------

The project uses a **dual error handling model** (see `MGIS error handling documentation <https://thelfer.github.io/mgis/web/error_handling.html>`_):

1. **Non-throwing mode**: For functions that accept a ``Context``.
2. **Throwing mode**: For constructors or functions marked with ``attributes::Throwing``.

Non-throwing Mode
~~~~~~~~~~~~~~~~~

- **Never throw exceptions** in functions that accept a ``Context``.
- Report errors using the ``Context``:

  .. code-block:: c++

     bool myFunction(Context& ctx, ...) noexcept {
         if (error_condition) {
             return ctx.registerErrorMessage("descriptive error message");
         }
         return true;
     }

- Use ``isInvalid`` to check for errors in returned values:

  .. code-block:: c++

     auto result = someFunction(ctx, ...);
     if (isInvalid(result)) {
         // Handle error
         return ctx.registerErrorMessage("failed to do something");
     }

- Return types for non-throwing functions:

  - ``bool``: For functions that return a success/failure status.
  - ``std::optional<T>``: For functions that may return a value or nothing.
  - ``OptionalReference<T>``: For functions that may return a reference or nothing.
  - ``InvalidResult``: A special type representing an invalid result.

Throwing Mode
~~~~~~~~~~~~~

- Reserved for **constructors** or functions marked with ``attributes::Throwing``.
- Use the ``raise`` function to throw exceptions:

  .. code-block:: c++

     MyClass::MyClass(...) : member(throwing, ...) {
         if (error_condition) {
             raise("descriptive error message");
         }
     }

- Use the ``throwing`` attribute to call functions that may throw:

  .. code-block:: c++

     const auto value = get<int>(throwing, params, "ParameterName");

Error Propagation
~~~~~~~~~~~~~~~~~

- To generate an exception for a function using ``Context`` to report
  errors, use a throwing handler, as follows:

  .. code-block:: c++

     auto or_raise = ctx.getThrowingFailureHandler();
     someFunction(ctx, ...) | or_raise;

Parameter Handling
-------------------

For a detailed guide on using ``Parameter``, ``Parameters``, and ``ParametersValidator``,
see the :ref:`parameter handling page <mfem_mgis_developer_guide_parameter>`.

1. **Use ``ParametersValidator``**: Prefer ``ParametersValidator`` over the deprecated
   ``checkParameters`` functions for new code.

2. **Validate Early**: Validate parameters at the beginning of functions or constructors.

3. **Type Safety**: Use template arguments to restrict parameter types when possible.

4. **Required Parameters**: Mark required parameters explicitly using the ``required`` option.

5. **Incompatible Parameters**: Use ``addIncompatibleParametersList`` to enforce mutual
   exclusivity between parameters.

6. **Accessing Parameters**:

   - Non-throwing mode (recommended):

     .. code-block:: c++

        auto ctx = Context{};
        const auto value = get<int>(ctx, params, "ParameterName");
        if (isInvalid(value)) {
            return ctx.registerErrorMessage("invalid parameter");
        }

   - Throwing mode (for constructors or functions with ``attributes::Throwing``):

     .. code-block:: c++

        const auto value = get<int>(throwing, params, "ParameterName");

7. **Checking Parameter Existence**: Use the ``contains`` function:

   .. code-block:: c++

      if (contains(params, "ParameterName")) {
          // Parameter exists
      }

Memory Management
-----------------

1. **Smart Pointers**: Prefer ``std::unique_ptr`` and ``std::shared_ptr`` over raw pointers.
   Use ``make_unique`` and ``make_shared`` helper functions.

2. **Ownership**: Be explicit about ownership. Use raw pointers or references for non-owning
   relationships.

3. **Move Semantics**: Use move semantics for efficient transfers of resources.
   Mark move constructors and move assignment operators as ``noexcept``.

4. **Copy Semantics**: If copying is expensive or not meaningful, delete the copy constructor
   and copy assignment operator.

Class Design
------------

1. **Use struct Over class**: Prefer ``struct`` for types where public members
   are declared first. This aligns with the project's convention of declaring public members
   before private or protected ones.

2. **Virtual Destructors**: All base classes must have a virtual destructor.
   Mark it as ``noexcept`` and ``= default`` if possible.

   .. code-block:: c++

      virtual ~MyBaseClass() noexcept = default;

2. **Override Specifier**: Always use the ``override`` specifier for virtual functions
   that override a base class method.

3. **Final Specifier**: Use the ``final`` specifier for classes or methods that should not
   be further derived or overridden.

4. **Default Member Functions**: Use ``= default`` for default constructors, destructors,
   copy constructors, and assignment operators when appropriate.

5. **Deleted Member Functions**: Use ``= delete`` to explicitly delete member functions
   that should not be used.

6. **Rule of Five**: If you define any of the copy constructor, copy assignment operator,
   move constructor, move assignment operator, or destructor, you should define all of them.

Testing
-------

1. **Test Framework**: Use the TFEL test framework for unit tests (see `TFELTests documentation <https://thelfer.github.io/tfel/web/TFELTests.html>`_).

2. **Test Naming**: Test files should be named ``<ClassOrFeature>Test.cxx``.

3. **Test Structure**: Each test case should be a class inheriting from ``tfel::tests::TestCase``.

4. **Test Assertions**: Use the ``TFEL_TESTS_CHECK``, ``TFEL_TESTS_CHECK_EQUAL``,
   ``TFEL_TESTS_ASSERT``, and similar macros for assertions.

5. **Error Handling in Tests**: Tests should verify both success and failure cases,
   including error messages.

Documentation
-------------

1. **Doxygen**: All public APIs must be documented using Doxygen comments.

2. **Examples**: Provide usage examples in the documentation where helpful.

3. **Cross-References**: Use ``\see``, ``\ref``, and ``\sa`` to link to related documentation.

4. **Code Comments**: Use inline comments sparingly. Prefer self-documenting code.
   When necessary, use ``//`` for short comments and ``/*! ... */`` for Doxygen comments.

Performance
-----------

1. **Avoid Copies**: Use references, pointers, or move semantics to avoid unnecessary copies.

2. **Reserve Capacity**: Reserve capacity for containers (e.g., ``std::vector``) when the size
   is known in advance.

3. **Algorithms**: Prefer standard library algorithms (e.g., ``std::sort``, ``std::find``) over
   hand-written loops when appropriate.

4. **Profiling**: Use the ``MGIS/Profiling.hxx`` utilities for profiling critical sections.

Modern C++ Features
-------------------

1. **Use Modern Features**: Leverage C++20 features such as:

   - Concepts
   - Ranges
   - ``std::span``
   - Designated initializers
   - ``[[nodiscard]]`` attribute
   - ``std::optional``
   - ``std::variant``

2. **Avoid Legacy Features**: Avoid C-style casts, raw arrays, and manual memory management.

3. **Type Safety**: Use ``enum class`` instead of plain ``enum`` for type safety.

4. **Constants**: Use ``constexpr`` for compile-time constants.

Miscellaneous
-------------

1. **Boolean Naming**: Use ``is``, ``has``, or ``can`` prefixes for boolean functions:

   .. code-block:: c++

      bool isValid(...) noexcept;
      bool hasFeature(...) const noexcept;
      bool canDoSomething(...) const noexcept;

2. **Getter/Setter Naming**: Use the noun form for getters and ``set`` prefix for setters:

   .. code-block:: c++

      auto getValue() const noexcept;
      void setValue(...) noexcept;

3. **Avoid Abbreviations**: Use descriptive names. Avoid abbreviations unless they are
   widely understood (e.g., ``ctx`` for context, ``params`` for parameters).

4. **Consistency**: Follow the existing naming and coding conventions in the codebase.

5. **Line Length**: Keep lines under 80 characters where possible.

6. **Trailing Commas**: Use trailing commas in lists, parameter lists, and similar constructs
   for easier diffs and version control.

7. **Initialization**: Prefer uniform initialization (braces ``{}``) over parentheses ``()``.
