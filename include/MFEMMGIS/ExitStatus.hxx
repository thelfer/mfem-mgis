/*!
 * \file   MFEM/MGIS/ExitStatus.hxx
 * \brief  This file declares the ExitStatus class.
 * \date   15/05/2023
 */

#ifndef LIB_MFEM_MGIS_EXIT_STATUS_HXX
#define LIB_MFEM_MGIS_EXIT_STATUS_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief a base class to describe the exit status of a function or a method.
   *
   * This class can be inherited to provide additional information to the
   * caller.
   */
  struct [[nodiscard]] ExitStatus {
    //! \brief this status indicates that the function succeeded
    static const ExitStatus success;
    /*!
     * \brief this status indicates that the function succeeded
     * but the results may be unreliable.
     *
     * For example, one may consider that if a phase transition
     * model predicts a total transition from one phase
     * to another during the time step, this results might be
     * unreliable and the time step shall be reduced
     * to capture all the relevant physical phenomena.
     *
     * \note unreliable results might occur during non
     * linear resolutions. So unreliable results shall only
     * be checked after convergence.
     */
    static const ExitStatus unreliableResults;
    /*!
     * \brief this status indicates that the function failed,
     * but that the caller could do something.
     *
     * For example, consider the case of a non linear behaviour.
     * The integration step may fail because the current estimate
     * of the strain increment is too large, which can happen
     * during the equilibrium iterations. In this case,
     * the non linear solver may use a linesearch-like
     * strategy and reduce the amplitude of the last
     * correction to the estimate of the displacement.
     */
    static const ExitStatus recoverableError;
    /*!
     * \brief this status indicates that the function failed,
     * and that the caller can not do something about it.
     *
     * For example, consider the case of a non linear behaviour
     * which depends on the temperature. If the temperature passed
     * to behaviour is negative (in Kelvin), then there is nothing
     * that the non linear solver can do.
     */
    static const ExitStatus unrecoverableError;
    //! \brief default constructor
    ExitStatus() noexcept;
    //! \brief move constructor
    ExitStatus(ExitStatus &&) noexcept = default;
    //! \brief copy constructor
    ExitStatus(const ExitStatus &) noexcept = default;
    //! \brief move assignment
    ExitStatus &operator=(ExitStatus &&) noexcept = default;
    //! \brief standard assignment
    ExitStatus &operator=(const ExitStatus &) noexcept = default;
    /*!
     * \brief constructor from a boolean
     * \param[in] s: if true, the status is `success`, otherwise
     * `unrecoverableError`
     */
    explicit ExitStatus(const bool s) noexcept;
    /*!
     * \brief default comparison operators
     * \return the result of the three-way comparison
     */
    auto operator<=>(ExitStatus const &) const = default;
    /*!
     * \brief change the current status. This is equivalent to the assignment
     * operator
     *
     * \param[in] s: new status
     */
    void setStatus(const bool s) noexcept;
    /*!
     * \brief change the current status. This is equivalent to the assignment
     * operator and the `setStatus` method
     *
     * \param[in] s: new status
     */
    void update(const bool s) noexcept;
    /*!
     * \return if the status is either `success` or `unreliableResults`
     *
     * This method is meant to test if the result of a function call
     * does not prevent from continuing the current computation process.
     * The user must however ensure that the distinction between
     * `success` and `unreliableResults` is propagated appropriately.
     */
    [[nodiscard]] bool shallContinue() const noexcept;
    /*!
     * \return if the status is either `recoverableError` or
     * `unrecoverableError`
     *
     * This method is meant to test if the result of a function call
     * is likely to stop the current computation and that the best thing
     * to do is to return to the parent function up to the point where
     * the error is treated.
     */
    [[nodiscard]] bool shallStop() const noexcept;
    /*!
     * \brief update the current status if the given status is worse than the
     * current one
     *
     * \param[in] o: status used to update the current one
     * \return the result of `shallContinue` after the update
     */
    bool update(const ExitStatus &o) noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \brief synchronize the object among MPI processes
     * \param[in] c: MPI communicator
     */
    void synchronize(const MPI_Comm c) noexcept;
#endif /* MFEM_USE_MPI */
    //! \brief destructor
    inline ~ExitStatus() noexcept = default;

   private:
    /*!
     * \brief internal data structure to avoid implicit conversion from integer
     * values
     */
    struct Status {
      //! \brief value of a success
      static constexpr size_type success = 0;
      //! \brief value of unreliable results
      static constexpr size_type unreliableResults = 1;
      //! \brief value of a recoverable error
      static constexpr size_type recoverableError = 2;
      //! \brief value of an unrecoverable error
      static constexpr size_type unrecoverableError = 3;
      //! \brief default comparison operators
      auto operator<=>(Status const &) const = default;
      //! \brief held value
      size_type value;
    };
    //! \brief constructor from a Status
    ExitStatus(const Status) noexcept;
    //! \brief exit status
    Status status = Status{Status::success};
  };  // end of class ExitStatus

  inline const ExitStatus ExitStatus::success =
      ExitStatus(Status{.value = Status::success});
  inline const ExitStatus ExitStatus::unreliableResults =
      ExitStatus(Status{.value = Status::unreliableResults});
  inline const ExitStatus ExitStatus::recoverableError =
      ExitStatus(Status{.value = Status::recoverableError});
  inline const ExitStatus ExitStatus::unrecoverableError =
      ExitStatus(Status{.value = Status::unrecoverableError});

#ifdef MFEM_USE_MPI

  /*!
   * \return an exit status common to all processes. This common status is taken
   * as the worst case value.
   *
   * \param[in] s: exit status of the current process
   * \param[in] c: MPI communicator
   */
  MFEM_MGIS_EXPORT ExitStatus synchronize(const ExitStatus s,
                                          const MPI_Comm c) noexcept;

#endif /* MFEM_USE_MPI */

}  // end of namespace mfem_mgis

namespace mgis::internal {

  //! \brief specialisation for exit status
  template <>
  struct InvalidValueTraits<mfem_mgis::ExitStatus> {
    //! \brief tag indicating that this class is properly specialized
    static constexpr bool isSpecialized = true;
    //! \return an invalid exit status
    static mfem_mgis::ExitStatus getValue() noexcept {
      return mfem_mgis::ExitStatus::unrecoverableError;
    }
  };

}  // end of namespace mgis::internal

#include "MFEMMGIS/ExitStatus.ixx"

#endif LIB_MFEM_MGIS_EXIT_STATUS_HXX
