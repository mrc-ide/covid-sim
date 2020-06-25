/** \file  Messages.h
 *  \brief Message printing
 *
 * This file provides message printing routines along with logging levels to
 * enable verbosity of the output to be controlled by the user.
 *
 * Intended usage is as follows:
 *
 * Messages::out(Messages::Info) << "Insert message here " << ... << "\n";
 */

#ifndef COVIDSIM_MESSAGES_H_INCLUDED_
#define COVIDSIM_MESSAGES_H_INCLUDED_

#include <cstdlib>
#include <memory>
#include <iostream>

/* To save typing we want to import all the Messages::Type enum values into the
 *  Messages namespace.  In C++20 this can be done by using enum class and
 *  using enum, earlier on we just use a traditional namespace type.
 */
#if defined(__cpp_using_enum)
# define ENUM enum class
# define USING_ENUM(x) using enum x
#else
# define ENUM enum
# define USING_ENUM(x)
#endif

namespace Messages
{
  /** Type of messages, Less severe messages have higher numbers. */
  ENUM Type
  {
    Error = 0,        ///< Error messages
    Warning = 1,      ///< Warning messages
    Info = 2,         ///< Info messages
    Progress = 3,     ///< Messages which show progress status
    Debug = 4,        ///< Debug messages

    DefaultAbort = 0, ///< Default highest type to abort after printing
    Default = 3       ///< Default highest type to show.
  };
  USING_ENUM(Type);

  /** \brief Optional stream
   *
   * Effectively acts as an optional<ostream> but allows output directly to
   * it.
   *
   * Once a message that is at a level <= than the abort level has been
   * printed we will exit the program.
   *
   * \internal Implemented as a shared pointer to a reference counted
   * implementation class.  When the shared pointer count goes down to zero we
   * check to see if we should abort and do the abort then.
   */
  class OptStream
  {
  private:
    /** Internal implementation class. */
    class Impl
    {
    public:
      /** \brief           Constructor
       *  \param do_output Output text
       *  \param do_abort  Abort on destruction
       *  \param os        Stream to output to (if relevant)
       */
      Impl(bool do_output, bool do_abort, std::ostream& os) noexcept
        : os_(os), do_output_(do_output), do_abort_(do_abort)
      { }

      /* Implementation class can't be copied/moved or assigned.  */
      Impl(Impl const&) noexcept = delete;
      Impl& operator=(Impl const&) noexcept = delete;
      Impl(Impl&&) noexcept = delete;
      Impl& operator=(Impl&&) noexcept = delete;

      /** \brief Destructor.
       *
       * Will abort program if do_abort_ is true.
       */
      ~Impl() noexcept
      {
        if (do_abort_)
        {
          std::cerr << "ERROR: Terminating program\n";
          std::exit(2);
        }
      }

      /// \brief  Are we doing output?
      bool do_output() const noexcept { return do_output_; }

      /// \brief  Get the stream
      std::ostream& os() noexcept { return os_; }

    private:
      std::ostream& os_;    ///< Output stream
      bool do_output_;      ///< Should we do output?
      bool do_abort_;       ///< Should we abort?
    };

  public:
    /** \brief           Constructor
     *  \param do_output Output messages?
     *  \param do_abort  Abort once output is complete?
     *  \param os        Stream to output to (default stderr)
     */
    explicit OptStream(bool do_output = true, bool do_abort = true, std::ostream& os = std::cerr) noexcept
      : impl_(std::make_shared<Impl>(do_output, do_abort, os))
    { }

    /// \brief Are we doing output?
    bool do_output() const { return impl_->do_output(); }

    /// \brief Get the stream
    std::ostream& os() { return impl_->os(); }

    /** \brief       Start outputting
     *  \param  type Message type
     *  \return      Stream to ouput to.
     */
    static OptStream get(Type type)
    {
      return OptStream(type <= level_, type <= abort_level_);
    }

    /** \brief       Set output level
     *  \param level Level
     *
     * After this call all messages of type \a level and lower will be output,
     * but none of a type higher.
     */
    static void set_level(Type level) { level_ = level; }

    /** \brief       Set abort level
     *  \param level Level
     *
     * After this call all messages of type \a level and lower will cause an
     * the program to abort after printing.  Higher message types will allow
     * the program to carry on running.
     */
    static void set_abort_level(Type level) { abort_level_ = level; }

  private:
    static Type level_;         ///< Print level
    static Type abort_level_;   ///< Abort level

    std::shared_ptr<Impl> impl_;  ///< Shared pointer to implementation
  };

  /** \brief     Output operator
   *  \param  os Stream to output to
   *  \param  t  Object to output
   *  \return    Output stream
   */
  template<typename T>
  OptStream operator<<(OptStream os, T t)
  {
    if (os.do_output())
    {
      os.os() << t;
    }
    return os;
  }

  /** \brief        Get the output stream
   *  \param  level Message level
   *  \retrun       Output stream
   */
  inline OptStream out(Type level)
  {
    OptStream os = OptStream::get(level);
    switch (level)
    {
      case Error: os << "ERROR: "; break;
      case Warning: os << "WARNING: "; break;
      case Debug: os << "DEBUG: "; break;
      default: break;
    }

    return os;
  }
} // namespace Messages

#endif // COVIDSIM_MESSAGES_H_INCLUDED_
