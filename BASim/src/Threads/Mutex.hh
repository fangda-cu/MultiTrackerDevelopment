/**
 * \file
 * \author  Luca Fascione <lukes@wetafx.co.nz>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Mutex
 *
 * \warning (C) Weta Digital, 2010
 */

#ifndef THREADS_MUTEX_HH
#define THREADS_MUTEX_HH

#ifndef _MSC_VER
#include <tr1/memory>
#else
#include <memory>
#endif

namespace BASim
{
namespace threads
{

/// Mutex class
class Mutex
{
public:
    Mutex();
    ~Mutex();

    void Lock();
    void Unlock();
    bool TryLock();
    long int Owner() const;

private:
    struct Impl;

    std::tr1::shared_ptr<Impl> m_impl;
};

/// Scoped lock
template<typename LockableT>
class ScopedLock
{
public:
    explicit ScopedLock(LockableT* lockable) :
        m_lockable(lockable)
    {
        m_lockable->Lock();
    }
    ~ScopedLock()
    {
        m_lockable->Unlock();
    }
    inline bool Acquired() const
    {
        return true;
    }
    inline long int Owner() const
    {
        return m_lockable->Owner();
    }

private:
    // Copy forbidden
    ScopedLock(const ScopedLock& dup) :
        m_lockable(NULL)
    {
    }
    ScopedLock& operator=(const ScopedLock& dup)
    {
        m_lockable = NULL;
    }

    LockableT* m_lockable;
};

/// Scoped TryLock
template<typename LockableT>
class ScopedTryLock
{
public:
    explicit ScopedTryLock(LockableT* lockable) :
        m_Acquired(false), m_lockable(lockable)
    {
        m_Acquired = m_lockable->TryLock();
    }
    ~ScopedTryLock()
    {
        if (m_Acquired)
            m_lockable->Unlock();
    }
    inline bool Acquired() const
    {
        return m_Acquired;
    }
    inline long int Owner() const
    {
        return m_lockable->Owner();
    }

private:
    // Copy forbidden
    ScopedTryLock(const ScopedTryLock& dup) :
        m_Acquired(false), m_lockable(NULL)
    {
    }
    ScopedTryLock& operator=(const ScopedTryLock& dup)
    {
        m_Acquired = false;
        m_lockable = NULL;
    }

    bool m_Acquired;
    LockableT* m_lockable;
};

} // namespace threads
} // namespace weta


#endif
