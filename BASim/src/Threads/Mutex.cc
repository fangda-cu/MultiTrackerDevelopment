/**
 * \file
 * \author  Luca Fascione <lukes@wetafx.co.nz>
 * \version 1.0
 *
 * \section DESCRIPTION
 *
 * Mutex implementation
 *
 * \warning (C) Weta Digital, 2010
 */

#include "Mutex.hh"
#include <pthread.h>
#include <sys/types.h>
#ifndef _MSC_VER
#include <sys/syscall.h>
//#include <linux/unistd.h> //OS_LINUX
#include <unistd.h> //OS_MACOSX // TODO:THREADING
#else
#include <Windows.h>
#endif

#include <errno.h>
#include <iostream>
#include <assert.h>

#ifdef WETA_DEBUG
#define CHECK_STATUS( status, func_name ) \
{ if( status != 0 ) { std::cerr << "ERROR: thread "<< syscall(__NR_gettid) << ", " #func_name " returned " << status << std::endl; } }
#else
#define CHECK_STATUS( status, func_name )
#endif

namespace BASim
{
namespace threads
{

/// Mutex class
struct Mutex::Impl
{
    Impl() :
        m_owner(0)
    {
// TODO:THREADING      
      
//#ifdef WETA_DEBUG
//        int status =
//#endif
//        pthread_mutexattr_init(&m_attr);
//        CHECK_STATUS( status, pthread_mutexattr_init );
//#ifdef WETA_DEBUG
//        status = pthread_mutexattr_settype( &m_attr, PTHREAD_MUTEX_ERRORCHECK_NP );
//        CHECK_STATUS( status, pthread_mutexattr_settype );
//#endif
//
//#ifdef WETA_DEBUG
//        status =
//#endif
//        pthread_mutex_init(&m_mutex, &m_attr);
//        CHECK_STATUS( status, pthread_mutex_init );
    }
    ~Impl()
    {
// TODO:THREADING
      
//#ifdef WETA_DEBUG
//        int status =
//#endif
//        pthread_mutex_destroy(&m_mutex);
//        CHECK_STATUS( status, pthread_mutex_destroy );
    }

    pthread_mutex_t m_mutex;
    pthread_mutexattr_t m_attr;
    long int m_owner;
};

Mutex::Mutex() :
    m_impl(new Impl)
{
}

Mutex::~Mutex()
{
}

long int Mutex::Owner() const
{
    return m_impl->m_owner;
}

void Mutex::Lock()
{
// TODO:THREADING
  
//#ifdef WETA_DEBUG
//    int status =
//#endif
//    pthread_mutex_lock(&m_impl->m_mutex);
//    CHECK_STATUS( status, pthread_mutex_lock );
//    //m_impl->m_owner = mach_thread_self(); //OS_MACOSX
//    m_impl->m_owner = syscall(__NR_gettid); //OS_LINUX
}

void Mutex::Unlock()
{  
// TODO:THREADING
  
//#ifdef WETA_DEBUG
//    int status =
//#endif
//    pthread_mutex_unlock(&m_impl->m_mutex);
//    CHECK_STATUS( status, pthread_mutex_unlock );
//    m_impl->m_owner = 0; // WARNING: assumes the mutex is not recursive
}

bool Mutex::TryLock()
{
  
// TODO:THREADING
  
//    int status = pthread_mutex_trylock(&m_impl->m_mutex);
//    if (status == 0)
//    {
//        //m_impl->m_owner = mach_thread_self(); //OS_MACOSX
//        m_impl->m_owner = syscall(__NR_gettid); //OS_LINUX
//        return true;
//    }
//    if (status == EBUSY)
//        return false; CHECK_STATUS( status, pthread_mutex_trylock );
//    return false;
   return false;
}

} // namespace threads
} // namespace weta
