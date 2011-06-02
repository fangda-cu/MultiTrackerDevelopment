/**
 * \file
 * \author  Mark Leone <mleone@wetafx.co.nz>
 *
 * \section DESCRIPTION
 *
 * \warning Copyright 2011 Weta Digital.  All rights reserved.
 */

#include "TextLog.hh"
#include <boost/thread/locks.hpp>
#include <cstdarg>
#include <iostream>

namespace BASim {

// Override the severity of the message with the specified id.
void 
TextLog::SetSeverity(const MsgInfo::Id& id, MsgInfo::Severity severity)
{
    boost::lock_guard<boost::mutex> lock(m_mutex);
    m_severities[id] = severity;
}


// Write message info, including severity, id, and optional timestamp.
void
TextLog::writeInfo(const MsgInfo& info) const
{
    // Write timestamp.
    if (m_useTimestamps)
        m_stream << "[" << m_timer.elapsed() << " sec] "; 

    // We're careful about formatting if the severity or message id are empty.
    const char* severityStr = MsgInfo::SeverityToString(info.GetSeverity());
    if (*severityStr && !info.GetId().empty()) {
        m_stream << severityStr << ' ' << info.GetId() << ": ";
    }
    else if (*severityStr) {
        m_stream << severityStr << ": ";
    }
    else if (!info.GetId().empty()) {
        m_stream << info.GetId() << ": ";
    }
}

// Prepare the output stream for a message with the given info.  Returns NULL if the message
// should be suppressed.  Otherwise the log is locked and the message info is written
// (i.e. severity and message id).
std::ostream* 
TextLog::getStream(const MsgInfo& info, boost::mutex::scoped_lock* lock, bool writeInfo)
{
    boost::mutex::scoped_lock localLock(m_mutex);

    // Get the message severity, which might have been overridden
    SeverityMap::const_iterator it = m_severities.find(info.GetId());
    MsgInfo::Severity severity = it == m_severities.end() ? info.GetSeverity() : it->second;

    // Suppress if the severity doesn't meet the specified level.
    if (severity < m_level)
        return NULL;

    // Supressed multiple occurrences of the same message id if requested.
    if (info.GetFrequency() == MsgInfo::kOncePerId) {
        bool isNew = m_prevIds.insert(info.GetId()).second;
        if (!isNew)
            return NULL;
    }

    // Write the message info.
    if (writeInfo) 
        this->writeInfo(info);

    // Transfer the lock to the caller and return the output stream.
    lock->swap(localLock);
    return &m_stream;
}

// Write a message, optionally suppressing duplicates.  The message is prefixed with the severity.
void
TextLog::Write(const MsgInfo& info, const char* format, ...)
{
    // Check whether the message should be suppressed.  If not, a lock is acquired
    // (the message info is not written).
    boost::mutex::scoped_lock lock;
    std::ostream* stream = getStream(info, &lock, false);
    if (stream == NULL)
        return;

    // Format the message.  For now we use vsnprintf with a fixed-sized local buffer.
    // TODO: release lock while formatting?
    char msg[kMaxMsgSize];
    va_list ap;
    va_start(ap, format);
    vsnprintf(msg, kMaxMsgSize, format, ap);
    va_end(ap);

    // Suppress identical duplicate messages if requested.
    if (info.GetFrequency() == MsgInfo::kOncePerMessage) {
        bool isNew = m_prevMsgs.insert(msg).second;
        if (!isNew)
            return;
    }

    writeInfo(info);
    *stream << msg;
    return;
}

// Write a string, without any locking or allocation (for signal handling safety).  The
// message is prefixed with the severity.
void
TextLog::WriteSafely(const MsgInfo& info, const char* msg) const
{
    // TODO: actually, it's not safe to write to a std::ostream from a signal handler (nor flush).
    // We could use a file-descriptor-based stream implementation:
    // http://www.josuttis.com/cppcode/fdstream.html 
    if (info.GetSeverity() >= m_level) {
        writeInfo(info);
        m_stream << msg << std::endl;
    }
}


// Global variable. Must be initialized somewhere before any logging.
TextLog* g_log;

}
