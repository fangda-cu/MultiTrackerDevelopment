/**
 * \file
 * \author  Mark Leone <mleone@wetafx.co.nz>
 *
 * \section DESCRIPTION
 *
 * \warning Copyright 2011 Weta Digital.  All rights reserved.
 */
#ifndef TEXT_LOG_HH
#define TEXT_LOG_HH

#include "MsgInfo.hh"
#include <boost/thread/mutex.hpp>
#include <boost/timer.hpp>
#include <iosfwd>
#include <memory>
#include <tr1/unordered_set>
#include <tr1/unordered_map>

namespace BASim {

/**
   TextLog provides message logging with configurable severity levels and duplicate suppression.
   Each message has a severity level and an identifer.  Messages can be formatted via printf-style
   format strings:

       log->Write(ErrorMsg("E001"), "Expected value %f\n", x)
   
   Messages can also be formatted via streaming, allowing any value supporing the << operator to
   be written to a log:

       ErrorStream(log, "E001") << "Expected value " << x << "\n";

   By default, messages with low severity are suppressed.  The severity level of a message can be
   changed (e.g. to suppress a warning, or turn a warning into an error).  A message can also
   specify an optional frequency (see MsgInfo::Frequency), allowing duplicates to be suppressed
   either based on the message id or its full text.
   
   TextLog is thread-safe, and care is taken to avoid commingling the text of messages from
   different threads.
*/
class TextLog {
public:
    /// Construct a log that writes to the given stream, suppressing messages below the specified
    /// level of severity.
    explicit TextLog(std::ostream& stream,
                     MsgInfo::Severity level = MsgInfo::kDefaultSeverity,
                     bool useTimestamps = false) :
        m_stream(stream),
        m_level(level),
        m_useTimestamps(useTimestamps)
    { }

    /// Disable the message with the given id.
    void Disable(const MsgInfo::Id& id) { SetSeverity(id, MsgInfo::kSilent); }

    /// Override the severity of the message with the specified id.
    void SetSeverity(const MsgInfo::Id& id, MsgInfo::Severity severity);

    /// A log stream wraps a log in scoped lock for the duration of a single message.  Any value
    /// that supports operator<< can be written to a log stream.  A log stream can optionally
    /// discard all messages.
    class Stream {
    public:
        /// Construct a log stream from the given log for a message with the given info.  The
        /// stream might all output if the message is suppressed.  Otherwise the log is locked via
        /// a scoped lock.
        Stream(TextLog* log, const MsgInfo& info)
        {
            m_out = log->getStream(info, &m_lock);
        }

        /// Any value that supports operator<< can be written to a log.
        template<typename T>
        Stream& operator<<(const T& value)
        {
            if (m_out) {
                *m_out << value;
            }
            return *this;
        }
    private:
        std::ostream* m_out;    // NULL if output is suppressed.
        boost::mutex::scoped_lock m_lock;

        // Copying and assignment are prohibited.
        Stream(const Stream&);
        Stream& operator=(const Stream&);
    };

    /// Write a message using printf-style formatting, prefixed with its severity and message id.
    /// The message might be suppressed based on its severity and message id (see MsgInfo.hh).
    void Write(const MsgInfo& info, const char* format, ...) 
        __attribute__ ((format (printf, 3, 4)));

    /// Write a string, without any locking or allocation (for signal handling safety).  The
    /// message is prefixed with the severity and message id.  Duplicate messages are not
    /// suppressed.
    void WriteSafely(const MsgInfo& info, const char* msg) const;

protected:
    // Streams need access to protected methods.
    friend class Stream;

    /// Prepare the output stream for a message with the given info.  Returns NULL if the message
    /// should be suppressed.  Otherwise the log is locked and the message info is optionally
    /// written (i.e. severity and message id).
    std::ostream* getStream(const MsgInfo& info, 
                            boost::mutex::scoped_lock* lock,
                            bool writeInfo = true);

private:
    // For now we format messages with vsnprintf using a fixed-sized local buffer.
    static const unsigned kMaxMsgSize = 4096;

    // Most operations lock the log for their duration.
    mutable boost::mutex m_mutex;

    // The output stream.
    std::ostream& m_stream;

    // Messages below this level of severity are suppressed.
    MsgInfo::Severity m_level;

    // Message severity can be overridden on a per-message basis (generalizing -Woff, -Werror).
    typedef std::tr1::unordered_map<MsgInfo::Id, MsgInfo::Severity> SeverityMap;
    SeverityMap m_severities;

    // When requested, the message id is saved and subsequent messages with the same id are
    // suppressed.
    std::tr1::unordered_set<MsgInfo::Id> m_prevIds;

    // When requested, the message text is saved and subsequent identical messages are suppressed.
    std::tr1::unordered_set<std::string> m_prevMsgs;

    // Timer for optional timestamps.
    boost::timer m_timer;
    bool m_useTimestamps;

    // Write message info, including severity, id, and optional timestamp.
    void writeInfo(const MsgInfo& info) const;
};

/// This template provides a convenient syntax for constructing a log stream for a message with a
/// fixed severity.
template<MsgInfo::Severity severity>
class LogStream : public TextLog::Stream {
public:
    LogStream(TextLog* log,
              const MsgInfo::Id& id, 
              MsgInfo::Frequency frequency = MsgInfo::kAlways) :
        TextLog::Stream(log, MsgInfo(severity, id, frequency))
    { }
};

typedef LogStream<MsgInfo::kSilent> SilentStream;
typedef LogStream<MsgInfo::kTrace> TraceStream;
typedef LogStream<MsgInfo::kDebug> DebugStream;
typedef LogStream<MsgInfo::kCopious> CopiousStream;
typedef LogStream<MsgInfo::kInfo> InfoStream;
typedef LogStream<MsgInfo::kNotice> NoticeStream;
typedef LogStream<MsgInfo::kWarning> WarningStream;
typedef LogStream<MsgInfo::kError> ErrorStream;
typedef LogStream<MsgInfo::kSevere> SevereStream;
typedef LogStream<MsgInfo::kForced> ForcedStream;

extern TextLog* g_log;

}

#endif // ndef TEXT_LOG_HH
