/**
 * \file
 * \author  Mark Leone <mleone@wetafx.co.nz>
 *
 * \section DESCRIPTION
 *
 * \warning Copyright 2011 Weta Digital.  All rights reserved.
 */
#ifndef MSG_INFO_HH
#define MSG_INFO_HH

#include <cassert>
#include <cstdarg>
#include <iosfwd>
#include <string>

namespace weta {
namespace logging {

class MsgInfo {
public:
    /// Message ids are represented as strings.
    typedef std::string Id;

    /// Message severity.
    enum Severity {
        kSilent,
        kTrace,
        kDebug,
        kCopious,
        kInfo,
        kNotice,
        kWarning,
        kError,
        kSevere,
        kForced,
    };

    /// Message frequency.
    enum Frequency {
        kOncePerId,
        kOncePerMessage,
        kAlways,
    };

    /// The default severity level.
    static const Severity kDefaultSeverity = kNotice;

    /// Convert severity level to string.
    inline static const char* SeverityToString(Severity severity)
    {
        switch (severity) {
          case kSilent: return "";
          case kTrace: return "Trace";
          case kDebug: return "Debug";
          case kCopious: return "Copious";
          case kInfo: return "Info";
          case kNotice: return "Notice";
          case kWarning: return "Warning";
          case kError: return "Error";
          case kSevere: return "Severe";
          case kForced: return "";
        }
        assert(false && "Unhandled severity");
        return "";
    }

    /// Construct a message identifier.
    MsgInfo(Severity severity,
            const Id& id, 
            Frequency frequency = kAlways) :
        m_id(id),
        m_severity(severity),
        m_frequency(frequency)
    { }

    /// Get the message id.
    const Id& GetId() const { return m_id; }

    /// Get message severity.
    Severity GetSeverity() const { return m_severity; }

    /// Get the message frequency.
    Frequency GetFrequency() const { return m_frequency; }

private:
    Id m_id;                    ///< The message identifier.
    Severity m_severity;        ///< Message severity.
    Frequency m_frequency;      ///< Message frequency 
};

template<MsgInfo::Severity severity>
class SeverityInfo : public MsgInfo {
public:
    explicit SeverityInfo(const Id& id, Frequency frequency = kAlways) :
        MsgInfo(severity, id, frequency)
    { }
};

typedef SeverityInfo<MsgInfo::kSilent> SilentMsg;
typedef SeverityInfo<MsgInfo::kTrace> TraceMsg;
typedef SeverityInfo<MsgInfo::kDebug> DebugMsg;
typedef SeverityInfo<MsgInfo::kCopious> CopiousMsg;
typedef SeverityInfo<MsgInfo::kInfo> InfoMsg;
typedef SeverityInfo<MsgInfo::kNotice> NoticeMsg;
typedef SeverityInfo<MsgInfo::kWarning> WarningMsg;
typedef SeverityInfo<MsgInfo::kError> ErrorMsg;
typedef SeverityInfo<MsgInfo::kSevere> SevereMsg;
typedef SeverityInfo<MsgInfo::kForced> ForcedMsg;

} // namespace logging
} // namespace weta

#endif // ndef MSG_INFO_HH
