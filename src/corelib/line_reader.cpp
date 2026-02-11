#include "line_reader.hpp"

#include <errno.h>

using namespace std;

/// Result codes for I/O operations.
/// @note
///     Exceptions (if any) thrown by IReader/IWriter interfaces should be
///     treated as unrecoverable errors (eRW_Error).
/// @sa
///   IReader, IWriter, IReaderWriter
enum ERW_Result {
    eRW_NotImplemented = -1,  ///< Action / information is not available
    eRW_Success = 0,          ///< Everything is okay, I/O completed
    eRW_Timeout,              ///< Timeout expired, try again later
    eRW_Error,                ///< Unrecoverable error, no retry possible
    eRW_Eof                   ///< End of data, should be considered permanent
};

CBufferLineReader::CBufferLineReader(const char* filename):
    m_Eof(false),
    m_UngetLine(false),
    //m_BufferSize(32*1024),
    m_BufferSize(4*1024*1024),
    m_Buffer(new char[m_BufferSize]),
    m_Pos(m_Buffer),
    m_End(m_Pos),
    m_InputPos(0),
    m_LineNumber(0) 
{
    hbn_gzopen(m_Reader, filename, "r");
    x_ReadBuffer();
}

CBufferLineReader::~CBufferLineReader()
{
    delete[] m_Buffer;
    hbn_gzclose(m_Reader);
}

bool CBufferLineReader::AtEOF() const 
{
    return m_Eof && !m_UngetLine && m_Pos >= m_End;
}

char CBufferLineReader::PeekChar() const 
{
    hbn_assert(!AtEOF());
    // if at EOF - undefined behavior
    if (AtEOF()) {
        return *m_Pos;
    }

    // if line was ungot - return its first symbol
    if (m_UngetLine) {
        // if line is empty - return 0
        if (!m_Line.second) {
            return 0;
        }
        return m_Line.first[0];
    }
    
    // if line is empty return 0
    if (*m_Pos == '\n' || *m_Pos == '\r') {
        return 0;
    }

    return *m_Pos;
}

void CBufferLineReader::UngetLine()
{
    hbn_assert(!m_UngetLine && m_Line.second);
    // if after Ungetline() or after constructor - noop
    if (m_UngetLine || m_Line.second == 0) {
        return;
    }
    --m_LineNumber;
    m_UngetLine = true;
}

CBufferLineReader& CBufferLineReader::operator++()
{
    // if at EOF - noop
    if (AtEOF()) {
        m_Line = std::pair<const char*, size_t>(nullptr, 0);
        return *this;
    }
    ++m_LineNumber;
    if (m_UngetLine) {
        hbn_assert(m_Line.second);
        m_UngetLine = false;
        return *this;
    }
    // check if we are at the buffer end
    const char* start = m_Pos;
    const char* end = m_End;
    for (const char* p = start; p < end; ++p) {
        if (*p == '\n') {
            m_Line = std::pair<const char*, size_t>(start, p - start);
            m_LastReadSize = p + 1 - start;
            m_Pos = ++p;
            if (p == end) {
                m_String.assign(m_Line.first, m_Line.second);
                m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
                x_ReadBuffer();
            }
            return *this;
        } else if (*p == '\r') {
            m_Line = std::pair<const char*, size_t>(start, p - start);
            m_LastReadSize = p + 1 - start;
            m_Pos = ++p;
            if (p == end) {
                m_String.assign(m_Line.first, m_Line.second);
                m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
                if (x_ReadBuffer()) {
                    p = m_Pos;
                    if (*p == '\n') {
                        m_Pos = p + 1;
                        ++m_LastReadSize;
                    }
                }
                return *this;
            }
            if (*p != '\n') {
                return *this;
            }
            ++m_LastReadSize;
            m_Pos = ++p;
            if (p == end) {
                m_String.assign(m_Line.first, m_Line.second);
                m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
                x_ReadBuffer();
            }
            return *this;
        }
    }
    x_LoadLong();
    return *this;
}

void CBufferLineReader::x_LoadLong()
{
    const char* start = m_Pos;
    const char* end = m_End;
    m_String.assign(start, end);
    while (x_ReadBuffer()) {
        start = m_Pos;
        end = m_End;
        for (const char* p = start; p < end; ++p) {
            char c = *p;
            if (c == '\r' || c == '\n') {
                m_String.append(start, p - start);
                m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
                m_LastReadSize = m_Line.second + 1;
                if (++p == end) {
                    m_String.assign(m_Line.first, m_Line.second);
                    m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
                    if (x_ReadBuffer()) {
                        p = m_Pos;
                        end = m_End;
                        if (p < end && c == '\r' && *p == '\n') {
                            ++p;
                            m_Pos = p;
                            ++m_LastReadSize;
                        }
                    }
                } else {
                    if (c == '\r' && *p == '\n') {
                        if (++p == end) {
                            x_ReadBuffer();
                            p = m_Pos;
                            ++m_LastReadSize;
                        }
                    }
                    m_Pos = p;
                }
                return;
            }
        }
        m_String.append(start, end - start);
    }
    m_Line = std::pair<const char*, size_t>(m_String.c_str(), m_String.size());
    m_LastReadSize = m_Line.second;
    return;
}


/// Read as many as "count" bytes into a buffer pointed to by the "buf"
/// argument.  Always store the number of bytes actually read (0 if read
/// none) via the pointer "bytes_read", if provided non-NULL.
/// Return non-eRW_Success code if EOF / error condition has been
/// encountered during the operation (some data may have been read,
/// nevertheless, and reflected in "*bytes_read").
/// Special case:  if "count" is passed as 0, then the value of "buf" must
/// be ignored, and no change should be made to the state of the input
/// device (but may return non-eRW_Success to indicate that the input
/// device has already been in an error condition).
static ERW_Result
x_load_data(gzFile stream, char* buffer, int count, int* bytes_read)
{
    ERW_Result result = eRW_Success;
    int n = gzread(stream, buffer, count);
    if (n < count) {
        if (gzeof(stream)) {
            result = eRW_Eof;
        } else {
            int e;
            const char* why = gzerror(stream, &e);
            if (e) HBN_ERR("%s", why);
        }
    }
    *bytes_read = n;
    return result;
}

bool CBufferLineReader::x_ReadBuffer() 
{
    m_InputPos += (m_End - m_Buffer);
    m_Pos = m_End = m_Buffer;
    while (1) {
        int size;
        ERW_Result result = x_load_data(m_Reader, m_Buffer, m_BufferSize, &size);
        switch (result) {
            case eRW_NotImplemented:
            case eRW_Error:
                HBN_ERR("Read Error");
                break;
            case eRW_Timeout:
                // keep spinning around
                break;
            case eRW_Eof:
                m_Eof = true;
                // fall through
            case eRW_Success:
                m_End = m_Pos + size;
                return eRW_Success || size > 0;
            default:
                hbn_assert(0);
        }
    }
    return false;
}

std::pair<const char*, size_t> CBufferLineReader::operator*() const 
{
    hbn_assert(!m_UngetLine);
    // after ungetline - undefined behavior
    if (m_UngetLine) {
        return std::pair<const char*, size_t>(nullptr, 0);
    }
    // Right ater contructor (m_LineNumber is 0 and UngetLine() was not run,
    // the latter was already checked) - return NULL
    if (m_Line.second == 0) {
        return std::pair<const char*, size_t>(nullptr, 0);
    }
    return m_Line;
}

size_t CBufferLineReader::GetPosition() const 
{
    size_t offset = m_Pos - m_Buffer;
    if (m_UngetLine) {
        offset -= m_LastReadSize;
    }
    return m_InputPos + offset;
}

size_t CBufferLineReader::GetLineNumber() const 
{
    /* Right after constructor (m_LineNumber is 0 and UngetLine() was not run)
       - 0 */
    /* If at EOF - returns the number of the last string */
    /* After UngetLine() - number of the previous string */
    /* Not at EOF, not after UngetLine() - number of the current string */
    return m_LineNumber;
}
