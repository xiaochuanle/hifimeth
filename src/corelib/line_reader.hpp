#ifndef __LINE_READER_HPP
#define __LINE_READER_HPP

#include "hbn_aux.hpp"

#include <string>

class CBufferLineReader
{
public:
    /*
    Read from the file, "-" (but not "./") means standard input
    An explicit call to operator++ or
    ReadLine() will bbe necessaray to fetch the first line
    */
    CBufferLineReader(const char* filename);

    ~CBufferLineReader();

    bool            AtEOF() const;
    char            PeekChar() const;
    CBufferLineReader&  operator++();
    void            UngetLine();
    std::pair<const char*, size_t>     operator*() const;
    size_t          GetPosition() const;
    size_t          GetLineNumber() const;
    const char*     GetFileName() const { return m_FileName.c_str(); }

private:
    CBufferLineReader(const CBufferLineReader&);
    CBufferLineReader& operator=(const CBufferLineReader&);

private:
    void x_LoadLong();
    bool x_ReadBuffer();
    
private:
    std::string     m_FileName;
    gzFile          m_Reader;
    bool            m_Eof;
    bool            m_UngetLine;
    size_t          m_LastReadSize;
    size_t          m_BufferSize;
    char*           m_Buffer;
    const char*     m_Pos;
    const char*     m_End;
    std::pair<const char*, size_t>     m_Line;
    std::string     m_String;
    size_t          m_InputPos;
    size_t          m_LineNumber;
};

typedef CBufferLineReader HbnLineReader;

#endif // __LINE_READER_HPP