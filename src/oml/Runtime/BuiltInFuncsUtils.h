/**
* @file BuiltInFuncsUtils.h
* @date November 2015
* Copyright (C) 2015-2022 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair's dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair's trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/

#ifndef __BUILTINFUNCSUTILS__
#define __BUILTINFUNCSUTILS__

#include <deque>

#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
//!
//! \brief Utility class for helper methods used by built-in functions
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsUtils
{
public:
    //!
    //! Constructor
    //!
    BuiltInFuncsUtils() {}
    //!
    //! Destructor
    //!
    ~BuiltInFuncsUtils() {}

    //!
    //! Returns true if file exists - works with wide characters
    //! \param name Given file name
    //
    static bool FileExists( const std::string& name);
    //
    //! Gets file names matching pattern in current directory
    //! \param pattern Pattern to match
    //!
    static std::vector<std::string> GetMatchingFiles( const std::string& pattern);
    //!
    //! Gets current working directory
    //!
    static std::string GetCurrentWorkingDir();
    //!
    //! Returns absolute path
    //! \param path Given path
    //!
    static std::string GetAbsolutePath( const std::string& path);
    //!
    //! Gets relative path
    //! \param path Given path
    //! \param dir  Parent directory
    //!
    static std::string GetRelativePath(const std::string& path,
                                       const std::string& dir);
    //!
    //! Adds trailing slash
    //! \param path Given path
    //!
    static void AddTrailingSlash(std::string& path);
    //!
    //! Strips trailing slash
    //! \param path Given path
    //!
    static void StripTrailingSlash(std::string& path);
    //!
    //! Strips trailing newline/endline
    //! \param str Given string
    //!
    static void StripTrailingNewline(std::string& str);
    //!
    //! Returns normalized path for the operating system
    //! \param path Given path
    static std::string Normpath( const std::string& path);
    //! Returns environment variable
    //! \param name Environment variable name
    //!
    static std::string GetEnv( const std::string& name);
    //!
    //! Returns true if given file path is absolute
    //! \param path Given path
    //!
    static bool IsAbsolutePath( const std::string& path);
    //! Returns true if this is the root dir (for windows)
    //! \param path Given path
    //!
    static bool IsRootDir( const std::string& path);
    //!
    //! Returns true if the given absolute path is a directory
    //! \param path Given absolute path
    //!
    static bool IsDir( const std::string& path);
    //!
    //! Returns true if the given absolute path is a file
    //! \param path Given absolute path
    //!
    static bool IsFile(const std::string& path);
    //! Returns base directory of the given path/file
    //! \param path Given path/file
    //!
    static std::string GetBaseDir( const std::string& path);
    //!
    //! Returns the base name if a file or the directory name for the given path
    //! \param path Given path/file
    //!
    static std::string GetBaseName( const std::string& path);
    //!
    //! Returns extension for a given file name
    //! \param path Given file name
    //!
    static std::string GetFileExtension( const std::string& filename);

    //!
    //! Returns true if the given double is an integer
    //! \param val Given value
    //!
    static bool IsInt( double val);
    //!
    //! Checks char bounds and returns valid char
    //! \param eval    Evaluator interface
    //! \param val     Given value
    //! \param showerr True if error needs to be thrown
    //!
    static unsigned char GetValidChar( EvaluatorInterface& eval, 
                                       double              val,
                                       bool                showerr = false);
    //!
    //! Gets ordered string
    //! \param cur Given currency
    //!
    static std::string GetOrderedStringVal( const Currency& cur);
    //!
    //! Checks math status and throws error/reports warning if not ok
    //! \param eval Evaluator interface
    //! \param stat Status 
    //!
    static void CheckMathStatus( EvaluatorInterface& eval, 
                                 hwMathStatus        stat);
    //!
    //! Returns currency after reading matrix row
    //! \param eval  Evaluator interface
    //! \param input Input currency
    //! \param row   Row to read
    //!
    static Currency ReadRow( EvaluatorInterface& eval,
                             const Currency&     input, 
                             int                 row);
    //!
    //! Sets warning message
    //! \param eval Evaluator interface
    //! \param str  String to print 
    //!
    static void SetWarning( EvaluatorInterface& eval, 
                            const std::string&  str);
    //!
    //! Returns currency which gives formatted output of strings. Strings will be 
    //! padded to width of largest value
    //! \param in             Vector of strings
    //! \param handleEmptyStr True if empty strings are set to NaN
    //! \param handleNaN      True if special padding is needed for NaN
    //! \param alignright     True if right align, else it is left aligned
    //! \param pad            Character to pad
    //!
    static Currency FormatOutput( const std::vector<std::string>& in,
                                  bool                            handleEmptyStr,
                                  bool                            handleNaN,
                                  bool                            alignright,
                                  char                            pad);
    //!
    //! Utility function which reads currency for sizes, returns true on success
    //! \param cur      Positive integer, inf or real matrix specifying size
    //! \param varindex Index of size currency, needed for throwing errors
    //! \param rows     Number of rows
    //! \param cols     Number of columns
    //! \param isMtx    True if size is specified in a matrix
    //!
    static void GetSizeSpecifications( const Currency& cur,   
                                       int             varindex,
                                       double&         rows, 
                                       double&         cols,
                                       bool&           isMtx);
    //!
    //! Parses given string input and gets a vector of format options
    //! \param String which contains format descriptions eg '%f %s'
    //! \param Valid formats
    //! \param True if format was valid
    //!
    static std::vector<std::string> GetFormats(const std::string&, const std::string&, bool&);
    //!
    //!
    //! Utility to set slice in an ND matrix
    //! \param in    Given matrix
    //! \param index Slice index (0-based)
    //! \param lhs   MatrixN to which slice is set
    //!
    static void SetMatrixNSlice( const hwMatrix* mtx,
                                 size_t          index,
                                 hwMatrixN*      lhs);
    //!
    //! True for double values like Nan/Inf/-Inf, which ignore standard prinf format
    //! \param val Value to print
    //!
    static bool NeedsSpecialPrintfFormat( double val);
    //!
    //! Gets string for doubles like Nan/Inf/-Inf, which ignore standard prinf format
    //! \param format Printf format
    //! \param val    Value to print
    //!
    static std::string GetSpecialPrintfFormat( const std::string& format,
                                               double             val);
    //!
    //! Returns true if there is a toolbox function of the given name
    //! \param eval Evaluator interface
    //! \param name Name of the function
    //!
    static bool IsToolboxFunction( EvaluatorInterface eval, 
                                   const std::string& name);
    //!
    //! Throws an error if the given file index is invalid
    //! \param eval            Evaluator interface
    //! \param fileid          File id
    //! \param argindex        Input argument index
    //! \param checkstdstreams True if std streams need to be checked
    //!
    static void CheckFileIndex( EvaluatorInterface eval,
                                int                fileid,
                                int                argindex,
                                bool               checkstdstreams);
    //!
    //! Gets file id from input currency
    //! \param eval Evaluator interface
    //! \param cur  Input currency
    //! \param idx  Input currency index
    //!
    static int GetFileId( EvaluatorInterface eval,
                          const Currency&    cur,
                          int                idx = 1);
    //!
    //! Returns true if std::cout buffer needs to be flushed
    //! \param Evaluator interface
    //! \param File descriptor
    static bool IsFlushCout(EvaluatorInterface,
                            int);
    //!
    //! Returns true after parsing input and gets formats
    //! \param in        Input string
    //! \param basefmts  Base format for the specific formats
    //! \param rawfmts   Raw format specified
    //! \param customFmt True if custom formats are allowed
    //! \param err       Error string
    //!
    bool GetFormats(const std::string&        in,
                    std::vector<std::string>& basefmts,
                    std::vector<std::string>& rawfmts,
                    std::string&              err);
    //!
    //! Returns true after parsing input and gets formats
    //! \param in        Input string
    //! \param basefmts  Base format for the specific formats
    //! \param rawfmts   Raw format specified
    //! \param usefmts   True if this format is used in output
    //! \param err       Error string
    //!
    bool GetFormats(const std::string&        in,
                    std::vector<std::string>& basefmts,
                    std::vector<std::string>& rawfmts,
                    std::vector<bool>&        usefmts,
                    std::string&              err);
    //!
    //! Returns string trimmed from leading characters
    //! \param in Input string
    //! \param trim Chars to trim
    //!
    static std::string LTrim(const std::string& in,
                             const std::string& trim = " \t");
    //!
    //! Trims trailing characters
    //! \param in Input string
    //! \param trim Chars to trim
    //!
    static std::string RTrim(const std::string& in,
                             const std::string& trim = " \t");

    //!
    //! Throws regex error
    //! \param code Error code
    //!
    void ThrowRegexError(std::regex_constants::error_type code);

    // utf8 support utilities
    //! 
    //! Returns normalized path for operating system, supports Unicode
    //! \param path Given path
    //!
    static std::wstring GetNormpathW(const std::wstring& path);
    //!
    //! Returns true if there are wide characters
    //! \param Input string
    //!
    static bool HasWideChars(const std::string&);
    //!
    //! Returns true if a real matrix has wide characters
    //! \param Matrix
    //!
    static bool HasWideChars(const hwMatrix*);
#ifdef OS_WIN
    //!
    //! Gets absolute path, supports unicode on Windows
    //! \param input Input string
    //!
    std::wstring GetAbsolutePathW(const std::wstring& input);
    //!
    //! Returns true if given path is absolute, supports unicode on Windows
    //! \param input Input string
    //!
    bool IsAbsolutePathW(const std::wstring& input);

    //!
    //! Gets current working directory on windows
    //!
    std::wstring GetCurrentWorkingDirW();
    //!
    //! Strips trailing slash for wide strings
    //! \param in Input string
    //!
    std::wstring StripTrailingSlashW(const std::wstring& in);
    //!
    //! Strips trailing newlines for wide strings
    //! \param in Wide string input
    //!
    std::wstring StripTrailingNewlineW(const std::wstring& in);
    //! 
    //! Returns true if normalized path exists on disk, supports Unicode
    //! \param path Given path
    //!
    bool DoesPathExistW(const std::wstring& path);
    //!
    //! Adds a trailing slash, supports unicode
    //! \param path Given path
    //!
    void AddTrailingSlashW(std::wstring& path);

    static std::wstring GetFileContentsW(const std::wstring&);

#endif
    //!
    //! Returns extension for a given file name
    //! \param Given file name
    //!
    std::wstring GetFileExtensionW(const std::wstring&);

    //!
    //! Converts std::string to std::wstring
    //! \param in Input string
    //!
    static std::wstring StdString2WString(const std::string& in);
    //!
    //! Converts std::wstring to std::string, supports Unicode
    //! \param in Wide string input
    //!
    static std::string WString2StdString(const std::wstring& in);

    //!
    //! Returns true if given path exists on disk, supports Unicode
    //! \param path Given path
    //!
    bool DoesPathExist(const std::string& path);

    //!
    //! Returns path after stripping multiple slashes and normalizing for OS
    //! \param path Given path
    //!
    std::string StripMultipleSlashesAndNormalize(const std::string& path);
    //!
    //! Create a chained display id so that all displays with the same id can be
    //! chained and deleted even if one of them in the chain is deleted. Used 
    //! for displaying multiline strings
    //!
    long CreateChainedDisplayId();
    //!
    //! Helper method to set environment variable
    //! \param name Env variable to set
    //! \param val  Value to set
    //!
    void SetEnvVariable(const std::string& name,
                        const std::string& val);
    //!
    //! Returns true if file is encoded
    //! \param eval Evaluator interface
    //! \param fid  File id
    //!
    bool IsFileEncoded(EvaluatorInterface eval,
                       int                fid);

    //!
    //! Gets size of widest character - works with Unicode
    //! \param str Input string
    //!
    size_t GetWidestCharSize(const std::string& str);
    //!
    //! Sets pagination mode with boolean for compatibility with previous versions
    //! \param Boolean value
    //!
    void SetPaginate(bool);
    //!
    //! Gets pagination mode as a string value
    //!
    std::string GetPaginate() const;
    //!
    //! Sets pagination mode with string value
    //! \param String value
    //!
    void SetPaginate(const std::string&);
    //!
    //! Creates a temporary folder and returns its name
    //!
    static std::string CreateTempFolder();
    //!
    //! Creates the directory, if it does not exist
    //! \param dir Directory name
    //!
    static void Mkdir(const std::string& dir);

    //!
    //! Returns file pointer after opening in given mode, supports Unicode
    //! \param File name
    //! \param Mode
    //!
    std::FILE* FileOpen(const std::string&,
                        const std::string&);
    //!
    //! Returns true if output log is open
    //!
    static bool IsOutputLogOpen();
    //!
    //! Opens output log in append mode
    //!
    static void OpenOutputLogForAppend();
    //!
    //! Closes output log
    //!
    static void CloseOutputLog();
    //!
    //! Saves string to output log
    //! \param Given string
    //!
    static void SaveToOutputLog(const std::string&);
    //!
    //! Returns true successful in converting string to scalar/complex
    //! \param in       Input string
    //! \param rval     Real part
    //! \param ival     Imaginary part, if it exists
    //! \param isscalar True if value is scalar
    //!
    static bool ReadNumber(const std::string& in,
                           double&            rval,
                           double&            ival,
                           bool&              isscalar);
    //!
    //! Helper method to get plus operator
    //!
    static int GetPlusOperator();
    //!
    //! Helper method to get minus operator
    //!
    static int GetMinusOperator();
    //!
    //! Reads file contents if it exists, supports unicode. No errors are thrown if file is invalid
    //! \param Filename
    //!
    static std::string GetFileContents(const std::string&);
    //!
    //! Clears variable with default name "ans" from the current scope in the evaluator
    //! \param Evaluator Interface
    //! \param Variable name
    //! 
    static void ClearAnsVariable(EvaluatorInterface, const std::string&);
    //!
    //! Converts deque to Currency - Replaces ContainerToMatrix
    //! \param Container of doubles
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Deque2Currency(const std::deque<double>&, bool = true);
    //!
    //! Converts deque to Currency - Replaces ContainerToMatrix
    //! \param Container of hwComplex
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Deque2Currency(const std::deque<hwComplex>&, bool = true);
    //!
    //! Converts deque to Currency
    //! \param Container of hwComplex
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Deque2Currency(const std::deque<std::string>&, bool = true);
    //!
    //! Converts vector to Currency - Replaces ContainerToMatrix
    //! \param Vector of doubles
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Vector2Currency(const std::vector<double>&, bool = true);
    //!
    //! Converts vector to Currency - Replaces ContainerToMatrix
    //! \param Vector of hwComplex
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Vector2Currency(const std::vector<hwComplex>&, bool = true);
    //!
    //! Converts vector to Currency - Replaces ContainerToMatrix
    //! \param Vector of ints
    //! \param True if resulting matrix needs to only have one row
    //! 
    static Currency Vector2Currency(const std::vector<int>&, bool = true);
    //!
    //! Converts string vector to cell array Currency - Replaces containerToCellArray
    //! \param Vector of strings
    //! \param True if resulting cell array needs to only have one row
    //! 
    static Currency Vector2Currency(const std::vector<std::string>&, bool = true);

private: 
    //!
    //! Returns true if strtod conversion is successful
    //! \param in  Input string
    //! \param end End string
    //! \param val Value converted by strtod
    //!
    static bool IsValidStrtodResult(const std::string& in,
                                    const std::string& end,
                                    double             val);
};

#endif // __BUILTINFUNCSUTILS__


