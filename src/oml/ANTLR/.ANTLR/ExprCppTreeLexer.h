/** \file
 *  This C header file was generated by $ANTLR version 3.5
 *
 *     -  From the grammar source file : ExprCppTree.g
 *     -                            On : 2018-12-03 09:59:20
 *     -                 for the lexer : ExprCppTreeLexerLexer
 *
 * Editing it, at least manually, is not wise.
 *
 * C language generator and runtime by Jim Idle, jimi|hereisanat|idle|dotgoeshere|ws.
 *
 *
 * The lexer 
ExprCppTreeLexer

has the callable functions (rules) shown below,
 * which will invoke the code for the associated rule in the source grammar
 * assuming that the input stream is pointing to a token/text stream that could begin
 * this rule.
 *
 * For instance if you call the first (topmost) rule in a parser grammar, you will
 * get the results of a full parse, but calling a rule half way through the grammar will
 * allow you to pass part of a full token stream to the parser, such as for syntax checking
 * in editors and so on.
 *
 * The parser entry points are called indirectly (by function pointer to function) via
 * a parser context typedef pExprCppTreeLexer, which is returned from a call to ExprCppTreeLexerNew().
 *
 * As this is a generated lexer, it is unlikely you will call it 'manually'. However
 * the methods are provided anyway.
 *
 * The methods in pExprCppTreeLexer are  as follows:
 *
 *  - 
 void
      pExprCppTreeLexer->T__109(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__110(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__111(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__112(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__113(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__114(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__115(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->T__116(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->IF(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->WHILE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->FOR(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->SWITCH(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->CASE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->OTHERWISE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->GLOBAL(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->PERSISTENT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->RETURN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->END(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->TRY(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->FUNCTION(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->CLEAR(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->ELSEIF(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->OTHER(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->CLASSDEF(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->METHODS(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->PROPERTIES(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->EQUAL(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->NEQUAL(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->PLUS(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->MINUS(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->TIMES(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->DIV(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->ETIMES(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->EDIV(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LDIV(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->ELDIV(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->QUOTE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->EQUOTE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->COLON(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->POW(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LTHAN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->GTHAN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LEQ(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->GEQ(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->DOT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->DOTPOW(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->AND(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->OR(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LAND(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LOR(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->NEGATE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->SEMIC(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->COMMA(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->DOTDOT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LBRACKET(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->RBRACKET(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LPAREN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->RPAREN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LCURLY(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->RCURLY(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->PERCENT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->AMP(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->PCTLCURLY(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->RCURLYPCT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->ASSIGN(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->LETTER(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->UNICHAR(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->DIGIT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->IDENT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->IDENT2(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->UNICHAR_IDENT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->NUMBERHACK(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->NUMBER(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->HEXVAL(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->HEXHACK(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->NEWLINE(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->WS(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->CONT(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->HACK(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->HACKB(pExprCppTreeLexer)
 *  - 
 void
      pExprCppTreeLexer->Tokens(pExprCppTreeLexer)
 *
 * The return type for any particular rule is of course determined by the source
 * grammar file.
 */
// [The "BSD license"]
// Copyright (c) 2005-2009 Jim Idle, Temporal Wave LLC
// http://www.temporal-wave.com
// http://www.linkedin.com/in/jimidle
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef	_ExprCppTreeLexer_H
#define _ExprCppTreeLexer_H
/* =============================================================================
 * Standard antlr3 C runtime definitions
 */
#include    <antlr3.h>

/* End of standard antlr 3 runtime definitions
 * =============================================================================
 */

#ifdef __cplusplus
extern "C" {
#endif

// Forward declare the context typedef so that we can use it before it is
// properly defined. Delegators and delegates (from import statements) are
// interdependent and their context structures contain pointers to each other
// C only allows such things to be declared if you pre-declare the typedef.
//
typedef struct ExprCppTreeLexer_Ctx_struct ExprCppTreeLexer, * pExprCppTreeLexer;



	#include "../../Runtime/ANTLRoverride.h"


#ifdef	ANTLR3_WINDOWS
// Disable: Unreferenced parameter,							- Rules with parameters that are not used
//          constant conditional,							- ANTLR realizes that a prediction is always true (synpred usually)
//          initialized but unused variable					- tree rewrite variables declared but not needed
//          Unreferenced local variable						- lexer rule declares but does not always use _type
//          potentially unitialized variable used			- retval always returned from a rule
//			unreferenced local function has been removed	- susually getTokenNames or freeScope, they can go without warnigns
//
// These are only really displayed at warning level /W4 but that is the code ideal I am aiming at
// and the codegen must generate some of these warnings by necessity, apart from 4100, which is
// usually generated when a parser rule is given a parameter that it does not use. Mostly though
// this is a matter of orthogonality hence I disable that one.
//
#pragma warning( disable : 4100 )
#pragma warning( disable : 4101 )
#pragma warning( disable : 4127 )
#pragma warning( disable : 4189 )
#pragma warning( disable : 4505 )
#pragma warning( disable : 4701 )
#endif

/** Context tracking structure for 
ExprCppTreeLexer

 */
struct ExprCppTreeLexer_Ctx_struct
{
    /** Built in ANTLR3 context tracker contains all the generic elements
     *  required for context tracking.
     */
    pANTLR3_LEXER    pLexer;

     void
     (*mT__109)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__110)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__111)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__112)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__113)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__114)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__115)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mT__116)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mIF)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mWHILE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mFOR)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mSWITCH)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCASE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mOTHERWISE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mGLOBAL)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPERSISTENT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mRETURN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mEND)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mTRY)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mFUNCTION)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCLEAR)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mELSEIF)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mOTHER)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCLASSDEF)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mMETHODS)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPROPERTIES)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mEQUAL)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mNEQUAL)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPLUS)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mMINUS)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mTIMES)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mDIV)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mETIMES)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mEDIV)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLDIV)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mELDIV)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mQUOTE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mEQUOTE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCOLON)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPOW)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLTHAN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mGTHAN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLEQ)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mGEQ)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mDOT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mDOTPOW)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mAND)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mOR)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLAND)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLOR)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mNEGATE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mSEMIC)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCOMMA)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mDOTDOT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLBRACKET)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mRBRACKET)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLPAREN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mRPAREN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLCURLY)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mRCURLY)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPERCENT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mAMP)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mPCTLCURLY)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mRCURLYPCT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mASSIGN)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mLETTER)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mUNICHAR)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mDIGIT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mIDENT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mIDENT2)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mUNICHAR_IDENT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mNUMBERHACK)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mNUMBER)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mHEXVAL)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mHEXHACK)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mNEWLINE)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mWS)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mCONT)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mHACK)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mHACKB)	(struct ExprCppTreeLexer_Ctx_struct * ctx);

     void
     (*mTokens)	(struct ExprCppTreeLexer_Ctx_struct * ctx);
    const char * (*getGrammarFileName)();
    void            (*reset)  (struct ExprCppTreeLexer_Ctx_struct * ctx);
    void	    (*free)   (struct ExprCppTreeLexer_Ctx_struct * ctx);
};

// Function protoypes for the constructor functions that external translation units
// such as delegators and delegates may wish to call.
//
ANTLR3_API pExprCppTreeLexer ExprCppTreeLexerNew         (
pANTLR3_INPUT_STREAM
 instream);
ANTLR3_API pExprCppTreeLexer ExprCppTreeLexerNewSSD      (
pANTLR3_INPUT_STREAM
 instream, pANTLR3_RECOGNIZER_SHARED_STATE state);

/** Symbolic definitions of all the tokens that the 
lexer
 will work with.
 * \{
 *
 * Antlr will define EOF, but we can't use that as it it is too common in
 * in C header files and that would be confusing. There is no way to filter this out at the moment
 * so we just undef it here for now. That isn't the value we get back from C recognizers
 * anyway. We are looking for ANTLR3_TOKEN_EOF.
 */
#ifdef	EOF
#undef	EOF
#endif
#ifdef	Tokens
#undef	Tokens
#endif
#define EOF      -1
#define T__109      109
#define T__110      110
#define T__111      111
#define T__112      112
#define T__113      113
#define T__114      114
#define T__115      115
#define T__116      116
#define AMP      4
#define AND      5
#define ASSIGN      6
#define CASE      7
#define CELL_ARRAY      8
#define CELL_PARAM_LIST      9
#define CELL_VAL      10
#define CLASSDEF      11
#define CLEAR      12
#define COLON      13
#define COMMA      14
#define CONDITIONAL      15
#define CONT      16
#define CTRANSP      17
#define DIGIT      18
#define DIV      19
#define DOT      20
#define DOTDOT      21
#define DOTPOW      22
#define DUMMY      23
#define DYN_FIELD      24
#define EDIV      25
#define ELDIV      26
#define ELSE      27
#define ELSEIF      28
#define END      29
#define EQUAL      30
#define EQUOTE      31
#define ETIMES      32
#define FIELD      33
#define FOR      34
#define FUNC      35
#define FUNCTION      36
#define FUNC_DEF      37
#define FUNC_HANDLE      38
#define FUNC_LIST      39
#define GEQ      40
#define GLOBAL      41
#define GTHAN      42
#define HACK      43
#define HACKB      44
#define HEXHACK      45
#define HEXVAL      46
#define HML_STRING      47
#define IDENT      48
#define IDENT2      49
#define ID_LIST      50
#define IF      51
#define IF2      52
#define INLINE_INDEX      53
#define INLINE_INDEX_CELL      54
#define INPLACE      55
#define LAND      56
#define LBRACKET      57
#define LCURLY      58
#define LDIV      59
#define LEQ      60
#define LETTER      61
#define LOR      62
#define LPAREN      63
#define LTHAN      64
#define MATRIX      65
#define METHODS      66
#define MINUS      67
#define MR_FUNC      68
#define NEGATE      69
#define NEQUAL      70
#define NEWLINE      71
#define NUMBER      72
#define NUMBERHACK      73
#define OR      74
#define OTHER      75
#define OTHERWISE      76
#define PARAM_LIST      77
#define PCTLCURLY      78
#define PERCENT      79
#define PERSISTENT      80
#define PLUS      81
#define POW      82
#define PROPERTIES      83
#define QUOTE      84
#define QUOTE2      85
#define RBRACKET      86
#define RCURLY      87
#define RCURLYPCT      88
#define RETURN      89
#define RPAREN      90
#define SEMIC      91
#define STATEMENT_LIST      92
#define STMT      93
#define STRUCT      94
#define SWITCH      95
#define TIMES      96
#define TRANSP      97
#define TRY      98
#define UGLY1      99
#define UGLY2      100
#define UGLY3      101
#define UGLY4      102
#define UMINUS      103
#define UNICHAR      104
#define UNICHAR_IDENT      105
#define VECTOR      106
#define WHILE      107
#define WS      108
#ifdef	EOF
#undef	EOF
#define	EOF	ANTLR3_TOKEN_EOF
#endif

#ifndef TOKENSOURCE
#define TOKENSOURCE(lxr) lxr->pLexer->rec->state->tokSource
#endif

/* End of token definitions for ExprCppTreeLexer
 * =============================================================================
 */
/** } */

#ifdef __cplusplus
}
#endif

#endif

/* END - Note:Keep extra line feed to satisfy UNIX systems */
