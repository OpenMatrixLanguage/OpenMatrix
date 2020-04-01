/**
* @file ANTLRData.cpp
* @date April 2015
* Copyright (C) 2015-2019 Altair Engineering, Inc.  
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

#include "ANTLRData.h"
#include <iostream>

ANTLRData::ANTLRData(pANTLR3_INPUT_STREAM in_input, bool delay_parser)
{
	input  = in_input; // we own this now
	lex    = ExprCppTreeLexerNew(input);
	tokens = antlr3CommonTokenStreamSourceNew(8192+1, TOKENSOURCE(lex));

	if (!delay_parser)
		parser = ExprCppTreeParserNew(tokens);
	else
		parser = NULL;
}

void ANTLRData::CreateParser(pANTLR3_COMMON_TOKEN_STREAM in_tokens)
{
	tokens = in_tokens;
	parser = ExprCppTreeParserNew(tokens);
}

ANTLRData::~ANTLRData()
{
	if (parser) // we won't always have a parser
		parser->free(parser);
	tokens->free(tokens);
	lex->free(lex);
	input->close(input);
}

pANTLR3_INPUT_STREAM ANTLRData::InputFromFilename(const std::string& filename)
{
	return antlr3FileStreamNew((pANTLR3_UINT8)filename.c_str(), ANTLR3_ENC_UTF8);
}

pANTLR3_INPUT_STREAM ANTLRData::InputFromExpression(const std::string& instring, const std::string& use_filename)
{
	const unsigned char* temp_str = reinterpret_cast<const unsigned char *> (instring.c_str());
	return antlr3StringStreamNew((pANTLR3_UINT8)temp_str, ANTLR3_ENC_UTF8, (ANTLR3_UINT32)instring.size(), (pANTLR3_UINT8)use_filename.c_str());
}

void ANTLRData::DumpTokenStream(bool include_hidden)
{
	pANTLR3_VECTOR vec = tokens->getTokens(tokens);
	int num_tokens = vec->size(vec);

	for (int j=0; j<num_tokens; ++j)
	{
		pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);
		pANTLR3_STRING str = tok->getText(tok);

		if ((tok->getChannel(tok) == HIDDEN) && !include_hidden)
			continue;

		if (tok->getType(tok) == NEWLINE)
			std::cout << "[[NEWLINE]]" << std::endl;
		else if (tok->getType(tok) == WS)
			std::cout << "[[" << str->chars << "]]" << std::endl;
		else
			std::cout << str->chars << "    " << tok->getType(tok) << std::endl;
	}
}

//! Utility to get string token
const char* GetStringToken(pANTLR3_VECTOR vec, int index)
{
	pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, index);
	pANTLR3_STRING       str = tok->getText(tok);
	return (const char*)str->chars;
}

//! Preprocesses token stream
void ANTLRData::PreprocessTokenStream(pANTLR3_COMMON_TOKEN_STREAM& tokens)
{
	pANTLR3_VECTOR vec        = tokens->getTokens(tokens);
	int            num_tokens = vec->size(vec);

	bool in_string     = false;
	bool in_comment    = false;
	bool in_ml_comment = false;
	
	// One pass for vector spacing
	for (int j=0; j<num_tokens; ++j)
	{
		pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);
		int token_type = tok->getType(tok);

		if (tok->getType(tok) == QUOTE)
		{
			if (in_string)
			{
				in_string = false;
			}
			else
			{
				if (j > 0)
				{
					pANTLR3_COMMON_TOKEN prev_tok      = (pANTLR3_COMMON_TOKEN)vec->get(vec, j-1);
					int                  prev_tok_type = prev_tok->getType(prev_tok);          

					if ((prev_tok_type != RBRACKET) && (prev_tok_type != RPAREN) && (prev_tok_type != IDENT) && (prev_tok_type != DOT))
						in_string = true;
				}
				else
				{
					in_string = true;
				}
			}
		}

		if (token_type == PERCENT && !in_string)
			in_comment = true;

		if (token_type == CONT && !in_string)
		{
			for (int k = j; k < num_tokens; ++k)
			{
				pANTLR3_COMMON_TOKEN next_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);
				next_tok->setChannel(next_tok, HIDDEN);

				if (next_tok->getType(next_tok) == NEWLINE)
				{
					j = k;
					break;
				}
			}
		}

		if (token_type == PCTLCURLY && !in_string)
		{
			for (int k=j; k<num_tokens; ++k)
			{
				pANTLR3_COMMON_TOKEN next_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);

				next_tok->setChannel(next_tok, HIDDEN);

				// If the first token is on a hidden channel, ANTLR doesn't know
				// what to do and generates a syntax error.  So don't let that happen.
				if (k == 0)
					next_tok->setType(next_tok, NEWLINE);

				if (next_tok->getType(next_tok) == RCURLYPCT)
				{
					j=k;
					break;
				}
			}

			continue;
		}

		if (token_type == NEWLINE)
		{
			in_string = false; // strings cannot span lines
			in_comment = false;
		}

		if (token_type == CONT)
		{
			if (in_comment && !in_string)
			{
				tok->setChannel(tok, 0);
				tok->setType(tok, NEWLINE);
			}
		}

		int exit_type = RBRACKET;

		if (((token_type == LBRACKET) || (token_type == LCURLY)) && !in_string && !in_comment)
		{
			if (token_type == LCURLY)
				exit_type = RCURLY;

			int brace_token = j;
			int start_count = 1;

			for (int k=j+1; k<num_tokens; k++)
			{
				pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);
				int token_type = tok->getType(tok);

				pANTLR3_COMMON_TOKEN prev_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-1);

				// There's a string in the cell or matrix, so skip that whole thing
				// but make sure we don't have a transpose operator
				bool prev_tok_check = true;
				
				if (prev_tok->getType(prev_tok) == RPAREN)
					prev_tok_check = false;

				if (prev_tok->getType(prev_tok) == IDENT)
					prev_tok_check = false;

				if ((token_type == QUOTE) && prev_tok_check)
				{
					int  q;
					bool true_string = true;

					for (q=k+1; q<num_tokens; q++)
					{
						pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);

						if (tok->getType(tok) == NEWLINE)
						{
							true_string = false;
							break;
						}

						if (tok->getType(tok) == QUOTE)
							break;
					}

					if (true_string)
					{
						for (q=k+1; q<num_tokens; q++)
						{
							pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);

							if (tok->getType(tok) == QUOTE)
								break;
						}

						k=q;
					}
				}

				if (token_type == LPAREN)
				{
					int q;
					int paren_count = 1;

					for (q=k+1; q<num_tokens; q++)
					{
						pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);
						int token_type = tok->getType(tok);

						if (token_type == LPAREN)
							paren_count++;
						else if (token_type == RPAREN)
							paren_count--;

						if (paren_count == 0)
							break;
					}

					k=q;
				}
				else if ((token_type == LBRACKET) && (exit_type == RBRACKET))
				{
					start_count++;
				}

				if ((token_type == WS)  && (k != brace_token+1))
				{
					pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k+1);

					if (!tok2)
						continue;
					
					int tok2_type = tok2->getType(tok2);

					pANTLR3_COMMON_TOKEN tok_prev = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-1);

					if (tok_prev->getType(tok_prev) == CONT)
					{
						tok_prev = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-2);

						if (tok_prev->getType(tok_prev) == WS)
							tok_prev = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-3);
					}
					
					int prev_type = tok_prev->getType(tok_prev);
					int next_type = tok2_type;

					bool prev_token_flag = false;
					bool next_token_flag = false;

					if ((prev_type != COMMA) && (prev_type != COLON) && (prev_type != SEMIC) && (prev_type != NEWLINE) && (prev_type != TIMES) && (prev_type != DIV) && (prev_type != ASSIGN) && (prev_type != EQUAL))
						prev_token_flag = true;

					if ((next_type == MINUS) || (next_type == PLUS)  || (next_type == QUOTE) || (next_type == EQUOTE) || (next_type == LPAREN) || (next_type == LCURLY) || (next_type == CONT))
						next_token_flag = true;

					if ((next_type == LPAREN) && ((prev_type == PLUS) || (prev_type == MINUS)))
						next_token_flag = false;

					if (((next_type == PLUS) || (next_type == MINUS)) && ((prev_type == PLUS) || (prev_type == MINUS)))
						next_token_flag = false;

					if ((next_type == LPAREN) && (prev_type == RPAREN))
					{
						// need to find the character before the previous (
						// if it's an @, then this is an anonymous function definition
						// and we don't want the comma

						for (int m = k-3; m > 0; m--)
						{
							pANTLR3_COMMON_TOKEN test_tok  = (pANTLR3_COMMON_TOKEN)vec->get(vec, m);
							int test_tok_type = test_tok->getType(test_tok);

							if (test_tok_type == LPAREN)
							{
								pANTLR3_COMMON_TOKEN test_tok2  = (pANTLR3_COMMON_TOKEN)vec->get(vec, m-1);
								int test_tok2_type = test_tok2->getType(test_tok2);

								if (test_tok2_type == AMP)
									next_token_flag = false;

								break;
							}
						}
					}

					if (prev_token_flag && next_token_flag)
					{
                        if (k + 2 < (int)vec->count)
                        {
                            pANTLR3_COMMON_TOKEN tok3 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k + 2);
                            int tok3_type = tok3->getType(tok3);

                            // the special case is for strings that start with a space
                            if ((tok3_type != WS) || ((tok3_type == WS) && (tok2_type == QUOTE)))
                            {
                                tok->setChannel(tok, 0);
                                tok->setType(tok, COMMA);
                            }
                        }
					}
				}
				else if (tok->getType(tok) == CONT)
				{
					tok->setChannel(tok, HIDDEN);

					pANTLR3_COMMON_TOKEN tok2 = NULL;

					for (int q = k+1; q<num_tokens; q++)
					{
						tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);
						tok2->setChannel(tok2, HIDDEN);

						if (tok2->getType(tok2) == NEWLINE)
							break;
					}

					pANTLR3_COMMON_TOKEN prev_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-1);

					bool skip = false;

					// check if previous token is already a comma
					if (prev_tok && prev_tok->getType(prev_tok) == COMMA)
					{
						skip = true;
					}
					else if (prev_tok && prev_tok->getType(prev_tok) == WS)
					{
						pANTLR3_COMMON_TOKEN prev_tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-2);

						if (prev_tok2 && prev_tok2->getType(prev_tok2) == COMMA)
							skip = true;
					}

					if (tok2 && (tok2->getType(tok2) == QUOTE) && !skip)
					{
						tok->setChannel(tok, 0);
						tok->setType(tok, COMMA);
					}
				}
				else if (token_type == exit_type)
				{
					--start_count;

					if (start_count == 0)
					{
						j = k;
						break;
					}
				}
			}
		}
	}

	// second pass for end
	in_string = false;

	for (int j=0; j<num_tokens; ++j)
	{
		pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);

		if (tok->getChannel(tok) == HIDDEN)
			continue;

		int token_type = tok->getType(tok);

		if (token_type == PERCENT)
		{
			// spin until newline
			for (int k=j+1; k<num_tokens; ++k)
			{
				pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);

				if (tok2->getType(tok2) == NEWLINE)
				{
					j = k;
					break;
				}
			}
		}

		if (token_type == QUOTE)
		{
			int q;

			pANTLR3_COMMON_TOKEN prev_tok = NULL;
					
			if (j>0)
				prev_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j-1);

			bool real_quote = true;

			if (prev_tok)
			{
				if (prev_tok->getType(prev_tok) == DOT)
					real_quote = false;
				else if (prev_tok->getType(prev_tok) == RBRACKET)
					real_quote = false;
				else if (prev_tok->getType(prev_tok) == RPAREN)
					real_quote = false;
				else if (prev_tok->getType(prev_tok) == RCURLY)
					real_quote = false;
				else if (prev_tok->getType(prev_tok) == IDENT)
					real_quote = false;
				else if (prev_tok->getType(prev_tok) == NUMBER)
					real_quote = false;
			}

			if (real_quote)
			{
				for (q=j+1; q<num_tokens; ++q)
				{
					pANTLR3_COMMON_TOKEN tok3 = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);
					pANTLR3_STRING str = tok3->getText(tok3);
					char* c_str3 = (char*)str->chars;

					if (tok3->getType(tok3) == QUOTE)
						break;
				}

				j=q;
			}
		}

		if (token_type == NUMBERHACK)
		{
			pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, j+1);

			if (tok2->getType(tok2) == TIMES)
			{
				tok->setType(tok, NUMBER);
				tok2->setType(tok2, ETIMES);
			}
			else if (tok2->getType(tok2) == DIV)
			{
				tok->setType(tok, NUMBER);
				tok2->setType(tok2, EDIV);
			}
			else if (tok2->getType(tok2) == POW)
			{
				tok->setType(tok, NUMBER);
				tok2->setType(tok2, DOTPOW);
			}
			else if (tok2->getType(tok2) == LDIV)
			{
				tok->setType(tok, NUMBER);
				tok2->setType(tok2, ELDIV);
			}
			else if (tok2->getType(tok2) == QUOTE)
			{
				tok->setType(tok, NUMBER);
				tok2->setType(tok2, QUOTE);
			}
		}
		
		if ((token_type == LPAREN) || (token_type == LCURLY))
		{
			int paren_count = 0;
			int brace_count = 0;

			if (token_type == LPAREN)
				paren_count = 1;
			else
				brace_count = 1;

			for (int k=j+1; k<num_tokens; ++k)
			{
				pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);
				int tok2_type = tok2->getType(tok2);

				if (tok2_type == LPAREN)
				{
					paren_count++;
				}
				else if (tok2_type == RPAREN)
				{
					paren_count--;

					if ((paren_count == 0) && (brace_count == 0))
					{
						j = k;
						break;
					}
				}
				else if (tok2_type == LCURLY)
				{
					brace_count++;
				}
				else if (tok2_type == RCURLY)
				{
					brace_count--;

					if ((paren_count == 0) && (brace_count == 0))
					{
						j = k;
						break;
					}
				}
				else if (tok2_type == NEWLINE)
				{
					if (tok2->getChannel(tok2) != HIDDEN)
						break;
				}
				else if (tok2_type == QUOTE)
				{
					int q;

					pANTLR3_COMMON_TOKEN prev_tok = NULL;
					
					if (k>0)
						prev_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, k-1);

					bool real_quote = true;

					if (prev_tok)
					{
						if (prev_tok->getType(prev_tok) == DOT)
							real_quote = false;
						else if (prev_tok->getType(prev_tok) == RBRACKET)
							real_quote = false;
						else if (prev_tok->getType(prev_tok) == RPAREN)
							real_quote = false;
						else if (prev_tok->getType(prev_tok) == IDENT)
							real_quote = false;
						else if (prev_tok->getType(prev_tok) == NUMBER)
							real_quote = false;
					}

					if (real_quote)
					{
						for (q=k+1; q<num_tokens; ++q)
						{
							pANTLR3_COMMON_TOKEN tok3 = (pANTLR3_COMMON_TOKEN)vec->get(vec, q);
							pANTLR3_STRING str = tok3->getText(tok3);
							char* c_str3 = (char*)str->chars;

							if (tok3->getType(tok3) == QUOTE)
								break;
						}

						k=q;
					}
				}
				
				if (tok2_type == END)
				{
					tok2->setType(tok2, IDENT);
					pANTLR3_STRING str = tok2->getText(tok2);
					str->set(str, "end");
					tok2->setText(tok2, str);
				}
			}
		}
		else if (token_type == IDENT)
		{
			if (j < (num_tokens-2))
			{
				pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, j+1);
				pANTLR3_COMMON_TOKEN tok3 = (pANTLR3_COMMON_TOKEN)vec->get(vec, j+2);

				// What I'd like to do here is replace tok3 with a LPAREN
				// and then insert an RPAREN after tokn.  However, there is no insert
				// method on the ANTLR3_VECTOR, so the best I can do for now is to remove
				// the quotes. This does have a functional impact on a case like:
				// a 'b c' unfortunately.
				if ((tok2->getType(tok2) == WS) && (tok3->getType(tok3) == QUOTE))
				{
					tok3->setChannel(tok3, HIDDEN);

					for (int k=j+3; k<num_tokens; ++k)
					{
						pANTLR3_COMMON_TOKEN tokn = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);		

						if (tokn->getType(tokn) == QUOTE)
						{
							tokn->setChannel(tokn, HIDDEN);
							j = k+1;
							break;
						}
					}
				}
			}
		}
	}

	// do this in separate loops because this one depends on the end process being completed
	for (int j=0; j<num_tokens; ++j)
	{
		pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);

		if (tok->getChannel(tok) == HIDDEN)
			continue;

		if (tok->getType(tok) == QUOTE)
		{
			in_string = !in_string;
		}

		if ((tok->getType(tok) == EQUOTE))
		{
			for (int k=j+1; k<num_tokens-1; k++)
			{
				pANTLR3_COMMON_TOKEN tok2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);

				if (tok2->getType(tok2) == QUOTE)
				{
					if (!in_string)
					{
						tok->setType(tok, QUOTE);
						pANTLR3_STRING str = tok->getText(tok);
						str->set(str, "'");
						tok->setText(tok, str);

						tok2->setType(tok2, EQUOTE);
						str = tok2->getText(tok2);
						str->set(str, "''");
						tok2->setText(tok2, str);

						in_string = true;
					}
				}
				else if (tok2->getType(tok2) != EQUOTE)
				{
					break;
				}
			}
		}

		if (tok->getType(tok) == NEWLINE)
			in_string = false;
	}
}