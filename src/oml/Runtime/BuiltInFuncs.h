/**
* @file BuiltInFuncs.h
* @date October 2013
* Copyright (C) 2013-2020 Altair Engineering, Inc.  
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

#ifndef __BuiltInFuncs_h
#define __BuiltInFuncs_h

#define HW_ERROR_INVALIDOPTION(str) "Error: invalid option '" + str + "'"
#define HW_ERROR_NOSUCHFILE(file) "Error: no such file " + file

#include <map>
#include <array>
#include <regex>
#include <deque>
#include <cmath>
#include <math.h>
#include "EvaluatorInt.h"
#include "Currency.h"
#include "ErrorInfo.h"

struct BuiltinFunc;
void mapBuiltInFuncs(std::map<std::string, BuiltinFunc>* std_functions);

// built-in functions

bool oml_ferror(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mlock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_munlock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mislocked(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_flintmax(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_subsasgn(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_subsref(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fileparts(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nargoutchk(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_narginchk(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_assert(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_builtin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_atan2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
//! Returns true if successful in joining path names to build complete filename
bool oml_fullfile(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_struct2cell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cart2sph(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sph2cart(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cart2pol(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_isvector(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_deblank(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pol2cart(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nargin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nargout(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_bitand(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_bitor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_bitxor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isdir(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rmfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ispc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ismac(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isunix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_blanks(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_assignin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_evalin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_struct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_iscell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isfinite(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_and(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_or(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_uplus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_not(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_uminus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rdivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ldivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mldivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mrdivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_plus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_times(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mtimes(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_minus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_power(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mpower(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_date(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_restorepath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_gt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_geq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_leq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_eq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ne(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_linsolve(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pinv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_normalize(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mkdir(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_setenv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cond(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_sort(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isnumeric(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_all(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_any(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isscalar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_error(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_warning(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isinf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isnan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isstr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ismatrix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_issymmetric(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_issquare(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ishermitian(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_genvarname(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sub2ind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_intersect(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_setdiff(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_setxor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_double(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_regexp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cell2struct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cplxpair(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ind2sub(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isprime(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_conv2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_vertcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isletter(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isspace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_horzcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
//! Recursively displays contents of cell array
bool oml_celldisp(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_gcd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pathsep(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_func2str(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_str2func(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strncmpi(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strcmpi(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cellfun(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_datenum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_issorted(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_getfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ismember(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_angle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_which(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_linspace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isglobal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fprintf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_printf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_str2num(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fieldnames(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sprintf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_feval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rmpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_addpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_addhiddenpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_addpath2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_path(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_repmat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_deg2rad(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rad2deg(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_acosd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_asind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cosd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cross(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isstruct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fgetl(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_realmax(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_realmin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cell2mat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_run(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_run2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fflush(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_frewind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fwrite(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_inf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fread(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_seek_end(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_seek_set(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_seek_cur(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_stdin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_stdout(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_stderr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fseek(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_feof(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_islogical(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ftell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fgets(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fclose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strsplit(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ischar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_char(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strjoin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strrep(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strncmp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strcmp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lower(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_upper(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_strfind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isvarname(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_unique(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_tic(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_toc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_factor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lcm(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_polyval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_triu(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_tril(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_iscellstr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
//! Concatenates n-dimensional array objects along given dimension.
bool oml_cat(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_poly(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_permute(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_fix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rem(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_reshape(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_permute(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_squeeze(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_norm(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_i(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
//! Returns true if successul, output is a vector of nonzero elements of matrix
bool oml_find(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_eval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_length(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_ndims(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_size(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_iscomplex(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isreal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isequal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_isempty(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_getenv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_clock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_balance(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_union(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rank(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_primes(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_hypot(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rank(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_qr(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lu(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_schur(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_svd(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_eps(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_eig(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_chol(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_complex(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_transpose(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_eye(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sign(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_dot(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_kron(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pow2(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_diff(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_diag(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sparse(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_issparse(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nnz(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_speye(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_spones(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_full(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_conv(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_e(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pi(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ones(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_zeros(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_trace(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_det(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rcond(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_mod(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_imag(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_real(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_exp(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_conj(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_prod(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_sum(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cumsum(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cumprod(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_accumarray(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_log10(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_log2(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_log(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_abs(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cos(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_acos(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cosh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_acosh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sin(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_asin(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sinh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_asinh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_tan(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_atan(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_tanh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_atanh(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_sqrt(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_max(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_min(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
HML2DLL_DECLS bool oml_ceil(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_floor(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_round(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_print(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_inv(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_pwd(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_memoryuse(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_logical(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_clear(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_end(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nargout(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_nargin(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_true(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_false(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_refcnt(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_continue(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_break(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lasterr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_lastwarn(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_secret(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_meshgrid(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ndgrid(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_exit(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_cputime(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_diary(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_numel(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_circshift(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_shiftdim(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_checksyntax(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_ast(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_writepfile(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_analyze(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_getmetadata(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_help(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_helptest(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_helpmodule(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
//! Throws an error without stack information [errormsgonly]
bool oml_erromsgonly(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
//! Sets a warning without stack information [warningmsgonly]
bool oml_warningmsgonly(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_getargc(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_getargv(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_now(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);

//! List variables in the current session matching given pattern
bool oml_who(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
bool oml_logspace(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rot90(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_rehash(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
bool oml_verbose(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);

// function replacement methods
HML2DLL_DECLS void _OML_Error(EvaluatorInterface& eval, std::vector<Currency>::const_iterator start, std::vector<Currency>::const_iterator end);
bool keywordFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, Currency val);
HML2DLL_DECLS bool noConditionFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*realFunc) (double), hwComplex (*cplxFunc)(const hwComplex&));
bool conditionFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*realFunc) (double), hwComplex (*cplxFunc)(const hwComplex&), hwComplex (*cnvrtFunc)(double), bool (*conditionFunc)(double));
bool roundingFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*roundFunc) (double));
bool createCommonMatrix(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, Currency newval);
bool createCommonNDMatrix(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, Currency newval);
bool limitFunc(EvaluatorInterface& eval, const std::vector<Currency> &inputs, std::vector<Currency> &outputs, float (*floatFunc)(), double (*doubleFunc)());
bool scalarToLogicalFunc(const std::vector<Currency> &inputs, std::vector<Currency> &outputs, double (*checker)(double));
bool sortBasedOperation(EvaluatorInterface &eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, bool allowRepeatIndices,
    std::deque<std::string> (*workerCell)(std::deque<std::string>&, std::deque<std::string>&), std::deque<hwMatrix> (*workerRow)(std::deque<hwMatrix>&, std::deque<hwMatrix>&),
    std::deque<double> (*workerReal)(std::deque<double>&, std::deque<double>&), std::deque<hwComplex> (*workerCplx)(std::deque<hwComplex>&, std::deque<hwComplex>&));
bool allowComplexToLogical(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, bool (*checker)(double), bool complexBoth);
bool ineqOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool eqOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool binOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool binOperatorFuncVararg(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool logOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool unOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op);
bool anyall(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*func)(EvaluatorInterface&, const hwMatrix*), bool defaultval);
bool polCartConverter(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
    std::pair<double,double> (*realfunc)(double,double),
    std::pair<hwComplex,hwComplex> (*cplxfunc)(const hwComplex&, const hwComplex&),
    std::vector<Currency> (*helperfunc)(EvaluatorInterface& eval, const std::vector<Currency>& inputs));
bool sphCartConverter(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
    void (*realfunc)(double, double, double, double&, double&, double&),
    void (*cplxfunc)(const hwComplex&, const hwComplex&, const hwComplex&, hwComplex&, hwComplex&, hwComplex&),
    std::vector<Currency> (*helperfunc)(EvaluatorInterface& eval, const std::vector<Currency>& inputs));

// condition methods
bool nonnegative(double inpt);
bool absAtMostOne(double inpt);
bool atLeastOne(double inpt);
bool isint(const hwComplex &c);
bool isnonnegint(double d);
bool iszero(double d);
bool iszero(const hwComplex &c);

// input-reading methods
HML2DLL_DECLS bool boolFromCur(const Currency& cur);
HML2DLL_DECLS hwMatrix* matrixCopyFromInput(const Currency &input, bool allowString);
HML2DLL_DECLS hwMatrix* safeMatrixCopyFromInput(const Currency& input, bool allowString);
HML2DLL_DECLS int getFileFromInput(EvaluatorInterface &eval, const Currency &input);
FILE* makeTempFile();
std::string readTmpFile(FILE* file);
HML2DLL_DECLS std::string readOption(EvaluatorInterface& eval, const Currency &input);
std::vector<std::array<double, 6> > datesFromInput(const std::vector<Currency> &inputs, int *m, int *n);
void get1DStringsFromInput(const Currency& cur, std::vector<std::string>& vec);

// output-related methods
int getNumOutputs(const EvaluatorInterface &eval);
hwMatrix* vectorToDiag(hwMatrix *mtx);
hwMatrix* dequeToMatrix(const std::deque<hwComplex> &deq, bool row = true);
hwMatrix* dequeToMatrix(const std::deque<hwMatrix> &deq);

//! Converts a container to matrix (both real/complex)
// if nrows or ncols are negative, treated as the size of deq
// if both are -1, a column vector is returned
// otherwise assumes nrows*ncols >= container.size()
template<typename T>
inline hwMatrix* containerToMatrix(const T &container, bool row = true)
{
    int containerSize = static_cast<int>(container.size());

    // Although the matrix is created as real, if there is a complex element in
    // the values, SetElement will flip matrix type to complex
    hwMatrix* ret = (row)? 
    EvaluatorInterface::allocateMatrix(1, containerSize, hwMatrix::REAL) :
    EvaluatorInterface::allocateMatrix(containerSize, 1, hwMatrix::REAL);

    int matrixSize = ret->Size();
    // Check both the size of the container and matrix as the matrix size could
    // be smaller than the container size requested.
    for (int i = 0; i < containerSize && i < matrixSize; ++i)
        ret->SetElement(i, container[i]);  // Don't assign values directly, see note above

    return ret;
}
//! Converts a container to matrix (both real/complex)
// pads any extra space with zeros
// if nrows or ncols are negative, treated as the size of deq
// if both are -1, a column vector is returned
// otherwise assumes nrows*ncols >= container.size()
template<typename T>
inline hwMatrix* containerToMatrix(const T &container, int nrows, int ncols)
{
    hwMatrix* ret     = NULL;
    int containerSize = static_cast<int>(container.size());

    // Although the matrix is created as real, if there is a complex element in
    // the values, SetElement will flip matrix type to complex
    if (nrows < 0)
        if (ncols < 0)
            ret = EvaluatorInterface::allocateMatrix(containerSize, 1, hwMatrix::REAL);
        else
            ret = EvaluatorInterface::allocateMatrix((int)(ceil(container.size() / (double) ncols)), ncols, hwMatrix::REAL);
    else if (ncols < 0)
        ret = EvaluatorInterface::allocateMatrix(nrows, (int)(ceil(container.size() / (double) nrows)), hwMatrix::REAL);
    else
        ret = EvaluatorInterface::allocateMatrix(nrows, ncols, hwMatrix::REAL);

    int matrixSize = ret->Size();

    int i = 0;  // Outside as the loop any elements not set by container need to 
                // be initialized to 0

    // Check both the size of the container and matrix as the matrix size could
    // be smaller than the container size requested.
    for (; i < matrixSize && i < containerSize; i++)
        ret->SetElement(i, (double)container[i]); // Don't assign values directly, see note above

    while (i < matrixSize) 
        (*ret)(i++) = 0.0;

    return ret;
}

template <typename T>
inline HML_CELLARRAY* containerToCellArray(const T &container, bool row)
{
    int size = (int)container.size();
    HML_CELLARRAY *ret;
    if (row)
        ret = EvaluatorInterface::allocateCellArray(1, size);
    else
        ret = EvaluatorInterface::allocateCellArray(size, 1);

    for (int i = 0; i < size; i++)
        (*ret)(i) = container[i];

    return ret;
}

Currency addStringMask(const hwMatrix *m);
inline Currency addStringMask(double d) { return addStringMask(EvaluatorInterface::allocateMatrix(1, 1, d)); }

// math methods
double mod(double x, double y);
double signum(double x);
double doubleMod(double a, double b);
int gcd(int x, int y, int *s, int *t);
hwComplex gcd(hwComplex x, hwComplex y, hwComplex *s, hwComplex *t);
hwComplex signum(const hwComplex &cplx);
Currency curMultElems(Currency &x, Currency &y);
inline int bitAnd(int x, int y) { return x & y; }
inline int bitOr (int x, int y) { return x | y; }
inline int bitXor(int x, int y) { return x ^ y; }

// sort-related methods

// uses a binary search to find the index of val -- v must be sorted
// if reverseReturn is true, this will return -1 if val is in v and
// the index val should be inserted at if it is not in v
// if allowDuplicate is true, the index value of an instance of val
// or where it should be inserted will always be returned
template<typename T>
inline int indexOf(const std::deque<T> &d, T val, bool reverseReturn = false, bool allowDuplicate = false)
{
    int size = static_cast<int>(d.size());
    int minIndex = 0;
    int midIndex = 0;
    int maxIndex = size;
    while (minIndex <= maxIndex && minIndex < size)
    {
        midIndex = (minIndex + maxIndex) / 2;
        T curVal = d[midIndex];

        if (curVal == val)
            return (!allowDuplicate && reverseReturn ? -1 : midIndex);
        else if (curVal < val)
            minIndex = midIndex + 1;
        else
            maxIndex = midIndex - 1;
    }
    return (allowDuplicate || reverseReturn ? minIndex : -1);
}

int indexOf(const std::deque<std::string> &v, const std::string &val, bool reverseReturn = false, bool allowDuplicate = false);
int indexOf(const std::deque<hwComplex> &v, const hwComplex &val, bool reverseReturn = false, bool allowDuplicate = false);
int indexOf(const std::deque<hwMatrix> &v, hwMatrix &val, bool reverseReturn = false, bool allowDuplicate = false);
template <typename T>
inline bool hasElement(std::deque<T> &v, double val) { return (indexOf(v, val, false) != -1); }
bool complexLessThan(const hwComplex &cplx1, const hwComplex &cplx2);
inline bool complexGreaterThan(const hwComplex &cplx1, const hwComplex &cplx2) { return complexLessThan(cplx2, cplx1); }
bool rowVecLessThan(const hwMatrix *mtx1, const hwMatrix *mtx2);
inline bool rowVecLessThan(const hwMatrix &mtx1, const hwMatrix &mtx2) { return rowVecLessThan(&mtx1, &mtx2); }
inline bool rowVecGreaterThan(const hwMatrix *mtx1, const hwMatrix *mtx2) { return rowVecLessThan(mtx2, mtx1); }

// string-related methods
std::string doEscapeSequences(const std::string&);
char escapeChar(char c);
bool isstr(HML_CELLARRAY *cell);
hwMatrix* tostr(EvaluatorInterface& eval, const hwMatrix *m, bool throwError);
hwMatrix* tostr(EvaluatorInterface& eval, HML_CELLARRAY *cell, bool throwError);
HML_CELLARRAY* tocellstr(EvaluatorInterface& eval, HML_CELLARRAY *c, bool throwError);
HML2DLL_DECLS Currency toCurrencyStr(EvaluatorInterface& eval, const Currency &c, bool throwError, bool keepcell);
std::string convertToLower(std::string str);
hwMatrix* convertToLower(const hwMatrix *str);
hwMatrix* convertToUpper(const hwMatrix *str);
HML_CELLARRAY* convertToUpper(const HML_CELLARRAY *cell);
HML_CELLARRAY* convertToLower(const HML_CELLARRAY *cell);
Currency convertToUpper(const Currency &cur);
Currency convertToLower(const Currency &cur);
std::vector<Currency> convertToLower(const std::vector<Currency>::const_iterator start, const std::vector<Currency>::const_iterator end);
void cellAppendStr(EvaluatorInterface& eval, HML_CELLARRAY *outcell, const hwMatrix *newstr, bool throwError);
hwMatrix* trimright(const hwMatrix *totrim);
bool getSingularStringOrCell(Currency input, hwMatrix *&outstr, HML_CELLARRAY *&outcell, bool &usecell);
void buildString(std::vector<Currency> vec, hwMatrix *str);
std::string orderedStringVal(const Currency &cur);
std::string orderedStringVal(const hwMatrix *strmtx);

// date-related methods
double daysInMonths(int monthsPassed, int year);
int getDaysInMonth(int month);
bool isLeapYear(int year);

// individual function runners
Currency deblankHelper(EvaluatorInterface& eval, const Currency& cur);
void removeFields(EvaluatorInterface& eval, StructData* strct, const hwMatrix* fields);
bool _strncompHelper(std::string s1, std::string s2, int n);
bool dostrcmp(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int n = 0);
HML2DLL_DECLS Currency tripleCurrencyFunc(EvaluatorInterface& eval, const Currency &cur1, const Currency &cur2, const Currency &cur3, void* overlap, Currency (*func)(EvaluatorInterface&, const Currency*, const Currency*, const Currency*, void*));
Currency doStrRep(EvaluatorInterface& eval, const Currency *searchcur, const Currency *patcur, const Currency *repcur, void *dooverlap);
void dostrcat(EvaluatorInterface& eval, Currency &out, HML_CELLARRAY *cell);
Currency dostrcat(EvaluatorInterface& eval, const hwMatrix *first, const hwMatrix *second);
void dostrcat(EvaluatorInterface& eval, Currency &out, const hwMatrix *newstr);
bool dofgets(EvaluatorInterface &eval, const std::vector<Currency> &inputs, std::string &readline);
Currency doScalarToLogical(const Currency &cur, double (*checker)(double));
void cd(std::string &dir, EvaluatorInterface &eval);
std::vector<Currency> doregexp(EvaluatorInterface& eval, const std::string &search, const std::string &pattern, std::vector<Currency>::const_iterator inputIter, const std::vector<Currency>::const_iterator inputEnd);
bool isempty(const Currency &input);
HML2DLL_DECLS std::string strrep(EvaluatorInterface& eval, const std::string &search, const std::string &pattern, const std::string &replace, bool overlap);
HML2DLL_DECLS std::vector<Currency> mtxFun(EvaluatorInterface& eval, const std::vector<Currency>& inputs, int numOutputs, std::vector<Currency> (*method)(EvaluatorInterface&, const std::vector<Currency>&), const std::vector<Currency>& extras = std::vector<Currency>());
bool twoMatrixCaller(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, hwMathStatus(hwMatrix::*func)(const hwMatrix&, const hwMatrix&));

HML2DLL_DECLS bool checkMakeComplex(EvaluatorInterface& eval, hwMatrix *mtx1, hwMatrix *mtx2);
inline bool checkMakeComplex(EvaluatorInterface& eval, HML_CELLARRAY* cell1, HML_CELLARRAY* cell2) { return false; }
HML2DLL_DECLS bool checkMakeComplex(EvaluatorInterface& eval, hwMatrix *mtx1, hwMatrix *mtx2, hwMatrix *mtx3);

// general common helpers
void getMatchingSizes(const hwMatrix*, const hwMatrix*, int* m, int* n);
bool isvarname(const std::string&);
void sleep(double sec);
hwMatrix* makeRealCopy(const hwMatrix* cplxMtx);
bool isActuallyComplex(const Currency& cur);
bool isActuallyReal(const Currency& cur);
template <typename T1, typename T2>
inline bool sameSize(const T1 *m1, const T2 *m2) { return (m1->M() == m2->M() && m1->N() == m2->N()); }
HML2DLL_DECLS std::string readString(const Currency& str, int index = 0);
hwMatrix* readRow(EvaluatorInterface& eval, const hwMatrix *mtx, int index = 0);
HML2DLL_DECLS Currency readRow(EvaluatorInterface& eval, const Currency &input, int index = 0);
bool checkIsDirectory(EvaluatorInterface& eval, std::string &str, bool throwError);
bool isDirectory(std::string& str, std::string* errmsg = nullptr);
std::string getAbsolutePath(EvaluatorInterface& eval, const Currency &cur);
std::string sprintf(EvaluatorInterface& eval, std::vector<Currency>::const_iterator iter, const std::vector<Currency>::const_iterator enditer);
std::string sprintf(EvaluatorInterface& eval, const std::string &tmplt, std::vector<Currency>::const_iterator rawiter, const std::vector<Currency>::const_iterator rawenditer);
HML2DLL_DECLS void getDimensionsFromInput(const std::vector<Currency> &inputs, int *m, int *n);
bool isequal(const Currency &observed, const Currency &expected, const Currency& tol);
bool isequal(const Currency &c1, const Currency &c2);
Currency getAtIndex(const Currency &cur, int index);
Currency getAtIndex(const Currency &cur, int i, int j);
std::string concatRowToString(EvaluatorInterface& eval, int row, const HML_CELLARRAY *cell);
HML2DLL_DECLS inline Currency getTrue() { Currency True(1.0); True.SetMask(Currency::MASK_LOGICAL); return True; }
HML2DLL_DECLS inline Currency getFalse() { Currency False(0.0); False.SetMask(Currency::MASK_LOGICAL); return False; }
Currency unnest(Currency nested, const std::string &errmsg = std::string("Error: input must be a string"));
void makeInt(hwComplex &c);
Currency transpose(EvaluatorInterface& eval, const Currency &cur);
hwMatrix *removePadding(const hwMatrix *mtx);
hwMatrix *getInnerMatrix(const hwMatrix *mtx, int top, int bottom, int left, int right);
std::deque<int> getprimes(int max);
HML2DLL_DECLS bool fileExists(const std::string &file_name);
bool isEscaped(const std::string &str, size_t index);
Currency findUnsorted(EvaluatorInterface& eval, std::deque<hwMatrix> vals, const hwMatrix *searchin, bool forward);
Currency findUnsorted(std::deque<double> vals, const hwMatrix *searchin, bool forward);
Currency findUnsorted(std::deque<hwComplex> vals, const hwMatrix *searchin, bool forward);
Currency findUnsorted(EvaluatorInterface& eval, std::deque<std::string> vals, const HML_CELLARRAY *searchin, bool forward);
Currency findSorted(EvaluatorInterface& eval, std::deque<hwMatrix> vals, const hwMatrix *searchin);
Currency findSorted(std::deque<double> vals, const hwMatrix *searchin);
Currency findSorted(std::deque<hwComplex> vals, const hwMatrix *searchin);
Currency findSorted(EvaluatorInterface& eval, std::deque<std::string> vals, const HML_CELLARRAY *searchin);
int posIntFromDouble(double d);
HML2DLL_DECLS void writeCol(EvaluatorInterface& eval, hwMatrix* mtx, hwMatrix* col, int index);
HML2DLL_DECLS void writeRow(EvaluatorInterface& eval, hwMatrix* mtx, hwMatrix* row, int index);
HML2DLL_DECLS void writeCol(EvaluatorInterface& eval, HML_CELLARRAY* cell, HML_CELLARRAY* col, int index);
HML2DLL_DECLS void writeRow(EvaluatorInterface& eval, HML_CELLARRAY* cell, HML_CELLARRAY* row, int index);

// data checkers
bool checkisfinite(double d);
bool isfinite(hwComplex c);
bool isfinite(hwMatrix *m);
HML2DLL_DECLS bool isint(double d);
HML2DLL_DECLS bool islonglong(double d);
HML2DLL_DECLS bool isposint(double d);
HML2DLL_DECLS bool isint(const hwMatrix *mtx);
HML2DLL_DECLS bool isinfinity(double d);
// HML2DLL_DECLS bool isfinite(double d);
HML2DLL_DECLS double realval(const hwMatrix*, int i, int j);
HML2DLL_DECLS double realval(const hwMatrix*, int index);
HML2DLL_DECLS double realvalorscalar(const hwMatrix* mtx, int i, int j);
HML2DLL_DECLS double realvalorscalar(const hwMatrix* mtx, int index);

// path-related helpers
//!
//! Checks if given string is a directory, strips trailing slash and adds to path
//! \param eval        Evaluator interface
//! \param str         Given string
//! \param appendToEnd True if path should be added to the end of the path list
//!
void checkAddPath(EvaluatorInterface& eval, 
                  std::string&        str, 
                  bool                appendToEnd);
//!
//! Returns path as string 
//! \param eval    Evaluator interface
//! \param pathSep Path separator
//!
std::string getPathString(EvaluatorInterface& eval, 
                          char                pathSep);
std::vector<std::string> separatePathNames(EvaluatorInterface& eval, Currency cur);

// unique helpers
// factor
void addValMultiplicity(long long value, bool useTwoVecs, std::vector<long long> &primes, std::vector<int> &multiplics);
// poly
hwMatrix* checkMatrixFinite(hwMatrix *m);
// isequal
bool cellArraysEqual(const HML_CELLARRAY *c1, const HML_CELLARRAY *c2, const Currency* tol = nullptr);
bool structsEqual(const StructData* s1, const StructData* s2, const Currency* tol = nullptr);
bool funcsEqual(FunctionInfo* f1, FunctionInfo* f2);
// cell2mat
void increasetemp(int *list, int &temp, int rowindex, int max);
// num2str
std::string makeString(EvaluatorInterface& eval, double d, int precision = -1);
std::string makeString(EvaluatorInterface& eval, double d, const std::string &format);
template <typename T>
std::string makeString(EvaluatorInterface& eval, const hwComplex &c, T precForm);
// isfield
inline bool isField(EvaluatorInterface& eval, const std::map<std::string, int> &fieldNames, const std::string& field) { return fieldNames.count(field) != 0ULL; }
inline bool isField(EvaluatorInterface& eval, const std::map<std::string, int> &fieldNames, const Currency &field) { return field.IsString() && fieldNames.count(readRow(eval, field).StringVal()); }
// ismember
void stringVecFromCurrency(EvaluatorInterface& eval, const Currency &input1, const Currency &input2, std::vector<Currency> &searchfor, std::vector<Currency> &searchin, int *m, int *n);
void stringVecFromCurrencyRows(EvaluatorInterface& eval, const Currency &input1, const Currency &input2, std::vector<Currency> &searchfor, std::vector<Currency> &searchin, int *m, int *n);
void tryPushBackString(const HML_CELLARRAY *cell, std::vector<Currency> &topush);
void tryPushBackString(const hwMatrix *mtx, std::vector<Currency> &topush);
void tryPushBackStringByRows(EvaluatorInterface& eval, const hwMatrix *mtx, std::vector<Currency> &topush);
// transpose
template <typename T>
Currency _transpose(EvaluatorInterface& eval, const T *source);
// horzcat and vertcat
Currency getFirstInputOrCell(const std::vector<Currency> &inputs);
// isprime
bool isprime(int p);
// cplxpair
bool complexLessThanTol(const hwComplex &cplx1, const hwComplex &cplx2, double tol);
int indexOfTol(const std::deque<hwComplex> &d, const hwComplex &val, double tol);
std::vector<Currency> makeAComplex(EvaluatorInterface& eval, const std::vector<Currency>& inputs);
// all
double all(EvaluatorInterface& eval, const hwMatrix* mtx);
// any
double any(EvaluatorInterface& eval, const hwMatrix* mtx);
// sort
template <bool ASCEND>
hwMatrix* sort(EvaluatorInterface& eval, const hwMatrix* vec, std::pair<int&, hwMatrix*>* index_data);
template <bool ASCEND>
HML_CELLARRAY* sort(EvaluatorInterface& eval, const HML_CELLARRAY* cell, std::pair<int&, hwMatrix*>* index_data);

// cputime helper
double getCPUTime();

// cart2sph
std::vector<Currency> cart2sphHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs);
void cart2sphHelper(double x, double y, double z, double& az, double& elev, double& r);

// sph2cart
std::vector<Currency> sph2cartHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs);
void sph2cartHelper(double az, double elev, double r, double& x, double& y, double& z);

// pol2cart
std::vector<Currency> pol2cartHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs);
std::pair<double, double> pol2cartHelper(double theta, double r);

// cart2pol
std::vector<Currency> cart2polHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs);
std::pair<double, double> cart2polHelper(double theta, double r);

// removePadding
void _countZeros(const hwMatrix *mtx, int *count, int stop_count, int stop_inner, int incr, bool horiz);
void _countZerosComplex(const hwMatrix *mtx, int *count, int stop_count, int stop_inner, int incr, bool horiz);

// set-related helpers
void dounique(std::deque<hwMatrix> &vals);
void dounique(std::deque<double> &vals);
void dounique(std::deque<hwComplex> &vals);
void dounique(std::deque<std::string> &vals);

// This allows all set-related functions to be simplified & combined
std::deque<hwMatrix> doSetFunc(std::deque<hwMatrix> &a, std::deque<hwMatrix> &b, std::back_insert_iterator<std::deque<hwMatrix> > (*stdmethod)(std::deque<hwMatrix>::iterator,
    std::deque<hwMatrix>::iterator, std::deque<hwMatrix>::iterator, std::deque<hwMatrix>::iterator, std::back_insert_iterator<std::deque<hwMatrix> >, bool (*)(const hwMatrix&, const hwMatrix&)));
std::deque<hwComplex> doSetFunc(std::deque<hwComplex> &a, std::deque<hwComplex> &b, std::back_insert_iterator<std::deque<hwComplex> > (*stdmethod)(std::deque<hwComplex>::iterator,
    std::deque<hwComplex>::iterator, std::deque<hwComplex>::iterator, std::deque<hwComplex>::iterator, std::back_insert_iterator<std::deque<hwComplex> >, bool (*)(const hwComplex&, const hwComplex&)));
std::deque<double> doSetFunc(std::deque<double> &a, std::deque<double> &b, std::back_insert_iterator<std::deque<double> > (*stdmethod)(std::deque<double>::iterator,
    std::deque<double>::iterator, std::deque<double>::iterator, std::deque<double>::iterator, std::back_insert_iterator<std::deque<double> >));
std::deque<std::string> doSetFunc(std::deque<std::string> &a, std::deque<std::string> &b, std::back_insert_iterator<std::deque<std::string> > (*stdmethod)(std::deque<std::string>::iterator,
    std::deque<std::string>::iterator, std::deque<std::string>::iterator, std::deque<std::string>::iterator, std::back_insert_iterator<std::deque<std::string> >));

template <typename T>
inline std::deque<T> dounion(std::deque<T> &a, std::deque<T> &b) { return doSetFunc(a, b, std::set_union); }
template <typename T>
inline std::deque<T> dointersect(std::deque<T> &a, std::deque<T> &b) { return doSetFunc(a, b, std::set_intersection); }
template <typename T>
inline std::deque<T> dosetxor(std::deque<T> &a, std::deque<T> &b) { dounique(a); dounique(b); return doSetFunc(a, b, std::set_symmetric_difference); }
template <typename T>
inline std::deque<T> dosetdiff(std::deque<T> &a, std::deque<T> &b) { dounique(a); dounique(b); return doSetFunc(a, b, std::set_difference); }
// override these cases because we need to make a copy of a and b - they are used later in sortBasedOperation
inline std::deque<std::string> dounion(std::deque<std::string> &a, std::deque<std::string> &b) { std::deque<std::string> aa(a); std::deque<std::string> bb(b); return dounion<std::string>(aa, bb); }
inline std::deque<std::string> dointersect(std::deque<std::string> &a, std::deque<std::string> &b) { std::deque<std::string> aa(a); std::deque<std::string> bb(b); return dointersect<std::string>(aa, bb); }
inline std::deque<std::string> dosetxor(std::deque<std::string> &a, std::deque<std::string> &b) { std::deque<std::string> aa(a); std::deque<std::string> bb(b); return dosetxor<std::string>(aa, bb); }
inline std::deque<std::string> dosetdiff(std::deque<std::string> &a, std::deque<std::string> &b) { std::deque<std::string> aa(a); std::deque<std::string> bb(b); return dosetdiff<std::string>(aa, bb); }

// Helper functions
//!
//! Helper function for oml_celldisp
//! \param eval     Evaluator interface
//! \param cell     Cell
//! \param cellname Cell name
//!
void celldisp(EvaluatorInterface& eval, 
              HML_CELLARRAY*      cell, 
              const std::string&  cellname);

//!
//! Helper function for unique command
//! \param eval      Evaluator interface
//! \param x         Input currency
//! \param cmprows   True if whole rows are compared
//! \param forward   True if forward search
//! \param outputIdx True if idx(i) of each element of output(y) in given 
//!                  curr(x) is returned such that y = x(i)
//! \param inputIdx  True if idx(j) of each element of input(x) in 
//!                  output(y) such that x = y(j)
//! \param outputs   Vector of outputs
//!
void UniqueHelperFunc( EvaluatorInterface&    eval,
                       const Currency&        x,
                       bool                   cmprows,
                       bool                   forward,
                       bool                   outputIdx,
                       bool                   inputIdx,
                       std::vector<Currency>& outputs);
//!
//! Helper function for unique command for 2D matrices
//! \param eval      Evaluator interface
//! \param x         Input currency
//! \param cmprows   True if whole rows are compared
//! \param forward   True if forward search
//! \param outputIdx True if idx(i) of each element of output(y) in given 
//!                  curr(x) is returned such that y = x(i)
//! \param inputIdx  True if idx(j) of each element of input(x) in 
//!                  output(y) such that x = y(j)
//! \param outputs   Vector of outputs
//!
void UniqueHelperFuncMtx( EvaluatorInterface&    eval,
                          const Currency&        x,
                          bool                   cmprows,
                          bool                   forward,
                          bool                   outputIdx,
                          bool                   inputIdx,
                          std::vector<Currency>& outputs);
//!
//! Helper function for unique command for cell arrays
//! \param eval      Evaluator interface
//! \param x         Input currency
//! \param forward   True if forward search
//! \param outputIdx True if idx(i) of each element of output(y) in given 
//!                  curr(x) is returned such that y = x(i)
//! \param inputIdx  True if idx(j) of each element of input(x) in 
//!                  output(y) such that x = y(j)
//! \param outputs   Vector of outputs
//!
void UniqueHelperFuncCell( EvaluatorInterface&    eval,
                           const Currency&        x,
                           bool                   forward,
                           bool                   outputIdx,
                           bool                   inputIdx,
                           std::vector<Currency>& outputs);

// Sparse matrix helper functions
// These functions reside here temporarily to avoid propogating an MKL header
// dependency to all /kernel clients. They may be moved to /math/sparse in
// a future release.

// Add two sparse matrices, sum = A + B
void SparseAdd(const hwMatrixS& A, const hwMatrixS& B, hwMatrixS& sum);

// Multiply a sparse matrix by a full matrix, prod = A * B
void SparseMult(const hwMatrixS& A, const hwMatrix& B, hwMatrix& prod);

// Multiply a full matrix by a sparse matrix, prod = A * B
void SparseMult(const hwMatrix& A, const hwMatrixS& B, hwMatrix& prod);

// Multiply a sparse matrix by a sparse matrix, prod = A * B
void SparseMult(const hwMatrixS& A, const hwMatrixS& B, hwMatrixS& prod);

// Divide a full matrix by a sparse matrix on the left side, Q = A \ B
void SparseDivideLeft(const hwMatrixS& A, const hwMatrix& B, hwMatrix& Q);

// Divide a full matrix by a sparse matrix on the right side, Q = A / B
void SparseDivideRight(const hwMatrix& A, const hwMatrixS& B, hwMatrix& Q);

// Set pivot threshold for MKL sparse matrix division
bool oml_mkl_sdpt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);

#endif
