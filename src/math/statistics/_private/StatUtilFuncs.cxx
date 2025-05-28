/**
* @file StatUtilFuncs.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "StatUtilFuncs.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
// Compute the variance of a strided vector
//------------------------------------------------------------------------------
double Variance(const double* x, int stride, int n, bool sampleVar)
{
    double data_zero;
    double sum = 0.0;
    double sumSq = 0.0;
    double value;
    double variance;

    if (n)
    {
        data_zero = (*x);
    }

    for (int i = 1; i < n; i++)
    {
        x += stride;
        value = (*x) - data_zero;        // shift mean to avoid overflow
        sum += value;
        sumSq += value * value;
    }

    if (sampleVar)
    {
        if (n != 1)
        {
            variance = (static_cast<double>(n) * sumSq - sum * sum) /
                static_cast<double>(n * (n - 1));
        }
        else // (n == 1)
        {
            variance = 0.0;
        }
    }
    else // population variance
    {
        variance = (static_cast<double>(n) * sumSq - sum * sum) /
            static_cast<double>(n * n);
    }

    return variance;
}
//------------------------------------------------------------------------------
// Returns the covariance of two vectors
//------------------------------------------------------------------------------
double Covariance(const double* x,
                  const double* y,
                  int           n,
                  bool          sampleCov)
{
    // vectors are stored as double* arrays
    double numer = 0.0;
    double denom = 0.0;

    if (x && y)
    {
        if (n == 1)
        {
            return 0.0;
        }
        double xx;
        double yy;
        double xSum = 0.0;
        double ySum = 0.0;
        double xySum = 0.0;

        for (int i = 0; i < n; i++)
        {
            xx = x[i];
            yy = y[i];
            xSum += xx;
            ySum += yy;
            xySum += xx * yy;
        }

        numer = n * xySum - xSum * ySum;

        denom = (sampleCov) ? n * (n - 1.0) : n * n;
    }

    return numer / denom;
}
//------------------------------------------------------------------------------
// Returns the status and the covariance of two vectors
//------------------------------------------------------------------------------
hwMathStatus Covariance(const hwMatrix& X,
                        const hwMatrix& Y,
                        double&         covar,
                        bool            sampleCov)
{
    // vectors are stored as hwMatrix objects
    if (!X.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!X.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!Y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!Y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (X.Size() != Y.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    int n = X.Size();
    const double* x = X.GetRealData();
    const double* y = Y.GetRealData();

    covar = Covariance(x, y, n, sampleCov);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the covariance of two strided vectors
//------------------------------------------------------------------------------
double Covariance(const double* pt1, int stride1,
                  const double* pt2, int stride2, int n, bool sampleVar)
{
    double data_zero1;
    double data_zero2;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sumSq = 0.0;
    double value1;
    double value2;
    double covar;

    if (n)
    {
        data_zero1 = (*pt1);
        data_zero2 = (*pt2);
    }

    for (int i = 1; i < n; i++)
    {
        pt1 += stride1;
        pt2 += stride2;
        value1 = (*pt1) - data_zero1;        // shift mean to avoid overflow
        value2 = (*pt2) - data_zero2;        // shift mean to avoid overflow
        sum1 += value1;
        sum2 += value2;
        sumSq += value1 * value2;
    }

    if (sampleVar)
    {
        if (n != 1)
        {
            covar = (static_cast<double>(n) * sumSq - sum1 * sum2) /
                    static_cast<double>(n * (n - 1));
        }
        else // (n == 1)
        {
            covar = 0.0;
        }
    }
    else // population variance
    {
        covar = (static_cast<double>(n) * sumSq - sum1 * sum2) /
                static_cast<double>(n * n);
    }

    return covar;
}
//------------------------------------------------------------------------------
// Compute the moving median of a strided vector
//------------------------------------------------------------------------------
void MovingMedianHelper(const double* data, int numPts, int stride,
                        hwMatrix& window, int nb, bool includeNaN,
                        MovingWindowEndType endtype, double userVal,
                        std::vector<int>& unsortedIdx,
                        std::vector<int>& sortedIdx, double* movmedian)
{
    // The moving window on which the median is computed is
    // sorted and circular. Old points are removed by rotating
    // in the direction that requires the smallest number of
    // shifts to make room for each new point to be inserted.
    // The numeric values are inserted/sorted on the index range
    // [lowIdx, highIdx] and NaNs are placed outside of that range.
    int winSize = window.Size();
    int lowIdx = 0;
    int highIdx = -1;
    int NaNcount = 0;

    // initialize window with insert sort
    if (endtype == SHRINK || endtype == DISCARD)
    {
        for (int i = 0; i < nb; ++i)
        {
            window(winSize - 1 - i) = std::numeric_limits<double>::quiet_NaN();
            unsortedIdx[winSize - 1 - i] = i;
        }

        NaNcount = nb;
    }
    else if (endtype == USERVAL)
    {
        if (!IsNaN_T(userVal))
        {
            for (int i = 0; i < nb; ++i)
            {
                window(i) = userVal;
                unsortedIdx[i] = i;
            }

            highIdx = nb - 1;
        }
        else
        {
            for (int i = 0; i < nb; ++i)
            {
                window(winSize - 1 - i) = userVal;
                unsortedIdx[winSize - 1 - i] = i;
            }

            NaNcount = nb;
        }
    }
    else if (endtype == SAME)
    {
        if (!IsNaN_T(data[0]))
        {
            for (int i = 0; i < nb; ++i)
            {
                window(i) = data[0];
                unsortedIdx[i] = i;
            }

            highIdx = nb - 1;
        }
        else
        {
            for (int i = 0; i < nb; ++i)
            {
                window(winSize - 1 - i) = data[0];
                unsortedIdx[winSize - 1 - i] = i;
            }

            NaNcount = nb;
        }
    }
    else if (endtype == PERIODIC)
    {
        int kk = numPts - nb;

        // sort while wrapping
        for (int i = 0; i < nb; ++i)
        {
            int k = kk + i;
            double value = data[k * stride];

            if (IsNaN_T(value))
            {
                window(winSize - 1 - NaNcount) = value;
                unsortedIdx[winSize - 1 - NaNcount] = i;
                ++NaNcount;
                continue;
            }

            ++highIdx;
            int j = highIdx;

            for (; j > 0; --j)
            {
                if (value < window(j - 1))
                {
                    window(j) = window(j - 1);
                    unsortedIdx[j] = unsortedIdx[j - 1];
                }
                else
                {
                    break;
                }
            }

            window(j) = value;
            unsortedIdx[j] = i;
        }
    }

    for (int i = nb; i < winSize; ++i)
    {
        int k = i - nb;
        double value = data[k * stride];

        if (IsNaN_T(value))
        {
            window(winSize - 1 - NaNcount) = value;
            unsortedIdx[winSize - 1 - NaNcount] = i;
            ++NaNcount;
            continue;
        }

        ++highIdx;
        int j = highIdx;

        for (; j > 0; --j)
        {
            if (value < window(j - 1))
            {
                window(j) = window(j - 1);
                unsortedIdx[j] = unsortedIdx[j - 1];
            }
            else
            {
                break;
            }
        }

        window(j) = value;
        unsortedIdx[j] = i;
    }

    for (int i = 0; i < winSize; ++i)
    {
        sortedIdx[unsortedIdx[i]] = i;
    }

    // filter
    int medWinIdx = lowIdx + (highIdx - lowIdx) / 2;
    double median;

    if ((winSize - NaNcount) % 2 == 1)
        median = window(medWinIdx);
    else
        median = 0.5 * (window(medWinIdx) + window((medWinIdx + 1) % winSize));

    if (endtype == SHRINK || endtype == DISCARD)
    {
        if (includeNaN && NaNcount > nb)
        {
            movmedian[0] = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            movmedian[0] = median;
        }
    }
    else
    {
        if (includeNaN && NaNcount)
        {
            movmedian[0] = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            movmedian[0] = median;
        }
    }

    for (int i = 1; i < numPts; ++i)
    {
        int oldestRawIdx = i - nb - 1;
        int newestRawIdx = oldestRawIdx + winSize;
        int oldestWinIdx = sortedIdx[0];
        double newestWinVal;
        double oldestWinVal;

        if (oldestRawIdx >= 0)
            oldestWinVal = data[oldestRawIdx * stride];
        else if (endtype == SHRINK || endtype == DISCARD)
            oldestWinVal = std::numeric_limits<double>::quiet_NaN();
        else if (endtype == USERVAL)
            oldestWinVal = userVal;
        else if (endtype == SAME)
            oldestWinVal = data[(numPts - 1) * stride];
        else if (endtype == PERIODIC)
            oldestWinVal = data[(numPts + oldestRawIdx) * stride];
        else    // should not happen
            oldestWinVal = 0.0;

        if (IsNaN_T(oldestWinVal))
        {
            --NaNcount;
        }

        if (newestRawIdx < numPts)
            newestWinVal = data[newestRawIdx * stride];
        else if (endtype == SHRINK || endtype == DISCARD)
            newestWinVal = std::numeric_limits<double>::quiet_NaN();
        else if (endtype == USERVAL)
            newestWinVal = userVal;
        else if (endtype == SAME)
            newestWinVal = data[(numPts - 1) * stride];
        else if (endtype == PERIODIC)
            newestWinVal = data[(newestRawIdx - numPts) * stride];
        else    // should not happen
            newestWinVal = 0.0;

        // step 1: search for location of newest value
        int numHigher = (highIdx - medWinIdx + winSize) % winSize;
        int numLower = (medWinIdx - lowIdx + winSize) % winSize;
        int newestWinIdx = -1;
        int newExtremeVal = 0;  // (low = -1, high = 1)

        if (IsNaN_T(newestWinVal))
        {
            ++NaNcount;
            newestWinIdx = (lowIdx - 1 + winSize) % winSize;
            newExtremeVal = 1;  // non-numeric, but placed as if a new high
        }
        else
        {
            if (newestWinVal >= median)
            {
                for (int j = medWinIdx; j < medWinIdx + numHigher + 1; ++j)
                {
                    int jmod = j % winSize;

                    if (window(jmod) > newestWinVal)
                    {
                        newestWinIdx = jmod;    // assume shiftCase == 1 for now
                        break;
                    }
                }

                if (newestWinIdx == -1)
                {
                    // new maximum
                    newExtremeVal = 1;
                    newestWinIdx = (highIdx == winSize - 1) ? 0 : (highIdx + 1);
                }
            }
            else // (newestWinVal < median)
            {
                for (int j = medWinIdx + winSize; j > medWinIdx - numLower - 1 + winSize; --j)
                {
                    int jmod = j % winSize;

                    if (window(jmod) < newestWinVal)
                    {
                        newestWinIdx = jmod;    // assume shiftCase == -1 for now
                        break;
                    }
                }

                if (newestWinIdx == -1)
                {
                    // new minimum
                    newExtremeVal = -1;
                    newestWinIdx = (lowIdx == 0) ? winSize - 1 : (lowIdx - 1);
                }
            }
        }

        // step 2: choose shift direction
        int shiftCase = 0;      // (left = -1, right = 1)

        if (oldestWinIdx < newestWinIdx)
        {
            if (newestWinIdx - oldestWinIdx - 1 >= winSize / 2)
                shiftCase = 1;
            else
                shiftCase = -1;
        }
        else if (oldestWinIdx > newestWinIdx)
        {
            if (oldestWinIdx - newestWinIdx - 1 < winSize / 2)
                shiftCase = 1;
            else
                shiftCase = -1;
        }

        // step 3: correct for earlier shiftCase assumption if needed
        if (newestWinVal >= median)
        {
            if (shiftCase == -1)
                newestWinIdx = (newestWinIdx == 0) ? winSize - 1 : (newestWinIdx - 1);
        }
        else
        {
            if (shiftCase == 1)
                newestWinIdx = (newestWinIdx == winSize - 1) ? 0 : (newestWinIdx + 1);
        }

        // step 4: remove oldest value, insert newest
        if (oldestWinIdx < newestWinIdx)
        {
            if (shiftCase == -1) // shift left
            {
                for (int j = oldestWinIdx; j < newestWinIdx; ++j)
                {
                    window(j) = window(j + 1);
                    unsortedIdx[j] = unsortedIdx[j + 1] - 1;
                    sortedIdx[unsortedIdx[j]] = j;
                }

                // wrap around
                for (int j = newestWinIdx + 1; j < oldestWinIdx + winSize; ++j)
                {
                    --unsortedIdx[j % winSize];
                    sortedIdx[unsortedIdx[j % winSize]] = j % winSize;
                }
            }
            else // (shiftCase == 1), shift right
            {
                // wrap around
                for (int j = oldestWinIdx + winSize; j > newestWinIdx; --j)
                {
                    window(j % winSize) = window((j - 1) % winSize);
                    unsortedIdx[j % winSize] = unsortedIdx[(j - 1) % winSize] - 1;
                    sortedIdx[unsortedIdx[j % winSize]] = j % winSize;
                }

                for (int j = oldestWinIdx + 1; j < newestWinIdx; ++j)
                {
                    --unsortedIdx[j];
                    sortedIdx[unsortedIdx[j]] = j;
                }
            }
        }
        else if (oldestWinIdx > newestWinIdx)
        {
            if (shiftCase == -1) // shift left
            {
                // wrap around
                for (int j = oldestWinIdx; j < newestWinIdx + winSize; ++j)
                {
                    window(j % winSize) = window((j + 1) % winSize);
                    unsortedIdx[j % winSize] = unsortedIdx[(j + 1) % winSize] - 1;
                    sortedIdx[unsortedIdx[j % winSize]] = j % winSize;
                }

                for (int j = newestWinIdx + 1; j < oldestWinIdx; ++j)
                {
                    --unsortedIdx[j];
                    sortedIdx[unsortedIdx[j]] = j;
                }
            }
            else // (shiftCase == 1), shift right
            {
                for (int j = oldestWinIdx; j > newestWinIdx; --j)
                {
                    window(j) = window(j - 1);
                    unsortedIdx[j] = unsortedIdx[j - 1] - 1;
                    sortedIdx[unsortedIdx[j]] = j;
                }

                // wrap around
                for (int j = oldestWinIdx + 1; j < newestWinIdx + winSize; ++j)
                {
                    --unsortedIdx[j % winSize];
                    sortedIdx[unsortedIdx[j % winSize]] = j % winSize;
                }
            }
        }
        else // (oldestWinIdx == newestWinIdx)
        {
            // shiftCase = 0;

            for (int j = 0; j < newestWinIdx; ++j)
            {
                --unsortedIdx[j];
                sortedIdx[unsortedIdx[j]] = j;
            }

            for (int j = newestWinIdx + 1; j < winSize; ++j)
            {
                --unsortedIdx[j];
                sortedIdx[unsortedIdx[j]] = j;
            }
        }

        window(newestWinIdx) = newestWinVal;
        unsortedIdx[newestWinIdx] = winSize - 1;
        sortedIdx[unsortedIdx[newestWinIdx]] = newestWinIdx;

        // step 5: update low / high numeric indices
        if (newExtremeVal == -1)            // new low
        {
            lowIdx = newestWinIdx;
        }
        else if (lowIdx == oldestWinIdx)     // low value removed
        {
            if (shiftCase == 1 || newExtremeVal == 1)
            {
                lowIdx = (lowIdx == winSize - 1) ? 0 : (lowIdx + 1);
            }
        }
        else if (shiftCase == -1)           // check if low value shifted left
        {
            if (oldestWinIdx < newestWinIdx)
            {
                if (lowIdx > oldestWinIdx && lowIdx <= newestWinIdx)
                {
                    lowIdx = (lowIdx == 0) ? winSize - 1 : (lowIdx - 1);
                }
            }
            else if (oldestWinIdx > newestWinIdx)
            {
                if (lowIdx > oldestWinIdx || lowIdx <= (newestWinIdx + winSize) % winSize)
                {
                    lowIdx = (lowIdx == 0) ? winSize - 1 : (lowIdx - 1);
                }
            }
        }
        else if (shiftCase == 1)            // check if low value shifted right
        {
            if (oldestWinIdx < newestWinIdx)
            {
                if (lowIdx >= newestWinIdx || lowIdx < (oldestWinIdx + winSize) % winSize)
                {
                    lowIdx = (lowIdx == winSize - 1) ? 0 : (lowIdx + 1);
                }
            }
            else if (oldestWinIdx > newestWinIdx)
            {
                if (lowIdx >= newestWinIdx && lowIdx < oldestWinIdx)
                {
                    lowIdx = (lowIdx == winSize - 1) ? 0 : (lowIdx + 1);
                }
            }
        }

        highIdx = (lowIdx + (winSize - NaNcount) - 1) % winSize;

        // step 6: update median
        medWinIdx = (lowIdx + (winSize - NaNcount - 1) / 2 + winSize) % winSize;

        if ((winSize - NaNcount) % 2 == 1)
            median = window(medWinIdx);
        else
            median = 0.5 * (window(medWinIdx) + window((medWinIdx + 1) % winSize));

        if (endtype == SHRINK || endtype == DISCARD)
        {
            if (i < nb)
            {
                if (includeNaN && NaNcount > (nb - i))
                {
                    movmedian[i * stride] = std::numeric_limits<double>::quiet_NaN();

                    if (NaNcount - (nb - i) == 1)
                    {
                        window(medWinIdx) = window(medWinIdx + 1);
                        median = window(medWinIdx);
                    }
                }
                else
                {
                    movmedian[i * stride] = median;
                }
            }
            else if (i > numPts - (winSize - nb))
            {
                if (includeNaN && NaNcount > (i - nb - numPts + winSize))
                {
                    movmedian[i * stride] = std::numeric_limits<double>::quiet_NaN();

                    if (NaNcount - (nb + numPts - i) == 1)
                    {
                        window(medWinIdx) = window(medWinIdx + 1);
                        median = window(medWinIdx);
                    }
                }
                else
                {
                    movmedian[i * stride] = median;
                }
            }
            else
            {
                if (includeNaN && NaNcount)
                {
                    movmedian[i * stride] = std::numeric_limits<double>::quiet_NaN();

                    if (winSize - NaNcount == 1)
                    {
                        window(medWinIdx) = window(medWinIdx + 1);
                        median = window(medWinIdx);
                    }
                }
                else
                {
                    movmedian[i * stride] = median;
                }
            }
        }
        else
        {
            if (includeNaN && NaNcount)
            {
                movmedian[i * stride] = std::numeric_limits<double>::quiet_NaN();

                if (winSize - NaNcount == 1)
                {
                    window(medWinIdx) = window(medWinIdx + 1);
                    median = window(medWinIdx);
                }
            }
            else
            {
                movmedian[i * stride] = median;
            }
        }
    }
}
//------------------------------------------------------------------------------
// Compute the squared Euclidean distance between two ND points
//------------------------------------------------------------------------------
double SqEuclidDistance(const double* pt1, int stride1,
                        const double* pt2, int stride2, int n)
{
    double diff = (*pt2) - (*pt1);
    double normSq = diff * diff;

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        diff = (*pt2) - (*pt1);
        normSq += diff * diff;
    }

    return normSq;
}
//------------------------------------------------------------------------------
// Compute the Euclidean distance between two ND points
//------------------------------------------------------------------------------
double EuclidDistance(const double* pt1, int stride1,
                      const double* pt2, int stride2, int n)
{
    return sqrt(SqEuclidDistance(pt1, stride1, pt2, stride2, n));
}
//------------------------------------------------------------------------------
// Compute the Standarized Euclidean distance between two ND points
//------------------------------------------------------------------------------
double StdEuclidDistance(const double* pt1, int stride1,
                         const double* pt2, int stride2, int n)
{
    double var1 = Variance(pt1, stride1, n, true);
    double var2 = Variance(pt2, stride2, n, true);

    return sqrt(SqEuclidDistance(pt1, stride1, pt2, stride2, n) / (var1 * var2));
}
//------------------------------------------------------------------------------
// Compute the City block distance between two ND points
//------------------------------------------------------------------------------
double CityBlockDistance(const double* pt1, int stride1,
                         const double* pt2, int stride2, int n)
{
    double diff = (*pt2) - (*pt1);
    double norm = fabs(diff);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        diff = (*pt2) - (*pt1);
        norm += fabs(diff);
    }

    return norm;
}
//------------------------------------------------------------------------------
// Compute the Cosine distance between two ND points
//------------------------------------------------------------------------------
double CosineDistance(const double* pt1, int stride1,
                      const double* pt2, int stride2, int n)
{
    double dp = (*pt2) * (*pt1);
    double s1 = (*pt1) * (*pt1);
    double s2 = (*pt2) * (*pt2);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        dp += (*pt2) * (*pt1);
        s1 += (*pt1) * (*pt1);
        s2 += (*pt2) * (*pt2);
    }

    return 1.0 - dp / (sqrt(s1) * sqrt(s2));
}
//------------------------------------------------------------------------------
// Compute the Correlation distance between two ND points
//------------------------------------------------------------------------------
double CorrelationDistance(const double* pt1, int stride1,
                           const double* pt2, int stride2, int n)
{
    double cov = Covariance(pt1, stride1, pt2, stride2, n, true);
    double var1 = Variance(pt1, stride1, n, true);
    double var2 = Variance(pt2, stride2, n, true);

    return 1.0 - cov / (sqrt(var1) * sqrt(var2));
}
//------------------------------------------------------------------------------
// Compute the Hamming distance between two ND points
//------------------------------------------------------------------------------
double HammingDistance(const double* pt1, int stride1,
                       const double* pt2, int stride2, int n)
{
    int d = (*pt2 != *pt1) ? 1 : 0;

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        d += (*pt2 != *pt1) ? 1 : 0;
    }

    return static_cast<double> (d) / n;
}
//------------------------------------------------------------------------------
// Compute the Jaccard distance between two ND points
//------------------------------------------------------------------------------
double JaccardDistance(const double* pt1, int stride1,
                       const double* pt2, int stride2, int n)
{
    int sn = 0;
    int sd = 0;

    if ((*pt2 != 0.0) | (*pt1 != 0.0))
    {
        if (*pt2 != *pt1)
            ++sn;

        ++sd;
    }

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;

        if ((*pt2 != 0.0) | (*pt1 != 0.0))
        {
            if (*pt2 != *pt1)
                ++sn;

            ++sd;
        }
    }

    return static_cast<double> (sn) / static_cast<double> (sd);
}
//------------------------------------------------------------------------------
// Compute the Chebyshev distance between two ND points
//------------------------------------------------------------------------------
double ChebyshevDistance(const double* pt1, int stride1,
                         const double* pt2, int stride2, int n)
{
    double diff = (*pt2) - (*pt1);
    double norm = fabs(diff);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        diff = fabs((*pt2) - (*pt1));

        if (diff > norm)
            norm = diff;
    }

    return norm;
}
//------------------------------------------------------------------------------
// Compute the Minkowski distance between two ND points
//------------------------------------------------------------------------------
double MinkowskiDistance(const double* pt1, int stride1,
                         const double* pt2, int stride2, int n, int p)
{

    if (p < 1)
        return 0.0;

    double diff = fabs((*pt2) - (*pt1));
    double normp = pow(diff, p);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        diff = fabs((*pt2) - (*pt1));
        normp += pow(diff, p);
    }

    return pow(normp, 1.0 / static_cast<double> (p));
}
//------------------------------------------------------------------------------
// Compute the L1 distance between two ND points
//------------------------------------------------------------------------------
double L1Distance(const double* pt1, int stride1,
                  const double* pt2, int stride2, int n)
{
    double diff = (*pt2) - (*pt1);
    double norm = fabs(diff);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        diff = (*pt2) - (*pt1);
        norm += fabs(diff);
    }

    return norm;
}
//------------------------------------------------------------------------------
// Compute the Chi-squared distance between two ND points
//------------------------------------------------------------------------------
double ChiSqDistance(const double* pt1, int stride1,
                     const double* pt2, int stride2, int n)
{
    double sum  = (*pt2) + (*pt1);
    double diff = (*pt2) - (*pt1);
    double norm = diff * diff / sum;

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        sum  = (*pt2) + (*pt1);
        diff = (*pt2) - (*pt1);
        norm += diff * diff / sum;
    }

    norm /= 2.0;

    return norm;
}
//------------------------------------------------------------------------------
// Compute the Earth Mover's distance between two ND points
//------------------------------------------------------------------------------
double EarthMoversDistance(const double* pt1, int stride1,
                           const double* pt2, int stride2, int n)
{
    double cumsum1 = (*pt1);
    double cumsum2 = (*pt2);
    double norm = fabs(cumsum1 - cumsum2);

    for (int i = 1; i < n; ++i)
    {
        pt1 += stride1;
        pt2 += stride2;
        cumsum1 += (*pt1);
        cumsum2 += (*pt2);
        norm += fabs(cumsum1 - cumsum2);
    }

    return norm;
}
//------------------------------------------------------------------------------
// Returns the status and the total sum of squares of a data vector
//------------------------------------------------------------------------------
hwMathStatus TotalSumOfSquares(const hwMatrix& data, double& sst)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    double sum = 0.0;
    double sumSq = 0.0;
    double value;
    double data_zero;

    if (n)
    {
        data_zero = data(0);
    }
    for (int i = 1; i < n; i++)
    {
        value = data(i) - data_zero;
        sum += value;
        sumSq += value * value;
    }

    sst = sumSq - sum * sum / (double) n;

    return hwMathStatus();
}
