/**
* @file GeometryFuncs.cxx
* @date December, 2018
* Copyright (C) 2015-2019 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair’s dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair’s trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/
#include "GeometryFuncs.h"
#include "hwMatrix.h"
#include "hwMathStatus.h"

extern "C"
{
    #include "libqhull_r/qhull_ra.h"
}

//------------------------------------------------------------------------------
// Computes the 2D convex hull
//------------------------------------------------------------------------------
hwMathStatus ConvexHull(const hwMatrix&    x,
                        const hwMatrix&    y,
                        const std::string& options,
                        hwMatrixI&         hull,
                        double&            area,
                        FILE*              errfile)
{
    // check inputs
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!y.IsEmptyOrVector())
    {
        throw hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

	int numPts = x.Size();

	if (y.Size() != numPts)
	{
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
	}

    if (numPts < 3)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL_PNTS, 1, 2);
    }

    // call Qhull
    qhT qh_qh;                // Qhull's data structure.  First argument of most calls
    qhT* qh = &qh_qh;
    QHULL_LIB_CHECK
    qh_zero(qh, errfile);
    boolT ismalloc = false;   // True if qhull should free points in qh_freeqhull() or reallocation
    char flags[250];          // option flags for qhull, see qh-quick.htm
    int exitcode;             // 0 if no error from qhull
    int k = 0;
    int dim = 2;
    hwMatrix points(numPts, 2, hwMatrix::REAL);

    for (int i = 0; i < numPts; ++i)
    {
        points(k++) = x(i);
        points(k++) = y(i);
    }

    sprintf(flags, "qhull %s", options.c_str());

    exitcode = qh_new_qhull(qh, dim, numPts, points.GetRealData(), ismalloc,
                            flags, NULL, errfile);

    if (exitcode)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL);
    }

    facetT*   facet;            // set by FORALLfacets
    hwMatrixI index(qh->num_facets, 2, hwMatrixI::REAL);
    k = 0;

    // 'qh->facet_list' contains the convex hull
    FORALLfacets
    {
        vertexT*  vertex;
        vertexT** vertexp;

        FOREACHvertex_(facet->vertices)
        {
            index(k++) = 1 + qh_pointid(qh, vertex->point);     // one based
        }
    }

    // get one vertex index per facet
    hwMathStatus status = hull.Dimension(qh->num_vertices + 1, 1, hwMatrixI::REAL);

    if (!status.IsOk())
    {
        return status;
    }

    hull(0)    = index(0);
    int nextid = index(1);

    for (int i = 1; i < qh->num_facets; ++i)
    {
        for (int j = 2; j < 2*qh->num_vertices; ++j)
        {
            if (index(j) == nextid)
            {
                hull(i) = nextid;

                if (j % 2 == 0)
                {
                    nextid = index(j + 1);
                    index(j + 1) = 0;
                }
                else
                {
                    nextid = index(j - 1);
                    index(j - 1) = 0;
                }

                break;
            }
        }
    }

    hull(qh->num_vertices) = hull(0);

    // positive area upon input is a request for the calculation
    if (area > 0.0)
    {
        FORALLfacets
        {
            double dist;

            if (!facet->normal)
                continue;

            if (facet->upperdelaunay && qh->ATinfinity)
                continue;

            facet->f.area = area = qh_facetarea(qh, facet);
            facet->isarea = True;
            qh->totarea += area;
            qh_distplane(qh, qh->interior_point, facet, &dist);
            qh->totvol += -dist * area / qh->hull_dim;
        }

        area = qh->totvol;
    }

    // clean up Qhull memory
    qh_freeqhull(qh, !qh_ALL);

    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);

    if (curlong || totlong)
    {
        status(HW_MATH_WARN_EOL);
    }
    
    return status;
}
//------------------------------------------------------------------------------
// Computes the ND convex hull
//------------------------------------------------------------------------------
hwMathStatus ConvexHulln(const hwMatrix&    P,
                         const std::string& options,
                         hwMatrixI&         hull,
                         double&            volume,
                         FILE*              errfile)
{
    // check inputs
    if (!P.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int numDim = P.M();
    int numPts = P.N();

    if (numDim < 2)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL_DIMS, 1);
    }

    if (numPts < 3)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL_PNTS, 1);
    }

    // call Qhull
    qhT qh_qh;                // Qhull's data structure.  First argument of most calls
    qhT* qh = &qh_qh;
    QHULL_LIB_CHECK
    qh_zero(qh, errfile);
    boolT ismalloc = false;   // True if qhull should free points in qh_freeqhull() or reallocation
    char flags[250];          // option flags for qhull, see qh-quick.htm
    int exitcode;             // 0 if no error from qhull
    int k = 0;

    sprintf(flags, "qhull %s", options.c_str());

    exitcode = qh_new_qhull(qh, numDim, numPts, const_cast<double*> (P.GetRealData()), ismalloc,
                            flags, NULL, errfile);

    if (exitcode)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL);
    }

    hwMathStatus status = hull.Dimension(qh->num_facets, numDim, hwMatrixI::REAL);

    if (!status.IsOk())
    {
        return status;
    }

    facetT* facet;            // set by FORALLfacets
    int     i = 0;

    // 'qh->facet_list' contains the convex hull
    FORALLfacets
    {
        vertexT*  vertex;
        vertexT** vertexp;
        k = 0;

        if (numDim == 3)
        {
            setT* vertices = qh_facet3vertex(qh, facet);

            FOREACHvertex_(vertices)
            {
                hull(i, k++) = 1 + qh_pointid(qh, vertex->point);     // one based
            }

            qh_settempfree(qh, &vertices);
        }
        else
        {
            if (facet->toporient ^ qh_ORIENTclock)
            {
                FOREACHvertex_(facet->vertices)
                    hull(i, k++) = 1 + qh_pointid(qh, vertex->point);
            }
            else
            {
                FOREACHvertexreverse12_(facet->vertices)
                    hull(i, k++) = 1 + qh_pointid(qh, vertex->point);
             }
        }

        ++i;
    }

    // positive volume upon input is a request for the calculation
    if (volume > 0.0)
    {
        double area;
        double dist;

        FORALLfacets
        {
            if (!facet->normal)
                continue;

            if (facet->upperdelaunay && qh->ATinfinity)
                 continue;

            facet->f.area = area = qh_facetarea(qh, facet);
            facet->isarea = True;
            qh->totarea  += area;
            qh_distplane(qh, qh->interior_point, facet, &dist);
            qh->totvol += -dist * area / qh->hull_dim;
        }

        volume = qh->totvol;
    }

    // clean up Qhull memory
    qh_freeqhull(qh, !qh_ALL);

    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);

    if (curlong || totlong)
    {
        status(HW_MATH_WARN_EOL);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the ND Delaunay triangulation
//------------------------------------------------------------------------------
hwMathStatus Delaunayn(const hwMatrix&    P,
                       const std::string& options,
                       hwMatrixI&         triang,
                       FILE*              errfile)
{
    // check inputs
    if (!P.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int numDim = P.M();
    int numPts = P.N();

    if (numDim < 2)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL_DIMS, 1);
    }

    if (numPts < 3)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL_PNTS, 1);
    }

    // call Qhull
    qhT qh_qh;                // Qhull's data structure.  First argument of most calls
    qhT* qh = &qh_qh;
    QHULL_LIB_CHECK
    qh_zero(qh, errfile);
    boolT ismalloc = false;   // True if qhull should free points in qh_freeqhull() or reallocation
    char flags[250];          // option flags for qhull, see qh-quick.htm
    int exitcode;             // 0 if no error from qhull
    int k = 0;

    sprintf(flags, "qhull d %s", options.c_str());

    exitcode = qh_new_qhull(qh, numDim, numPts, const_cast<double*> (P.GetRealData()), ismalloc,
                            flags, NULL, errfile);

    if (exitcode)
    {
        return hwMathStatus(HW_MATH_ERR_QHULL);
    }

    qh_triangulate(qh);

    facetT*   facet;            // set by FORALLfacets
    vertexT*  vertex;
    vertexT** vertexp;

    int nf = 0;
    int i = 0;

    FORALLfacets
    {
        if (!facet->upperdelaunay)
        nf++;

        // Double check.  Non-simplicial facets will cause segfault below
        if (!facet->simplicial)
        {
            return hwMathStatus(HW_MATH_ERR_QHULL_NS_FACET);
        }
    }

    hwMathStatus status = triang.Dimension(nf, numDim + 1, hwMatrixI::REAL);

    if (!exitcode)
    {
        FORALLfacets
        {
            if (!facet->upperdelaunay)
            {
                int j = 0;

                FOREACHvertex_(facet->vertices)
                {
                    triang(i, j++) = 1 + qh_pointid(qh, vertex->point);
                }
                i++;
            }
        }
    }

    // clean up Qhull memory
    qh_freeqhull(qh, !qh_ALL);

    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);

    if (curlong || totlong)
    {
        status(HW_MATH_WARN_EOL);
    }

    return status;
}
