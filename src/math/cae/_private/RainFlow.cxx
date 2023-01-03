/**
* @file RainFlow.cxx
* @date February 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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
#include "RainFlow.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "hwMatrix.h"

/*
 * macro definations
 */
#define READDATA  maxActiveIndex++;\
                  if (maxActiveIndex == numPeaks) break;\
                  numActivePeaks++;\
                  if (activePeaks.size() < numActivePeaks)\
                      activePeaks.push_back(maxActiveIndex);\
                  else\
                      activePeaks[numActivePeaks - 1] = maxActiveIndex;\
                  continue;

RainFlow::RainFlow(void)
{
	_xchanneldata = nullptr;
	_ychanneldata = nullptr;
	_rawrainrange = nullptr;
	_rawrainmean = nullptr;
	minpoint = std::numeric_limits<double>::quiet_NaN();
	maxpoint = std::numeric_limits<double>::quiet_NaN();
	rangebins = nullptr;
	meanbins = nullptr;
	meanrangematrix = nullptr;
	countmatrix = nullptr;

	m_HysteresisValue = 0;
	user_nbins = 64;

	m_debug_flag = 0;
	m_PrefferedBin = NOTSET;
	m_IsRangePair = false;

	char* env_prefferedbin = getenv("HW_RAINFLOW_PREFERBIN");

	if (env_prefferedbin)
	{
		char* pEnd;
		long prefferedbin_flag = strtol(env_prefferedbin, &pEnd, 10);
		if (prefferedbin_flag == 1)
		{
			m_PrefferedBin = HIGH;
		}
		else if (prefferedbin_flag == 0)
		{
			m_PrefferedBin = LOW;
		}
	}
}

RainFlow::~RainFlow(void)
{
	if (_rawrainrange)
		delete[] _rawrainrange;

	if (_rawrainmean)
		delete[] _rawrainmean;

	if (rangebins)
		delete rangebins;

	if (meanbins)
		delete meanbins;

	if (meanrangematrix)
		delete meanrangematrix;

	if (countmatrix)
		delete countmatrix;

	if (m_debug_flag)
		m_OutStream_debug.close();
}

void RainFlow::setDebugFlag(STD::string strInputFilePath)
{
	char* env_debug = getenv("HW_RAINFLOW_DEBUG");
	long debug_flag = 0;
	if (env_debug)
	{
		char* pEnd;
		debug_flag = strtol(env_debug, &pEnd, 10);
	}

	if (debug_flag > 0)
	{
		m_debug_flag = 1;
		STD::string strdebugfilename = strInputFilePath;
		size_t endindex = strdebugfilename.rfind(".");
		if (endindex != STD::string::npos)
			strdebugfilename = strdebugfilename.substr(0, endindex);

		strdebugfilename = strdebugfilename + ".csv";
		m_OutStream_debug.open(strdebugfilename.c_str());
	}
}

void RainFlow::setUser_NBins(STD::string bins)
{
	char * pEnd;
	user_nbins = strtol (bins.c_str(),&pEnd,10);
}

void RainFlow::setHysteresis(char* hysteresis)
{
	char * pEnd;
	m_HysteresisValue = strtod (hysteresis,&pEnd);
}

void RainFlow::setHysteresis(int hysteresisValue)
{
	m_HysteresisValue = hysteresisValue;
}

void RainFlow::ExtractPeaksFromSignalData()
{
	int             j;
	//double          m1, m2;

	_xpeaks.push_back(_xchanneldata[0]);
	_ypeaks.push_back(_ychanneldata[0]);

	//No of points
	int i = _sigAttribs.nVals;

	STD::list<double>::iterator itr = _ypeaks.begin();
	for (j = 1; j < i - 1; j++) {

		if ((*itr > _ychanneldata[j] && _ychanneldata[j] < _ychanneldata[j + 1])
		 || (*itr < _ychanneldata[j] && _ychanneldata[j] > _ychanneldata[j + 1]))
		{
			_xpeaks.push_back(_xchanneldata[j]);
			_ypeaks.push_back(_ychanneldata[j]);
			itr++;
		}
	}
	_xpeaks.push_back(_xchanneldata[i - 1]);
	_ypeaks.push_back(_ychanneldata[i - 1]);

	if (m_debug_flag)
	{
		m_OutStream_debug << "PeaksX, PeaksY" << STD::endl;
		/*for (int x = 0; x < _npeaks; x++)
		{
			m_OutStream_debug<<_xpeaks[x]<<","<<_ypeaks[x]<<endl;
		}*/
		STD::list <double>::iterator itr_xpeaks = _xpeaks.begin();
		STD::list <double>::iterator itr_ypeaks = _ypeaks.begin();
		for (; itr_ypeaks != _ypeaks.end(); itr_xpeaks++, itr_ypeaks++)
		{
			m_OutStream_debug << *itr_xpeaks << "," << *itr_ypeaks << STD::endl;
		}
		m_OutStream_debug << STD::endl;
	}
}

void RainFlow::HysteresisGateCheck()
{
	if (m_HysteresisValue == 0)
		return;

	STD::list<double>::iterator itr_prev;
	STD::list<double>::iterator itr = _ypeaks.begin();
	itr++;

	//bool islast_op_erase = false;
	for (; itr != _ypeaks.end(); )
	{
		itr_prev = itr;
		itr_prev--;

		//if(islast_op_erase)
		//{
		//	STD::list<double>::iterator itr_next = itr;
		//	itr_next++;
		//	if(itr_next != _ypeaks.end())
		//	{
		//	}
		//}

		double currentRange = fabs(*itr_prev - *itr);

		if (currentRange < m_HysteresisValue)
		{
			//if(*itr_prev > 0 && *itr > 0 && *itr > *itr_prev)
			//{
			//	*itr_prev = *itr;
			//}
			//else if(*itr_prev < 0 && *itr < 0 && *itr < *itr_prev)
			//{
			//	*itr_prev = *itr;
			//}
			itr = _ypeaks.erase(itr);
		}
		else
			itr++;
	}

	if (m_debug_flag)
	{
		m_OutStream_debug << "After Hysteresis removal" << STD::endl;
		/*for (int x = 0; x < _npeaks; x++)
		{
			m_OutStream_debug<<_xpeaks[x]<<","<<_ypeaks[x]<<endl;
		}*/
		STD::list <double>::iterator itr_xpeaks = _xpeaks.begin();
		STD::list <double>::iterator itr_ypeaks = _ypeaks.begin();
		for (; itr_ypeaks != _ypeaks.end(); itr_xpeaks++, itr_ypeaks++)
		{
			m_OutStream_debug << *itr_xpeaks << "," << *itr_ypeaks << STD::endl;
		}
		m_OutStream_debug << STD::endl;
	}

	if (_ypeaks.size() < 3)
		return;

	STD::list<double>::iterator itr_next;
	itr = _ypeaks.begin();
	itr++;
	itr_next = itr;
	itr_next++;

	for (; itr_next != _ypeaks.end(); itr_next++)
	{
		itr_prev = itr;
		itr_prev--;

		if (!((*itr_prev > *itr && *itr < *itr_next) || (*itr_prev < *itr && *itr > *itr_next)))
		{
			itr = _ypeaks.erase(itr);
		}
		else
		{
			itr++;
		}
		itr_next = itr;

		if (itr == _ypeaks.end())
			break;
	}
	if (m_debug_flag)
	{
		m_OutStream_debug << "After Hysteresis removal, Peaks only" << STD::endl;
		/*for (int x = 0; x < _npeaks; x++)
		{
			m_OutStream_debug<<_xpeaks[x]<<","<<_ypeaks[x]<<endl;
		}*/
		STD::list <double>::iterator itr_xpeaks = _xpeaks.begin();
		STD::list <double>::iterator itr_ypeaks = _ypeaks.begin();
		for (; itr_ypeaks != _ypeaks.end(); itr_xpeaks++, itr_ypeaks++)
		{
			m_OutStream_debug << *itr_xpeaks << "," << *itr_ypeaks << STD::endl;
		}
		m_OutStream_debug << STD::endl;
	}
}

void RainFlow::OrganisePeaks()
{
	if (IsNaN_T(minpoint))
	{
		minpoint = *(_ypeaks.begin());
	}

	if (IsNaN_T(maxpoint))
	{
		maxpoint = *(_ypeaks.begin());
	}

	STD::list<double>::iterator itr;
	STD::list<double>::iterator itr_peak = _ypeaks.begin();
	int splicepntIndex = 0;
	int i = 0;
	for (itr = _ypeaks.begin(); itr != _ypeaks.end(); itr++)
	{
		//if( fabs(pInputVect[*itr]) > fabs(pInputVect[*itr_peak]) )
		//	itr_peak = itr;
		if (*itr > *itr_peak)
		{
			itr_peak = itr;
			splicepntIndex = i;
		}
		i++;

		if (minpoint > *itr)
			minpoint = *itr;

		if (maxpoint < *itr)
			maxpoint = *itr;
	}

	STD::list<double>::iterator itr_x = _xpeaks.begin();
	advance(itr_x, splicepntIndex);

	double max_x = *itr_x;
	double maxpeak = *itr_peak;
	
	if (m_IsRangePair == true)
	{
		_xorgpeaks.assign(_xpeaks.begin(), _xpeaks.end());
		_yorgpeaks.assign(_ypeaks.begin(), _ypeaks.end());

		_xpeaks.clear();
		_ypeaks.clear();
		return;
	}

	if (itr_peak == _ypeaks.begin())
	{
		_xpeaks.push_back(max_x);
		_ypeaks.push_back(maxpeak);
		_xorgpeaks.assign(_xpeaks.begin(), _xpeaks.end());
		_yorgpeaks.assign(_ypeaks.begin(), _ypeaks.end());

		_xpeaks.clear();
		_ypeaks.clear();
		return;
	}

	_xpeaks.splice(_xpeaks.end(), _xpeaks, _xpeaks.begin(), itr_x);
	_ypeaks.splice(_ypeaks.end(), _ypeaks, _ypeaks.begin(), itr_peak);
	_xpeaks.push_back(max_x);
	_ypeaks.push_back(maxpeak);

	int joinIndex = (int)_ypeaks.size() - 1 - splicepntIndex;

	//The first point of the orginal curve now follows the last point of the orginal curve.
	//So, check for peaks has to done for the first and last point of the orginal curve. 

	//Moving to last point of orginal curve
	int orginallastpntindex = joinIndex - 1;

	if (orginallastpntindex > 0 /*&& joinIndex < (_ypeaks.size() - 2)*/)
	{
		STD::list<double>::iterator itr_next;
		STD::list<double>::iterator itr_prev;
		itr_x = _xpeaks.begin();
		itr = _ypeaks.begin();
		advance(itr_x, orginallastpntindex);
		advance(itr, orginallastpntindex);
		for (int j = orginallastpntindex; j <= joinIndex; j++)
		{
			if (itr == _ypeaks.end())
				break;

			itr_prev = itr_next = itr;
			itr_prev--;
			itr_next++;
			if (itr_next == _ypeaks.end())
				break;

			if (!((*itr_prev > *itr && *itr < *itr_next) || (*itr_prev < *itr && *itr > *itr_next)))
			{
				itr_x = _xpeaks.erase(itr_x);
				itr = _ypeaks.erase(itr);
			}
			else
			{
				itr_x++;
				itr++;
			}
		}
	}

	_xorgpeaks.assign(_xpeaks.begin(), _xpeaks.end());
	_yorgpeaks.assign(_ypeaks.begin(), _ypeaks.end());

	_xpeaks.clear();
	_ypeaks.clear();

	//_norgpeaks = _yorgpeaks.size();

	//if(m_debug_flag)
	//{
	//	//m_OutStream_debug<<"OrganisedPeaksX, OrganisedPeaksY"<<endl;

	//	for (size_t j = 0; j < _yorgpeaks.size(); j++) {
	//		//_xorgpeaks.push_back(j);
	//		//printf("%d %lf\n",j,_yorgpeaks[j]);
	//		m_OutStream_debug<<j<<","<<_yorgpeaks[j]<<endl;
	//	}
	//	//printf("\n\n");
	//	m_OutStream_debug<<endl;
	//}
}

void RainFlow::FindExtremePeaks()
{
	if (IsNaN_T(minpoint))
	{
		minpoint = *(_ypeaks.begin());
	}

	if (IsNaN_T(maxpoint))
	{
		maxpoint = *(_ypeaks.begin());
	}

	STD::list<double>::iterator itr;

	for (itr = _ypeaks.begin(); itr != _ypeaks.end(); itr++)
	{
		if (*itr < minpoint)
			minpoint = *itr;

		if (*itr > maxpoint)
			maxpoint = *itr;
	}
}

void RainFlow::ExtractRainflowCounts()
{
	// Extraction of Rainflow Counts algorithm was provided by Narayan Rangarajan
	int              numActivePeaks, maxActiveIndex;
	STD::vector<int> activePeaks;
	double           x, y, yrange, ymean;

	int numPeaks = static_cast<int> (_yorgpeaks.size());

	if (numPeaks < 3)
		return;

	maxActiveIndex = -1;
	numActivePeaks = 0;

	while (true)
	{
		if (numActivePeaks < 3)
		{
			READDATA
		}

		x = fabs(_yorgpeaks[activePeaks[numActivePeaks - 1]] - _yorgpeaks[activePeaks[numActivePeaks - 2]]);
		y = fabs(_yorgpeaks[activePeaks[numActivePeaks - 2]] - _yorgpeaks[activePeaks[numActivePeaks - 3]]);

		if (x < y)
		{
			READDATA
		}

		yrange = y;
		ymean = (_yorgpeaks[activePeaks[numActivePeaks - 2]] + _yorgpeaks[activePeaks[numActivePeaks - 3]]) / 2.;

		if (m_debug_flag)
		{
			//m_OutStream_debug<<xrange<<","<<xmean<<endl;
		    // double slope = e[numActivePeaks - 3] / e[numActivePeaks - 2];
			if (activePeaks[numActivePeaks - 2] > activePeaks[numActivePeaks - 1])
				m_OutStream_debug << activePeaks[numActivePeaks - 2] << "," << activePeaks[numActivePeaks - 1]/*<<","<<slope*/ << STD::endl;
			else
				m_OutStream_debug << activePeaks[numActivePeaks - 1] << "," << activePeaks[numActivePeaks - 2]/*<<","<<slope*/ << STD::endl;

		}

		RainflowResult rainflowres;
		rainflowres.mean = ymean;
		rainflowres.range = yrange;
		rainflowres.start_id = _xorgpeaks[activePeaks[numActivePeaks - 3]];
		rainflowres.end_id = _xorgpeaks[activePeaks[numActivePeaks - 2]];
		m_vRainflowResult.push_back(rainflowres);

		activePeaks[numActivePeaks - 3] = activePeaks[numActivePeaks - 1];
		numActivePeaks -= 2;
	}

	int numCycles = static_cast<int> (m_vRainflowResult.size());

	countmatrix = new hwMatrix(numCycles, 5, hwMatrix::REAL);

	for (int i = 0; i < numCycles; ++i)
	{
		(*countmatrix)(i, 0) = 1.0;
		(*countmatrix)(i, 1) = m_vRainflowResult[i].range;
		(*countmatrix)(i, 2) = m_vRainflowResult[i].mean;
		(*countmatrix)(i, 3) = m_vRainflowResult[i].start_id + 1;
		(*countmatrix)(i, 4) = m_vRainflowResult[i].end_id + 1;
	}
}

void RainFlow::ComputeRangeBins()
{
	if (IsNaN_T(minpoint) || IsNaN_T(maxpoint))
	{
		return;
	}

	double increment = (maxpoint - minpoint) / (user_nbins);
	double firstbin = 0 + (increment / 2);

	rangebins = new hwMatrix(1, user_nbins, hwMatrix::REAL);

	(*rangebins)(0) = firstbin;

	for (int i = 1; i < user_nbins; i++)
	{
		(*rangebins)(i) = (*rangebins)(i - 1) + increment;
	}
}

void RainFlow::ComputeMeanBins()
{
	if (IsNaN_T(minpoint) || IsNaN_T(maxpoint))
	{
		return;
	}

	double increment = (maxpoint - minpoint) / (user_nbins);
	double firstbin = minpoint + (increment / 2);	// centered bin

	meanbins = new hwMatrix(user_nbins, 1, hwMatrix::REAL);

	(*meanbins)(0) = firstbin;

	for (int i = 1; i < user_nbins; i++)
	{
		(*meanbins)(i) = (*meanbins)(i - 1) + increment;
	}
}

void RainFlow::ComputeMeanRangeMatrix()
{
	if (IsNaN_T(minpoint) || IsNaN_T(maxpoint))
	{
		return;
	}

	meanrangematrix = new hwMatrix(user_nbins, user_nbins, hwMatrix::REAL);
	meanrangematrix->SetElements(0.0);

	double increment_range = (maxpoint - minpoint)/(user_nbins);
	double increment_mean = (maxpoint - minpoint)/(user_nbins);
	double abs_minpoint = - minpoint;
	double mean;
	double range;
	int size = countmatrix->M();

	for (int i = 0; i < size; i++)
	{
		mean = (*countmatrix)(i, 2);
		range = (*countmatrix)(i, 1);

		if (range < m_HysteresisValue)
			continue;

		int rbin = (int)floor((range / increment_range)/*+0.5*/);
		int mbin = (int)floor(((mean + abs_minpoint) / increment_mean)/*+0.5*/);

		if (rbin >= user_nbins)
			rbin = user_nbins - 1;
		else if (rbin <= 0)
			rbin = 0;
		else if (rbin == (user_nbins - 1))
		{
		}
		else
		{
			double prev = (*rangebins)(rbin - 1);
			double current = (*rangebins)(rbin);
			double next = (*rangebins)(rbin + 1);
			double mid_prevcurrent = (prev + current) / 2;
			double mid_currentnext = (current + next) / 2;
			if (range < mid_prevcurrent)
				rbin = rbin - 1;
			else if (range > mid_currentnext)
				rbin = rbin + 1;
			//if(range == mid_prevcurrent || range == mid_currentnext)
			//	cout<<"Rainflow range:"<<range<<"is at intersection of bins\n";

			if (range == mid_prevcurrent && m_PrefferedBin == LOW)
			{
				rbin = rbin - 1;
			}
			if (range == mid_currentnext && m_PrefferedBin == HIGH)
			{
				rbin = rbin + 1;
			}
		}

		if (mbin >= user_nbins)
			mbin = user_nbins - 1;
		else if (mbin <= 0)
			mbin = 0;
		else if (mbin == (user_nbins - 1))
		{
		}
		else
		{
			double prev = (*meanbins)(mbin - 1);
			double current = (*meanbins)(mbin);
			double next = (*meanbins)(mbin + 1);
			double mid_prevcurrent = (prev + current) / 2;
			double mid_currentnext = (current + next) / 2;
			if (mean < mid_prevcurrent)
				mbin = mbin - 1;
			else if (mean > mid_currentnext)
				mbin = mbin + 1;
			// if(mean == mid_prevcurrent || mean == mid_currentnext)
			//	 cout<<"Rainflow mean:"<<mean<<"is at intersection of bins\n";

			if (mean == mid_prevcurrent && m_PrefferedBin == LOW)
			{
				mbin = mbin - 1;
			}
			if (mean == mid_currentnext && m_PrefferedBin == HIGH)
			{
				mbin = mbin + 1;
			}
		}

		(*meanrangematrix)(mbin, rbin) += (*countmatrix)(i, 0);
	}
}

int RainFlow::Evaluate(const hwMatrix* _countMatrix)
{
	if (_countMatrix)
	{
		// _countMatrix is computed externally
		countmatrix = const_cast<hwMatrix*>(_countMatrix);
		ExtractPeaksFromSignalData();
		FindExtremePeaks();
	}
	else
	{
		// countmatrix is computed internally
		ExtractPeaksFromSignalData();

		if (_ypeaks.size() <= 1)
			return 1;

		// HysteresisGateCheck();
		OrganisePeaks();
		ExtractRainflowCounts();
	}
	
	ComputeRangeBins();
	ComputeMeanBins();
	ComputeMeanRangeMatrix();

	if (_countMatrix)
		countmatrix = nullptr;

	return 0;
}
