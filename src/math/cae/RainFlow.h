/**
* @file RainFlow.h
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
#ifndef RAINFLOW_H
#define RAINFLOW_H

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#define STD std

template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

enum PrefferedBin{
	LOW,
	HIGH,
	NOTSET
};

// Signal Attribute Struct
struct SignalAttributes {

	double start;
	double end;
	double step;

	int nVals;
    int nChans;

};

struct RainflowResult {
	double mean;
	double range;
	double start_id;
	double end_id;
};

class RainFlow
{
	// virtual void DisplayBins() = 0;

protected:
	int  _id;
	char _channelformat[10];

	// Channel Characterstic
	double _scale;
	char   _polarity[10];
	double _offset;

	// Direction
	char   _dir[10];   // Not applicable for motsub and sfosub

	// Channel Signal Attributes
	struct SignalAttributes _sigAttribs;

	// Channel data
	const double* _xchanneldata ;
	const double* _ychanneldata ;

	STD::vector <double> _yorgpeaks;
	STD::vector <double> _xorgpeaks;

	// Peak data
	STD::list <double> _xpeaks;
	STD::list <double> _ypeaks;

	// Raw rainflow
	double* _rawrainrange ;
	double* _rawrainmean ;

	double m_HysteresisValue;
	int user_nbins;
	double minpoint;
	double maxpoint;

	bool m_IsRangePair;

	bool m_debug_flag;
	STD::ofstream m_OutStream_debug;
	PrefferedBin m_PrefferedBin;

	STD::vector <RainflowResult> m_vRainflowResult;
	hwMatrix* rangebins;
	hwMatrix* meanbins;
	hwMatrix* meanrangematrix;
	hwMatrix* countmatrix;

public:
	RainFlow(void);
	virtual ~RainFlow(void);

    void setXchanneldata(const double* x) {_xchanneldata = x;}
    void setYchanneldata(const double* y) {_ychanneldata = y;}
    void setsigstart(double start) {_sigAttribs.start = start;}
    void setsigNC(int nchan) {_sigAttribs.nChans = nchan;}
    void setsigNV(int nvals) {_sigAttribs.nVals = nvals;}
    void setUser_NBins(int bins) {user_nbins = bins;}
    void ComputeRangeBins();
	void ComputeMeanBins();
	const hwMatrix* getrangebins() const { return rangebins; }
	const hwMatrix* getmeanbins() const { return meanbins; }
	const hwMatrix* getmeanrangematrix()const { return meanrangematrix; }
	const hwMatrix* getcountmatrix() const { return countmatrix; }

    void setDebugFlag(STD::string strInputFilePath);
	void setUser_NBins(STD::string bins);
	void setHysteresis(int hysteresisValue);
	void setHysteresis(char* hysteresis);

	int Evaluate(const hwMatrix* _countMatrix = nullptr);

protected:
	void ExtractPeaksFromSignalData();
	void HysteresisGateCheck();
	void OrganisePeaks();
	void FindExtremePeaks();
	void ExtractRainflowCounts();
	void ComputeMeanRangeMatrix();
};

#endif // RAINFLOW_H


