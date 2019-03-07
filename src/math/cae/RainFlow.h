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

enum PrefferedBin{
	LOW,
	HIGH,
	NOTSET
};

//Signal Attribute Struct
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
};

class RainFlow
{
	// virtual void DisplayBins() = 0;

protected:
	int  _id;
	char _channelformat[10];

	//Channel Characterstic
	double _scale;
	char   _polarity[10];
	double _offset;

	//Direction
	char   _dir[10];   // Not applicable for motsub and sfosub


	//Channel Signal Attributes
	struct SignalAttributes _sigAttribs;

	//Channel data
	double* _xchanneldata ;
	double* _ychanneldata ;

	STD::vector <double> _yorgpeaks;

	//Peak data
	STD::list <double> _xpeaks;
	STD::list <double> _ypeaks;

	//Raw rainflow
	double* _rawrainrange ;
	double* _rawrainmean ;

	double m_HysteresisValue;
	int user_nbins;
	double *minpoint;
	double *maxpoint;

	bool m_IsRangePair;

	bool m_debug_flag;
	STD::ofstream m_OutStream_debug;
	PrefferedBin m_PrefferedBin;

	STD::vector <double> meanbins;
	STD::vector <double> rangebins;
	STD::vector <RainflowResult> m_vRainflowResult;
	STD::vector< STD::vector<double> > meanrangematrix;

public:
	RainFlow(void);
	virtual ~RainFlow(void);

    void setXchanneldata(double* x) {_xchanneldata = x;}
    void setYchanneldata(double* y) {_ychanneldata = y;}
    void setminpoint(double* min) {minpoint = min;}
    void setmaxpoint(double* max) {maxpoint = max;}
    void setsigstart(double start) {_sigAttribs.start = start;}
    void setsigNC(int nchan) {_sigAttribs.nChans = nchan;}
    void setsigNV(int nvals) {_sigAttribs.nVals = nvals;}
    void setUser_NBins(int bins) {user_nbins = bins;}
    void ComputeRangeBins();
	void ComputeMeanBins();
    const STD::vector<double>& getmeanbins() { return meanbins; }
    const STD::vector<double>& getrangebins() { return rangebins; }
    const STD::vector< STD::vector<double> >& getmeanrangematrix() { return meanrangematrix; }

    void setDebugFlag(STD::string strInputFilePath);

	void setUser_NBins(STD::string bins);

	void setHysteresis(int hysteresisValue);

	void setHysteresis(char* hysteresis);

	int Evaluate(); 

protected:

	void ExtractPeaksFromSignalData ();

	void HysteresisGateCheck();

	int ScanMax() ;

	int OrganisePeaks();

	void ExtractRainflowCounts();

    void BinExtractedFromRainflowCounts();

	// Utility function
	const char* StringTrim(const char*);

	// void ComputeRangeBins();
	// void ComputeMeanBins();
	void ComputeMeanRangeMatrix();

};

#endif // RAINFLOW_H


