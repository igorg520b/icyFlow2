#include "cz.h"

icy::CZ::CZ()
{

}

void icy::CZ::AcceptTentativeValues()
{
    for(int i=0;i<3;i++) {
        pmax[i] = _pmax;
        tmax[i] = _tmax;
    }
    if(_failed) failed = true;
    avgDn = _avgDn;
    avgDt = _avgDt;
    avgTn = _avgTn;
    avgTt = _avgTt;
    if (maxAvgDn < avgDn) maxAvgDn = avgDn;
    if (maxAvgDt < avgDt) maxAvgDt = avgDt;
}
