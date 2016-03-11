#include <MSCreate.h>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/OS/Path.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace ulastai;
using namespace casa;
using namespace std;

// Define the variables shared between the functions.
Table itsAntTab("ANTENNA", TableLock(TableLock::AutoNoReadLocking));
vector<double> itsRa;
vector<double> itsDec;
Array<double>  itsAntPos;
bool   itsWriteAutoCorr;
//bool   itsWriteImagerCol;
int    itsNPart;
int    itsNBand;
int    itsNFreq[]={16,32,64};
int    itsTileSizeFreq=-1;
int    itsTileSizeRest=-1;
int    itsNFlags=8;
vector<double> itsStartFreq;
vector<double> itsStepFreq;
double itsStartTime;
double itsEndTime;
string itsMsName="obs.MS";
string itsFlagColumn="flag";




class FakeDataSource
  :public RawDataSource
{
private:
  double t0;
  double t;
  std::vector<std::pair<int,int> > bl;
public:
  FakeDataSource()
    :RawDataSource(780),t0(0),t(0)
  {
    Quantity qn;
    MVTime::read (qn, "2016/02/03 00:00:00", true);
    t=qn.getValue ("s");
    t0=t;
    for(int i=0;i<39;++i)
      {
	for(int j=i+1;j<40;++j)
	  {
	    bl.push_back(std::make_pair(i,j));
	  }
      }
    
  }

  bool doFetchOne()
  {
    t+=3;
    if(t-t0>3600)
      {
	return false;
      }
    else
      {
	return true;
      }
  }

  std::pair<int,int> doAntennaPair(int b)const
  {
    return bl[b];
  }

  casa::Array<casa::Complex> doData(int field,int band,int bl)const
  {
    assert(band>=0&&band<itsNBand);
    int n=itsNFreq[band];
    IPosition shape(2,1,n);
    Array<Complex> d(shape);
    d=Complex(1,1);
    return d;
  }

  casa::Array<casa::Bool> doFlags(int field,int band,int bl)const
  {
    assert(band>=0&&band<itsNBand);
    int n=itsNFreq[band];
    IPosition shape(2,1,n);
    Array<Bool> d(shape);
    d=false;
    return d;
  }

  double doTime()const
  {
    return t;
  }

  

}fds;




void createMS (int nband, int bandnr, const string& msName)
{
  //int nfpb = itsNFreq/itsNBand;
  MSCreate msmaker(msName, itsStartTime, 1,
                   itsAntTab, itsWriteAutoCorr,
		   itsFlagColumn, itsNFlags);
  for (int i=0; i<nband; ++i) {
    // Determine middle of band.
    
    double freqRef = itsStartFreq[bandnr] + itsNFreq[i]*itsStepFreq[bandnr]/2;
    msmaker.addBand (itsNFreq[i], freqRef, itsStepFreq[bandnr]);
    ++bandnr;
  }
  for (uint i=0; i<itsRa.size(); ++i) {
    msmaker.addField (itsRa[i], itsDec[i]);
  }
  //for (int i=0; i<itsNTime; ++i) {
  for(;fds.fetchOne();)
    {
      std::cerr<<setprecision(20);
      //std::cerr<<t<<std::endl;
      msmaker.writeTimeStep(fds);
    }
}

void doOne (int seqnr, const string& msName)
{
  int nbpp = itsNBand / itsNPart;
  // Form the MS name.
  // If it contains %d, use that to fill in the seqnr.
  // Otherwise append _seqnr to the name.
  string name=msName;
  //if (msName.find ("%d") != string::npos) {
  //name = formatString (msName.c_str(), seqnr);
  //} else {
  //name = msName + toString (seqnr, "_p%d");
    //}
  // Create the MS.
  createMS (itsNBand, seqnr*nbpp, name);
}

void doAll()
{
  for (int i=0; i<itsNPart; ++i) {
    doOne (i, itsMsName);
  }

}


int main (int argc, char** argv)
{
  //doAll();
  const double pi=3.1415926;
  itsRa.push_back(0);
  itsDec.push_back(pi/2.);
  
  
  //ROArrayColumn<double> posCol(tab, "POSITION");
  //itsAntPos = posCol.getColumn();

  //cout<<itsAntPos.shape()[0]<<endl;
  //cout<<itsAntPos.shape()[1]<<endl;
  itsWriteAutoCorr=false;
  

  itsNPart=1;
  itsNBand=3;
  
    
  itsNPart=1;
  //doOne(1,"aa.MS");


  for(int i=0;i<itsNBand;++i)
    {
      itsStartFreq.push_back(50E6+i*10E6);
      itsStepFreq.push_back(200e6/8192);
    }

  
  Quantity qn;
  assert(MVTime::read (qn, "2016/02/03 00:00:00", true));
  itsStartTime = qn.getValue ("s");
  
  doOne(0,itsMsName);
  
  return 0;
}
