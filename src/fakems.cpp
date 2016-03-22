#include <mscreate.hpp>
#include <casa/Quanta/MVTime.h>
#include <casa/Quanta/MVAngle.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Vector.h>
#include <casa/OS/Path.h>
#include <casa/Quanta/MVPosition.h>
#include <measures/Measures/MPosition.h>
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
Table its_ant_tab("ANTENNA", TableLock(TableLock::AutoNoReadLocking));
vector<double> its_ra;
vector<double> its_dec;
Array<double>  its_ant_pos;
bool   its_write_auto_corr;
//bool   itsWriteImagerCol;
int    its_nparts;
int    its_nbands;
int    its_nfreqs[]={16,32,64};
int    its_tile_size_freq=-1;
int    its_tile_size_rest=-1;
int    its_nflags=8;
vector<double> its_start_freq;
vector<double> its_step_freq;
double its_start_time;
double its_end_time;
string its_ms_name="obs.MS";





class fake_data_source
  :public raw_data_source
{
private:
  double t0;
  double t;
  std::vector<std::pair<int,int> > bl;
public:
  fake_data_source()
    :raw_data_source(780),t0(0),t(0)
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

  bool do_fetch_one()
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

  std::pair<int,int> do_antenna_pair(int b)const
  {
    return bl[b];
  }

  casa::Array<casa::Complex> do_data(int field,int band,int bl)const
  {
    //cout<<"bl"<<endl;
    assert(band>=0&&band<its_nbands);
    int n=its_nfreqs[band];
    IPosition shape(2,1,n);
    Array<Complex> d(shape);
    d=Complex(1,1);
    return d;
  }

  casa::Array<casa::Bool> do_flags(int field,int band,int bl)const
  {
    assert(band>=0&&band<its_nbands);
    int n=its_nfreqs[band];
    IPosition shape(2,1,n);
    Array<Bool> d(shape);
    d=false;
    return d;
  }

  double do_time()const
  {
    return t;
  }

  

}fds;




void createms (int nband, int bandnr, const string& ms_name)
{
  //int nfpb = its_nfreqs/its_nbands;
  
  ROArrayColumn<double> pos_col(its_ant_tab, "POSITION");
  Array<double> its_ant_pos(pos_col.getColumn());
  Matrix<double> ant_pos(its_ant_pos);
  double mx(0),my(0),mz(0);
  for(int i=0;i<its_ant_pos.shape()[1];++i)
    {
      mx+=ant_pos(0,i);
      my+=ant_pos(1,i);
      mz+=ant_pos(2,i);
    }
  mx/=its_ant_pos.shape()[1];
  my/=its_ant_pos.shape()[1];
  mz/=its_ant_pos.shape()[1];
  mx=ant_pos(0,0);
  my=ant_pos(1,0);
  my=ant_pos(2,0);
  
  mscreate msmaker(ms_name, its_start_time, 1,
                   its_ant_tab,
		   casa::MPosition(casa::MVPosition(mx,my,mz),MPosition::ITRF),
		   its_write_auto_corr);
  for (int i=0; i<nband; ++i) {
    // Determine middle of band.
    
    double freq_ref = its_start_freq[bandnr] + its_nfreqs[i]*its_step_freq[bandnr]/2;
    msmaker.add_band (its_nfreqs[i], freq_ref, its_step_freq[bandnr]);
    ++bandnr;
  }
  for (uint i=0; i<its_ra.size(); ++i) {
    msmaker.add_field (its_ra[i], its_dec[i]);
  }
  //for (int i=0; i<itsNTime; ++i) {
  for(;fds.fetch_one();)
    {
      std::cerr<<setprecision(20);
      //std::cerr<<t<<std::endl;
      msmaker.write_time_step(fds);
    }
}

void do_one (int seqnr, const string& ms_name)
{
  int nbpp = its_nbands / its_nparts;
  // Form the MS name.
  // If it contains %d, use that to fill in the seqnr.
  // Otherwise append _seqnr to the name.
  string name=ms_name;
  //if (ms_name.find ("%d") != string::npos) {
  //name = formatString (ms_name.c_str(), seqnr);
  //} else {
  //name = ms_name + toString (seqnr, "_p%d");
    //}
  // Create the MS.
  createms (its_nbands, seqnr*nbpp, name);
}

void doAll()
{
  for (int i=0; i<its_nparts; ++i) {
    do_one (i, its_ms_name);
  }

}


int main (int argc, char** argv)
{
  //doAll();
  const double pi=3.1415926;
  its_ra.push_back(0);
  its_dec.push_back(pi/2.);
  
  
  //ROArrayColumn<double> posCol(tab, "POSITION");
  //its_ant_pos = posCol.getColumn();

  //cout<<its_ant_pos.shape()[0]<<endl;
  //cout<<its_ant_pos.shape()[1]<<endl;
  its_write_auto_corr=false;
  

  its_nparts=1;
  its_nbands=3;
  
    
  its_nparts=1;
  //do_one(1,"aa.MS");


  for(int i=0;i<its_nbands;++i)
    {
      its_start_freq.push_back(50E6+i*10E6);
      its_step_freq.push_back(200e6/8192);
    }

  
  Quantity qn;
  assert(MVTime::read (qn, "2016/02/03 00:00:00", true));
  its_start_time = qn.getValue ("s");
  
  do_one(0,its_ms_name);
  
  return 0;
}
