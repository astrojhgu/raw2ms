//#define AIPS_ARRAY_INDEX_CHECK 1
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
//#include <ms/MeasurementSets/MSMainColumns.h>
#include <ms/MeasurementSets.h>
#include <measures/Measures/MBaseline.h>

#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <date_time.hpp>
#include <map>
#include <string>
#include <vector>
#include <fio.h>

//using namespace blitz;
using namespace std;
using namespace casa;
using namespace ulastai;
const double pi=atan(1)*4;
const int nch=8192;
const int img_size=1024;
const double max_bl=2640;
const double max_freq=200E6;
const double c=2.99792458E8;
const double max_uv=max_bl/(c/max_freq);
int main(int argc,char* argv[])
{
  if(argc!=3)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <ms table> <outfile>"<<endl;
      return -1;
    }
  
  MeasurementSet mstab(argv[1],TableLock(TableLock::AutoNoReadLocking));
  
  cout<<mstab.nrow()<<std::endl;

  ROMSColumns columns(mstab);

  const ROArrayColumn< Double > & uvw_column(columns.uvw());

  const ROMSDataDescColumns& desc(columns.dataDescription());

  const ROScalarColumn< Int >& spwId(desc.spectralWindowId());

  const ROMSSpWindowColumns& spw(columns.spectralWindow());

  const ROArrayColumn< Double >& chan_freq_column(spw.chanFreq());

  const ROScalarColumn<Int>& descID(columns.dataDescId());

  const ROArrayColumn< Bool >& flag_column(columns.flag());
  //cout<<descID.nrow()<<endl;

  //const Vector<Double> v(uvw.get(100000));
  //const ROArrayColumn< Complex >& data_column(columns.data());
  const ROArrayColumn< Complex >& data_column(columns.correctedData());
  //const ROArrayColumn< Complex >& data_column(columns.modelData());
  //const ArrayColumn<Complex> data_column(mstab,"MODEL_DATA");
  //std::cout<<v<<std::endl;

  blitz::Array<double,2> mxr(img_size,img_size);
  blitz::Array<double,2> mxi(img_size,img_size);
  blitz::Array<long,2> cnt(img_size,img_size);
  mxr=0;
  mxi=0;
  cnt=0;
  

  for(int i=0;i<columns.nrow();++i)
    {
      int did=descID.get(i);
      int spwid=spwId.get(did);
      const Vector<double> chan_freq(chan_freq_column.get(spwid));
      const Vector<Double> uvw(uvw_column.get(i));
      const Vector<Complex> data(data_column.get(i));
      const Vector<Bool> flag(flag_column.get(i));
      //cout<<i<<endl;
      if(i%10000==0)
	{
	  cout<<i/(double)columns.nrow()<<" "<<chan_freq.size()<<endl;
	}

      for(int ch=0;ch<data.size();++ch)
	{
	  //cout<<d.real()<<endl;
	  if(!flag[ch])
	    {
	      Complex d(data[ch]);
	      double freq=chan_freq[ch];
	      double lambda=c/freq;
	      double u_l=uvw[0]/lambda;
	      double v_l=uvw[1]/lambda;
	      int iu=u_l/max_uv*(img_size/2)+img_size/2;
	      int iv=v_l/max_uv*(img_size/2)+img_size/2;

	      if(iu>=0&&iu<img_size&&iv>=0&&iv<img_size)
		{
		  mxr(iu,iv)+=d.real();
		  mxi(iu,iv)+=d.imag();
		  cnt(iu,iv)+=1;
		}
	    }
	} 
    }

  for(int i=0;i<img_size;++i)
    {
      for(int j=0;j<img_size;++j)
	{
	  if(cnt(i,j)>0)
	    {
	      mxr(i,j)/=cnt(i,j);
	      mxi(i,j)/=cnt(i,j);
	    }
	}
    }
  std::string prefix(argv[2]);
  cfitsfile ff;
  ff.create((prefix+"_r.fits").c_str());
  ff<<mxr;
  ff.close();
  ff.create((prefix+"_i.fits").c_str());
  ff<<mxi;
  ff.close();
}
