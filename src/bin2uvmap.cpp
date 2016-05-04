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
  if(argc!=8)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <prefix> <antenna table> <date> <ant1> <ant2> <delay> <outfile>>"<<endl;
      return -1;
    }
  std::string prefix(argv[1]);
  
  Table ant_tab(argv[2],TableLock(TableLock::AutoNoReadLocking));

  ROArrayColumn<Double> antPosCol(ant_tab,"POSITION");
  Matrix<double> ant_poses(Array<Double>(antPosCol.getColumn()));

  int nantennas=ant_poses.shape()[1];
  ROScalarColumn<String> antNameCol(ant_tab,"NAME");
  Array<String> antNames(antNameCol.getColumn());
  
  std::map<std::string,std::vector<double> > ant_pos_map;

  for(int i=0;i<nantennas;++i)
    {
      double x=ant_poses(0,i);
      double y=ant_poses(1,i);
      double z=ant_poses(2,i);

      string name=antNames(IPosition(2,i,0));
      cerr<<name<<" "<<x<<" "<<y<<" "<<z<<endl;

      ant_pos_map[name]=std::vector<double>{x,y,z};
    }

  std::string ant1_name(argv[4]);
  std::string ant2_name(argv[5]);
  cerr<<ant1_name<<" -- "<<ant2_name<<endl;
  auto iter1=ant_pos_map.find(ant1_name);
  auto iter2=ant_pos_map.find(ant2_name);

  assert(iter1!=ant_pos_map.end()&&
	 iter2!=ant_pos_map.end());

  auto ant1_pos=iter1->second;
  auto ant2_pos=iter2->second;

  std::vector<double> bl_coord;
  for(int i=0;i<3;++i)
    {
      bl_coord.push_back(ant2_pos[i]-ant1_pos[i]);
    }
  MBaseline bl(MVBaseline(MVPosition(bl_coord[0],
				     bl_coord[1],
				     bl_coord[2])),
	       MBaseline::ITRF);
  std::string date_str(argv[3]);
  std::string time_file_name(prefix+"/time-0-"+date_str+".txt");
  std::string bin_file_name(prefix+"/"+ant1_name+ant2_name+"-0-"+date_str+".bin");
  ifstream ifs_time(time_file_name.c_str());
  ifstream ifs_bin(bin_file_name.c_str());

  assert(ifs_bin.is_open());

  std::vector<std::complex<float> > data(nch);
  blitz::Array<double,2> mx(img_size,img_size);
  mx=0;
  blitz::Array<long,2> cnt(img_size,img_size);
  cnt=0;
  double delay=std::stod(argv[6]);
  for(;;)
    {
      std::string time_line;
      std::getline(ifs_time,time_line);
      if(!ifs_time.good())
	{
	  break;
	}
      ifs_bin.read((char*)(data.data()),nch*2*sizeof(float));
      assert(ifs_bin.good());
      
      double current_time=parse_21cma_date(time_line);
      casa::Vector<double> uvw(mscreate::calc_uvw(bl,current_time,0,pi/2));
      double u=uvw(0);
      double v=uvw(1);
      double w=uvw(2);
      //cout<<u<<" "<<v<<endl;

      for(int i=0;i<data.size();++i)
	{
	  if (i<6144||i>=6400)
	    {
	      //continue;
	    }
	  double freq=i/(double)nch*max_freq;
	  double lbd=c/freq;
	  double u1=u/lbd;
	  double v1=v/lbd;
	  int iu=u1/max_uv*(img_size/2)+img_size/2;
	  int iv=v1/max_uv*(img_size/2)+img_size/2;
	  if(iu>=0&&iu<img_size&&
	     iv>=0&&iv<img_size)
	    {
	      auto d=data[i];
	      double phase=delay*1.47/lbd*2*pi;
	      d*=std::exp(std::complex<float>(0,1)*phase);
	      mx(iu,iv)=d.real();
	      cnt(iu,iv)=1;
	    }
	}
    }

  for(int i=0;i<img_size;++i)
    {
      for(int j=0;j<img_size;++j)
	{
	  if(cnt(i,j)>0)
	    {
	      mx(i,j)/=cnt(i,j);
	    }
	}
    }

  cfitsfile ff;
  ff.create(argv[7]);
  ff<<mx;
  
}
