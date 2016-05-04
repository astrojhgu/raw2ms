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

using namespace std;
using namespace casa;
using namespace ulastai;
const double pi=atan(1)*4;
int main(int argc,char* argv[])
{
  if(argc!=5)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <antenna table> <time file> <ant1> <ant2>"<<endl;
      return -1;
    }

  Table ant_tab(argv[1],TableLock(TableLock::AutoNoReadLocking));

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

  cerr<<argv[3]<<" -- "<<argv[4]<<endl;
  auto iter1=ant_pos_map.find(std::string(argv[3]));
  auto iter2=ant_pos_map.find(std::string(argv[4]));

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


  ifstream ifs_time(argv[2]);
  for(;;)
    {
      std::string time_line;
      std::getline(ifs_time,time_line);
      if(!ifs_time.good())
	{
	  break;
	}
      double current_time=parse_21cma_date(time_line);
      casa::Vector<double> uvw(mscreate::calc_uvw(bl,current_time,0,pi/2));
      double u=uvw(0);
      double v=uvw(1);
      double w=uvw(2);
      cout<<u<<" "<<v<<endl;
    }
  
}
