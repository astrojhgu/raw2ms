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
  std::vector<double> ant1_pos{2225079.88002 , -5440041.37753 , -2481724.59803};
  std::vector<double> ant2_pos{2224981.09778 , -5440131.25039 , -2481621.06637};
  
  
  std::vector<double> bl_coord;
  for(int i=0;i<3;++i)
    {
      bl_coord.push_back(ant2_pos[i]-ant1_pos[i]);
    }
  MBaseline bl(MVBaseline(MVPosition(bl_coord[0],
				     bl_coord[1],
				     bl_coord[2])),
	       MBaseline::ITRF);
  
  
  double current_time=4860027420.0;
  double ra=1.40920681045;
  double dec=-0.636322083578;
  casa::Vector<double> uvw(mscreate::calc_uvw(bl,current_time,ra,dec));
  double u=uvw(0);
  double v=uvw(1);
  double w=uvw(2);
  cout<<u<<" "<<v<<" "<<w<<endl;
}
