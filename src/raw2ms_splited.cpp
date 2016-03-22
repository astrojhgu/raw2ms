#include <mscreate.hpp>
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
#include <casa/Quanta/MVPosition.h>
#include <measures/Measures/MPosition.h>
#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <date_time.hpp>

using namespace ulastai;
using namespace casa;
using namespace std;

const double c=2.99792458E8;//m/s
const double freq_per_ch=200E6/8192.0;
const double pi=atan(1)*4;
const int ch_max=8191;
const int ch_min=0;
const double RA=0;
const double Dec=pi/2;
const bool writeAutoCorr=false;
const int nchannels=8192;

struct data_flag_pair
{
  casa::Array<casa::Complex> data;
  casa::Array<casa::Bool> flag;
};

class visb_by_baseline_source
  :public raw_data_source
{
private:
  std::vector<shared_ptr<ifstream> > files;
  std::vector<std::pair<int,int> > baseline_nodes;
  std::vector<std::vector<data_flag_pair> > data_buffer;
  std::vector<std::pair<int,int> > chlimits;
  std::ifstream time_file;
  double current_time;
  double start_time;
public:
  visb_by_baseline_source(const std::vector<std::string>& antenna_names,
		       const std::string& data_path,
		       const std::string& date,
		       const std::vector<std::pair<int,int> >& _chlimits)
    :raw_data_source(antenna_names.size()*(antenna_names.size()+1)/2),chlimits(_chlimits),time_file(data_path+"/time-0-"+date+".txt"),current_time(0)
  {
    assert(time_file.is_open());
    ifstream ifs_time_tmp(data_path+"/time-0-"+date+".txt");
    std::string time_line;
    getline(ifs_time_tmp,time_line);
    start_time=parse_21cma_date(time_line)-3.2;
    ifs_time_tmp.close();
    
    int nantennas=antenna_names.size();
    int nbl=nantennas*(nantennas+1)/2;
    for(int i=0;i<nantennas;++i)
      {
	for(int j=i;j<nantennas;++j)
	  {
	    baseline_nodes.push_back({i,j});
	    files.push_back(std::shared_ptr<std::ifstream>(new ifstream(data_path+"/"+antenna_names[i]+antenna_names[j]+"-0-"+date+".bin")));
	    assert(files.back()->is_open());
	  }
      }
    for(int i=0;i<chlimits.size();++i)
      {
	IPosition data_shape(2,1,chlimits[i].second-chlimits[i].first);
	std::vector<data_flag_pair> d(nbl);
	for(int j=0;j<nbl;++j)
	  {
	    d[j].data.resize(data_shape);
	    d[j].flag.resize(data_shape);
	  }
	data_buffer.push_back(d);
      }
  }

  double get_start_time()const
  {
    return start_time;
  }

  bool do_fetch_one()override
  {
    std::vector<float> df(nchannels*2);
    std::string time_line;
    std::getline(time_file,time_line);
    if(!time_file.good())
      {
	return false;
      }

    current_time=parse_21cma_date(time_line);
    
    for(int i=0;i<files.size();++i)
      {
	files[i]->read((char*)(df.data()),nchannels*2*sizeof(float));
	assert(files[i]->good());
	for(int j=0;j<chlimits.size();++j)
	  {
	    int ch_lower,ch_upper;
	    ch_lower=chlimits[j].first;
	    ch_upper=chlimits[j].second;
	    for(int k=ch_lower;k!=ch_upper;++k)
	      {
		data_buffer[j][i].data(IPosition(2,0,k-ch_lower))=Complex(df[k*2],df[k*2+1]);
		if(std::isnan(df[k*2])||
		   std::isnan(df[k*2+1]))
		  {
		    data_buffer[j][i].flag(IPosition(2,0,k-ch_lower))=true;
		  }
		else
		  {
		    data_buffer[j][i].flag(IPosition(2,0,k-ch_lower))=false;
		  }
		
	      }
	  }
      }
    return true;
  }

  std::pair<int,int> do_antenna_pair(int bl)const
  {
    return baseline_nodes.at(bl);
  }

  casa::Array<casa::Complex> do_data(int field,int band,int bl)const
  {
    return data_buffer.at(band).at(bl).data;
  }

  casa::Array<casa::Bool> do_flags(int field,int band,int bl)const
  {
    return data_buffer.at(band).at(bl).flag;
  }

  double do_time()const
  {
    return current_time;
  }
};

std::pair<int,int> parse_ch(const std::string& s)
{
  std::pair<int,int> result;
  size_t n=s.find_first_of(":");
  assert(n!=std::string::npos);
  std::string s1=s.substr(0,n);
  std::string s2=s.substr(n+1,s.size());
  istringstream iss;
  iss.str(s1);
  iss>>result.first;
  iss.clear();
  iss.str(s2);
  iss>>result.second;
  assert(result.second>result.first);
  return result;
}

int main (int argc, char** argv)
{
  if(argc<6)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <antenna table> <out name> <date> <input path> <ch1:ch2> [ch3:ch4] [ch5:ch6]..."<<std::endl;
      return -1;
    }
  std::string antenna_tab_name(argv[1]);
  std::string out_prefix(argv[2]);
  std::string date(argv[3]);
  std::string input_path(argv[4]);

  std::vector<std::pair<int,int> > chlimits;
  for(int i=5;i<argc;++i)
    {
      std::string chs(argv[i]);
      chlimits.push_back(parse_ch(chs));
      assert(chlimits.back().first<chlimits.back().second);
      assert(chlimits.back().first>=2048&&chlimits.back().second<=8192);
    }
  
  Table ant_tab(antenna_tab_name, TableLock(TableLock::AutoNoReadLocking));

  ROArrayColumn<double> pos_col(ant_tab, "POSITION");
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
  
  casa::MPosition array_pos(casa::MVPosition(mx,my,mz),MPosition::ITRF);

  
  ROScalarColumn<String> antNameCol(ant_tab,"NAME");
  Array<String> antNames(antNameCol.getColumn());
  std::vector<std::string> antNameVec;
  
  for(int i=0;i<40;++i)
    {
      antNameVec.push_back(antNames(IPosition(2,i,0)));
      std::cout<<antNameVec.back()<<std::endl;
    }

  visb_by_baseline_source vbs(antNameVec,
			   input_path,
			   date,
			   chlimits);

  
  //mscreate msmaker(out_prefix, vbs.get_start_time(), 1,  ant_tab, true, "flag", 8);
  std::vector<std::shared_ptr<mscreate> > msmakers;
  for(int i=0;i<chlimits.size();++i)
    {
      int ch_lower=chlimits[i].first;
      int ch_upper=chlimits[i].second;

      std::string out_name=out_prefix+std::to_string(ch_lower)+":"+std::to_string(ch_upper)+".MS";
      auto p=std::shared_ptr<mscreate>(new mscreate(out_name,vbs.get_start_time(),1,ant_tab,array_pos,true));
      
      msmakers.push_back(p);
      int nch=ch_upper-ch_lower;
      //double fref=ch_lower*freq_per_ch+freq_per_ch/2.0;
      double fref=(ch_lower+ch_upper)/2.0*freq_per_ch;
      msmakers.back()->add_band(nch,fref,freq_per_ch);
      msmakers.back()->add_field(0.0,pi/2);
    }

  
  for(int i=0;;++i)
    {
      cout<<"reading..."<<endl;
      if(!vbs.fetch_one())
	{
	  break;
	}
      cout<<i<<endl;
      cout<<"writing"<<endl;
      for(auto i:msmakers)
	{
	  i->write_time_step(vbs);
	}
    }
  
  return 0;
}
