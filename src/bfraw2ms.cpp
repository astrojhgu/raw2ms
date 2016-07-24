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
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <date_time.hpp>

using namespace ulastai;
using namespace casa;
using namespace std;

const double c=2.99792458E8;//m/s
const double freq_per_ch=200E6/8192.0;
const double dt=11;
const double pi=atan(1)*4;
double RA=0;
double Dec=pi/2;
int nchannels=0;
int ch_beg;
int ch_end;
double delay;

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
  std::vector<data_flag_pair> data_buffer;
  std::pair<int,int> chlimit;
  std::ifstream time_file;
  double current_time;
  double start_time;
public:
  visb_by_baseline_source(const std::vector<std::string>& antenna_names,
			  const std::string& data_path,
			  const std::pair<int,int>& _chlimit)
    :raw_data_source(antenna_names.size()*(antenna_names.size()+1)/2),chlimit(_chlimit),time_file(data_path+"/time.txt"),current_time(0)
  {
    assert(time_file.is_open());
    ifstream ifs_time_tmp(data_path+"/time.txt");
    std::string time_line;
    getline(ifs_time_tmp,time_line);
    start_time=parse_21cma_date(time_line)-dt;
    ifs_time_tmp.close();
    
    int nantennas=antenna_names.size();
    int nbl=nantennas*(nantennas+1)/2;
    for(int i=0;i<nantennas;++i)
      {
	for(int j=i;j<nantennas;++j)
	  {
	    baseline_nodes.push_back({i,j});
	    files.push_back(std::shared_ptr<std::ifstream>(new ifstream(data_path+"/"+antenna_names[i]+antenna_names[j]+".bin")));
	    std::cerr<<data_path+"/"+antenna_names[i]+antenna_names[j]+".bin"<<std::endl;
	    assert(files.back()->is_open());
	  }
      }
    //exit(0);
    IPosition data_shape(2,1,chlimit.second-chlimit.first);
    std::vector<data_flag_pair> d(nbl);
    for(int j=0;j<nbl;++j)
      {
	d[j].data.resize(data_shape);
	d[j].flag.resize(data_shape);
      }
    data_buffer=d;
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
	auto antenna_pair=do_antenna_pair(i);
	
	files[i]->read((char*)(df.data()),nchannels*2*sizeof(float));
	assert(files[i]->good());
	int ch_lower,ch_upper;
	ch_lower=chlimit.first;
	ch_upper=chlimit.second;
	for(int k=0;k!=nchannels;++k)
	  {
	    double freq=(ch_beg+k)*freq_per_ch;
	    double dphase=(antenna_pair.first==antenna_pair.second)?0:delay/(c/freq)*2*pi;
	    data_buffer[i].data(IPosition(2,0,k))=Complex(df[k*2],df[k*2+1])*exp(Complex(0,1)*dphase);
	    if(std::isnan(df[k*2])||
	       std::isnan(df[k*2+1]))
	      {
		data_buffer[i].flag(IPosition(2,0,k))=true;
	      }
	    else
	      {
		data_buffer[i].flag(IPosition(2,0,k))=false;
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
    assert(band==0);
    return data_buffer.at(bl).data;
  }

  casa::Array<casa::Float> do_sigma(int field,int band,int bl)const
  {
    IPosition data_shape(1,1);
    casa::Array<Float> sigma(data_shape);
    auto p=do_antenna_pair(bl);
    if(p.first==p.second)
      {
	sigma=1/std::sqrt(freq_per_ch*dt);
      }
    else
      {
	sigma=1/std::sqrt(2*freq_per_ch*dt);
      }
	  
    return sigma;
  }

  casa::Array<casa::Bool> do_flags(int field,int band,int bl)const
  {
    return data_buffer.at(bl).flag;
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
  if(argc<9)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <antenna table> <out name> <input path> <ra> <dec> <chbeg> <chend> <delay>"<<std::endl;
      return -1;
    }
  std::string antenna_tab_name(argv[1]);
  std::string out_name(argv[2]);
  std::string input_path(argv[3]);
  RA=std::stod(std::string(argv[4]));
  Dec=std::stod(std::string(argv[5]));
  ch_beg=std::stoi(std::string(argv[6]));
  ch_end=std::stoi(std::string(argv[7]));
  nchannels=ch_end-ch_beg;
  std::cerr<<"nch="<<nchannels<<std::endl;
  delay=std::stod(std::string(argv[8]));
  std::cerr<<"delay="<<delay<<std::endl;
  //exit(0);
  
  std::pair<int,int> chlimit{ch_beg,ch_end};

  Table ant_tab(antenna_tab_name, TableLock(TableLock::AutoNoReadLocking));
  ROScalarColumn<String> antNameCol(ant_tab,"NAME");
  Array<String> antNames(antNameCol.getColumn());
  std::vector<std::string> antNameVec;
  
  for(int i=0;i<2;++i)
    {
      antNameVec.push_back(antNames(IPosition(2,i,0)));
      std::cout<<antNameVec.back()<<std::endl;
    }

  visb_by_baseline_source vbs(antNameVec,
			      input_path,
			      chlimit);
  
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
  
  
  mscreate msmaker(out_name, vbs.get_start_time(), 1,  ant_tab, casa::MPosition(casa::MVPosition(mx,my,mz),MPosition::ITRF),true);
  msmaker.set_correct_w(true);
  int nch=ch_end-ch_beg;
  double fref=(ch_beg+ch_end)/2.0*freq_per_ch;
  msmaker.add_band(nch,fref,freq_per_ch);
  msmaker.add_field(RA,Dec);

  for(int i=0;;++i)
    {
      cout<<"reading..."<<endl;
      if(!vbs.fetch_one())
	{
	  break;
	}
      cout<<i<<endl;
      cout<<"writing"<<endl;
      msmaker.write_time_step(vbs);
    }
  
  return 0;
}
