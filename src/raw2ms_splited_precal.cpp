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
#include <measures/Measures/MBaseline.h>

#include <cassert>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <date_time.hpp>
#include <fio.h>
using namespace ulastai;
using namespace casa;
using namespace std;

const double c=2.99792458E8;//m/s
const double freq_per_ch=200E6/8192.0;
const double pi=atan(1)*4;
const double dt=3.2;
const int ch_max=8191;
const int ch_min=0;
const double RA=0;
const double Dec=pi/2;
const bool writeAutoCorr=false;
const int nchannels=8192;
const int nantennas=40;
const int img_size=1024;
const double max_uv=2640.00/(c/200E6);
blitz::Array<double,2> mx(img_size,img_size);
blitz::Array<long,2> wgt(img_size,img_size);

std::vector<vector<double> > ant_pos_tab;

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
  const std::vector<std::vector<std::vector<std::complex<float> > > >& gain_vec;
  std::vector<int> time_step;
  std::ifstream time_file;
  double current_time;
  double start_time;
  long current_idx;
  long current_time_step_id;
public:
  visb_by_baseline_source(const std::vector<std::string>& antenna_names,
			  const std::vector<std::vector<std::vector<std::complex<float> > > >& _gain_vec,
			  const std::vector<int>& _time_step,
			  const std::string& data_path,
			  const std::string& date,
			  const std::vector<std::pair<int,int> >& _chlimits)
    :gain_vec(_gain_vec),time_step(_time_step),raw_data_source(antenna_names.size()*(antenna_names.size()+1)/2),chlimits(_chlimits),time_file(data_path+"/time-0-"+date+".txt"),current_time(0),current_idx(0),current_time_step_id(0)
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
    std::vector<std::complex<float> > df(nchannels);
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

	std::vector<double>& ant1_pos=ant_pos_tab.at(antenna_pair.first);
	std::vector<double>& ant2_pos=ant_pos_tab.at(antenna_pair.second);

	double blx=ant2_pos[0]-ant1_pos[0];
	double bly=ant2_pos[1]-ant1_pos[1];
	double blz=ant2_pos[2]-ant1_pos[2];
	
	MBaseline bl(MVBaseline(MVPosition(blx,bly,blz)),
		     MBaseline::ITRF);
	
	casa::Vector<double> uvw(mscreate::calc_uvw(bl,current_time,0,pi/2));
	double u=uvw[0];
	double v=uvw[1];

	//std::cerr<<u<<" "<<v<<std::endl;
	//double delay=delay_vec[antenna_pair.second]-delay_vec[antenna_pair.first];
	files[i]->read((char*)(df.data()),nchannels*2*sizeof(float));
	assert(files[i]->good());
	for(int j=0;j<chlimits.size();++j)
	  {
	    int ch_lower,ch_upper;
	    ch_lower=chlimits[j].first;
	    ch_upper=chlimits[j].second;
	    for(int k=ch_lower;k!=ch_upper;++k)
	      {
		std::complex<float> gg=gain_vec.at(k).at(current_time_step_id).at(antenna_pair.second)*std::conj(gain_vec.at(k).at(current_time_step_id).at(antenna_pair.first));
		
		double freq=k/(double)nchannels*200E6;
		double lambda=c/freq;
		double u_lambda=u/lambda;
		double v_lambda=v/lambda;
		//double dphase=delay*1.47/(c/freq)*2*pi;
		auto corrected_data=df[k]*gg;
		data_buffer[j][i].data(IPosition(2,0,k-ch_lower))=corrected_data;

		
		if(std::isnan(corrected_data.real())||
		   std::isnan(corrected_data.imag()))
		  {
		    data_buffer[j][i].flag(IPosition(2,0,k-ch_lower))=true;
		  }
		else
		  {
		    data_buffer[j][i].flag(IPosition(2,0,k-ch_lower))=false;

		    int iu=u_lambda/(max_uv)*(img_size/2)+img_size/2;
		    int iv=v_lambda/(max_uv)*(img_size/2)+img_size/2;
		    
		    if(iu>=0&&iu<img_size&&
		       iv>=0&&iv<img_size&&
		       antenna_pair.first!=antenna_pair.second)
		      {
			//std::cerr<<iu<<" "<<iv<<std::endl;
			//assert(std::isfinite(corrected_data.real()));
			mx(iu,iv)+=corrected_data.real();
			wgt(iu,iv)+=1;
		      }   
		  }
	      }
	  }
      }
    ++current_idx;
    if(time_step.at(current_time_step_id)<=current_idx)
      {
	++current_time_step_id;

	blitz::Array<double,2> img(img_size,img_size);
	img=0;
	for(int i=0;i<img_size;++i)
	  {
	    for(int j=0;j<img_size;++j)
	      {
		if(wgt(i,j)>0)
		  {
		    img(i,j)=mx(i,j)/wgt(i,j);
		  }
		else
		  {
		    img(i,j)=0;
		  }
	      }
	  }
	cfitsfile ff;
	ff.create("uv.fits");
	ff<<img;
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
  if(argc<8)
    {
      std::cerr<<"Usage:"<<argv[0]<<" <antenna table> <gain_prefix> <step file> <out name> <date> <input path> <ch1:ch2> [ch3:ch4] [ch5:ch6]..."<<std::endl;
      return -1;
    }
  mx=0;
  wgt=0;
  std::string antenna_tab_name(argv[1]);
  std::string gain_prefix(argv[2]);
  std::string time_step_name(argv[3]);
  std::string out_prefix(argv[4]);
  std::string date(argv[5]);
  std::string input_path(argv[6]);


  std::vector<std::vector<std::vector<std::complex<float> > > > gain_vec(nchannels);
  for(int i=2048;i<nchannels;++i)
    {
      std::string ch_str=std::to_string(i);
      std::ifstream ifs_gain((gain_prefix+ch_str+".dat").c_str());
      std::cerr<<gain_prefix+ch_str+".dat"<<std::endl;
      assert(ifs_gain.is_open());
      for(int j=0;;++j)
	{
	  std::vector<float> gain_buffer(nantennas*2);
	  ifs_gain.read((char*)gain_buffer.data(),sizeof(float)*2*nantennas);
	  if(!ifs_gain.good())
	    {
	      break;
	    }
	  
	  std::vector<std::complex<float> > gain(nantennas);
	  for(int k=0;k<nantennas;++k)
	    {
	      float ampl=gain_buffer[k+nantennas];
	      float phase=gain_buffer[k];
	      gain.at(k)=std::complex<float>(ampl*std::cos(phase),ampl*std::sin(phase));
	      if(k==10&&i==5000)
		{
		  std::cerr<<gain[k]<<std::endl;
		}
	    }
	  gain_vec[i].push_back(gain);
	}
    }

  std::vector<int> time_step_vec;
  std::ifstream time_step_file(time_step_name.c_str());
  for(;;)
    {
      int s=0;
      int f=0;
      time_step_file>>s>>f;
      if(!time_step_file.good())
	{
	  break;
	}
      time_step_vec.push_back(s);
    }


  std::vector<std::pair<int,int> > chlimits;
  for(int i=7;i<argc;++i)
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
  ofstream ant_pos_dump("ant_pos.qdp");
  ant_pos_dump<<std::setprecision(20)<<endl;
  for(int i=0;i<its_ant_pos.shape()[1];++i)
    {
      mx+=ant_pos(0,i);
      my+=ant_pos(1,i);
      mz+=ant_pos(2,i);
      std::vector<double> pos(3);
      pos[0]=ant_pos(0,i);
      pos[1]=ant_pos(1,i);
      pos[2]=ant_pos(2,i);
      ant_pos_tab.push_back(pos);
      ant_pos_dump<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
      ant_pos_dump.flush();
    }
  ant_pos_dump.close();
  mx/=its_ant_pos.shape()[1];
  my/=its_ant_pos.shape()[1];
  mz/=its_ant_pos.shape()[1];
  
  casa::MPosition array_pos(casa::MVPosition(mx,my,mz),MPosition::ITRF);

  
  ROScalarColumn<String> antNameCol(ant_tab,"NAME");
  Array<String> antNames(antNameCol.getColumn());
  std::vector<std::string> antNameVec;
  std::map<std::string,double> ant_delay;
  
  for(int i=0;i<nantennas;++i)
    {
      antNameVec.push_back(antNames(IPosition(2,i,0)));
    }

  visb_by_baseline_source vbs(antNameVec,
			      gain_vec,
			      time_step_vec,
			      input_path,
			      date,
			      chlimits);

  
  //mscreate msmaker(out_prefix, vbs.get_start_time(), 1,  ant_tab, true, "flag", 8);
  std::vector<std::shared_ptr<mscreate> > msmakers;
  for(int i=0;i<chlimits.size();++i)
    {
      int ch_lower=chlimits[i].first;
      int ch_upper=chlimits[i].second;

      std::string out_name=out_prefix+std::to_string(ch_lower)+"-"+std::to_string(ch_upper)+".MS";
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
      cerr<<"reading..."<<endl;
      if(!vbs.fetch_one())
	{
	  break;
	}
      cerr<<i<<endl;
      cerr<<"writing"<<endl;
      for(auto i:msmakers)
	{
	  i->write_time_step(vbs);
	}
    }
  
  return 0;
}
