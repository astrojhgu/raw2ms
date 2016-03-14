#include <date_time.hpp>
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
  cout<<setprecision(20);
  cout<<parse_21cma_date("Thu Aug 15 00:00:09 2013")<<endl;
  cout<<parse_21cma_date("Thu Jan 15 00:01:09 2013")<<endl;
}
