#include "../include/util.h"


/********* helper methods ******************/

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
Error(const std::string & s)
{
  std::cerr<<s<<std::endl;
  exit(1);
};

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
bool
util::
has_only_spaces(const std::string &str)
{
  for (std::string::const_iterator it = str.begin(); it != str.end(); ++it)
    {
      if (*it != ' ') return false;
    }
  return true;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<std::string>
util::
strip_white_space(std::vector<std::string> &s)
{

  std::vector<std::string> my_s ;

  for (auto &i : s){
    //my_s.push_back(std::regex_replace(i,std::regex("//s+"), "")) ;
    i.erase(i.find_last_not_of(" \n\r\t")+1);;
    i.erase(0, i.find_first_not_of(" \n\r\t"));;
    my_s.push_back(i);
  }

  return my_s ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
pp(const std::string &s)
{

  std::cout << std::endl << std::endl << s << std::endl << std::endl ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
print(std::vector<std::string> &s)
{
  std::cout << std::endl;
  for (auto i : s) std::cout << i << " " ;
  std::cout << std::endl ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
print(std::vector<float> &v)
{
  std::cout << std::endl;
  for (auto i : v) std::cout << i << " " ;
  std::cout << std::endl ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
print(std::vector<int> &v)
{
  std::cout << std::endl;
  for (auto i : v) std::cout << i << " " ;
  std::cout << std::endl ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
print_run_details()
{
  std::string time_and_user_string =  time_and_user();

  std::cout << "execution details: " << time_and_user_string << std::endl ;

  std::string my_path(getcwd(NULL,0));
  std::cout << "working directory : " << my_path << std::endl ;

  std::cout << std::endl;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
compile_time(const char *my_file, const char *my_func, const char *my_date, const char *my_time)
{
  std::cout << std::endl;
  std::cout << "file = " << my_file << " : function = " << my_func << std::endl ;
  std::cout << "compiled on " << my_date << " at " << my_time << " " << std::endl ;
  print_run_details() ;

}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::string
util::
time_and_user()
{

  std::string time_and_user_string ;

  time_t rawtime;
  struct tm * timeinfo;
  char *lgn;
  struct passwd *pw;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  time_and_user_string += asctime(timeinfo) ;
  time_and_user_string.pop_back() ;


  std::string username = " by user: " ;

  //sh if ((lgn = getlogin()) == NULL || (pw = getpwnam(lgn)) == NULL) {
  //sh   fprintf(stderr, "Get of user information failed.\n"); exit(1); // ??? WHY EXIT IF THIS FAILS ??? it breaks on nohup
  //sh }
  if ((lgn = getlogin()) == NULL) {
    // std::cout << "getlogin failed: NULL\n" << std::endl ;
    return time_and_user_string + username + "nohup" ;
  }
  if ((pw = getpwnam(lgn)) == NULL) {
    // std::cout << "getpwnam failed: NULL\n" << std::endl ;
    return time_and_user_string + username + "nohup" ;
  }

  return time_and_user_string + username + pw->pw_name ;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<int>
util::
unique_vector_int(const std::vector<int> & in)
{
  std::vector<int> out(in.size());
  std::vector<int>::iterator it;
  it = std::unique_copy(in.begin(), in.end(), out.begin());
  std::sort(out.begin(), it);
  it = std::unique_copy(out.begin(), it, out.begin());
  out.resize(std::distance(out.begin(),it));
  return out;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<std::string>
util::
unique_vector_string(const std::vector<std::string> & in)
{
  std::vector<std::string> out = strip_white_space(const_cast<std::vector<std::string>&>(in));
  std::vector<std::string>::iterator it;
  it = std::unique_copy(in.begin(), in.end(), out.begin());
  std::sort(out.begin(), it);
  it = std::unique_copy(out.begin(), it, out.begin());
  out.resize(std::distance(out.begin(),it));
  return out;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
unsigned long int
util::
random_seed()
{
  unsigned long int seed;
  struct timeval tv;

  gettimeofday(&tv,0);
  seed = tv.tv_sec + tv.tv_usec;

  return seed;
  //return 0;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
util::
get_random_float(const float a, const float b)
{
  return ((b-a)*((float)rand()/RAND_MAX))+a;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
long long
util::
diff_ms(struct timeval t1, struct timeval t2)
{
  return (((t1.tv_sec - t2.tv_sec) * 1000000) +
          (t1.tv_usec - t2.tv_usec));
  //(t1.tv_usec - t2.tv_usec))/1000;
}

///
/// @par Detailed description
/// ...
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
util::
profiling(struct timeval & t1, struct timeval & t2, const std::string & moduleName)
{
  gettimeofday(&t2, NULL);
  //printf("PROFILING INFO >> module: %20s starts: %20ld takes: %20lld microseconds\n", moduleName.c_str(), long((t1.tv_sec)*1000000L+(t1.tv_usec)), util::diff_ms(t2, t1));
  printf("PROFILING INFO >> module: %20s takes: %20lld microseconds\n", moduleName.c_str(), util::diff_ms(t2, t1));
  gettimeofday(&t1, NULL);
}
