// Code that exercises potential C++ regex library bug
#include <regex>
#include <iostream>
#include <string>
using namespace std;

int main(int argc,
	 char *argv[])
{
  string input = "CAT[23]";
  std::regex e("^(.*)\\[([0-9]*)\\]$", std::regex::extended);
  cout <<  std::regex_replace(input, e, "\\1;\\2",
			      std::regex_constants::match_default | std::regex_constants::format_sed)
			      << endl;
  exit(0);
}
