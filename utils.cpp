using namespace std;
#include "utils.h"

string strOfFloat(const float& f){
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4) << f;
    std::string mystring = ss.str();
    return mystring;
}
