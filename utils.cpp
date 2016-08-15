using namespace std;
#include "utils.h"

//Converts a float in a string that can be used to write correctly floats to describe colors in OFF files
string strOfFloat(const float& f){
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4) << f;
    std::string mystring = ss.str();
    return mystring;
}
