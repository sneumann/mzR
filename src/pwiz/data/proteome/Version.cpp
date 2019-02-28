// This file was automatically generated.
// You should not modify it manually, as it may be re-generated.
//
// Licensed under the Apache License, Version 2.0 (the License); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an AS IS BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#include "Version.hpp"
#include <sstream>

#ifdef PWIZ_USER_VERSION_INFO_H // in case you need to add any info version of your own
#include PWIZ_USER_VERSION_INFO_H  // must define PWIZ_USER_VERSION_INFO_H_STR for use below
#endif

namespace pwiz {
namespace proteome {


int Version::Major()                {return 3;}
int Version::Minor()                {return 0;}
int Version::Revision()             {return 19059;}
std::string Version::LastModified()   {return "c161365";}
std::string Version::Branch()   {return "master";}
std::string Version::str()
{
	std::ostringstream v;
	v << Major() << '.' << Minor() << '.' << Revision();
#ifdef PWIZ_USER_VERSION_INFO_H
	v << " (" << PWIZ_USER_VERSION_INFO_H_STR << ")";
#endif
	return v.str();
}


} // namespace proteome
} // namespace pwiz
