// $Id: FailTest.cpp 6141 2014-05-05 21:03:47Z chambm $
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//


#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    cerr << "This is a test of the emergency test failure system." << endl;
    cout << "DON'T PANIC!" << endl;
    cerr << "This is only a test." << endl;
    return 42;
}