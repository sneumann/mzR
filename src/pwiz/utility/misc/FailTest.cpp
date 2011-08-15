// $Id: FailTest.cpp 1814 2010-02-16 22:52:44Z chambm $

#include <iostream>
using namespace std;

int main()
{
    cerr << "This is a test of the emergency test failure system." << endl;
    cout << "DON'T PANIC!" << endl;
    cerr << "This is only a test." << endl;
    return 42;
}
