#include <boost/thread.hpp>

extern "C"
{

#include "H5pubconf.h"
#include "H5private.h"
#include "H5TSprivate.h"
#include <stdio.h>

void boost_first_thread_init_really()
{
    H5TS_win32_first_thread_init(NULL, NULL, NULL);
}

int boost_first_thread_init()
{
    static struct {boost::once_flag flag;} flag = {BOOST_ONCE_INIT};
    boost::call_once(flag.flag, boost_first_thread_init_really);
    return 1;
}

}
