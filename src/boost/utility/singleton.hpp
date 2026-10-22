// 
// Copyright (c) 2003-2005 David Abrahams
// Copyright (c) 2004 Daniel Weinig
// Distributed under the Boost Software License, Version 1.0.
// See http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_SINGLETON_HPP_INCLUDED
#define BOOST_SINGLETON_HPP_INCLUDED

#include <boost/noncopyable.hpp>

namespace boost {

// Tag type for restricted constructor access
struct restricted {};

namespace detail {

template<typename T>
struct singleton_holder
{
    static T& instance()
    {
        static T obj(restricted());
        return obj;
    }
};

} // namespace detail

template<typename T>
class singleton : private noncopyable
{
public:
    static T& instance()
    {
        return detail::singleton_holder<T>::instance();
    }
};

} // namespace boost

#endif // BOOST_SINGLETON_HPP_INCLUDED
