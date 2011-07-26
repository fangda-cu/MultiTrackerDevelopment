/*
 * Typelist.hh
 *
 *  Created on: 20/07/2011
 *      Author: jaubry
 *
 *  Ideas and code from Modern C++ Design by Andrei Alexandrescu
 *
 */

#ifndef TYPELIST_HH_
#define TYPELIST_HH_

namespace strandsim
{

class NullType
{
};

template<typename T, typename U> struct Typelist
{
    typedef T Head;
    typedef U Tail;
};

template<typename TList> struct Length;

template<> struct Length<NullType>
{
    enum
    {
        value = 0
    };
};

template<typename T, typename U> struct Length<Typelist<T, U> >
{
    enum
    {
        value = Length<U>::value + 1
    };
};

}

#endif /* TYPELIST_HH_ */
