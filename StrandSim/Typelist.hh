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


template<typename TList, template<typename > class FunctionT> struct TypelistIterator;

template<template<typename > class FunctionT> struct TypelistIterator<NullType, FunctionT>
{
    static void execute()
    {
    }
};

template<typename T, typename U, template<typename > class FunctionT> struct TypelistIterator<
        Typelist<T, U> , FunctionT>
{
    static void execute()
    {
        FunctionT<T>::execute();
        TypeListIterator<U, FunctionT>::execute();
    }
};

}

#endif /* TYPELIST_HH_ */
