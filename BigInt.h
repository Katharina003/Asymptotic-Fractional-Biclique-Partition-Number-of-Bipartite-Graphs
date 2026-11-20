#pragma once

#include <boost/multiprecision/cpp_int.hpp>


/*
Define the alias BigInt for multiprecision integers. 
BigInts can get arbritary large and are only bounded by memory limits, which makes
them perfect for our use case.
*/
using BigInt = boost::multiprecision::cpp_int;