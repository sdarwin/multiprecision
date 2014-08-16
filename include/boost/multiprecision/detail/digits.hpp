///////////////////////////////////////////////////////////////
//  Copyright 2012 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_

#ifndef BOOST_MP_DIGITS_HPP
#define BOOST_MP_DIGITS_HPP

namespace boost{ namespace multiprecision{ namespace detail{

inline unsigned long digits10_2_2(unsigned long d10)
{
   const unsigned long long d10_uLL = static_cast<unsigned long long>(d10);

   return static_cast<unsigned long>((d10_uLL * 1000uLL) / 301u + ((d10_uLL * 1000uLL) % 301u ? 2u : 1u));
}

inline unsigned long digits2_2_10(unsigned long d2)
{
   const unsigned long long d2_uLL = static_cast<unsigned long long>(d2);

   return static_cast<unsigned long>((d2_uLL * 301uLL) / 1000u);
}

}}} // namespaces

#endif
