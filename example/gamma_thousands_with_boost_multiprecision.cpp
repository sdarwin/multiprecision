
//      Copyright Christopher Kormanyos 2014 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/cstdfloat.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace gamma { namespace thousands {

namespace detail
{
  struct outfile_parameters
  {
    static const std::size_t number_of_digits_extra_trunc = 20U;
    static const std::size_t number_of_digits_per_word    = 10U;
    static const std::size_t number_of_words_per_line     = 5U;
    static const std::size_t number_of_digits_per_line    = (number_of_digits_per_word * number_of_words_per_line);
    static const std::size_t number_of_lines_per_group    = 10U;
  };

  template<typename float_type>
  std::ostream& report_tgamma_timing(std::ostream& os, const boost::float_least64_t elapsed)
  {
    return os << "============================================================" << '\n'
              << "Computed "
              << static_cast<std::uint64_t>(std::numeric_limits<float_type>::digits10 - 1)
              << " digits of tgamma(1/2).\n"
              << "Total computation time : "
              << std::fixed
              << std::setprecision(2)
              << elapsed
              << " seconds"
              << '\n'
              << "============================================================"
              << '\n';
  }
} //  namespace detail

template<typename float_type>
bool print_tgamma()
{
  // Calculate the value of tgamma(1/2). Use the clock() function
  // to obtain the total time of the tgamma(1/2) calculation.

  std::cout << "Computing tgamma(1/2)..." << '\n';

  const std::clock_t start = std::clock();
  const float_type tgamma_of_one_half = boost::math::tgamma(float_type(1) / 2);
  const std::clock_t stop = std::clock();

  // Evaluate the time that was required for the tgamma(1/2) calculation.
  const boost::float_least64_t elapsed = (  static_cast<boost::float_least64_t>(stop)
                                          - static_cast<boost::float_least64_t>(start)
                                         )
                                         /
                                         static_cast<boost::float_least64_t>(CLOCKS_PER_SEC);

  // Report the time of the tgamma(1/2) calculation to the console.
  static_cast<void>(detail::report_tgamma_timing<float_type>(std::cout, elapsed));

  using std::fabs;
  using std::sqrt;
  const float_type delta = fabs(tgamma_of_one_half - sqrt(boost::math::constants::pi<float_type>()));

  static const float_type tolerance = std::numeric_limits<float_type>::epsilon() * 1.0E+06F;

  const bool tgamma_result_is_ok = (delta < tolerance);

  std::cout << "Check tgamma tolerance : "
            << std::setprecision(4)
            << std::scientific
            << delta
            << ". The result seems to be : "
            << (tgamma_result_is_ok ? "OK." : "NOT OK!")
            << '\n';

  // Report that we are writing the output file.
  std::cout << "Writing the output file.\n";

  // Pipe the value of tgamma(1/2) into a stringstream object.
  std::stringstream ss;

  const std::streamsize precision_to_print =   static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10)
                                             + static_cast<std::streamsize>(detail::outfile_parameters::number_of_digits_extra_trunc);

  ss << std::setprecision(precision_to_print) << tgamma_of_one_half;

  // Extract the string value of gamma(1/2).
  std::string str = ss.str();

  while(str.length() < static_cast<std::string::size_type>(precision_to_print + 1))
  {
    str.push_back(static_cast<char>('0'));
  }

  // Print tgamma(1/2) in the following format.
  //
  // tgamma(1/2) = 1.
  // 7724538509 0551602729 8167483341 1451827975 4945612238  : 50
  // 7128213807 7898529112 8459103218 1374950656 7385446654  : 100
  // ...
  //
  // Here, the digits after the decimal point are grouped
  // in sets of digits per line. The running digit number
  // is reported at the end of each line. The digit grouping
  // is defined with parameters listed above.

  bool result_is_ok = false;

  const std::string::size_type pos_find_one = str.find('1');
  const std::string::size_type pos_find_dot = str.find('.');

  if(  (pos_find_one == static_cast<std::string::size_type>(0U))
     &&(pos_find_dot == static_cast<std::string::size_type>(1U)))
  {
    std::string::size_type pos = 2U;

    // Create the output file and open it.
    std::ofstream output_file("tgamma.out");

    if(output_file.is_open())
    {
      // Report the time of the tgamma(1/2) calculation to the output file.
      static_cast<void>(detail::report_tgamma_timing<float_type>(output_file, elapsed));

      // Print the first line of tgamma(1/2) in the first file.
      output_file << "tgamma(1/2) = " << str.substr(0, pos) << '\n';

      // Extract the digits after the decimal point in a loop.
      // Insert spaces and newlines in an easy-to-read format.

      do
      {
        const std::size_t number_of_digits_remaining = str.length() - pos;

        const std::size_t number_of_digits_in_substring = (std::min)(number_of_digits_remaining,
                                                                     detail::outfile_parameters::number_of_digits_per_word);

        output_file << str.substr(pos, number_of_digits_in_substring) << " ";

        pos += number_of_digits_in_substring;

        const std::string::size_type pos2 = pos - 2U;

        if((pos2 % detail::outfile_parameters::number_of_digits_per_line) == 0U)
        {
          // A single line has ended. Print the running digit count
          // and a newline.
          output_file << " : " << pos2 << '\n';

          if((pos2 % (detail::outfile_parameters::number_of_lines_per_group * detail::outfile_parameters::number_of_digits_per_line)) == 0U)
          {
            // A group of lines has ended.

            if(pos >= (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc))
            {
              // The output is finished. Do nothing and break from the loop condition below.
              ;
            }
            else
            {
              // The group of lines is full, but there is still space in the
              // current file. Simply print a standalone newline character.
              output_file << '\n';
            }
          }
        }
      }
      while(pos < (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc));

      // Close the output file.
      output_file.close();

      result_is_ok = true;
    }
  }

  return result_is_ok;
}

} } // namespace gamma::thousands

namespace
{
  struct my_digits_of_gamma
  {
    static const std::uint_least32_t digits10 = UINT32_C(10000) + UINT32_C(1);
  };

  typedef boost::multiprecision::number<boost::multiprecision::gmp_float<my_digits_of_gamma::digits10>,
                                        boost::multiprecision::et_off>
  float_type;
}

int main()
{
  volatile bool print_tgamma_is_ok = gamma::thousands::print_tgamma<float_type>();

  static_cast<void>(print_tgamma_is_ok);
}
