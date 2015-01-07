
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

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace gamma { namespace thousands { namespace detail {

struct outfile_parameters
{
  static const std::size_t  number_of_digits_extra_trunc = 20U;
  static const std::size_t  number_of_digits_per_word    = 10U;
  static const std::size_t  number_of_words_per_line     = 5U;
  static const std::size_t  number_of_digits_per_line    = (number_of_digits_per_word * number_of_words_per_line);
  static const std::size_t  number_of_lines_per_group    = 10U;
};

template<typename float_type>
std::ostream& report_gamma_timing(std::ostream& os, const double elapsed)
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

} } } // namespace gamma::thousands::detail

namespace gamma { namespace thousands {

template<typename float_type>
bool print_gamma()
{
  // Calculate the value of tgamma(1/2). Use the clock function to
  // obtain the total time of the tgamma(1/2) calculation.

  std::cout << "Computing tgamma(1/2)..." << '\n';

  const std::clock_t start = std::clock();

  const float_type gamma_of_one_half = boost::math::tgamma(float_type(1) / 2);

  const std::clock_t stop = std::clock();

  // Evaluate the time that was required for the tgamma(1/2) calculation.
  const double elapsed = (static_cast<double>(stop) - static_cast<double>(start) / static_cast<double>(CLOCKS_PER_SEC)) / 1000;

  // Report the time of the tgamma(1/2) calculation to the console.
  static_cast<void>(detail::report_gamma_timing<float_type>(std::cout, elapsed));

  using std::fabs;
  using std::sqrt;
  const float_type delta_check_of_tgamma_result = fabs(gamma_of_one_half - sqrt(boost::math::constants::pi<float_type>()));

  static const float_type tolerance = UINT32_C(1000000) * std::numeric_limits<float_type>::epsilon();

  const bool verification_of_tgamma_result_is_ok = (delta_check_of_tgamma_result < tolerance);

  std::cout << "Check tgamma tolerance: "
            << std::setprecision(4)
            << std::scientific
            << delta_check_of_tgamma_result
            << ". This seems to be "
            << (verification_of_tgamma_result_is_ok ? "OK." : "NOT OK!")
            << '\n';

  // Report that we are writing the output file.
  std::cout << "Writing the output file." << '\n';

  // Pipe the value of tgamma(1/2) into a stringstream object.
  std::stringstream ss;

  const std::streamsize precision_to_print =   static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10)
                                             + static_cast<std::streamsize>(detail::outfile_parameters::number_of_digits_extra_trunc);

  ss << std::setprecision(precision_to_print) << gamma_of_one_half;

  // Extract the string value of gamma(1/2).
  std::string str = ss.str();
  std::string::size_type p;

  while(str.length() < static_cast<std::size_t>(precision_to_print + 1))
  {
    str.push_back(static_cast<char>('0'));
  }

  // Print tgamma(1/2) in the following format.
  //
  // tgamma(1/2) = 1.
  // 7724538.....
  //
  // Here, the digits after the decimal point are grouped
  // in sets digits per line, and the running digit number
  // is reported at the end of each line. The digit grouping
  // is defined with parameters listed above.

  if(((p = str.find('1')) != 0U) || ((p = str.find('.')) != 1U))
  {
    return false;
  }
  else
  {
    p = 2U;
  }

  // Create the output file and open it.
  std::ofstream output_file("tgamma.out");

  const bool the_output_file_is_open = output_file.is_open();

  if(the_output_file_is_open == false)
  {
    return false;
  }

  // Report the time of the tgamma(1/2) calculation to the output file.
  static_cast<void>(detail::report_gamma_timing<float_type>(output_file, elapsed));

  // Print the first line of tgamma(1/2) in the first file.
  output_file << "tgamma(1/2) = " << str.substr(0, p) << '\n';

  // Extract the digits after the decimal point in a loop.
  // Insert spaces and newlines in an easy-to-read format.

  do
  {
    const std::size_t number_of_digits_remaining = str.length() - p;

    const std::size_t number_of_digits_in_substring = (std::min)(number_of_digits_remaining,
                                                                 detail::outfile_parameters::number_of_digits_per_word);

    output_file << str.substr(p, number_of_digits_in_substring) << " ";

    p += number_of_digits_in_substring;

    const std::string::size_type p2 = p - 2U;

    if((p2 % detail::outfile_parameters::number_of_digits_per_line) == 0U)
    {
      // A single line has ended. Print the running digit count
      // and a newline.
      output_file << " : " << p2 << '\n';

      if((p2 % (detail::outfile_parameters::number_of_lines_per_group * detail::outfile_parameters::number_of_digits_per_line)) == 0U)
      {
        // A group of lines has ended.

        if(p >= (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc))
        {
          // The output is finished. Do nothing and break from the loop below.
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
  while(p < (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc));

  // Close the output file.
  output_file.close();

  return true;
}

} } // namespace gamma::thousands

namespace
{
  struct my_digits_of_gamma
  {
    static const unsigned digits10 = 10000U + 1U;
  };

  typedef boost::multiprecision::number<boost::multiprecision::gmp_float<my_digits_of_gamma::digits10>,
                                        boost::multiprecision::et_off>
  float_type;

//  typedef float float_type;
}

int main()
{
  volatile bool print_gamma_is_ok = gamma::thousands::print_gamma<float_type>();

  static_cast<void>(print_gamma_is_ok);
}
