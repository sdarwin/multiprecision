
//          Copyright Christopher Kormanyos 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// This program computes millions of digits of pi, in fact
// up to one billion digits of pi, using Boost.Multiprecision
// combined with GMP (or MPIR).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace pi { namespace millions { namespace detail {

struct outfile_parameters
{
  static const std::size_t  number_of_digits_extra_trunc = static_cast<std::size_t>(10U);
  static const std::size_t  number_of_digits_per_column  = static_cast<std::size_t>(10U);
  static const std::size_t  number_of_columns_per_line   = static_cast<std::size_t>( 5U);
  static const std::size_t  number_of_digits_per_line    = static_cast<std::size_t>(number_of_digits_per_column * number_of_columns_per_line);
  static const std::size_t  number_of_lines_per_group    = static_cast<std::size_t>(10U);
};

// *****************************************************************************
// Function    : template<typename float_type>
//               const float_type& calculate_pi_template(const bool progress_is_printed_to_cout)
//
// Description : Compute pi using a quadratically convergent Gauss AGM,
//               in the Schoenhage variant. For a description of the algorithm,
//               see the book "Pi Unleashed", as shown in the comments below.
//               If the input parameter progress_is_printed_to_cout = true,
//               then the calculation progress will be printed to std::cout.
//
//               Book reference for "Pi Unleashed":
//               http://www.jjj.de/pibook/pibook.html
//
// *****************************************************************************
template<typename float_type>
const float_type& calculate_pi(const bool progress_is_printed_to_cout)
{
  using std::fabs;
  using std::sqrt;

  static bool is_init = false;

  static float_type val_pi;

  if(!is_init)
  {
    is_init = true;

    const std::regex rx("^[^e]+[e0+-]+([0-9]+)$");

    std::match_results<std::string::const_iterator> mr;

    std::stringstream ss;
    ss.setf(std::ios::scientific);
    ss.precision(static_cast<std::streamsize>(4));

    float_type a (1);
    float_type bB(0.5F);
    float_type s (0.5F);
    float_type t (0.375F);

    // This loop is designed for computing a maximum of a few billion
    // decimal digits of pi. The number of digits roughly doubles
    // with each iteration of the loop. After 20 iterations,
    // the precision is about 2.8 million decimal digits.
    // After 29 iterations, the precision is more than one
    // billion decimal digits.

    for(unsigned k = 1U; k < 64U; ++k)
    {
      a      += sqrt(bB);
      a      /= 2U;
      val_pi  = (a * a);
      bB      = (val_pi - t) * 2U;

      const float_type iterate_term((bB - val_pi) * (1ULL << k));

      s += iterate_term;

      // Extract the base-10 order of magnitude for a rough estimate of
      // the number of base-10 digits that have been obtained in this
      // iteration. Here, we produce a short printout of the iteration
      // term that is subsequently parsed with a regular expression
      // for extracting the base-10 order.

      // Note: We are only extracting a few digits from iterate_term.
      // So piping to stringstream is not exorbitantly costly here.
      ss << iterate_term;

      const std::uint64_t approximate_digits10_of_iteration_term =
        (std::regex_match(ss.str(), mr, rx) ? boost::lexical_cast<std::uint64_t>(mr[1U])
                                            : UINT64_C(0));

      if(progress_is_printed_to_cout)
      {
          std::cout << "Approximate base-10 digits of this iteration : "
                    << std::right
                    << std::setw(12)
                    << approximate_digits10_of_iteration_term
                    << '\n';
      }

      // Test the approximate base-10 digits of this iteration term. If we have
      // attained about half of the total desired digits with this iteration term,
      // then the calculation is finished.
      if(approximate_digits10_of_iteration_term > static_cast<std::uint64_t>((std::numeric_limits<float_type>::digits10 / 2) + 16))
      {
        break;
      }

      t = (val_pi + bB) / 4U;

      ss.str(std::string());
    }

    if(progress_is_printed_to_cout)
    {
      std::cout << "Iteration loop done, compute inverse" << '\n';
    }

    val_pi += bB;
    val_pi /= s;

    if(progress_is_printed_to_cout)
    {
      std::cout << "Pi calculation is done." << '\n';
    }
  }

  return val_pi;
}

template<typename float_type>
std::ostream& report_pi_timing(std::ostream& os, const double elapsed)
{
  return os << "============================================================" << '\n'
            << "Computed "
            << static_cast<std::uint64_t>(std::numeric_limits<float_type>::digits10 - 1)
            << " digits of pi.\n"
            << "Total computation time : "
            << std::fixed
            << std::setprecision(2)
            << elapsed
            << " seconds"
            << '\n'
            << "============================================================"
            << '\n';
}

} } } // namespace pi::millions::detail

namespace pi { namespace millions {

template<typename float_type>
bool print_pi()
{
  // Calculate the value of pi. When doing so, print the calculation
  // messages to the console. Use the clock function to obtain the
  // total time of the pi calculation.

  const std::clock_t start = std::clock();

  detail::calculate_pi<float_type>(true);

  const std::clock_t stop = std::clock();

  // Evaluate the time that was required for the pi calculation.
  const double elapsed = (static_cast<double>(stop) - static_cast<double>(start) / static_cast<double>(CLOCKS_PER_SEC)) / 1000;

  // Report the time of the pi calculation to the console.
  static_cast<void>(detail::report_pi_timing<float_type>(std::cout, elapsed));

  // Report that we are writing the output file.
  std::cout << "Writing the output file." << '\n';

  // Pipe the value of pi into a stringstream object.
  std::stringstream ss;

  const std::streamsize precision_to_print =   static_cast<std::streamsize>(std::numeric_limits<float_type>::digits10)
                                             + static_cast<std::streamsize>(detail::outfile_parameters::number_of_digits_extra_trunc);

  ss << std::setprecision(precision_to_print) << detail::calculate_pi<float_type>(false);

  // Extract the string value of pi.
  std::string str = ss.str();

  while(str.length() < static_cast<std::size_t>(precision_to_print + 1))
  {
    str.push_back(static_cast<char>('0'));
  }

  // Print pi in the following format.
  //
  // pi = 3.
  // 14159.....
  //
  // Here, the digits after the decimal point are grouped
  // in sets of digits per line, and the running digit number
  // is reported at the end of each line. The digit grouping
  // is defined with parameters listed above.

  std::string::size_type p;

  if(((p = str.find('3')) != 0U) || ((p = str.find('.')) != 1U))
  {
    return false;
  }
  else
  {
    p = 2U;
  }

  // Create a vector of output files and open them.
  std::ofstream output_file;

  // Create and open all of the output files.
  output_file.open("pi.out");

  // Verify that all of the output files are open.

  if(output_file.is_open() == false)
  {
    return false;
  }

  // Report the time of the pi calculation to the output file.
  static_cast<void>(detail::report_pi_timing<float_type>(output_file, elapsed));

  // Print the first line of pi in the first file.
  output_file << "Pi = " << str.substr(0, p) << '\n';

  // Extract the digits after the decimal point in a loop.
  // Insert spaces and newlines in an easy-to-read format.

  do
  {
    const std::size_t number_of_digits_remaining = str.length() - p;

    const std::size_t number_of_digits_in_substring = (std::min)(number_of_digits_remaining,
                                                                 detail::outfile_parameters::number_of_digits_per_column);

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

  // Close all of the output files.
  output_file.close();

  return true;
}

} } // namespace pi::millions

namespace
{
  struct my_digits_of_pi
  {
    static const unsigned digits10 = 1000000U + 1U;
  };

  typedef boost::multiprecision::number<boost::multiprecision::gmp_float<my_digits_of_pi::digits10>,
                                        boost::multiprecision::et_on>
  float_type;
}

int main()
{
  const bool print_pi_is_ok = pi::millions::print_pi<float_type>();

  static_cast<void>(print_pi_is_ok);
}
