
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
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/gmp.hpp>
//#include <boost/multiprecision/mpfr.hpp>

namespace pi { namespace millions { namespace detail {

struct outfile_parameters
{
  static const std::size_t  number_of_digits_extra_trunc = 10U;
  static const std::size_t  number_of_digits_per_word    = 10U;
  static const std::size_t  number_of_words_per_line     = 5U;
  static const std::size_t  number_of_digits_per_line    = (number_of_digits_per_word * number_of_words_per_line);
  static const std::size_t  number_of_lines_per_group    = 10U;
  static const std::int64_t number_of_digits_per_file    = 2000000000LL;
};

template<typename float_type>
bool get_output_files(std::vector<std::ofstream>& output_files)
{
  std::size_t number_of_files =   static_cast<std::size_t>( ((std::numeric_limits<float_type>::digits10 - 1) + outfile_parameters::number_of_digits_extra_trunc) / outfile_parameters::number_of_digits_per_file)
                                + static_cast<std::size_t>((((std::numeric_limits<float_type>::digits10 - 1) + outfile_parameters::number_of_digits_extra_trunc) % outfile_parameters::number_of_digits_per_file) != 0);

  output_files.resize(number_of_files);

  unsigned file_name_number = 0U;

  // Create and open all of the output files.
  std::for_each(output_files.begin(),
                output_files.end(),
                [&file_name_number](std::ofstream& of)
                {
                  std::stringstream ss;

                  ss << std::setw(6) << std::setfill('0') << file_name_number;

                  of.open(std::string(("pi_" + ss.str()) + ".out").c_str());

                  ++file_name_number;
                });

  // Verify that all of the output files are open.
  const bool all_output_files_are_open = std::all_of(output_files.begin(),
                                                     output_files.end(),
                                                     [](const std::ofstream& of) -> bool
                                                     {
                                                       return of.is_open();
                                                     });

  return all_output_files_are_open;
}

// *****************************************************************************
// Function    : template<typename float_type>
//               const float_type& calculate_pi_template(const bool b_trace)
//
// Description : Compute pi using a quadratically convergent Gauss AGM,
//               in the Schoenhage variant. For a description of the algorithm,
//               see the book "Pi Unleashed", as shown in the comments below.
//               If the input b_trace = true, then the calculation progress
//               will be output to std::cout.
//
//               Book reference for "Pi Unleashed":
//               http://www.jjj.de/pibook/pibook.html
//
// *****************************************************************************
template<typename float_type>
const float_type& calculate_pi(const bool b_trace)
{
  using std::fabs;
  using std::sqrt;

  static bool is_init = false;

  static float_type val_pi;

  if(!is_init)
  {
    is_init = true;

    static const std::regex rx("^[^e]+[e0+-]+([0-9]+)$");
    std::match_results<std::string::const_iterator> mr;
    std::stringstream ss;

    float_type a (1);
    float_type bB(0.5F);
    float_type s (0.5F);
    float_type t (0.375F);

    // This loop is designed for a maximum of a few billion
    // decimal digits of pi. The index k should reach no higher
    // than about 25 or 30. After 20 iterations, the precision
    // is about 1.4 million decimal digits.

    for(unsigned k = 1U; k < 64U; ++k)
    {
      a      += sqrt(bB);
      a      /= 2U;
      val_pi  = (a * a);
      bB      = (val_pi - t) * 2U;

      const float_type iterate_term((bB - val_pi) * (1ULL << k));

      s += iterate_term;

      ss.str(std::string());

      ss << std::scientific << std::setprecision(8) << iterate_term;

      const std::string str(ss.str());

      std::uint64_t approximate_digits_of_pi = (std::regex_match(str, mr, rx) ? boost::lexical_cast<std::uint64_t>(mr[1U])
                                                                              : UINT64_C(0));

      if(b_trace)
      {
        // Extract the base-10 order of magnitude for a rough estimate of
        // the digits in this iteration of the calculation. Here, we produce
        // a short printout of the iteration term that is subsequently
        // parsed with a regular expression for extracting the base-10 order.

          std::cout << "Approximate digits of pi : "
                    << std::right
                    << std::setw(12)
                    << approximate_digits_of_pi
                    << '\n';
      }

      // Test the significant digits of the last iteration change.
      // If it is large enough, then the calculation is finished.
      if(approximate_digits_of_pi > ((std::numeric_limits<float_type>::digits10 / 2) + 16))
      {
        break;
      }

      t = (val_pi + bB) / 4U;
    }

    if(b_trace) { std::cout << "Iteration loop done, compute inverse" << '\n'; }

    val_pi += bB;
    val_pi /= s;

    if(b_trace) { std::cout << "Pi calculation is done." << '\n'; }
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
  std::vector<std::ofstream> output_files;

  const bool all_output_files_are_open = detail::get_output_files<float_type>(output_files);

  if(all_output_files_are_open == false)
  {
    return false;
  }

  std::vector<std::ofstream>::size_type of_index = 0U;

  // Report the time of the pi calculation to the first output file.
  static_cast<void>(detail::report_pi_timing<float_type>(output_files[of_index], elapsed));

  // Print the first line of pi in the first file.
  output_files[of_index] << "Pi = " << str.substr(0, p) << '\n';

  // Extract the digits after the decimal point in a loop.
  // Insert spaces and newlines in an easy-to-read format.

  do
  {
    const std::size_t number_of_digits_remaining = str.length() - p;

    const std::size_t number_of_digits_in_substring = (std::min)(number_of_digits_remaining,
                                                                 detail::outfile_parameters::number_of_digits_per_word);

    output_files[of_index] << str.substr(p, number_of_digits_in_substring) << " ";

    p += number_of_digits_in_substring;

    const std::string::size_type p2 = p - 2U;

    if((p2 % detail::outfile_parameters::number_of_digits_per_line) == 0U)
    {
      // A single line has ended. Print the running digit count
      // and a newline.
      output_files[of_index] << " : " << p2 << '\n';

      if((p2 % (detail::outfile_parameters::number_of_lines_per_group * detail::outfile_parameters::number_of_digits_per_line)) == 0U)
      {
        // A group of lines has ended.

        const bool the_current_file_is_full = ((p2 % detail::outfile_parameters::number_of_digits_per_file) == 0U);

        if(the_current_file_is_full)
        {
          // The current file is full. Switch to a new file.
          ++of_index;
        }
        else if(p >= (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc))
        {
          // The output is finished. Do nothing and break from the loop below.
          ;
        }
        else
        {
          // The group of lines is full, but there is still space in the
          // current file. Simply print a standalone newline character.
          output_files[of_index] << '\n';
        }
      }
    }
  }
  while(p < (str.length() - detail::outfile_parameters::number_of_digits_extra_trunc));

  // Close all of the output files.
  std::for_each(output_files.begin(),
                output_files.end(),
                [](std::ofstream& of)
                {
                  of.close();
                });

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
  volatile bool print_pi_is_ok = pi::millions::print_pi<float_type>();

  static_cast<void>(print_pi_is_ok);
}
