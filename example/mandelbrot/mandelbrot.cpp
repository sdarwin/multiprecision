///////////////////////////////////////////////////////////////////////////////
//      Copyright Christopher Kormanyos 2015 - 2017, 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

// This example uses Boost.Multiprecision to implement
// a high-precision Mandelbrot iteration and visualization.
// Graphic file creation uses Boost.Gil-old to wrap JPEG.
// Color-strething in combination with the histogram method
// is used for creating vivid landscapes.

// TBD: The color stretching and histogram methods
// should be investigated and possibly refactored.
// At the moment, they are programmed in a non-intuitive
// way that is difficult to understand.

// TBD: At the moment, a single color scheme is used.
// It is implemented in three individual "color"
// lambda functions. This is something that would
// benefit from some kind of configuration-ability.
// Can we use user-supplied RGB color functions?

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <boost/gil/extension/io/jpeg/old.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

//#define BOOST_MANDELBROT_01_FULL
//#define BOOST_MANDELBROT_03_TOP
//#define BOOST_MANDELBROT_04_SWIRL
#define BOOST_MANDELBROT_05_SEAHORSES
//#define BOOST_MANDELBROT_06_BRANCHES
//#define BOOST_MANDELBROT_07_SEAHORSE_VALLEY
//#define BOOST_MANDELBROT_08_DEEP_DIVE_01
//#define BOOST_MANDELBROT_09_DEEP_DIVE_02

namespace my_concurrency
{
  template<typename index_type,
            typename callable_function_type>
  void parallel_for(index_type             start,
                    index_type             end,
                    callable_function_type parallel_function)
  {
    // Estimate the number of threads available.
    static const unsigned int number_of_threads_hint =
      std::thread::hardware_concurrency();

    static const unsigned int number_of_threads =
      ((number_of_threads_hint == 0U) ? 4U : number_of_threads_hint);

    // Set the size of a slice for the range functions.
    index_type n = index_type(end - start) + index_type(1);

    index_type slice =
      static_cast<index_type>(std::round(n / static_cast<double> (number_of_threads)));

    slice = (std::max)(slice, index_type(1));

    // Inner loop.
    auto launch_range =
      [&parallel_function](index_type index_lo, index_type index_hi)
      {
        for(index_type i = index_lo; i < index_hi; ++i)
        {
          parallel_function(i);
        }
      };

    // Create the thread pool and launch the jobs.
    std::vector<std::thread> pool;

    pool.reserve(number_of_threads);

    index_type i1 = start;
    index_type i2 = (std::min)(index_type(start + slice), end);

    for(index_type i = 0U; ((index_type(i + index_type(1U)) < number_of_threads) && (i1 < end)); ++i)
    {
      pool.emplace_back(launch_range, i1, i2);

      i1 = i2;

      i2 = (std::min)(index_type(i2 + slice), end);
    }

    if(i1 < end)
    {
      pool.emplace_back(launch_range, i1, end);
    }

    // Wait for the jobs to finish.
    for(std::thread& thread_in_pool : pool)
    {
      if(thread_in_pool.joinable())
      {
        thread_in_pool.join();
      }
    }
  }
} // namespace my_concurrency

// Declare a base class for the Mandelbrot configuration.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations>
class mandelbrot_config_base
{
public:
  static const std::uint_fast32_t max_iterations = MaxIterations;

  typedef NumericType mandelbrot_config_numeric_type;

  virtual ~mandelbrot_config_base() { }

  const mandelbrot_config_numeric_type& x_lo() const { return my_x_lo; }
  const mandelbrot_config_numeric_type& x_hi() const { return my_x_hi; }
  const mandelbrot_config_numeric_type& y_lo() const { return my_y_lo; }
  const mandelbrot_config_numeric_type& y_hi() const { return my_y_hi; }

  virtual int mandelbrot_fractional_resolution() const = 0;

  virtual const mandelbrot_config_numeric_type& step() const = 0;

  std::uint_fast32_t width() const
  {
    const std::uint_fast32_t non_justified_width =
      static_cast<std::uint_fast32_t>((my_x_hi - my_x_lo) / this->step());

    return non_justified_width;
  }

  std::uint_fast32_t height() const
  {
    const std::uint_fast32_t non_justified_height =
      static_cast<std::uint_fast32_t>((my_y_hi - my_y_lo) / this->step());

    return non_justified_height;
  }

protected:
  const mandelbrot_config_numeric_type my_x_lo;
  const mandelbrot_config_numeric_type my_x_hi;
  const mandelbrot_config_numeric_type my_y_lo;
  const mandelbrot_config_numeric_type my_y_hi;

  mandelbrot_config_base(const mandelbrot_config_numeric_type& xl,
                         const mandelbrot_config_numeric_type& xh,
                         const mandelbrot_config_numeric_type& yl,
                         const mandelbrot_config_numeric_type& yh)
    : my_x_lo(xl),
      my_x_hi(xh),
      my_y_lo(yl),
      my_y_hi(yh) { }

private:
  mandelbrot_config_base() : my_x_lo(),
                             my_x_hi(),
                             my_y_lo(),
                             my_y_hi() { }
};

// Make a template class that represents the Mandelbrot configuration.
// This class automatically creates sensible parameters based on
// the resolution of the fixed-point type supplied in the template
// parameter. If a custom pixel count is required, the step()
// method can be modified accordingly.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations,
         const int MandelbrotFractionalResolution>
class mandelbrot_config : public mandelbrot_config_base<NumericType, MaxIterations>
{
private:
  typedef mandelbrot_config_base<NumericType, MaxIterations> base_class_type;

public:
  static_assert(MandelbrotFractionalResolution < -1,
                "The Mandelbrot fractional resolution should be less than -1");

  mandelbrot_config(const typename base_class_type::mandelbrot_config_numeric_type& xl,
                    const typename base_class_type::mandelbrot_config_numeric_type& xh,
                    const typename base_class_type::mandelbrot_config_numeric_type& yl,
                    const typename base_class_type::mandelbrot_config_numeric_type& yh)
    : base_class_type(xl, xh, yl, yh),
      my_step()
  {
    using std::ldexp;

    my_step = ldexp(typename base_class_type::mandelbrot_config_numeric_type(1U), MandelbrotFractionalResolution);
  }

  mandelbrot_config(const std::string& str_xl,
                    const std::string& str_xh,
                    const std::string& str_yl,
                    const std::string& str_yh)
    : base_class_type(boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_xl),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_xh),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_yl),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_yh)),
      my_step()
  {
    using std::ldexp;

    my_step = ldexp(typename base_class_type::mandelbrot_config_numeric_type(1U), MandelbrotFractionalResolution);
  }

  mandelbrot_config(const char* pc_xl,
                    const char* pc_xh,
                    const char* pc_yl,
                    const char* pc_yh)
    : base_class_type(boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_xl)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_xh)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_yl)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_yh))),
      my_step()
  {
    using std::ldexp;

    my_step = ldexp(typename base_class_type::mandelbrot_config_numeric_type(1U), MandelbrotFractionalResolution);
  }

  virtual ~mandelbrot_config() { }

private:
  typename base_class_type::mandelbrot_config_numeric_type my_step;

  virtual int mandelbrot_fractional_resolution() const { return MandelbrotFractionalResolution; }

  virtual const typename base_class_type::mandelbrot_config_numeric_type& step() const { return my_step; }
};

// This class generates the rows of the mandelbrot iteration.
// The coordinates are set up according to the Mandelbrot configuration.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations>
class mandelbrot_generator
{
public:
  static const std::uint_fast32_t max_iterations = MaxIterations;

  using mandelbrot_config_type = mandelbrot_config_base<NumericType, max_iterations>;

  mandelbrot_generator(const mandelbrot_config_type& config)
    : mandelbrot_config_object   (config),
      mandelbrot_image           (config.width(), config.height()),
      mandelbrot_view            (boost::gil::rgb8_view_t()),
      mandelbrot_iteration_matrix(mandelbrot_config_object.width(),
                                  std::vector<std::uint_fast32_t>(mandelbrot_config_object.height())),
      mandelbrot_color_histogram (max_iterations + 1U, UINT32_C(0))
  {
    mandelbrot_view = boost::gil::view(mandelbrot_image);
  }

  ~mandelbrot_generator() = default;

  void generate_mandelbrot_image()
  {
    // Setup the x-axis coordinates.
    std::vector<NumericType> x_values(mandelbrot_config_object.width());
    std::vector<NumericType> y_values(mandelbrot_config_object.height());

    const NumericType xy_step(mandelbrot_config_object.step());

    // Initialize the x-axis and y-axis coordinates one time only.
    NumericType x_coord(mandelbrot_config_object.x_lo());
    NumericType y_coord(mandelbrot_config_object.y_hi());

    for(std::size_t j_col = 0U; j_col < x_values.size(); ++j_col)
    {
      x_values[j_col] = x_coord;

      x_coord += xy_step;
    }

    for(std::size_t i_row = 0U; i_row < y_values.size(); ++i_row)
    {
      y_values[i_row] = y_coord;

      y_coord -= xy_step;
    }

    std::atomic_flag mandelbrot_iteration_lock = ATOMIC_FLAG_INIT;

    std::size_t unordered_parallel_row_count = 0U;

    static const NumericType four(4U);

    my_concurrency::parallel_for
    (
      std::size_t(0U),
      y_values.size(),
      [&mandelbrot_iteration_lock, &unordered_parallel_row_count, &x_values, &y_values, this](std::size_t i_row)
      {
        while(mandelbrot_iteration_lock.test_and_set()) { ; }
        ++unordered_parallel_row_count;
        std::cout << "Calculating Mandelbrot image at row "
                  << std::setw(6)
                  << unordered_parallel_row_count
                  << " of "
                  << std::setw(6)
                  << mandelbrot_config_object.height()
                  << " total. Have patience."
                  << "\r";
        mandelbrot_iteration_lock.clear();

        for(std::size_t j_col = 0U; j_col < x_values.size(); ++j_col)
        {
          NumericType zr (0U);
          NumericType zi (0U);
          NumericType zr2(0U);
          NumericType zi2(0U);

          // Use an optimized complex-numbered multiplication scheme.
          // Thereby reduce the main work of the Mandelbrot iteration to
          // three real-valued multiplications and several real-valued
          // addition/subtraction operations.

          std::uint_fast32_t iteration_result = UINT32_C(0);

          // Perform the iteration sequence for generating the Mandelbrot set.
          // Here is the main work of the program.

          while((iteration_result < max_iterations) && ((zr2 + zi2) < four))
          {
            // Optimized complex multiply and add.
            zi *= zr;

            zi  = (zi  + zi)  + y_values[i_row];
            zr  = (zr2 - zi2) + x_values[j_col];

            zr2 = zr * zr;
            zi2 = zi * zi;

            ++iteration_result;
          }

          while(mandelbrot_iteration_lock.test_and_set()) { ; }
          mandelbrot_iteration_matrix[j_col][i_row] = iteration_result;
          ++mandelbrot_color_histogram[iteration_result];
          mandelbrot_iteration_lock.clear();
        }
      }
    );

    const std::uint_fast32_t total_pixels =
      (  static_cast<std::uint_fast32_t>(mandelbrot_config_object.width ())
       * static_cast<std::uint_fast32_t>(mandelbrot_config_object.height()));

    // Perform color-stretching using the histogram approach.
    // Convert the histogram entries such that a given entry contains
    // the sum of its own entries plus all previous entries. This provides
    // a set of scale factors for the color. The histogram approach
    // automatically scales to the distribution of pixels in the image.

    const std::uint_fast32_t mandelbrot_sum =
      std::accumulate(mandelbrot_color_histogram.begin(),
                      mandelbrot_color_histogram.end(),
                      std::uint_fast32_t(0U),
      [&total_pixels](std::uint_fast32_t& sum, std::uint_fast32_t& histogram_entry) -> std::uint_fast32_t
      {
        sum += histogram_entry;

        const double sum_div_total_pixels =
          static_cast<double>(sum) / static_cast<double>(total_pixels);

        const double histogram_scale = std::pow(sum_div_total_pixels, 1.2);

        std::uint_fast32_t scaled_histogram_value =
          static_cast<std::uint_fast32_t>(histogram_scale * static_cast<double>(0xFFU));

        if(scaled_histogram_value < 0xFFU)
        {
          scaled_histogram_value =
            static_cast<std::uint_fast32_t>(std::pow(double(scaled_histogram_value), 1.00));
        }

        histogram_entry = UINT32_C(0xFF) - scaled_histogram_value;

        return sum;
      });

    static_cast<void>(mandelbrot_sum);

    for(std::uint_fast32_t i_row = UINT32_C(0); i_row < mandelbrot_config_object.height(); ++i_row)
    {
      for(std::uint_fast32_t j_col = UINT32_C(0); j_col < mandelbrot_config_object.width(); ++j_col)
      {
        const std::uint_fast32_t color = mandelbrot_color_histogram[mandelbrot_iteration_matrix[j_col][i_row]];

        const std::array<std::uint_fast32_t (*)(const std::uint_fast32_t&), 3U> color_functions =
        {
          [](const std::uint_fast32_t& the_color) -> std::uint_fast32_t
          {
            const double color_phase = (double(the_color) / 255.0) * (3.1415926535897932385 * 8.0);

            const double my_color_red = (std::sin(color_phase) / 2.0) + 0.5;

            return static_cast<std::uint_fast32_t>(my_color_red * 255.0);
          },

          [](const std::uint_fast32_t& the_color) -> std::uint_fast32_t
          {
            return the_color;
          },

          [](const std::uint_fast32_t& the_color) -> std::uint_fast32_t
          {
            return (the_color * the_color) / 255;
          }
        };

        const std::uint_fast32_t color_r = ((color <= 4U) ? color : color_functions[0U](color));
        const std::uint_fast32_t color_g = ((color <= 4U) ? color : color_functions[1U](color));
        const std::uint_fast32_t color_b = ((color <= 4U) ? color : color_functions[2U](color));

        // Mix the color supplied in the template hue parameters.
        const std::uint8_t rh = static_cast<std::uint8_t>((255U * color_r) / UINT32_C(255));
        const std::uint8_t gh = static_cast<std::uint8_t>((255U * color_g) / UINT32_C(255));
        const std::uint8_t bh = static_cast<std::uint8_t>((255U * color_b) / UINT32_C(255));

        const boost::gil::rgb8_pixel_t the_color  = boost::gil::rgb8_pixel_t(rh, gh, bh);

        mandelbrot_view(j_col, i_row) = boost::gil::rgb8_pixel_t(the_color);
      }
    }

    boost::gil::jpeg_write_view("mandelbrot.jpg", mandelbrot_view);

    std::cout << std::endl
              << "The ouptput file mandelbrot.jpg has been written"
              << std::endl;
  }

private:
  const mandelbrot_config_type&                mandelbrot_config_object;

  boost::gil::rgb8_image_t                     mandelbrot_image;
  boost::gil::rgb8_view_t                      mandelbrot_view;

  std::vector<std::vector<std::uint_fast32_t>> mandelbrot_iteration_matrix;
  std::vector<std::uint_fast32_t>              mandelbrot_color_histogram;
};

int main()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>,
                                        boost::multiprecision::et_off>
  numeric_type;

  #if defined BOOST_MANDELBROT_01_FULL

    // This is the classic full immage from (-2.0, -1.0) ... (05, 1.0).
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(2000), -10>;

    const mandelbrot_config_type mandelbrot_config_object(-2.000L, +0.500L,
                                                          -1.000L, +1.000L);

  #elif defined BOOST_MANDELBROT_03_TOP

    // This is a view of an upper part of the image (near the top of the classic full view).
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(1000), -13>;

    const mandelbrot_config_type mandelbrot_config_object(-0.130L - 0.282L, -0.130L + 0.282L,
                                                          +0.856L - 0.282L, +0.856L + 0.282L);

  #elif defined BOOST_MANDELBROT_04_SWIRL

    // This is a fanning swirl image.
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(2000), -23>;

    const mandelbrot_config_type mandelbrot_config_object(-0.749730L - 0.0002315L, -0.749730L + 0.0002315L,
                                                          -0.046608L - 0.0002315L, -0.046608L + 0.0002315L);

  #elif defined BOOST_MANDELBROT_05_SEAHORSES

    // This is a swirly seahorse image.
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(2000), -50>;

    const mandelbrot_config_type
      mandelbrot_config_object(-0.7453983606668L - 1.72E-12L, -0.7453983606668L + 1.72E-12L,
                               +0.1125046349960L - 1.72E-12L, +0.1125046349960L + 1.72E-12L);

  #elif defined BOOST_MANDELBROT_06_BRANCHES

    // This is a spiral image of branches.
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(2000), -49>;

    const mandelbrot_config_type mandelbrot_config_object(+0.3369844464873L - 4.2E-12L, +0.3369844464873L + 4.2E-12L,
                                                          +0.0487782196791L - 4.2E-12L, +0.0487782196791L + 4.2E-12L);

  #elif defined BOOST_MANDELBROT_07_SEAHORSE_VALLEY

    // This is an image from the seahorse valley.
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(1000), -16>;

    const mandelbrot_config_type
      mandelbrot_config_object("-0.748", "-0.700",
                               "+0.222", "+0.270");

  #elif defined BOOST_MANDELBROT_08_DEEP_DIVE_01

    // This is a deep zoom image.
    // Note: Use 128 decimal digits for this iteration.

    static_assert(std::numeric_limits<numeric_type>::digits10 >= 128,
                  "Error: Please use 128 or more decimal digits for deep dive 01.");

    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(800), -366>;

    const numeric_type delta("+1.16E-107");
    const numeric_type cx   ("-1.9999999991382701187582747629086949883168091366368209595068022727154702772791898403544767055386190962248152412805947511E+00");
    const numeric_type cy   ("+1.3148954435076375751362475668065050021517005209120957095294493435305489940275245944710958864319980774657032331030784899E-14");

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #elif defined BOOST_MANDELBROT_09_DEEP_DIVE_02

    // This is a deep zoom image.
    // Note: Use 80 decimal digits for this iteration.

    static_assert(std::numeric_limits<numeric_type>::digits10 >= 80,
                  "Error: Please use 80 or more decimal digits for deep dive 02.");

    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(8000), -192>;

    const numeric_type delta("+1.78E-55");
    const numeric_type cx   (numeric_type("-1.295189082147777457017064177185681926706566460884888469217455"));
    const numeric_type cy   (numeric_type("+0.440936982678320138880903678356262612113214627431396203682665"));

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #else

    #error: Mandelbrot image type is not defined!

  #endif

  using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

  using mandelbrot_generator_type = mandelbrot_generator<mandelbrot_numeric_type, mandelbrot_config_type::max_iterations>;

  const std::clock_t start = std::clock();

  mandelbrot_generator_type mandelbrot_generator(mandelbrot_config_object);

  mandelbrot_generator.generate_mandelbrot_image();

  const float elapsed = (float(std::clock()) - float(start)) / float(CLOCKS_PER_SEC);

  std::cout << "Time for calculation: "
            << elapsed
            << "s"
            << std::endl;
}
