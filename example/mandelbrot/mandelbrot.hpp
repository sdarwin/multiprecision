///////////////////////////////////////////////////////////////////////////////
//      Copyright Christopher Kormanyos 2015 - 2017, 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MULTIPRECISION_MANDELBROT_2020_03_12_HPP_
#define BOOST_MULTIPRECISION_MANDELBROT_2020_03_12_HPP_

// This example uses Boost.Multiprecision to implement
// a high-precision Mandelbrot iteration and visualization.
// Graphic file creation uses Boost.Gil-old to wrap JPEG.
// Color-strething in combination with the histogram method
// is used for creating vivid landscapes.

// TBD: The color stretching and histogram methods
// should be investigated and possibly refactored.
// At the moment, they are programmed in a less intuitive
// way that might be difficult to understand.

// TBD: At the moment, a single color scheme is used.
// It is implemented in three individual "color"
// lambda functions. This is something that would
// benefit from some kind of configuration-ability.
// Can we use user-supplied RGB color functions?

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <string>
#include <thread>
#include <vector>

#include <boost/gil/extension/io/jpeg/old.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/lexical_cast.hpp>

namespace boost { namespace multiprecision { namespace mandelbrot {

namespace detail {

namespace my_concurrency {
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

} // namespace detail

// Declare a base class for the Mandelbrot configuration.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations>
class mandelbrot_config_base
{
public:
  static const std::uint_fast32_t max_iterations = MaxIterations;

  using mandelbrot_config_numeric_type = NumericType;

  virtual ~mandelbrot_config_base() { }

  const mandelbrot_config_numeric_type& x_lo() const { return my_x_lo; }
  const mandelbrot_config_numeric_type& x_hi() const { return my_x_hi; }
  const mandelbrot_config_numeric_type& y_lo() const { return my_y_lo; }
  const mandelbrot_config_numeric_type& y_hi() const { return my_y_hi; }

  virtual int mandelbrot_fractional_resolution() const = 0;

  virtual const mandelbrot_config_numeric_type& step() const = 0;

  std::uint_fast32_t integral_width() const
  {
    const std::uint_fast32_t non_justified_width =
      static_cast<std::uint_fast32_t>(my_width / this->step());

    return non_justified_width;
  }

  std::uint_fast32_t integral_height() const
  {
    const std::uint_fast32_t non_justified_height =
      static_cast<std::uint_fast32_t>(my_height / this->step());

    return non_justified_height;
  }

protected:
  const mandelbrot_config_numeric_type my_x_lo;
  const mandelbrot_config_numeric_type my_x_hi;
  const mandelbrot_config_numeric_type my_y_lo;
  const mandelbrot_config_numeric_type my_y_hi;
  const mandelbrot_config_numeric_type my_width;
  const mandelbrot_config_numeric_type my_height;

  mandelbrot_config_base(const mandelbrot_config_numeric_type& xl,
                         const mandelbrot_config_numeric_type& xh,
                         const mandelbrot_config_numeric_type& yl,
                         const mandelbrot_config_numeric_type& yh)
    : my_x_lo  (xl),
      my_x_hi  (xh),
      my_y_lo  (yl),
      my_y_hi  (yh),
      my_width (my_x_hi - my_x_lo),
      my_height(my_y_hi - my_y_lo){ }

private:
  mandelbrot_config_base() = default;
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
  using base_class_type = mandelbrot_config_base<NumericType, MaxIterations>;

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
      mandelbrot_image           (config.integral_width(), config.integral_height()),
      mandelbrot_view            (boost::gil::rgb8_view_t()),
      mandelbrot_iteration_matrix(config.integral_width(),
                                  std::vector<std::uint_fast32_t>(config.integral_height())),
      mandelbrot_color_histogram (max_iterations + 1U, UINT32_C(0))
  {
    mandelbrot_view = boost::gil::view(mandelbrot_image);
  }

  ~mandelbrot_generator() = default;

  void generate_mandelbrot_image(const char* pstr_filename, std::ostream& os)
  {
    // Setup the x-axis and y-axis coordinates.
    // TBD: Consider potential program redesign with unequal x/y ranges.

    const NumericType xy_step(mandelbrot_config_object.step());

    std::vector<NumericType> x_values(mandelbrot_config_object.integral_width());
    std::vector<NumericType> y_values(mandelbrot_config_object.integral_height());

    NumericType x_coord(mandelbrot_config_object.x_lo());
    NumericType y_coord(mandelbrot_config_object.y_hi());

    for(auto& x : x_values) { x = x_coord; x_coord += xy_step; }
    for(auto& y : y_values) { y = y_coord; y_coord -= xy_step; }

    std::atomic_flag mandelbrot_iteration_lock = ATOMIC_FLAG_INIT;

    std::size_t unordered_parallel_row_count = 0U;

    static const NumericType four(4U);

    detail::my_concurrency::parallel_for
    (
      std::size_t(0U),
      y_values.size(),
      [&mandelbrot_iteration_lock, &unordered_parallel_row_count, &os, &x_values, &y_values, this](std::size_t j_row)
      {
        while(mandelbrot_iteration_lock.test_and_set()) { ; }
        ++unordered_parallel_row_count;
        os << "Calculating Mandelbrot image at row "
           << unordered_parallel_row_count
           << " of "
           << y_values.size()
           << " total: "
           << std::fixed
           << std::setprecision(1)
           << (100.0F * float(unordered_parallel_row_count)) / float(y_values.size())
           << "%. Have patience."
           << "\r";
        mandelbrot_iteration_lock.clear();

        for(std::size_t i_col = 0U; i_col < x_values.size(); ++i_col)
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

            zi  = (zi  + zi)  + y_values[j_row];
            zr  = (zr2 - zi2) + x_values[i_col];

            zr2 = zr; zr2 *= zr;
            zi2 = zi; zi2 *= zi;

            ++iteration_result;
          }

          while(mandelbrot_iteration_lock.test_and_set()) { ; }
          mandelbrot_iteration_matrix[i_col][j_row] = iteration_result;
          ++mandelbrot_color_histogram[iteration_result];
          mandelbrot_iteration_lock.clear();
        }
      }
    );

    const std::uint_fast32_t total_pixels = static_cast<std::uint_fast32_t>(x_values.size() * y_values.size());

    // Perform color-stretching using the histogram approach.
    // Convert the histogram entries such that a given entry contains
    // the sum of its own entries plus all previous entries. This provides
    // a set of scale factors for the color. The histogram approach
    // automatically scales to the distribution of pixels in the image.

    os << std::endl;
    os << "Perform color-stretching using the histogram approach." << std::endl;

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

    for(std::uint_fast32_t j_row = UINT32_C(0); j_row < y_values.size(); ++j_row)
    {
      for(std::uint_fast32_t i_col = UINT32_C(0); i_col < x_values.size(); ++i_col)
      {
        const std::uint_fast32_t color = mandelbrot_color_histogram[mandelbrot_iteration_matrix[i_col][j_row]];

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

        mandelbrot_view(i_col, j_row) = boost::gil::rgb8_pixel_t(the_color);
      }
    }

    const std::string str_filename(pstr_filename);

    boost::gil::jpeg_write_view(str_filename, mandelbrot_view);

    os << std::endl
       << "The ouptput file "
       << str_filename
       << " has been written"
       << std::endl;
  }

private:
  const mandelbrot_config_type&                mandelbrot_config_object;

  boost::gil::rgb8_image_t                     mandelbrot_image;
  boost::gil::rgb8_view_t                      mandelbrot_view;

  std::vector<std::vector<std::uint_fast32_t>> mandelbrot_iteration_matrix;
  std::vector<std::uint_fast32_t>              mandelbrot_color_histogram;
};

} } } // namespace boost::multiprecision::mandelbrot

#endif // BOOST_MULTIPRECISION_MANDELBROT_2020_03_12_HPP_
