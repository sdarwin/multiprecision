///////////////////////////////////////////////////////////////////////////////
//      Copyright Christopher Kormanyos 2015 - 2017, 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//

// A well-known example of fractal is the Mandelbrot set,
// which is based upon the function z_{n+1} = z_{n}^2 + c.
// A common way of coloring Mandelbrot images is by taking
// the number of iterations required to reach non-bounded
// divergence from c and then assigning that value a color.

// This example uses Boost.Multiprecision to implement
// a high-precision Mandelbrot iteration and visualization.
// Graphic file creation uses Boost.Gil (old) to wrap JPEG.
// Color-strething in combination with the histogram method
// is used for creating vivid images. The default color
// scheme uses stretched, amplified and modulated black
// and white coloring. The Mandelbrot iteration is carried
// out with hardware concurrency with multiple threads.
// The multithreading dispatcher used 3/4 of the available
// CPU cores that can be found using hardware concurrency.

// The Mandelbrot set consists of those points c in the
// complex plane for which the iteration
//   z_{n+1} = z_{n}^2 + c,
// with z_{0} = 0 stays bounded.

// Interesting iteration points could be points having an orbit.
// An orbit of length n is a sequence of z_{n} with
//   z_{1} = c, z_{2}, ..., z{n},
// such that z_{n} = z_{1} and z_{n} != z_{k} with (1 < k < n).
// In order to find these, numerical methods are needed.
// The equation z_{n} = z_{1} can only be solved in closed form
// by hand for small n. A point c of order n will also show up
// as a point of order n_{m}, for some m > 1. Mark these points
// in your set.

// Any point that is inside the Mandelbrot set and close to the
// boundary between the set and its complement as well as any point
// outside the Mandelbrot set that is close to this boundary is an
// interesting point. The closer you are to the boundary, the more
// you need to zoom in to see the interesting parts. In particular,
// all points on the x-axis between -2 and 1/4 are in the Mandelbrot set.
// Especially close to x = -2 (from the right), the point (x, 0)
// is arbitrarily close to the boundary. So try the point (eps - 2, 0)
// for a small (eps > 0). Some Mandelbrot softwares use a strategy that
// zooms in, continually trying to find a point close to the boundary
// while zooming, and uses that as the zoom point.

// Points on the boundary of the Mandelbrot set can potentially
// have the most interesting orbits.  The easiest boundary points
// to compute are:
//   * The spike along the negative real axis
//   * The boundary of the main cardioid:
//       r = (1 - cos(theta))/2,
//       x = r*cos(theta)+0.25,
//       y = r*sin(theta)
//   * The boundary of the period 2 disk:
//       r = 0.25,
//       x = r*cos(theta)-1,
//       y = r*sin(theta)

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <thread>
#include <vector>

#include <boost/gil/extension/io/jpeg/old.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/lexical_cast.hpp>

#define BOOST_MULTIPRECISION_MANDELBROT_01_FULL                  1
#define BOOST_MULTIPRECISION_MANDELBROT_03_TOP                   3
#define BOOST_MULTIPRECISION_MANDELBROT_04_SWIRL                 4
#define BOOST_MULTIPRECISION_MANDELBROT_05_SEAHORSES             5
#define BOOST_MULTIPRECISION_MANDELBROT_06_BRANCHES              6
#define BOOST_MULTIPRECISION_MANDELBROT_07_SEAHORSE_VALLEY       7
#define BOOST_MULTIPRECISION_MANDELBROT_08_DEEP_DIVE_01          8
#define BOOST_MULTIPRECISION_MANDELBROT_09_DEEP_DIVE_02          9
#define BOOST_MULTIPRECISION_MANDELBROT_10_ZOOM_WIKI_01         10
#define BOOST_MULTIPRECISION_MANDELBROT_11_ZOOM_WIKI_02         11
#define BOOST_MULTIPRECISION_MANDELBROT_12_ZOOM_WIKI_03         12
#define BOOST_MULTIPRECISION_MANDELBROT_20_ZOOM_VERY_DEEP_00    20
#define BOOST_MULTIPRECISION_MANDELBROT_21_ZOOM_VERY_DEEP_01    21
#define BOOST_MULTIPRECISION_MANDELBROT_22_ZOOM_VERY_DEEP_02    22
#define BOOST_MULTIPRECISION_MANDELBROT_30_ZOOM_ANOTHER_00      30
#define BOOST_MULTIPRECISION_MANDELBROT_31_ZOOM_ANOTHER_01      31
#define BOOST_MULTIPRECISION_MANDELBROT_32_ZOOM_ANOTHER_02      32

#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_CPP_DEC_FLOAT       101
#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_GMP_FLOAT           102

#if !defined(BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX)
#define BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX BOOST_MULTIPRECISION_MANDELBROT_05_SEAHORSES
#endif

#if !defined(BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE)
#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_CPP_DEC_FLOAT
//#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_GMP_FLOAT
#endif

#if  (BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE == BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_CPP_DEC_FLOAT)
#include <boost/multiprecision/cpp_dec_float.hpp>
#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME cpp_dec_float
#elif(BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE == BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_GMP_FLOAT)
#include <boost/multiprecision/gmp.hpp>
#define BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME gmp_float
#else
#error BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_TYPE is undefined.
#endif

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

  static const unsigned int number_of_threads_total =
    ((number_of_threads_hint == 0U) ? 4U : number_of_threads_hint);

  // Use only 3/4 of the available cores.
  static const unsigned int number_of_threads = number_of_threads_total - (number_of_threads_total / 8U);

  // Set the size of a slice for the range functions.
  index_type n = index_type(end - start) + index_type(1);

  index_type slice =
    static_cast<index_type>(std::round(n / static_cast<double>(number_of_threads)));

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

class color_functions_base
{
public:
  virtual ~color_functions_base() = default;

  virtual std::uint_fast32_t color_function_r(const std::uint_fast32_t&) const = 0;
  virtual std::uint_fast32_t color_function_g(const std::uint_fast32_t&) const = 0;
  virtual std::uint_fast32_t color_function_b(const std::uint_fast32_t&) const = 0;

protected:
  color_functions_base() = default;

  static std::uint_fast32_t color_phaser_01(const std::uint_fast32_t& c)
  {
    const double color_phase = (double(c) / 255.0) * (3.1415926535897932385 * 8.0);

    const double my_color = (std::sin(color_phase) / 2.0) + 0.5;

    return static_cast<std::uint_fast32_t>(my_color * 255.0);
  }
};

class color_functions_bw final : public color_functions_base
{
public:
  color_functions_bw() = default;

  virtual ~color_functions_bw() = default;

private:
  virtual std::uint_fast32_t color_function_r(const std::uint_fast32_t& c) const { return color_phaser_01(c); }
  virtual std::uint_fast32_t color_function_g(const std::uint_fast32_t& c) const { return color_phaser_01(c); }
  virtual std::uint_fast32_t color_function_b(const std::uint_fast32_t& c) const { return color_phaser_01(c); }
};

class color_functions_pretty final : public color_functions_base
{
public:
  color_functions_pretty() = default;

  virtual ~color_functions_pretty() = default;

private:
  virtual std::uint_fast32_t color_function_r(const std::uint_fast32_t& c) const
  {
    return color_phaser_01(c);
  }

  virtual std::uint_fast32_t color_function_g(const std::uint_fast32_t& c) const
  {
    return c;
  }

  virtual std::uint_fast32_t color_function_b(const std::uint_fast32_t& c) const
  {
    return static_cast<std::uint_fast32_t>((double(c) * double(c)) / 255.0);
  }
};

class color_stretches_base
{
public:
  virtual ~color_stretches_base() = default;

  void init(const std::uint_fast32_t total_pixels)
  {
    my_total_pixels = total_pixels;
    my_sum          = 0U;
  }

  virtual void color_stretch(std::uint_fast32_t&) = 0;

protected:
  std::uint_fast32_t my_total_pixels;
  std::uint_fast32_t my_sum;

  color_stretches_base() : my_total_pixels(0U),
                           my_sum         (0U) { }

  color_stretches_base(const color_stretches_base& other)
    : my_total_pixels(other.my_total_pixels),
      my_sum         (other.my_sum) { }

  color_stretches_base& operator=(const color_stretches_base& other)
  {
    my_total_pixels = other.my_total_pixels;
    my_sum          = other.my_sum;

    return *this;
  }
};

class color_stretches_default final : public color_stretches_base
{
public:
  color_stretches_default() = default;

  virtual ~color_stretches_default() = default;

  color_stretches_default& operator=(const color_stretches_default& other)
  {
    if(this != &other)
    {
      static_cast<void>(color_stretches_base::operator=(other));
    }

    return *this;
  }

  virtual void color_stretch(std::uint_fast32_t& histogram_entry)
  {
    // Perform color stretching using the histogram approach.
    // Convert the histogram entries such that a given entry contains
    // the sum of its own entries plus all previous entries. This provides
    // a set of scale factors for the color. The histogram approach
    // automatically scales to the distribution of pixels in the image.

    my_sum += histogram_entry;

    const double sum_div_total_pixels =
      static_cast<double>(my_sum) / static_cast<double>(my_total_pixels);

    const double histogram_scale = std::pow(sum_div_total_pixels, 1.2);

    const std::uint_fast32_t scaled_histogram_value =
      static_cast<std::uint_fast32_t>(histogram_scale * static_cast<double>(0xFFU));

    histogram_entry = UINT32_C(0xFF) - scaled_histogram_value;
  }
};

} // namespace detail

// Declare a base class for the Mandelbrot configuration.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations>
class mandelbrot_config_base
{
public:
  static const std::uint_fast32_t max_iterations = MaxIterations;

  using mandelbrot_config_numeric_type = NumericType;

  virtual ~mandelbrot_config_base() = default;

  const mandelbrot_config_numeric_type& x_lo() const { return my_x_lo; }
  const mandelbrot_config_numeric_type& x_hi() const { return my_x_hi; }
  const mandelbrot_config_numeric_type& y_lo() const { return my_y_lo; }
  const mandelbrot_config_numeric_type& y_hi() const { return my_y_hi; }

  virtual const mandelbrot_config_numeric_type& step() const = 0;

  std::uint_fast32_t integral_width() const
  {
    const std::uint_fast32_t non_rounded_width2 = static_cast<std::uint_fast32_t>((my_width * 2) / this->step());

    return static_cast<std::uint_fast32_t>(non_rounded_width2 + 1U) / 2U;
  }

  std::uint_fast32_t integral_height() const
  {
    const std::uint_fast32_t non_rounded_height2 = static_cast<std::uint_fast32_t>((my_height * 2) / this->step());

    return static_cast<std::uint_fast32_t>(non_rounded_height2 + 1U) / 2U;
  }

protected:
  const mandelbrot_config_numeric_type my_x_lo;
  const mandelbrot_config_numeric_type my_x_hi;
  const mandelbrot_config_numeric_type my_y_lo;
  const mandelbrot_config_numeric_type my_y_hi;
  const mandelbrot_config_numeric_type my_width;
  const mandelbrot_config_numeric_type my_height;

  mandelbrot_config_base(const mandelbrot_config_base& other)
    : my_x_lo  (other.my_x_lo),
      my_x_hi  (other.my_x_hi),
      my_y_lo  (other.my_y_lo),
      my_y_hi  (other.my_y_hi),
      my_width (other.my_width),
      my_height(other.my_height) { }

  mandelbrot_config_base(const mandelbrot_config_numeric_type& xl,
                         const mandelbrot_config_numeric_type& xh,
                         const mandelbrot_config_numeric_type& yl,
                         const mandelbrot_config_numeric_type& yh)
    : my_x_lo  (xl),
      my_x_hi  (xh),
      my_y_lo  (yl),
      my_y_hi  (yh),
      my_width (my_x_hi - my_x_lo),
      my_height(my_y_hi - my_y_lo) { }

  mandelbrot_config_base& operator=(const mandelbrot_config_base& other)
  {
    my_x_lo   = other.my_x_lo;
    my_x_hi   = other.my_x_hi;
    my_y_lo   = other.my_y_lo;
    my_y_hi   = other.my_y_hi;
    my_width  = other.my_width;
    my_height = other.my_height;

    return *this;
  }

private:
  mandelbrot_config_base() = delete;
};

// Make a template class that represents the Mandelbrot configuration.
// This class automatically creates sensible parameters based on
// the resolution of the fixed-point type supplied in the template
// parameter. If a custom pixel count is required, the step()
// method can be modified accordingly.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations,
         const std::uint_fast32_t PixelCountX>
class mandelbrot_config final : public mandelbrot_config_base<NumericType, MaxIterations>
{
private:
  using base_class_type = mandelbrot_config_base<NumericType, MaxIterations>;

public:
  mandelbrot_config(const typename base_class_type::mandelbrot_config_numeric_type& xl,
                    const typename base_class_type::mandelbrot_config_numeric_type& xh,
                    const typename base_class_type::mandelbrot_config_numeric_type& yl,
                    const typename base_class_type::mandelbrot_config_numeric_type& yh)
    : base_class_type(xl, xh, yl, yh),
      my_step(base_class_type::my_width / PixelCountX) { }

  mandelbrot_config(const std::string& str_xl,
                    const std::string& str_xh,
                    const std::string& str_yl,
                    const std::string& str_yh)
    : base_class_type(boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_xl),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_xh),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_yl),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(str_yh)),
      my_step(base_class_type::my_width / PixelCountX) { }

  mandelbrot_config(const char* pc_xl,
                    const char* pc_xh,
                    const char* pc_yl,
                    const char* pc_yh)
    : base_class_type(boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_xl)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_xh)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_yl)),
                      boost::lexical_cast<typename base_class_type::mandelbrot_config_numeric_type>(std::string(pc_yh))),
      my_step(base_class_type::my_width / PixelCountX) { }

  mandelbrot_config(const mandelbrot_config& other)
    : base_class_type(other),
      my_step(other.my_step) { }

  virtual ~mandelbrot_config() = default;

  mandelbrot_config& operator=(const mandelbrot_config& other)
  {
    if(this != other)
    {
      static_cast<void>(base_class_type::operator=(other));

      my_step = other.my_step;
    }

    return *this;
  }

private:
  typename base_class_type::mandelbrot_config_numeric_type my_step;

  virtual const typename base_class_type::mandelbrot_config_numeric_type& step() const { return my_step; }
};

// This class generates the rows of the mandelbrot iteration.
// The coordinates are set up according to the Mandelbrot configuration.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations>
class mandelbrot_generator final
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

  void generate_mandelbrot_image(const std::string&                  str_filename,
                                 const detail::color_functions_base& color_functions = detail::color_functions_bw(),
                                       detail::color_stretches_base& color_stretches = detail::color_stretches_default(),
                                       std::ostream&                 output_stream   = std::cout)
  {
    // Setup the x-axis and y-axis coordinates.

    std::vector<NumericType> x_values(mandelbrot_config_object.integral_width());
    std::vector<NumericType> y_values(mandelbrot_config_object.integral_height());

    {
      const NumericType local_step(mandelbrot_config_object.step());

      NumericType x_coord(mandelbrot_config_object.x_lo());
      NumericType y_coord(mandelbrot_config_object.y_hi());

      for(auto& x : x_values) { x = x_coord; x_coord += local_step; }
      for(auto& y : y_values) { y = y_coord; y_coord -= local_step; }
    }

    static const NumericType four(4U);

    std::atomic_flag mandelbrot_iteration_lock = ATOMIC_FLAG_INIT;

    std::size_t unordered_parallel_row_count = 0U;

    detail::my_concurrency::parallel_for
    (
      std::size_t(0U),
      y_values.size(),
      [&mandelbrot_iteration_lock, &unordered_parallel_row_count, &output_stream, &x_values, &y_values, this](std::size_t j_row)
      {
        while(mandelbrot_iteration_lock.test_and_set()) { ; }
        ++unordered_parallel_row_count;
        output_stream << "Calculating Mandelbrot image at row "
                      << unordered_parallel_row_count
                      << " of "
                      << y_values.size()
                      << ". Total processed so far: "
                      << std::fixed
                      << std::setprecision(1)
                      << (100.0 * double(unordered_parallel_row_count)) / double(y_values.size())
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

    output_stream << std::endl;

    output_stream << "Perform color stretching." << std::endl;
    apply_color_stretches(x_values, y_values, color_stretches);

    output_stream << "Apply color functions." << std::endl;
    apply_color_functions(x_values, y_values, color_functions);

    output_stream << "Write output JPEG file " << str_filename << "." << std::endl;
    boost::gil::jpeg_write_view(str_filename, mandelbrot_view);
  }

private:
  const mandelbrot_config_type&                mandelbrot_config_object;

  boost::gil::rgb8_image_t                     mandelbrot_image;
  boost::gil::rgb8_view_t                      mandelbrot_view;

  std::vector<std::vector<std::uint_fast32_t>> mandelbrot_iteration_matrix;
  std::vector<std::uint_fast32_t>              mandelbrot_color_histogram;

  void apply_color_stretches(const std::vector<NumericType>& x_values,
                             const std::vector<NumericType>& y_values,
                             detail::color_stretches_base& color_stretches)
  {
    color_stretches.init(static_cast<std::uint_fast32_t>(x_values.size() * y_values.size()));

    for(auto& histogram_entry : mandelbrot_color_histogram)
    {
      color_stretches.color_stretch(histogram_entry);
    }
  }

  void apply_color_functions(const std::vector<NumericType>& x_values,
                             const std::vector<NumericType>& y_values,
                             const detail::color_functions_base& color_functions)
  {
    for(std::uint_fast32_t j_row = UINT32_C(0); j_row < y_values.size(); ++j_row)
    {
      for(std::uint_fast32_t i_col = UINT32_C(0); i_col < x_values.size(); ++i_col)
      {
        const std::uint_fast32_t color = mandelbrot_color_histogram[mandelbrot_iteration_matrix[i_col][j_row]];

        // Get the three hue values.
        const std::uint_fast32_t color_r = ((color <= 4U) ? color : color_functions.color_function_r(color));
        const std::uint_fast32_t color_g = ((color <= 4U) ? color : color_functions.color_function_g(color));
        const std::uint_fast32_t color_b = ((color <= 4U) ? color : color_functions.color_function_b(color));

        // Mix the color from the hue values.
        const std::uint8_t rh = static_cast<std::uint8_t>((255U * color_r) / UINT32_C(255));
        const std::uint8_t gh = static_cast<std::uint8_t>((255U * color_g) / UINT32_C(255));
        const std::uint8_t bh = static_cast<std::uint8_t>((255U * color_b) / UINT32_C(255));

        const boost::gil::rgb8_pixel_t the_color  = boost::gil::rgb8_pixel_t(rh, gh, bh);

        mandelbrot_view(i_col, j_row) = boost::gil::rgb8_pixel_t(the_color);
      }
    }
  }
};

} } } // namespace boost::multiprecision::mandelbrot

int main()
{
  #if (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_01_FULL)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_01_FULL") + ".jpg";

    // This is the classic full immage.
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 4096U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;


    const mandelbrot_numeric_type delta_half(+1.25L);
    const mandelbrot_numeric_type cx        (-0.75L);
    const mandelbrot_numeric_type cy        (+0.0L);

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_03_TOP)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_03_TOP") + ".jpg";

    // This is a view of an upper part of the image (near the top of the classic full view).
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("0.282");
    const mandelbrot_numeric_type cx        ("-0.130");
    const mandelbrot_numeric_type cy        ("0.856");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_04_SWIRL)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_04_SWIRL") + ".jpg";

    // This is a fanning swirl image.
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("0.0002315");
    const mandelbrot_numeric_type cx        ("-0.749730");
    const mandelbrot_numeric_type cy        ("-0.046608");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_05_SEAHORSES)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_05_SEAHORSES") + ".jpg";

    // This is a swirly seahorse image.
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("1.76E-12");
    const mandelbrot_numeric_type cx        ("-0.7453983606667815");
    const mandelbrot_numeric_type cy        ("0.1125046349959942");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_06_BRANCHES)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_06_BRANCHES") + ".jpg";

    // This is a spiral image of branches.
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("4.2E-12");
    const mandelbrot_numeric_type cx        ("0.3369844464873");
    const mandelbrot_numeric_type cy        ("0.0487782196791");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_07_SEAHORSE_VALLEY)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_07_SEAHORSE_VALLEY") + ".jpg";

    // This is an image from the seahorse valley.
    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<31>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("0.024");
    const mandelbrot_numeric_type cx        ("-0.748");
    const mandelbrot_numeric_type cy        ("0.222");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_08_DEEP_DIVE_01)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_08_DEEP_DIVE_01") + ".jpg";

    // This is a deep zoom image.
    // Note: Use 143 or more decimal digits for this iteration.

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<143>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 2000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("1.25E-107");
    const mandelbrot_numeric_type cx        ("-1.99999999913827011875827476290869498831680913663682095950680227271547027727918984035447670553861909622481524124");
    const mandelbrot_numeric_type cy        ("0.00000000000001314895443507637575136247566806505002151700520912095709529449343530548994027524594471095886432006");

    static_assert(std::numeric_limits<mandelbrot_numeric_type>::digits10 >= 143,
                  "Error: Please use 143 or more decimal digits for BOOST_MULTIPRECISION_MANDELBROT_08_DEEP_DIVE_01.");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX == BOOST_MULTIPRECISION_MANDELBROT_09_DEEP_DIVE_02)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_09_DEEP_DIVE_02") + ".jpg";

    // This is a deep zoom image.
    // Note: Use 79 or more decimal digits for this iteration.

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<95>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 15000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("2.55E-55");
    const mandelbrot_numeric_type cx        ("-1.295189082147777457017064177185681926706566460884888469217456");
    const mandelbrot_numeric_type cy        ("0.440936982678320138880903678356262612113214627431396203682661");

    static_assert(std::numeric_limits<mandelbrot_numeric_type>::digits10 >= 95,
                  "Error: Please use 95 or more decimal digits for BOOST_MULTIPRECISION_MANDELBROT_09_DEEP_DIVE_02.");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_10_ZOOM_WIKI_00)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_10_ZOOM_WIKI_00") + ".jpg";

    // This is a medium zoom image from the zoom coordinates of:
    // http://en.wikipedia.org/wiki/File:Mandelbrot_sequence_new.gif
    // Note: Use 55 or more decimal digits for this iteration.

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<55>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 20000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("3.1E-25");
    const mandelbrot_numeric_type cx        ("-0.743643887037158704752191506114774");
    const mandelbrot_numeric_type cy        ("0.131825904205311970493132056385139");

    static_assert(std::numeric_limits<mandelbrot_numeric_type>::digits10 >= 55,
                  "Error: Please use 55 or more decimal digits for BOOST_MULTIPRECISION_MANDELBROT_10_ZOOM_WIKI_00.");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_11_ZOOM_WIKI_01)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_11_ZOOM_WIKI_01") + ".jpg";

    // This is a medium zoom image from the zoom coordinates of:
    // http://en.wikipedia.org/wiki/File:Mandelbrot_sequence_new.gif
    // Note: Use 55 or more decimal digits for this iteration.

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<55>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 30000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("3.3E-27");
    const mandelbrot_numeric_type cx        ("-0.743643887037158704752191506114774");
    const mandelbrot_numeric_type cy        ("0.131825904205311970493132056385139");

    static_assert(std::numeric_limits<mandelbrot_numeric_type>::digits10 >= 55,
                  "Error: Please use 55 or more decimal digits for BOOST_MULTIPRECISION_MANDELBROT_11_ZOOM_WIKI_01.");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_12_ZOOM_WIKI_02)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_12_ZOOM_WIKI_02") + ".jpg";

    // This is a medium zoom image from the zoom coordinates of:
    // http://en.wikipedia.org/wiki/File:Mandelbrot_sequence_new.gif
    // Note: Use 55 or more decimal digits for this iteration.

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<55>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 40000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("6.5E-28");
    const mandelbrot_numeric_type cx        ("-0.743643887037158704752191506114774");
    const mandelbrot_numeric_type cy        ("0.131825904205311970493132056385139");

    static_assert(std::numeric_limits<mandelbrot_numeric_type>::digits10 >= 55,
                  "Error: Please use 55 or more decimal digits for BOOST_MULTIPRECISION_MANDELBROT_12_ZOOM_WIKI_02.");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_20_ZOOM_VERY_DEEP_00)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_20_ZOOM_VERY_DEEP_00") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<365>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 60000U, 1536U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    // At the moment, this is my personal best deep dive with relatively
    // high precision. This deep dive features magnification of
    // approximately 1e310, which is 1 and 310 zeros. Generating this
    // image took about 15 hours on my machine.

    // TBD: Try to zoom even deeper with higher precision and more iterations.
    // The remarkable video here: https://www.youtube.com/watch?v=pCpLWbHVNhk
    // makes use of the same high-precision coordinate (shown below).
    // The author of that video, for instance, reports and beautifully
    // depicts achieving a deep zoom of approximately 1e1091 at this
    // staggeringly high-precision coordinate.

    const mandelbrot_numeric_type delta_half("4.4E-311");
    const mandelbrot_numeric_type cx        ("0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928");
    const mandelbrot_numeric_type cy        ("-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_21_ZOOM_VERY_DEEP_01)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_21_ZOOM_VERY_DEEP_01") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<567>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 160000U, 512U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("1.1E-501");
    const mandelbrot_numeric_type cx        ("0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928");
    const mandelbrot_numeric_type cy        ("-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_22_ZOOM_VERY_DEEP_02)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_22_ZOOM_VERY_DEEP_02") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<1048>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 700000U, 128U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("1.1E-990");
    const mandelbrot_numeric_type cx        ("0.360240443437614363236125244449545308482607807958585750488375814740195346059218100311752936722773426396233731729724987737320035372683285317664532401218521579554288661726564324134702299962817029213329980895208036363104546639698106204384566555001322985619004717862781192694046362748742863016467354574422779443226982622356594130430232458472420816652623492974891730419252651127672782407292315574480207005828774566475024380960675386215814315654794021855269375824443853463117354448779647099224311848192893972572398662626725254769950976527431277402440752868498588785436705371093442460696090720654908973712759963732914849861213100695402602927267843779747314419332179148608587129105289166676461292845685734536033692577618496925170576714796693411776794742904333484665301628662532967079174729170714156810530598764525260869731233845987202037712637770582084286587072766838497865108477149114659838883818795374195150936369987302574377608649625020864292915913378927790344097552591919409137354459097560040374880346637533711271919419723135538377394364882968994646845930838049998854075817859391340445151448381853615103761584177161812057928");
    const mandelbrot_numeric_type cy        ("-0.6413130610648031748603750151793020665794949522823052595561775430644485741727536902556370230689681162370740565537072149790106973211105273740851993394803287437606238596262287731075999483940467161288840614581091294325709988992269165007394305732683208318834672366947550710920088501655704252385244481168836426277052232593412981472237968353661477793530336607247738951625817755401065045362273039788332245567345061665756708689359294516668271440525273653083717877701237756144214394870245598590883973716531691124286669552803640414068523325276808909040317617092683826521501539932397262012011082098721944643118695001226048977430038509470101715555439047884752058334804891389685530946112621573416582482926221804767466258346014417934356149837352092608891639072745930639364693513216719114523328990690069588676087923656657656023794484324797546024248328156586471662631008741349069961493817600100133439721557969263221185095951241491408756751582471307537382827924073746760884081704887902040036056611401378785952452105099242499241003208013460878442953408648178692353788153787229940221611731034405203519945313911627314900851851072122990492499999999999999999991");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_30_ZOOM_ANOTHER_00)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_30_ZOOM_ANOTHER_00") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<287>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 60000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("9.5E-191");
    const mandelbrot_numeric_type cx        ("-1.740062382579339905220844167065825638296641720436171866879862418461182919644153056054840718339483225743450008259172138785492983677893366503417299549623738838303346465461290768441055486136870719850559269507357211790243666940134793753068611574745943820712885258222629105433648695946003865");
    const mandelbrot_numeric_type cy        ("0.0281753397792110489924115211443195096875390767429906085704013095958801743240920186385400814658560553615695084486774077000669037710191665338060418999324320867147028768983704831316527873719459264592084600433150333362859318102017032958074799966721030307082150171994798478089798638258639934");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_31_ZOOM_ANOTHER_01)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_31_ZOOM_ANOTHER_01") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<287>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 60000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("6.4E-201");
    const mandelbrot_numeric_type cx        ("-1.740062382579339905220844167065825638296641720436171866879862418461182919644153056054840718339483225743450008259172138785492983677893366503417299549623738838303346465461290768441055486136870719850559269507357211790243666940134793753068611574745943820712885258222629105433648695946003865");
    const mandelbrot_numeric_type cy        ("0.0281753397792110489924115211443195096875390767429906085704013095958801743240920186385400814658560553615695084486774077000669037710191665338060418999324320867147028768983704831316527873719459264592084600433150333362859318102017032958074799966721030307082150171994798478089798638258639934");

  #elif (BOOST_MULTIPRECISION_MANDELBROT_IMAGE_INDEX ==  BOOST_MULTIPRECISION_MANDELBROT_32_ZOOM_ANOTHER_02)

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MULTIPRECISION_MANDELBROT_32_ZOOM_ANOTHER_02") + ".jpg";

    using local_numeric_type      = boost::multiprecision::number<boost::multiprecision::BOOST_MULTIPRECISION_MANDELBROT_NUMBER_BACKEND_NAME<287>, boost::multiprecision::et_off>;
    using mandelbrot_config_type  = boost::multiprecision::mandelbrot::mandelbrot_config<local_numeric_type, 60000U, 2048U>;
    using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

    const mandelbrot_numeric_type delta_half("4.9E-211");
    const mandelbrot_numeric_type cx        ("-1.740062382579339905220844167065825638296641720436171866879862418461182919644153056054840718339483225743450008259172138785492983677893366503417299549623738838303346465461290768441055486136870719850559269507357211790243666940134793753068611574745943820712885258222629105433648695946003865");
    const mandelbrot_numeric_type cy        ("0.0281753397792110489924115211443195096875390767429906085704013095958801743240920186385400814658560553615695084486774077000669037710191665338060418999324320867147028768983704831316527873719459264592084600433150333362859318102017032958074799966721030307082150171994798478089798638258639934");

  #else

    #error: Mandelbrot image type is not defined!

  #endif

  const mandelbrot_config_type mandelbrot_config_object(cx - delta_half, cx + delta_half,
                                                        cy - delta_half, cy + delta_half);

  using mandelbrot_generator_type =
    boost::multiprecision::mandelbrot::mandelbrot_generator<mandelbrot_numeric_type,
                                                            mandelbrot_config_type::max_iterations>;

        boost::multiprecision::mandelbrot::detail::color_stretches_default local_color_stretches;
  const boost::multiprecision::mandelbrot::detail::color_functions_bw      local_color_functions;

  mandelbrot_generator_type mandelbrot_generator(mandelbrot_config_object);

  const std::clock_t start = std::clock();

  mandelbrot_generator.generate_mandelbrot_image(str_filename,
                                                 local_color_functions,
                                                 local_color_stretches);

  const double elapsed = (double(std::clock()) - double(start)) / double(CLOCKS_PER_SEC);

  std::cout << "Time for calculation: "
            << elapsed
            << "s"
            << std::endl;
}
