
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

// Declare a base class for the Mandelbrot configuration.
template<typename NumericType>
class mandelbrot_config_base
{
public:
  typedef NumericType mandelbrot_config_numeric_type;

  virtual ~mandelbrot_config_base() { }

  const mandelbrot_config_numeric_type& x_lo() const { return my_x_lo; }
  const mandelbrot_config_numeric_type& x_hi() const { return my_x_hi; }
  const mandelbrot_config_numeric_type& y_lo() const { return my_y_lo; }
  const mandelbrot_config_numeric_type& y_hi() const { return my_y_hi; }

  virtual std::uint_fast32_t max_iterations() const = 0;

  virtual int mandelbrot_fractional_resolution() const = 0;

  virtual const mandelbrot_config_numeric_type& step() const = 0;

  virtual std::uint_fast32_t width () const = 0;
  virtual std::uint_fast32_t height() const = 0;

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
  mandelbrot_config_base()
    : my_x_lo(mandelbrot_config_numeric_type()),
      my_x_hi(mandelbrot_config_numeric_type()),
      my_y_lo(mandelbrot_config_numeric_type()),
      my_y_hi(mandelbrot_config_numeric_type()) { }
};

// Make a template class that represents the Mandelbrot configuration.
// This class automatically creates sensible parameters based on
// the resolution of the fixed-point type supplied in the template
// parameter. If a custom pixel count is required, the step()
// method can be modified accordingly.
template<typename NumericType,
         const std::uint_fast32_t MaxIterations,
         const int MandelbrotFractionalResolution>
class mandelbrot_config : public mandelbrot_config_base<NumericType>
{
private:
  typedef mandelbrot_config_base<NumericType> base_class_type;

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

    my_step = ldexp(NumericType(1U), MandelbrotFractionalResolution);
  }

  virtual ~mandelbrot_config() { }

private:
  NumericType my_step;

  virtual std::uint_fast32_t max_iterations() const { return MaxIterations; }

  virtual int mandelbrot_fractional_resolution() const { return MandelbrotFractionalResolution; }

  virtual const NumericType& step() const { return my_step; }

  virtual std::uint_fast32_t width() const
  {
    const std::uint_fast32_t non_justified_width =
      static_cast<std::uint_fast32_t>((this->x_hi() - this->x_lo()) / step());

    const std::uint_fast32_t justified_width =
      non_justified_width + (4U - (non_justified_width % 4));

    return justified_width;
  }

  virtual std::uint_fast32_t height() const
  {
    const std::uint_fast32_t non_justified_height =
      static_cast<std::uint_fast32_t>((this->y_hi() - this->y_lo()) / step());

    const std::uint_fast32_t justified_height =
      non_justified_height + (4U - (non_justified_height % 4));

    return justified_height;
  }
};

class mandelbrot_inner_loop_object_base
{
public:
  virtual void loop() = 0;

  mandelbrot_inner_loop_object_base(const mandelbrot_inner_loop_object_base&) = default;

  virtual ~mandelbrot_inner_loop_object_base() = default;

  mandelbrot_inner_loop_object_base& operator=(const mandelbrot_inner_loop_object_base&) = default;

  std::uint_fast32_t get_iteration_result() const { return iteration_result; }

protected:
  const std::uint_fast32_t maximum_iterations;
  std::uint_fast32_t iteration_result;

  mandelbrot_inner_loop_object_base(const std::uint_fast32_t max_iter)
    : maximum_iterations(max_iter),
      iteration_result  (0U) { }

private:
  mandelbrot_inner_loop_object_base() = delete;
};

template<typename NumericType>
class mandelbrot_inner_loop_object final : public mandelbrot_inner_loop_object_base
{
public:
  mandelbrot_inner_loop_object(const std::uint_fast32_t max_iter,
                               const NumericType& x,
                               const NumericType& y)
    : mandelbrot_inner_loop_object_base(max_iter),
      cr (x),
      ci (y),
      zr (0U),
      zi (0U),
      zr2(0U),
      zi2(0U) { }

  mandelbrot_inner_loop_object(const mandelbrot_inner_loop_object&) = default;

  virtual ~mandelbrot_inner_loop_object() = default;

  mandelbrot_inner_loop_object& operator=(const mandelbrot_inner_loop_object&) = default;

  virtual void loop()
  {
    // Use an optimized complex-numbered multiplication scheme.
    // Thereby reduce the main work of the Mandelbrot iteration to
    // three real-valued multiplications and several real-valued
    // addition/subtraction operations.

    iteration_result = UINT32_C(0);

    // Perform the iteration sequence for generating the Mandelbrot set.
    // Here is the main work of the program.

    while((iteration_result < maximum_iterations) && ((zr2 + zi2) < 4))
    {
      // Optimized complex multiply and add.
      zi *= zr;

      zi  = (zi  + zi)  + ci;
      zr  = (zr2 - zi2) + cr;

      zr2 = zr * zr;
      zi2 = zi * zi;

      ++iteration_result;
    }
  }

private:
  const NumericType cr;
  const NumericType ci;

  NumericType zr;
  NumericType zi;

  NumericType zr2;
  NumericType zi2;

  mandelbrot_inner_loop_object() = delete;
};

namespace local
{
  std::vector<std::uint_fast32_t> mandelbrot_thread_result_vector_00;
  std::vector<std::uint_fast32_t> mandelbrot_thread_result_vector_01;
  std::vector<std::uint_fast32_t> mandelbrot_thread_result_vector_02;
  std::vector<std::uint_fast32_t> mandelbrot_thread_result_vector_03;

  std::atomic_bool mandelbrot_thread_exit_flag_00;
  std::atomic_bool mandelbrot_thread_exit_flag_01;
  std::atomic_bool mandelbrot_thread_exit_flag_02;
  std::atomic_bool mandelbrot_thread_exit_flag_03;

  std::atomic_bool mandelbrot_thread_do_loop_00;
  std::atomic_bool mandelbrot_thread_do_loop_01;
  std::atomic_bool mandelbrot_thread_do_loop_02;
  std::atomic_bool mandelbrot_thread_do_loop_03;

  template<typename NumericInputIteratorType>
  static void thread_routine_00(NumericInputIteratorType                                            x_first,
                                NumericInputIteratorType                                            x_last,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_val_00,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_step_00,
                                std::uint_fast32_t                                                  max_iter)
  {
    typedef typename std::iterator_traits<NumericInputIteratorType>::value_type local_value_type;

    const std::ptrdiff_t x_values_size = std::distance(x_first, x_last);

    mandelbrot_thread_result_vector_00.resize(std::size_t(x_values_size));

    while(std::atomic_load(&mandelbrot_thread_exit_flag_00) == false)
    {
      if(std::atomic_load(&mandelbrot_thread_do_loop_00) == true)
      {
        std::fill(mandelbrot_thread_result_vector_00.begin(),
                  mandelbrot_thread_result_vector_00.end(),
                  std::uint_fast32_t(0U));

        for(std::size_t col = 0U; col < std::size_t(x_values_size); ++col)
        {
          mandelbrot_inner_loop_object<local_value_type>
            my_mandelbrot_inner_loop_object_00(max_iter,
                                               *(x_first + col),
                                               y_val_00);

          my_mandelbrot_inner_loop_object_00.loop();

          mandelbrot_thread_result_vector_00[col] = my_mandelbrot_inner_loop_object_00.get_iteration_result();
        }

        y_val_00 -= 4 * y_step_00;

        std::atomic_store(&mandelbrot_thread_do_loop_00, false);
      }
    }
  }

  template<typename NumericInputIteratorType>
  static void thread_routine_01(NumericInputIteratorType                                            x_first,
                                NumericInputIteratorType                                            x_last,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_val_01,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_step_01,
                                std::uint_fast32_t                                                  max_iter)
  {
    typedef typename std::iterator_traits<NumericInputIteratorType>::value_type local_value_type;

    const std::ptrdiff_t x_values_size = std::distance(x_first, x_last);

    mandelbrot_thread_result_vector_01.resize(std::size_t(x_values_size));

    while(std::atomic_load(&mandelbrot_thread_exit_flag_01) == false)
    {
      if(std::atomic_load(&mandelbrot_thread_do_loop_01) == true)
      {
        std::fill(mandelbrot_thread_result_vector_01.begin(),
                  mandelbrot_thread_result_vector_01.end(),
                  std::uint_fast32_t(0U));

        for(std::size_t col = 0U; col < std::size_t(x_values_size); ++col)
        {
          mandelbrot_inner_loop_object<local_value_type>
            my_mandelbrot_inner_loop_object_01(max_iter,
                                               *(x_first + col),
                                               y_val_01);

          my_mandelbrot_inner_loop_object_01.loop();

          mandelbrot_thread_result_vector_01[col] = my_mandelbrot_inner_loop_object_01.get_iteration_result();
        }

        y_val_01 -= 4 * y_step_01;

        std::atomic_store(&mandelbrot_thread_do_loop_01, false);
      }
    }
  }

  template<typename NumericInputIteratorType>
  static void thread_routine_02(NumericInputIteratorType                                            x_first,
                                NumericInputIteratorType                                            x_last,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_val_02,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_step_02,
                                std::uint_fast32_t                                                  max_iter)
  {
    typedef typename std::iterator_traits<NumericInputIteratorType>::value_type local_value_type;

    const std::ptrdiff_t x_values_size = std::distance(x_first, x_last);

    mandelbrot_thread_result_vector_02.resize(std::size_t(x_values_size));

    while(std::atomic_load(&mandelbrot_thread_exit_flag_02) == false)
    {
      if(std::atomic_load(&mandelbrot_thread_do_loop_02) == true)
      {
        std::fill(mandelbrot_thread_result_vector_02.begin(),
                  mandelbrot_thread_result_vector_02.end(),
                  std::uint_fast32_t(0U));

        for(std::size_t col = 0U; col < std::size_t(x_values_size); ++col)
        {
          mandelbrot_inner_loop_object<local_value_type>
            my_mandelbrot_inner_loop_object_02(max_iter,
                                               *(x_first + col),
                                               y_val_02);

          my_mandelbrot_inner_loop_object_02.loop();

          mandelbrot_thread_result_vector_02[col] = my_mandelbrot_inner_loop_object_02.get_iteration_result();
        }

        y_val_02 -= 4 * y_step_02;

        std::atomic_store(&mandelbrot_thread_do_loop_02, false);
      }
    }
  }

  template<typename NumericInputIteratorType>
  static void thread_routine_03(NumericInputIteratorType                                            x_first,
                                NumericInputIteratorType                                            x_last,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_val_03,
                                typename std::iterator_traits<NumericInputIteratorType>::value_type y_step_03,
                                std::uint_fast32_t                                                  max_iter)
  {
    typedef typename std::iterator_traits<NumericInputIteratorType>::value_type local_value_type;

    const std::ptrdiff_t x_values_size = std::distance(x_first, x_last);

    mandelbrot_thread_result_vector_03.resize(std::size_t(x_values_size));

    while(std::atomic_load(&mandelbrot_thread_exit_flag_03) == false)
    {
      if(std::atomic_load(&mandelbrot_thread_do_loop_03) == true)
      {
        std::fill(mandelbrot_thread_result_vector_03.begin(),
                  mandelbrot_thread_result_vector_03.end(),
                  std::uint_fast32_t(0U));

        for(std::size_t col = 0U; col < std::size_t(x_values_size); ++col)
        {
          mandelbrot_inner_loop_object<local_value_type>
            my_mandelbrot_inner_loop_object_03(max_iter,
                                               *(x_first + col),
                                               y_val_03);

          my_mandelbrot_inner_loop_object_03.loop();

          mandelbrot_thread_result_vector_03[col] = my_mandelbrot_inner_loop_object_03.get_iteration_result();
        }

        y_val_03 -= 4 * y_step_03;

        std::atomic_store(&mandelbrot_thread_do_loop_03, false);
      }
    }
  }
} // namespace local

// This class generates the rows of the mandelbrot iteration.
// The coordinates are set up according to the Mandelbrot configuration.
template<typename NumericType>
class mandelbrot_generator
{
public:
  mandelbrot_generator(const mandelbrot_config_base<NumericType>& config)
    : mandelbrot_config_object(config),
      mandelbrot_image               (config.width(), config.height()),
      mandelbrot_view                (boost::gil::rgb8_view_t()),
      mandelbrot_iteration_matrix    (mandelbrot_config_object.width(),
                                      std::vector<std::uint_fast32_t>(mandelbrot_config_object.height())),
      mandelbrot_color_histogram     (static_cast<std::size_t>(config.max_iterations()) + 1U, UINT32_C(0))
  {
    mandelbrot_view = boost::gil::view(mandelbrot_image);
  }

  ~mandelbrot_generator() { }

  void generate_mandelbrot_image()
  {
    // Setup the x-axis coordinates.
    std::vector<NumericType> x_values(mandelbrot_config_object.width());

    const NumericType xy_step(mandelbrot_config_object.step());

    {
      // Initialize the x-axis coordinates (one time only).
      NumericType x_coord(mandelbrot_config_object.x_lo());

      for(NumericType& x : x_values)
      {
        x = x_coord;

        x_coord += xy_step;
      }
    }

    // Initialize the y-axis coordinate.
    const NumericType y_hi(mandelbrot_config_object.y_hi());
    const NumericType delta_y0 = y_hi - (0 * xy_step);
    const NumericType delta_y1 = y_hi - (1 * xy_step);
    const NumericType delta_y2 = y_hi - (2 * xy_step);
    const NumericType delta_y3 = y_hi - (3 * xy_step);
    const volatile std::uint_fast32_t max_iter = mandelbrot_config_object.max_iterations();

    std::thread thread_00(local::thread_routine_00<typename std::vector<NumericType>::const_iterator>,
                          x_values.cbegin(),
                          x_values.cend(),
                          delta_y0,
                          xy_step,
                          max_iter);

    std::thread thread_01(local::thread_routine_01<typename std::vector<NumericType>::const_iterator>,
                          x_values.cbegin(),
                          x_values.cend(),
                          delta_y1,
                          xy_step,
                          max_iter);

    std::thread thread_02(local::thread_routine_02<typename std::vector<NumericType>::const_iterator>,
                          x_values.cbegin(),
                          x_values.cend(),
                          delta_y2,
                          xy_step,
                          max_iter);

    std::thread thread_03(local::thread_routine_03<typename std::vector<NumericType>::const_iterator>,
                          x_values.cbegin(),
                          x_values.cend(),
                          delta_y3,
                          xy_step,
                          max_iter);

    // Loop through all the rows of pixels on the vertical
    // y-axis in the direction of decreasing y-value.

    std::uint_fast32_t row = UINT32_C(0);

    for(;;)
    {
      // Loop through the next columns of pixels on the horizontal
      // x-axis in the direction of increasing x-value.
      // Threads are used with one individual loop on the x-axis
      // localized in a given thread.

      std::atomic_store(&local::mandelbrot_thread_do_loop_00, true);
      std::atomic_store(&local::mandelbrot_thread_do_loop_01, true);
      std::atomic_store(&local::mandelbrot_thread_do_loop_02, true);
      std::atomic_store(&local::mandelbrot_thread_do_loop_03, true);

      for(;;)
      {
        if(   std::atomic_load(&local::mandelbrot_thread_do_loop_00) == false
           && std::atomic_load(&local::mandelbrot_thread_do_loop_01) == false
           && std::atomic_load(&local::mandelbrot_thread_do_loop_02) == false
           && std::atomic_load(&local::mandelbrot_thread_do_loop_03) == false)
        {
          break;
        }

        std::this_thread::sleep_for(std::chrono::microseconds(3U));
      }

      for(std::size_t col = 0U; col < x_values.size(); ++col)
      {
        const std::uint_fast32_t iteration_result_00 = local::mandelbrot_thread_result_vector_00[col];
        const std::uint_fast32_t iteration_result_01 = local::mandelbrot_thread_result_vector_01[col];
        const std::uint_fast32_t iteration_result_02 = local::mandelbrot_thread_result_vector_02[col];
        const std::uint_fast32_t iteration_result_03 = local::mandelbrot_thread_result_vector_03[col];

        mandelbrot_iteration_matrix[col][static_cast<std::size_t>(row) + 0U] = iteration_result_00;
        mandelbrot_iteration_matrix[col][static_cast<std::size_t>(row) + 1U] = iteration_result_01;
        mandelbrot_iteration_matrix[col][static_cast<std::size_t>(row) + 2U] = iteration_result_02;
        mandelbrot_iteration_matrix[col][static_cast<std::size_t>(row) + 3U] = iteration_result_03;

        ++mandelbrot_color_histogram[iteration_result_00];
        ++mandelbrot_color_histogram[iteration_result_01];
        ++mandelbrot_color_histogram[iteration_result_02];
        ++mandelbrot_color_histogram[iteration_result_03];
      }

      row += 4U;

      std::cout << "Calculating Mandelbrot image at row "
                << std::setw(6)
                << row
                << " of "
                << std::setw(6)
                << mandelbrot_config_object.height()
                << " total. Have patience."
                << "\r";

      if(row >= mandelbrot_config_object.height())
      {
        break;
      }
    }

    std::atomic_store(&local::mandelbrot_thread_exit_flag_00, true);
    std::atomic_store(&local::mandelbrot_thread_exit_flag_01, true);
    std::atomic_store(&local::mandelbrot_thread_exit_flag_02, true);
    std::atomic_store(&local::mandelbrot_thread_exit_flag_03, true);

    thread_00.join();
    thread_01.join();
    thread_02.join();
    thread_03.join();

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

    for(row = UINT32_C(0); row < mandelbrot_config_object.height(); ++row)
    {
      for(std::uint_fast32_t col = UINT32_C(0); col < mandelbrot_config_object.width(); ++col)
      {
        const std::uint_fast32_t color = mandelbrot_color_histogram[mandelbrot_iteration_matrix[col][row]];

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

        mandelbrot_view(col, row) = boost::gil::rgb8_pixel_t(the_color);
      }
    }

    boost::gil::jpeg_write_view("mandelbrot.jpg", mandelbrot_view);

    std::cout << std::endl
              << "The ouptput file mandelbrot.jpg has been written"
              << std::endl;
  }

private:
  const mandelbrot_config_base<NumericType>&   mandelbrot_config_object;
  boost::gil::rgb8_image_t                     mandelbrot_image;
  boost::gil::rgb8_view_t                      mandelbrot_view;
  std::vector<std::vector<std::uint_fast32_t>> mandelbrot_iteration_matrix;
  std::vector<std::uint_fast32_t>              mandelbrot_color_histogram;
};

int main()
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<32>,
                                        boost::multiprecision::et_off>
  numeric_type;

  #if defined BOOST_MANDELBROT_01_FULL

    // This is the classic full immage rendered in aqua tones (and black).
    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(2000), -9>;

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

//      mandelbrot_config_object("-0.74539836067085", "-0.74539836066275",
//                               "+0.11250463499195", "+0.11250463500005");

//      mandelbrot_config_object(-0.7453983606668L - 1.42E-12L, -0.7453983606668L + 1.42E-12L,
//                               +0.1125046349960L - 1.42E-12L, +0.1125046349960L + 1.42E-12L);

// OLD:
//     mandelbrot_config_object("-0.7453983606795", "-0.7453983606545",
//                              "+0.1125046349835", "+0.1125046350085");

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
                  "Error: Please use 80 or more decimal digits for deep dive 01.");

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

    using mandelbrot_config_type = mandelbrot_config<numeric_type, UINT32_C(10000), -192>;

    const numeric_type delta("+1.78E-55");
    const numeric_type cx   (numeric_type("-1.295189082147777457017064177185681926706566460884888469217455"));
    const numeric_type cy   (numeric_type("+0.440936982678320138880903678356262612113214627431396203682665"));

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #else

    #error: Mandelbrot image type is not defined!

  #endif

  typedef mandelbrot_config_type::mandelbrot_config_numeric_type mandelbrot_numeric_type;

  typedef mandelbrot_generator<mandelbrot_numeric_type> mandelbrot_generator_type;

  const std::clock_t start = std::clock();

  mandelbrot_generator_type* the_mandelbrot_generator = new mandelbrot_generator_type(mandelbrot_config_object);

  the_mandelbrot_generator->generate_mandelbrot_image();

  const float elapsed = (float(std::clock()) - float(start)) / float(CLOCKS_PER_SEC);

  std::cout << "Time for calculation: "
            << elapsed
            << "s"
            << std::endl;

  delete the_mandelbrot_generator;
}
