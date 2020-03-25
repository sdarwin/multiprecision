#include <ctime>
#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <mandelbrot.hpp>

#define BOOST_MANDELBROT_IMAGE_INDEX_01_FULL                  1
#define BOOST_MANDELBROT_IMAGE_INDEX_03_TOP                   3
#define BOOST_MANDELBROT_IMAGE_INDEX_04_SWIRL                 4
#define BOOST_MANDELBROT_IMAGE_INDEX_05_SEAHORSES             5
#define BOOST_MANDELBROT_IMAGE_INDEX_06_BRANCHES              6
#define BOOST_MANDELBROT_IMAGE_INDEX_07_SEAHORSE_VALLEY       7
#define BOOST_MANDELBROT_IMAGE_INDEX_08_DEEP_DIVE_01          8
#define BOOST_MANDELBROT_IMAGE_INDEX_09_DEEP_DIVE_02          9
#define BOOST_MANDELBROT_IMAGE_INDEX_10_ZOOM_WIKI_01         10

#if !defined(BOOST_MANDELBROT_IMAGE_INDEX)
#define BOOST_MANDELBROT_IMAGE_INDEX BOOST_MANDELBROT_IMAGE_INDEX_05_SEAHORSES
#endif

int main()
{
  #if (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_01_FULL)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_01_FULL") + ".jpg";

    // This is the classic full immage.
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(2000), -10>;

    const mandelbrot_config_type mandelbrot_config_object(-2.000L, +0.500L,
                                                          -1.000L, +1.000L);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_03_TOP)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_03_TOP") + ".jpg";

    // This is a view of an upper part of the image (near the top of the classic full view).
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(1000), -12>;

    const mandelbrot_config_type mandelbrot_config_object(-0.130L - 0.282L, -0.130L + 0.282L,
                                                          +0.856L - 0.282L, +0.856L + 0.282L);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_04_SWIRL)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_04_SWIRL") + ".jpg";

    // This is a fanning swirl image.
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(2000), -22>;

    const mandelbrot_config_type mandelbrot_config_object(-0.749730L - 0.0002315L, -0.749730L + 0.0002315L,
                                                          -0.046608L - 0.0002315L, -0.046608L + 0.0002315L);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_05_SEAHORSES)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_05_SEAHORSES") + ".jpg";

    // This is a swirly seahorse image.
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(2000), -48>;

    const mandelbrot_config_type
      mandelbrot_config_object(-0.7453983606667815L - 1.76E-12L, -0.7453983606667815L + 1.76E-12L,
                               +0.1125046349959942L - 1.76E-12L, +0.1125046349959942L + 1.76E-12L);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_06_BRANCHES)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_06_BRANCHES") + ".jpg";

    // This is a spiral image of branches.
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(2000), -47>;

    const mandelbrot_config_type mandelbrot_config_object(+0.3369844464873L - 4.2E-12L, +0.3369844464873L + 4.2E-12L,
                                                          +0.0487782196791L - 4.2E-12L, +0.0487782196791L + 4.2E-12L);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_07_SEAHORSE_VALLEY)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<31>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_07_SEAHORSE_VALLEY") + ".jpg";

    // This is an image from the seahorse valley.
    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(1000), -15>;

    const mandelbrot_config_type
      mandelbrot_config_object("-0.748", "-0.700",
                               "+0.222", "+0.270");

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_08_DEEP_DIVE_01)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<127>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_08_DEEP_DIVE_01") + ".jpg";

    // This is a deep zoom image.
    // Note: Use 128 or more decimal digits for this iteration.

    static_assert(std::numeric_limits<numeric_type>::digits10 >= 127,
                  "Error: Please use 127 or more decimal digits for BOOST_MANDELBROT_08_DEEP_DIVE_01.");

    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(2000), -365>;

    const numeric_type delta("+1.25E-107");
    const numeric_type cx   ("-1.99999999913827011875827476290869498831680913663682095950680227271547027727918984035447670553861909622481524124");
    const numeric_type cy   ("+0.00000000000001314895443507637575136247566806505002151700520912095709529449343530548994027524594471095886432006");

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX == BOOST_MANDELBROT_IMAGE_INDEX_09_DEEP_DIVE_02)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<79>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_09_DEEP_DIVE_02") + ".jpg";

    // This is a deep zoom image.
    // Note: Use 79 or more decimal digits for this iteration.

    static_assert(std::numeric_limits<numeric_type>::digits10 >= 79,
                  "Error: Please use 79 or more decimal digits for BOOST_MANDELBROT_09_DEEP_DIVE_02.");

    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(10000), -191>;

    const numeric_type delta("+2.15E-55");
    const numeric_type cx   (numeric_type("-1.295189082147777457017064177185681926706566460884888469217456"));
    const numeric_type cy   (numeric_type("+0.440936982678320138880903678356262612113214627431396203682661"));

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #elif (BOOST_MANDELBROT_IMAGE_INDEX ==  BOOST_MANDELBROT_IMAGE_INDEX_10_ZOOM_WIKI_01)

    using numeric_type = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<47>>;

    const std::string str_filename = "images/mandelbrot_" + std::string("BOOST_MANDELBROT_10_ZOOM_WIKI_01") + ".jpg";

    // This is a medium zoom image from the zoom coordinates of:
    // http://en.wikipedia.org/wiki/File:Mandelbrot_sequence_new.gif
    // Note: Use 39 or more decimal digits for this iteration.

    static_assert(std::numeric_limits<numeric_type>::digits10 >= 47,
                  "Error: Please use 47 or more decimal digits for BOOST_MANDELBROT_10_ZOOM_WIKI_01.");

    using mandelbrot_config_type = boost::multiprecision::mandelbrot::mandelbrot_config<numeric_type, UINT32_C(20000), -91>;

    const numeric_type delta("+3.0E-25");
    const numeric_type cx   (numeric_type("-0.743643887037158704752191506114774"));
    const numeric_type cy   (numeric_type("+0.131825904205311970493132056385139"));

    const mandelbrot_config_type
      mandelbrot_config_object(cx - delta, cx + delta,
                               cy - delta, cy + delta);

  #else

    #error: Mandelbrot image type is not defined!

  #endif

  using mandelbrot_numeric_type = mandelbrot_config_type::mandelbrot_config_numeric_type;

  using mandelbrot_generator_type =
    boost::multiprecision::mandelbrot::mandelbrot_generator<mandelbrot_numeric_type,
                                                            mandelbrot_config_type::max_iterations>;

  const std::clock_t start = std::clock();

        boost::multiprecision::mandelbrot::detail::color_stretches_default local_color_stretches;
  const boost::multiprecision::mandelbrot::detail::color_functions_bw      local_color_functions;

  mandelbrot_generator_type mandelbrot_generator(mandelbrot_config_object);

  mandelbrot_generator.generate_mandelbrot_image(str_filename,
                                                 local_color_functions,
                                                 local_color_stretches);

  const float elapsed = (float(std::clock()) - float(start)) / float(CLOCKS_PER_SEC);

  std::cout << "Time for calculation: "
            << elapsed
            << "s"
            << std::endl;
}
