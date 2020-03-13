#include <ctime>

#include <mandelbrot.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

//#define BOOST_MANDELBROT_01_FULL
//#define BOOST_MANDELBROT_03_TOP
//#define BOOST_MANDELBROT_04_SWIRL
#define BOOST_MANDELBROT_05_SEAHORSES
//#define BOOST_MANDELBROT_06_BRANCHES
//#define BOOST_MANDELBROT_07_SEAHORSE_VALLEY
//#define BOOST_MANDELBROT_08_DEEP_DIVE_01
//#define BOOST_MANDELBROT_09_DEEP_DIVE_02

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
