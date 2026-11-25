// Stand-alone driver for putidaFoam's microbeModel.C, to compare the
// solver's own rate and advance functions with the Python 0D model.
#include "microbeModel.H"

#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv)
{
    if (argc != 17)
    {
        std::fprintf(
            stderr,
            "usage: driver do_ml cap_n y_xs y_ms y_os x_max q_max K_o K_s "
            "X O2 G M t_adv dt n_updates\n");
        return 2;
    }
    using namespace microbemodel;
    int a = 1;
    const bool do_ml = std::atoi(argv[a++]) != 0;
    std::vector<double> p(10);
    const char* names[] = {"cap_n", "y_xs", "y_ms", "y_os",
                           "x_max", "q_max", "K_o", "K_s"};
    for (const char* name : names)
    {
        p[MM_param_keys[name]] = std::atof(argv[a++]);
    }
    p[MM_param_keys["o2_max"]] = 0.25906;
    p[MM_param_keys["kLa"]] = 0.0;
    std::vector<double> y(nvars, 0.0);
    y[X] = std::atof(argv[a++]);
    y[O2] = std::atof(argv[a++]);
    y[G] = std::atof(argv[a++]);
    y[M] = std::atof(argv[a++]);
    const double t_adv = std::atof(argv[a++]);
    const double dt = std::atof(argv[a++]);
    const int n_updates = std::atoi(argv[a++]);
    double now = 0.0;
    for (int i = 0; i < n_updates; i++)
    {
        const double our0 = get_our(y, nvars, p);
        advance(y, nvars, now, t_adv, dt, p, do_ml);
        const double our1 = get_our(y, nvars, p);
        now += t_adv;
        std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", now, y[X],
                    y[O2], y[G], y[M], y[H], 0.5 * (our0 + our1));
    }
    return 0;
}
