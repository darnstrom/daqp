#include <cmath>
#include "api.h"

static double automatic_rho(const double* H, int n) {
    double min_diag = DAQP_INF;
    double max_row_sum = 0.0;
    double fro_norm_sq = 0.0;

    for (int i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            const double hij = H[i * n + j];
            const double sym = 0.5 * (hij + H[j * n + i]);
            row_sum += std::abs(sym);
            fro_norm_sq += hij * hij;
            if (i == j && sym < min_diag)
                min_diag = sym;
        }
        if (row_sum > max_row_sum)
            max_row_sum = row_sum;
    }

    if (min_diag > 0.0 && max_row_sum > 0.0)
        return std::sqrt(min_diag * max_row_sum);
    return std::sqrt(fro_norm_sq) / 2.0;
}

int main() {
    constexpr int n = 2;
    constexpr int m = 2;
    constexpr int ms = 0;

    double H[4] = {1.0, 1.75, 0.0, 1.0};
    double f[2] = {2.0, 2.0};
    double A[4] = {1.0, 0.0, 0.0, 1.0};
    double bupper[2] = {1.0, 1.0};
    double blower[2] = {-1.0, -1.0};
    int sense[2] = {0, 0};
    DAQPProblem qp = {n, m, ms, H, f, A, bupper, blower, sense, nullptr, 0, 1};

    DAQPWorkspace work = {0};
    int exitflag = setup_daqp(&qp, &work, nullptr);
    if (exitflag < 0)
        return 1;

    const double expected_auto = automatic_rho(H, n);
    const double auto_tol = 1e-12 * (expected_auto > 1.0 ? expected_auto : 1.0);
    const bool auto_ok = std::abs(work.avi->rho - expected_auto) <= auto_tol;
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);

    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.rho_avi = 0.25;

    DAQPWorkspace work_override = {0};
    work_override.settings = &settings;
    exitflag = setup_daqp(&qp, &work_override, nullptr);
    if (exitflag < 0)
        return 1;

    const bool override_ok = std::abs(work_override.avi->rho - settings.rho_avi) <= 1e-15;
    work_override.settings = nullptr;
    free_daqp_workspace(&work_override);
    free_daqp_ldp(&work_override);

    return (auto_ok && override_ok) ? 0 : 1;
}
