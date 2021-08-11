#include "metric_Kerr_BL.hpp"

// Kerr spacetime
// Boyer-Lindquist coordinates (r, theta, phi)


/////////////////////////////////////////////

// Beta_i = 0 for all i

autodiff::var alphaBL(autodiff::var r, const Params_v p)
{
    return sqrt(1.0 - p.rs/r);
}

/////////////////////////////////////////////

autodiff::var grr(autodiff::var r, const Params_v p)
{
    return 1.-p.rs/r;
}

autodiff::var gthth(autodiff::var r)
{
    return 1./(r*r);
}

autodiff::var gphiphi(autodiff::var r, autodiff::var theta)
{
    return 1./(r*r*sin(theta)*sin(theta));
}




/////////////////////////////////////////////

// Beta_i = 0 for all i

autodiff::dual alphaBL_d(autodiff::dual r,const Params_d p)
{
    return sqrt(1.0 - p.rs/r);
}

/////////////////////////////////////////////

autodiff::dual grrd(autodiff::dual r, const Params_d p)
{
    return 1.-p.rs/r;
}

autodiff::dual gththd(autodiff::dual r)
{
    return 1./(r*r);
}

autodiff::dual gphiphid(autodiff::dual r, autodiff::dual theta)
{
    return 1./(r*r*sin(theta)*sin(theta));
}
