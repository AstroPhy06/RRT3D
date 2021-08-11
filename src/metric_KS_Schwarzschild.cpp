#include "metric_KS_Schwarzschild.hpp"


autodiff::dual KerrSchild_Schwarzschild::r2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3)
{
    return x1*x1+x2*x2+x3*x3;
}

autodiff::dual KerrSchild_Schwarzschild::r(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return sqrt(r2(x1,x2,x3));
}

autodiff::dual KerrSchild_Schwarzschild::H(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 2.0*_p.M/r(x1,x2,x3,x0);
}

/////////////////////////////////////////////

autodiff::dual KerrSchild_Schwarzschild::lx1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return x1/r(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::lx2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return x2/r(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::lx3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return x3/r(x1,x2,x3,x0);
}

/////////////////////////////////////////////


autodiff::dual KerrSchild_Schwarzschild::alpha(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    autodiff::dual alph = 1.0 / sqrt(1. + /*2. */ H(x1,x2,x3,x0));
    return alph;
}

autodiff::dual KerrSchild_Schwarzschild::B1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx1(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::B2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::B3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

/////////////////////////////////////////////

autodiff::dual KerrSchild_Schwarzschild::g11(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx1(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::g22(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx2(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::g33(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx3(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::g12(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::g13(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

autodiff::dual KerrSchild_Schwarzschild::g23(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx2(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

params KerrSchild_Schwarzschild::get_params()
{
    return _p;
}

