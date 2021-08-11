#include "metric_Kerr_Schild.hpp"



autodiff::dual KerrSchild::rho2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3)
{
    return x1*x1+x2*x2+x3*x3;
}

autodiff::dual KerrSchild::r2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 0.5*((rho2(x1,x2,x3)-_p.a*_p.a) + sqrt(pow(rho2(x1,x2,x3)-_p.a*_p.a,2) + 4*_p.a*_p.a*x3*x3));
}

autodiff::dual KerrSchild::r(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return sqrt(r2(x1,x2,x3,x0));
}

autodiff::dual KerrSchild::H(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 2.0*_p.M*pow(r(x1,x2,x3,x0),3.)/(pow(r2(x1,x2,x3,x0),2) + _p.a*_p.a*x3*x3);
}

/////////////////////////////////////////////

autodiff::dual KerrSchild::lx1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return (r(x1,x2,x3,x0)*x1+_p.a*x2)/(r2(x1,x2,x3,x0) + _p.a*_p.a);
}

autodiff::dual KerrSchild::lx2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return (r(x1,x2,x3,x0)*x2-_p.a*x1)/(r2(x1,x2,x3,x0) + _p.a*_p.a);
}

autodiff::dual KerrSchild::lx3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return x3/r(x1,x2,x3,x0);
}

/////////////////////////////////////////////


autodiff::dual KerrSchild::alpha(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    autodiff::dual alph = 1.0 / sqrt(1. + /*2. */ H(x1,x2,x3,x0));
    return alph;
}

autodiff::dual KerrSchild::B1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx1(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::B2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::B3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return /*2.*/H(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*alpha(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

/////////////////////////////////////////////

autodiff::dual KerrSchild::g11(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx1(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::g22(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx2(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::g33(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return 1.-/*2.*/H(x1,x2,x3,x0)*lx3(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::g12(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx2(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::g13(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx1(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

autodiff::dual KerrSchild::g23(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0)
{
    return -/*2.*/H(x1,x2,x3,x0)*lx2(x1,x2,x3,x0)*lx3(x1,x2,x3,x0);
}

params KerrSchild::get_params()
{
    return _p;
}

