#include <stdio.h>
#include <autodiff/forward.hpp>
#include <autodiff/reverse.hpp>


struct Params_d
{
    autodiff::dual rs;
};


struct Params_v
{
    autodiff::var rs;
};

/////////////////////////////////////////////

autodiff::var alphaBL(autodiff::var r, const Params_v p);

/////////////////////////////////////////////

autodiff::var grr(autodiff::var r, const Params_v p);
autodiff::var gthth(autodiff::var r);
autodiff::var gphiphi(autodiff::var r, autodiff::var theta);

///////////////////////////////////////////////////////////////
/*************************************************************/
///////////////////////////////////////////////////////////////


autodiff::dual alphaBL_d(autodiff::dual r,const Params_d p);

/////////////////////////////////////////////

autodiff::dual grrd(autodiff::dual r, const Params_d p);
autodiff::dual gththd(autodiff::dual r);
autodiff::dual gphiphid(autodiff::dual r, autodiff::dual theta);
