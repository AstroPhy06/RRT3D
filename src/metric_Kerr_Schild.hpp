#ifndef KS_H_
#define KS_H_


#include <stdio.h>
#include <autodiff/forward.hpp>
#include "metric.hpp"


class KerrSchild : public Metric
{

public:

	KerrSchild(params &p): Metric(p){};
	/// Metric specific

	autodiff::dual rho2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3);

	autodiff::dual r2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual r(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual H(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual lx1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual lx2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual lx3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	/// ADM formalism

	autodiff::dual alpha(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual B1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual B2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual B3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	/// Metric tensor

	autodiff::dual g11(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual g22(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual g33(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual g12(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual g13(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	autodiff::dual g23(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0);

	params get_params();

};

#endif
