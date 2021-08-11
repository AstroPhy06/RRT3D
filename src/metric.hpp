#ifndef METRIC_H_
#define METRIC_H_


// Each metric should provide
// metric components gxx
// Bi
// a
// distance to BH criteria?


typedef struct Params
{
    autodiff::dual a;
    autodiff::dual M;

    Params(double A,double m): a(A), M(m) {};

} params;


class Metric
{
protected:
	params _p;

public:

	Metric(params &p) : _p(p){}

	// Define components of the metric
	virtual autodiff::dual g11(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual g12(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual g13(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual g22(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual g23(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual g33(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	
	// required by 3+1 formalism
	virtual autodiff::dual alpha(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual B1(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual B2(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;
	virtual autodiff::dual B3(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;

        // used to compute distance from BH
        virtual autodiff::dual r(autodiff::dual x1, autodiff::dual x2, autodiff::dual x3, autodiff::dual x0) = 0;

	virtual params get_params() = 0;	
};

#endif
