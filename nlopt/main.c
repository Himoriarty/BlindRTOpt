#include <stdio.h>
#include <math.h>
#include "nlopt.h"
#define INF (1.7976931348623158e+308)
#define PI 3.141592653589793


typedef struct
{
    int    len;
    float *decay_data;
    double a, b;

}func_data;

double alpha_opt(int n, const double *x, double *grad, void *data) {
	func_data *ptr = (func_data*)data;
	/* FIXME ADD THE STUFF IN HERE TO DO WITH READING nld */
	int i;
	double sigma = 0.0;
	for (i = 0; i < ptr->len; i++) {
		sigma += pow(ptr->decay_data[i], 2) / pow(x[0] * pow(ptr->a, i) + (1 - x[0])*pow(ptr->b, i), 2);
	}
	sigma = sqrt(sigma / ptr->len);

	double sum = 0.0;
	double f1 = 0.0;
	double f3 = 0.0;

	for (i = 0; i < ptr->len; i++) {
		f1 += 0.5*pow(ptr->decay_data[i], 2)*pow(sigma*(x[0] * pow(ptr->a, i) + (1 - x[0])*pow(ptr->b, i)), -2);
		f3 += log(x[0] * pow(ptr->a, i) + (1 - x[0])*pow(ptr->b, i));
	}
	double fc = f1 + 0.5* ptr->len*log(2 * PI * pow(sigma, 2)) + f3;

	return fc;
}

double par3_opt(int n, const double *x, double *grad, void *data) {
	func_data *ptr = (func_data*)data;
	/* FIXME ADD THE STUFF IN HERE TO DO WITH READING nld */
	int i;
	double sigma = 0.0;
	for (i = 0; i < ptr->len; i++) {
		sigma += pow(ptr->decay_data[i], 2) / pow(x[2] * pow(x[0], i) + (1 - x[2])*pow(x[1], i), 2);
	}
	sigma = sqrt(sigma / ptr->len);

	double sum = 0.0;
	double f1 = 0.0;
	double f3 = 0.0;

	for (i = 0; i < ptr->len; i++) {
		f1 += 0.5*pow(ptr->decay_data[i], 2)*pow(sigma*(x[2] * pow(x[0], i) + (1 - x[2])*pow(x[1], i)), -2);
		f3 += log(x[2] * pow(x[0], i) + (1 - x[2])*pow(x[1], i));
	}
	double fc = f1 + 0.5* ptr->len*log(2 * PI * pow(sigma, 2)) + f3;

	return fc;
}


int main(void)
{
	double alpha_f_min = INF;
	double par3_f_min = INF;
	double lb[3] = { 0 };
	double ub[3] = {1};

	double x[3] = {0.5};
	
    float win[100];
    func_data data;
    int i;
    for (i = 100; i>=1; i--)
    {
        win[i-1] = i;
    }
    data.a = 0.992;
    data.b = 0.995;
    data.decay_data = win;
    data.len = 100;
   // alpha objective function
    // set up optimizer
    nlopt_opt opter = nlopt_create(NLOPT_LN_COBYLA, 1);
    // lower and upper bound
    nlopt_set_lower_bounds1(opter, lb[0]);
    nlopt_set_upper_bounds1(opter, ub[0]);
	//nlopt_set_maxeval(opter, 100);
 
	nlopt_set_min_objective(opter, (nlopt_func)alpha_opt, (void *) &data);
	nlopt_result result = nlopt_optimize(opter, &x, &alpha_f_min);
	printf("Minimum utility=%f, x=(%f)\n", alpha_f_min, x[0]);
    if (result)
        printf("Minimum utility=%f, x=(%f)\n", alpha_f_min, x[0]);
    // free
    nlopt_destroy(opter);

// a b alpha
	nlopt_opt opter2 = nlopt_create(NLOPT_LN_COBYLA, 3);
	// lower and upper bound
	lb[0] = lb[1] = 0.9954;  //0.9954 < a,b < 0.9999
	ub[0] = ub[1] = 0.9999;  
	lb[2] = 0; ub[2] = 1; // 0.0 < alpha < 1.0
	nlopt_set_lower_bounds(opter2, lb);
	nlopt_set_upper_bounds(opter2, ub);
	//nlopt_set_maxeval(opter2, 300);
	nlopt_set_min_objective(opter2, (nlopt_func)par3_opt, (void *)&data);
	x[2] = 0.5;//alpha
	x[0] = x[1] = 0.9976;//a, b
	nlopt_result result = nlopt_optimize(opter2, &x, &par3_f_min);
	printf("Minimum utility=%f, a=%f  b=%f  alpha=%f\n", par3_f_min, x[0], x[1], x[2]);


	if (result)
		printf("Minimum utility=%f, a=%f  b=%f  alpha=%f\n", par3_f_min, x[0],x[1],x[2]);
	// free
	nlopt_destroy(opter2);

    return 0;
}
