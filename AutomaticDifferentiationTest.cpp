#include "AutomaticDifferentiation.hpp"
#include <iostream>

template<typename T>
T f(std::vector<T> x){
    return x[0]*x[0]+x[1]*x[2]+1.0;
}

int main(int argc, char** argv)
{
    try{
        {
            using namespace AutomaticDifferentiation;
            FuncPtr<double> x0(new Variable<double>(0));
            FuncPtr<double> x1(new Variable<double>(1));
            FuncPtr<double> c(new Constant<double>(3.0));
            FuncPtr<double> y=x0*x0+x0*x1+x0/x1+c;

            std::cout << toString(y) << std::endl;

            const std::vector<double> val{10.0,2.0};
            std::cout << "y(val)=" << (*y)(val) << std::endl;

            FuncPtr<double> dy_0=(*y).derivative(0);
            std::cout << "dy_0(val)=" << (*dy_0)(val) << std::endl;

            std::cout << "original:" << std::endl << toString(dy_0) << std::endl;
            simplification(dy_0);
            std::cout << "simplified:" << std::endl << toString(dy_0) << std::endl;

            FuncPtr<double> dy_0_1=(*dy_0).derivative(1);
            std::cout << "dy_0_1(val)=" << (*dy_0_1)(val) << std::endl;

            std::cout << "original:" << std::endl << toString(dy_0_1) << std::endl;
            simplification(dy_0_1);
            std::cout << "simplified:" << std::endl << toString(dy_0_1) << std::endl;
        }
        std::cout << std::endl;
        {
            using namespace AutomaticDifferentiation;
            FuncPtr<double> x(new Variable<double>(0));
            FuncPtr<double> y=cosh(x);
            const std::vector<double> val{1.0};
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            FuncPtr<double> dy=(*y).derivative();
            std::cout << "dy(val)=" << (*dy)(val) << std::endl;
        }
        std::cout << std::endl;
        {
            using namespace AutomaticDifferentiation;
            FuncPtr<double> y=f(createVariables<double>(3));
            MatFuncPtr<double> jac=jacobian(y,3);
            MatFuncPtr<double> hes=hessian(y,3);
            const std::vector<double> val{10.0, 2.0, 5.0};
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            std::cout << "jac(y)(val)=" << std::endl << jac(val) << std::endl;
            std::cout << "hes(y)(val)=" << std::endl << hes(val) << std::endl;
        }
    }catch(std::string message){
        std::cerr << message << std::endl;
    }
}


