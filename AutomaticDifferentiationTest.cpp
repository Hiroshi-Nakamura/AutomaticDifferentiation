#include "AutomaticDifferentiation.hpp"
#include <iostream>

template<typename T>
T f(std::array<T,3> x){
    return x[0]*x[0]+x[1]+sin(x[2])+1.0;
}

int main(int argc, char** argv)
{
    try{
        {
            AutomaticDifferentiation::FuncPtr<double,2> x0(new AutomaticDifferentiation::Variable<double,2>(0));
            AutomaticDifferentiation::FuncPtr<double,2> x1(new AutomaticDifferentiation::Variable<double,2>(1));
            AutomaticDifferentiation::FuncPtr<double,2> c(new AutomaticDifferentiation::Constant<double,2>(3.0));
            AutomaticDifferentiation::FuncPtr<double,2> y=x0*x0+x0*x1+x0/x1+c;

            const std::array<double,2> val{10.0,2.0};
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,2> dy_0=(*y).derivative(0);
            std::cout << "dy_0(val)" << (*dy_0)(val) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,2> dy_0_1=(*dy_0).derivative(1);
            std::cout << "dy_0_1(val)" << (*dy_0_1)(val) << std::endl;
        }

        {
            AutomaticDifferentiation::FuncPtr<double,1> x(new AutomaticDifferentiation::Variable<double,1>(0));
            AutomaticDifferentiation::FuncPtr<double,1> y=cos(x);
            const std::array<double,1> alpha{M_PI/3.0};
            std::cout << "y(alpha)=" << (*y)(alpha) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,1> dy=(*y).derivative();
            std::cout << "dy(alpha)=" << (*dy)(alpha) << std::endl;
        }

        {
            auto y=f(AutomaticDifferentiation::createVariables<double,3>());
            auto jac=AutomaticDifferentiation::jacobian(y);
            auto hessian=AutomaticDifferentiation::hessian(y);
            const std::array<double,3> val{10.0, 2.0, 5.0};
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            std::cout << "jac(y)(val)=" << AutomaticDifferentiation::to_string<double,3>(jac(val)) << std::endl;
            std::cout << "hessian(y)(val)=" << AutomaticDifferentiation::to_string<double,3>(hessian(val)) << std::endl;
        }
    }catch(std::string message){
        std::cerr << message << std::endl;
    }
}


