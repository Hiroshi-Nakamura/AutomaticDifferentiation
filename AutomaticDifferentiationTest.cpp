#include "AutomaticDifferentiation.hpp"
#include <iostream>

template<typename T>
T f(std::array<T,3> x){
    return x[0]*x[0]+x[1]*x[2]+1.0;
}

int main(int argc, char** argv)
{
    try{
        {
            AutomaticDifferentiation::FuncPtr<double,2> x0(new AutomaticDifferentiation::Variable<double,2>(0));
            AutomaticDifferentiation::FuncPtr<double,2> x1(new AutomaticDifferentiation::Variable<double,2>(1));
            AutomaticDifferentiation::FuncPtr<double,2> c(new AutomaticDifferentiation::Constant<double,2>(3.0));
            AutomaticDifferentiation::FuncPtr<double,2> y=x0*x0+x0*x1+x0/x1+c;

#ifdef WITHOUT_EIGEN
            const std::array<double,2> val{10.0,2.0};
#else
            Eigen::Matrix<double,2,1> val{10.0,2.0};
#endif // WITHOUT_EIGEN
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,2> dy_0=(*y).derivative(0);
            std::cout << "dy_0(val)" << (*dy_0)(val) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,2> dy_0_1=(*dy_0).derivative(1);
            std::cout << "dy_0_1(val)" << (*dy_0_1)(val) << std::endl;
        }
        std::cout << std::endl;
        {
            AutomaticDifferentiation::FuncPtr<double,1> x(new AutomaticDifferentiation::Variable<double,1>(0));
            AutomaticDifferentiation::FuncPtr<double,1> y=exp(x);
#ifdef WITHOUT_EIGEN
            const std::array<double,1> val{M_PI/3.0};
#else
            const Eigen::Matrix<double,1,1> val{1.0};
#endif // WITHOUT_EIGEN
            std::cout << "y(val)=" << (*y)(val) << std::endl;
            AutomaticDifferentiation::FuncPtr<double,1> dy=(*y).derivative();
            std::cout << "dy(val)=" << (*dy)(val) << std::endl;
        }
        std::cout << std::endl;
        {
            auto y=f(AutomaticDifferentiation::createVariables<double,3>());
            auto jac=AutomaticDifferentiation::jacobian(y);
            auto hes=AutomaticDifferentiation::hessian(y);
#ifdef WITHOUT_EIGEN
            const std::array<double,3> val{10.0, 2.0, 5.0};
#else
            const Eigen::Matrix<double,3,1> val{10.0, 2.0, 5.0};
#endif // WITHOUT_EIGEN
            std::cout << "y(val)=" << (*y)(val) << std::endl;
#ifdef WITHOUT_EIGEN
            std::cout << "jac(y)(val)=" << std::endl << AutomaticDifferentiation::to_string<double,3>(jac(val)) << std::endl;
            std::cout << "hes(y)(val)=" << std::endl << AutomaticDifferentiation::to_string<double,3>(hes(val)) << std::endl;
#else
            std::cout << "jac(y)(val)=" << std::endl << jac(val) << std::endl;
            std::cout << "hes(y)(val)=" << std::endl << hes(val) << std::endl;
#endif // WITHOUT_EIGEN
        }
    }catch(std::string message){
        std::cerr << message << std::endl;
    }
}


