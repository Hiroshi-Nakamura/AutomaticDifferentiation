# AutomaticDifferentiation
Implementation of automatic differentiation by using of a Functor class and its operator overload.
You can easily run automatic differentiation program by including "AutomaticDifferentiation.hpp" to your sorce code.
This library also allows you higher-order automatic differenciation, such as Hessian.
The sample of usages are shown in "AutomaticDifferentiationTest.cpp".

Generally a Functor has an operator().
Of course, this AutomaticDifferentiation::Functor has the operator(),
whose argument is std::array standing for a mathmatical variable vector x.
AutomaticDifferentiation::Functor derives 3 types classes, Constant, Variable and Operator.
`AutomaticDifferentiation::Constant::operator(x)` always returns a constant value for any x.
`AutomaticDifferentiation::Variable::operator(x)` returns a value due to vector x.
`AutomaticDifferentiation::Operator::operator(x)` returns a value due to its operation type, SUM or PRODUCT or so on,
and due to the operated Functor(s).
For example, in the SUM operator case, the operator(x) returns the sum of the left and the right Functor's operator(x) outputs,
i.e. returns **LEFT(x)+RIGHT(x)**, where **LEFT** is the left operated functor, **RIGHT** is the right one.

The above are the basics of AutomaticDifferentiation::Functor.
After this preperation, we can easily calculate the derivative Functor.
The derivative of AutomaticDifferentiation::Constant is always AutomaticDifferentiation::Constant(0.0).
The derivative of AutomaticDifferentiation::Variable will be AutomaticDifferentiation::Constant(1.0) when derived by itself,
otherwise (derived by other varable) it will be AutomaticDifferentiation::Constant(0.0).
The derivative of AutomaticDifferentiation::Operator will be calculated by derivative chain rule.
For example, in the PRODUCT operator case, it will be **LEFT.derivertive()\*RIGHT+LEFT\*RIGHT.derivertive()**.
In the hpp code, this calculation is implemented by derivative(), whose argument is index of vector x.
If you want derivative by the first variable of vector x, type *derivative(0)*.

Note that the output of derivative() is not a value, but a Functor (exactly a pointer of Functor).
That's why we can calculate higher-order derivative by repeating derivative(). For example, y.derivative(1).derivative(0) will be ddy/dx_1 dx_0.

Because of class morphorism, FuncPtr, which is a std::shared_ptr of Functor, is prepared.
In your source code, don't use Functor itself but FuncPtr.
When you want to create Variable, type
    AutomaticDifferentiation::FuncPtr x0(new AutomaticDifferentiation::Variable(0));
In the case of Constant, type
    AutomaticDifferentiation::FuncPtr c(new AutomaticDifferentiation::Constant(3.0));
The argument of Constructor of Variable and Constant are different,
the formar (size_t) means the index of vector x, the latter (usually "double") means the constant value itself.
Usually the Operator is not created explicitly,
but by typing the equation, i.e.
    AutomaticDifferentiation::FuncPtr y=x0*x0+x0*x1+x0/x1+c
, you can take Operator instance y.

The Above are the fundamental part of this library. Adding this, some utilities are prepared.

It is nanutal that the argument of a mathmatical function is much long vector.
In such case, because we have to prepare
    AutomaticDifferentiation::FuncPtr x_1(new AutomaticDifferentiation::Variable(0));*,...,*AutomaticDifferentiation::FuncPtr x_n(new AutomaticDifferentiation::Variable(n-1));*, the function *AutomaticDifferentiation::createVariables()
is prepared.
By typing
    AutomaticDifferentiation::FuncPtr y=f(AutomaticDifferentiation::createVariables());
, you can take a Functor y standing for function f().

Jacobian and Hessian are popular.
The functions `AutomaticDifferentiation::jacobian()` and `AutomaticDifferentiation::hessian()` are prepared.
Then you can obtain them by typing
    auto jac=AutomaticDifferentiation::jacobian(y);* and *auto hes=AutomaticDifferentiation::hessian(y);

