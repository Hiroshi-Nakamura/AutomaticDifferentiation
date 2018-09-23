# AutomaticDifferentiation
Implementation of automatic differentiation by using of Functor class and its operator overload.
You can easily run automatic differentiation program by including "AutomaticDifferentiation.hpp" to your sorce code.
This library also allows you higher-order automatic differenciation, such as Hessian.
The sample of usages are shown in "AutomaticDifferentiationTest.cpp".

Generally a Functor has operator().
Of course, this AutomaticDifferentiation::Functor has operator(),
whose argument is std::array<T,DIM> standing for a mathmatical variable vector x.
AutomaticDifferentiation::Functor derives 3 types, Constant, Variable and Operator.
These are implemeted by each class.
AutomaticDifferentiation::Constant operator() always returns a constant value for any x.
AutomaticDifferentiation::Variable operator() returns a value due to vector x.
AutomaticDifferentiation::Operator operator() returns a value due to its operation type, sum or product or so on,
and due to the operated Functor(s).

The above are the basics of AutomaticDifferentiation::Functor.
After this preperation, we can easily derive the derivative of AutomaticDifferentiation::Functor.
The derivative of AutomaticDifferentiation::Constant is always AutomaticDifferentiation::Constant(0.0).
The derivative of AutomaticDifferentiation::Variable will be AutomaticDifferentiation::Constant(1.0) when derived by itself,
otherwise (derived by other varable) it will be AutomaticDifferentiation::Constant(0.0).


"AutomaticDifferentiation.cbp" is a project manage file for Code::Blocks.


