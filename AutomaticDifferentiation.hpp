#ifndef AUTOMATICDIFFERENTIATION_HPP_INCLUDED
#define AUTOMATICDIFFERENTIATION_HPP_INCLUDED

#include <eigen3/Eigen/Core>
#include <vector>
#include <memory>
#include <cmath>

#define OPERATOR_FUNC_DEFINE(operant,operant_name) \
    template<typename T> \
    FuncPtr<T> operator operant (const FuncPtr<T>& left, const FuncPtr<T>& right){ \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, left, right)); \
    } \
    template<typename T> \
    FuncPtr<T> operator operant (const FuncPtr<T>& left, const T& right_val){ \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, left, FuncPtr<T>(new Constant<T>(right_val)))); \
    } \
    template<typename T> \
    FuncPtr<T> operator operant (const T& left_val, const FuncPtr<T>& right){ \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, FuncPtr<T>(new Constant<T>(left_val)), right)); \
    }

#define MATH_FUNC_DEFINE(operant,operant_name) \
    template<typename T> \
    FuncPtr<T> operant(const FuncPtr<T>& functor){ \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, functor)); \
    }


namespace AutomaticDifferentiation {

    template<typename T>
    class Functor {
    public:
        virtual T operator()(const T* x) const =0; /// The same usage as a mathmatical function.
        T operator()(const std::vector<T>& x) const { return (*this)(x.data()); }; /// The same usage as a mathmatical function. NOT virtual.
        T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1>& x) const { return (*this)(x.data()); }; /// The same usage as a mathmatical function. NOT virtual.
        virtual std::shared_ptr<Functor<T>> derivative(size_t idx=0) const =0; /// Get derivative of this Functor. The return is also the derived class of Functor.
    };

    template<typename T>
    using FuncPtr=std::shared_ptr<Functor<T>>;

    template<typename T>
    class Constant : public Functor<T> {
    private:
        const T val;
    public:
        Constant(const T& _val) : val(_val) {}
        virtual T operator()(const T* x) const
        {
            return val;
        }
        virtual FuncPtr<T> derivative(size_t idx) const
        {
            return FuncPtr<T>(new Constant(0.0));
        }
    };

    enum FuncType { SUM, DIFFERENCE, PRODUCT, QUOTIENT, MINUS, COS, SIN, TAN, ACOS, ASIN, ATAN, EXP, LOG, SQRT };

    template<typename T>
    class Operator : public Functor<T> {
        const FuncType func_type;
        const FuncPtr<T> left;
        const FuncPtr<T> right;
    public:
        Operator(FuncType _func_type, const FuncPtr<T>& _left, const FuncPtr<T>& _right=nullptr) : func_type(_func_type), left(_left), right(_right) {}
        virtual T operator()(const T* x) const
        {
            switch(func_type){
            case FuncType::SUM:
                return (*left)(x) + (*right)(x);
                break;
            case FuncType::DIFFERENCE:
                return (*left)(x) - (*right)(x);
                break;
            case FuncType::PRODUCT:
                return (*left)(x) * (*right)(x);
                break;
            case FuncType::QUOTIENT:
                return (*left)(x) / (*right)(x);
                break;
            case FuncType::MINUS:
                return -(*left)(x);
                break;
            case FuncType::COS:
                return std::cos( (*left)(x) );
                break;
            case FuncType::SIN:
                return std::sin( (*left)(x) );
                break;
            case FuncType::TAN:
                return std::tan( (*left)(x) );
                break;
            case FuncType::ACOS:
                return std::acos( (*left)(x) );
                break;
            case FuncType::ASIN:
                return std::asin( (*left)(x) );
                break;
            case FuncType::ATAN:
                return std::atan( (*left)(x) );
                break;
            case FuncType::EXP:
                return std::exp( (*left)(x) );
                break;
            case FuncType::LOG:
                return std::log( (*left)(x) );
                break;
            case FuncType::SQRT:
                return std::sqrt( (*left)(x) );
                break;
            default:
                throw std::string("Not defined operator in ")+__func__;
            }
        }
        virtual FuncPtr<T> derivative(size_t idx) const
        {
            switch(func_type){
            case FuncType::SUM:
                return (*left).derivative(idx) + (*right).derivative(idx);
                break;
            case FuncType::DIFFERENCE:
                return (*left).derivative(idx) - (*right).derivative(idx);
                break;
            case FuncType::PRODUCT:
                return (*left).derivative(idx) * right + left * (*right).derivative(idx);
                break;
            case FuncType::QUOTIENT:
                return (*left).derivative(idx) / right - left * (*right).derivative(idx) / right / right;
                break;
            case FuncType::MINUS:
                return -(*left).derivative(idx);
                break;
            case FuncType::COS:
                return  - (*left).derivative(idx) * sin(left);
                break;
            case FuncType::SIN:
                return (*left).derivative(idx) * cos(left);
                break;
            case FuncType::TAN:
                return (*left).derivative(idx) / cos(left) / cos(left);
                break;
            case FuncType::ACOS:
                return  - (*left).derivative(idx) / sqrt(1.0 - left*left);
                break;
            case FuncType::ASIN:
                return (*left).derivative(idx)  / sqrt(1.0 - left*left);
                break;
            case FuncType::ATAN:
                return (*left).derivative(idx) / (left*left + 1.0);
                break;
            case FuncType::EXP:
                return (*left).derivative(idx) * exp(left);
                break;
            case FuncType::LOG:
                return (*left).derivative(idx) / (left);
                break;
            case FuncType::SQRT:
                return (*left).derivative(idx) / 2.0*sqrt(left);
                break;
            default:
                throw std::string("Not defined operator in ")+__func__;
            }
        }
    };

    template<typename T>
    class Variable : public Functor<T> {
    private:
        const size_t index;
    public:
        Variable(const size_t _index) : index(_index) {}
        virtual T operator()(const T* x) const
        {
            return x[index];
        }
        FuncPtr<T> derivative(size_t idx) const
        {
            if(idx==index) return FuncPtr<T>(new Constant<T>(1.0));
            else           return FuncPtr<T>(new Constant<T>(0.0));
        };
    };

    OPERATOR_FUNC_DEFINE(+,SUM);
    OPERATOR_FUNC_DEFINE(-,DIFFERENCE);
    OPERATOR_FUNC_DEFINE(*,PRODUCT);
    OPERATOR_FUNC_DEFINE(/,QUOTIENT);

    MATH_FUNC_DEFINE(cos,COS);
    MATH_FUNC_DEFINE(sin,SIN);
    MATH_FUNC_DEFINE(tan,TAN);
    MATH_FUNC_DEFINE(acos,ACOS);
    MATH_FUNC_DEFINE(asin,ASIN);
    MATH_FUNC_DEFINE(atan,ATAN);
    MATH_FUNC_DEFINE(exp,EXP);
    MATH_FUNC_DEFINE(log,LOG);
    MATH_FUNC_DEFINE(sqrt,SQRT);

    template<typename T>
    FuncPtr<T> operator-(const FuncPtr<T>& functor){
        return FuncPtr<T>(new Operator<T>(FuncType::MINUS, functor));
    }

    ///
    /// utulity-- create variables
    ///
    template<typename T>
    std::vector<FuncPtr<T>> createVariables(size_t dim)
    {
        std::vector<FuncPtr<T>> rtn;
        for(size_t i=0; i<dim; i++){
            rtn.emplace_back(new Variable<T>(i));
        }
        return rtn;
    }

    ///
    /// Matrix define for FuncPnr of Jacobian and Hessian
    ///
    template<typename T>
    class MatFuncPtr {
    private:
        const size_t nRows;
        const size_t nCols;
        FuncPtr<T>* func_ptr;
    public:
        MatFuncPtr(const size_t _nRows, const size_t _nCols=1) : nRows(_nRows), nCols(_nCols), func_ptr(new FuncPtr<T>[nRows*nCols]) {}
        ~MatFuncPtr(){ delete[] func_ptr; }
        FuncPtr<T>& operator()(const size_t row, const size_t col=0){ return func_ptr[nCols*row+col]; }
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> operator()(const T* x) const
        {
            Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rtn(nRows,nCols);
            for(size_t i=0; i<nRows; i++){
                for(size_t j=0; j<nCols; j++){
                    rtn(i,j)=(*func_ptr[nCols*i+j])(x);
                }
            }
            return rtn;
        }
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> operator()(const std::vector<T> x) const { return (*this)(x.data()); }
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> x) const { return (*this)(x.data()); }
    };

    ///
    /// utulity-- calculate Jacobian Functor
    ///
    template<typename T>
    MatFuncPtr<T> jacobian(const FuncPtr<T>& f, const size_t dim)
    {
        MatFuncPtr<T> rtn(dim);
        for(size_t i=0; i<dim; i++){
            rtn(i,0)=(*f).derivative(i);
        }
        return rtn;
    }

    ///
    /// utulity-- calculate Hessian Functor
    ///
    template<typename T>
    MatFuncPtr<T> hessian(const FuncPtr<T>& f, const size_t dim)
    {
        MatFuncPtr<T> rtn(dim,dim);
        MatFuncPtr<T> jac=jacobian(f,dim);
        for(size_t i=0; i<dim; i++){
            MatFuncPtr<T> jac_jac=jacobian(jac(i),dim);
            for(size_t j=0; j<dim; j++){
                rtn(i,j)=std::move(jac_jac(j));
            }
        }
        return rtn;
    }
}

#endif // AUTOMATICDIFFERENTIATION_HPP_INCLUDED
