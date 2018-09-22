#ifndef AUTOMATICDIFFERENTIATION_HPP_INCLUDED
#define AUTOMATICDIFFERENTIATION_HPP_INCLUDED

#ifndef WITHOUT_EIGEN
#include <eigen3/Eigen/Core>
#endif // WITHOUT_EIGEN
#include <array>
#include <memory>
#include <cmath>
#include <sstream>
#include <iomanip>

#define OPERATOR_FUNC_DEFINE(operant,operant_name) \
    template<typename T, int DIM> \
    FuncPtr<T,DIM> operator operant (const FuncPtr<T,DIM>& left, const FuncPtr<T,DIM>& right){ \
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::operant_name, left, right)); \
    } \
    template<typename T, int DIM> \
    FuncPtr<T,DIM> operator operant (const FuncPtr<T,DIM>& left, const T& right_val){ \
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::operant_name, left, FuncPtr<T,DIM>(new Constant<T,DIM>(right_val)))); \
    } \
    template<typename T, int DIM> \
    FuncPtr<T,DIM> operator operant (const T& left_val, const FuncPtr<T,DIM>& right){ \
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::operant_name, FuncPtr<T,DIM>(new Constant<T,DIM>(left_val)), right)); \
    }


namespace AutomaticDifferentiation {

    template<typename T, int DIM>
    class Functor {
    public:
#ifdef WITHOUT_EIGEN
        virtual T operator()(const std::array<T,DIM>& x) const =0; /// The same usage as a mathmatical function.
#else
        virtual T operator()(const Eigen::Matrix<T,DIM,1>& x) const =0; /// The same usage as a mathmatical function.
#endif // WITHOUT_EIGEN
        virtual std::shared_ptr<Functor<T,DIM>> derivative(size_t idx=0) const =0; /// Get derivative of this Functor. The return is also the derived class of Functor.
    };

    template<typename T, int DIM>
    using FuncPtr=std::shared_ptr<Functor<T,DIM>>;

    template<typename T, int DIM>
    class Constant : public Functor<T,DIM> {
    private:
        const T val;
    public:
        Constant(const T& _val) : val(_val) {}
#ifdef WITHOUT_EIGEN
        virtual T operator()(const std::array<T,DIM>& x) const
#else
        virtual T operator()(const Eigen::Matrix<T,DIM,1>& x) const
#endif // WITHOUT_EIGEN
        {
            return val;
        }
        virtual FuncPtr<T,DIM> derivative(size_t idx) const
        {
            return FuncPtr<T,DIM>(new Constant(0.0));
        }
    };

    enum FuncType { SUM, DIFFERENCE, PRODUCT, QUOTIENT, MINUS, SIN, COS };

    template<typename T, int DIM>
    class Operator : public Functor<T,DIM> {
        const FuncType func_type;
        const FuncPtr<T,DIM> left;
        const FuncPtr<T,DIM> right;
    public:
        Operator(FuncType _func_type, const FuncPtr<T,DIM>& _left, const FuncPtr<T,DIM>& _right=nullptr) : func_type(_func_type), left(_left), right(_right) {}
#ifdef WITHOUT_EIGEN
        virtual T operator()(const std::array<T,DIM>& x) const
#else
        virtual T operator()(const Eigen::Matrix<T,DIM,1>& x) const
#endif // WITHOUT_EIGEN
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
            case FuncType::SIN:
                return std::sin( (*left)(x) );
                break;
            case FuncType::COS:
                return std::cos( (*left)(x) );
                break;
            default:
                throw std::string("Not defined operator in ")+__func__;
            }
        }
        virtual FuncPtr<T,DIM> derivative(size_t idx) const
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
            case FuncType::SIN:
                return (*left).derivative(idx) * cos(left);
                break;
            case FuncType::COS:
                return  - (*left).derivative(idx) * sin(left);
                break;
            default:
                throw std::string("Not defined operator in ")+__func__;
            }
        }
    };

    template<typename T, int DIM>
    class Variable : public Functor<T,DIM> {
    private:
        const size_t index;
    public:
        Variable(const size_t _index) : index(_index) {}
#ifdef WITHOUT_EIGEN
        virtual T operator()(const std::array<T,DIM>& x) const
#else
        virtual T operator()(const Eigen::Matrix<T,DIM,1>& x) const
#endif // WITHOUT_EIGEN
        {
            return x[index];
        }
        FuncPtr<T,DIM> derivative(size_t idx) const
        {
            if(idx==index) return FuncPtr<T,DIM>(new Constant<T,DIM>(1.0));
            else           return FuncPtr<T,DIM>(new Constant<T,DIM>(0.0));
        };
    };

    OPERATOR_FUNC_DEFINE(+,SUM);
    OPERATOR_FUNC_DEFINE(-,DIFFERENCE);
    OPERATOR_FUNC_DEFINE(*,PRODUCT);
    OPERATOR_FUNC_DEFINE(/,QUOTIENT);

    template<typename T, int DIM>
    FuncPtr<T,DIM> operator-(const FuncPtr<T,DIM>& functor){
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::MINUS, functor));
    }

    template<typename T, int DIM>
    FuncPtr<T,DIM> sin(const FuncPtr<T,DIM>& functor){
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::SIN, functor));
    }

    template<typename T, int DIM>
    FuncPtr<T,DIM> cos(const FuncPtr<T,DIM>& functor){
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::COS, functor));
    }

    ///
    /// utulity-- create variables
    ///
    template<typename T, int DIM>
    std::array<FuncPtr<T,DIM>,DIM> createVariables()
    {
        std::array<FuncPtr<T,DIM>,DIM> rtn;
        for(size_t i=0; i<DIM; i++){
            rtn[i]=FuncPtr<T,DIM>(new Variable<T,DIM>(i));
        }
        return rtn;
    }

    ///
    /// utulity-- calculate Jacobian Functor
    ///
    template<typename T, int DIM>
    class VecFuncPtr {
    private:
        std::array<FuncPtr<T,DIM>,DIM> func_ptr;
    public:
        VecFuncPtr(const std::array<FuncPtr<T,DIM>,DIM>& _func_ptr) : func_ptr(_func_ptr) {}
#ifdef WITHOUT_EIGEN
        virtual std::array<T,DIM> operator()(const std::array<T,DIM>& x) const
        {
            std::array<T,DIM> rtn;
            for(size_t i=0; i<DIM; i++){
                rtn[i]=(*func_ptr[i])(x);
            }
            return rtn;
        }
#else
        virtual Eigen::Matrix<T,DIM,1> operator()(const Eigen::Matrix<T,DIM,1>& x) const
        {
            Eigen::Matrix<T,DIM,1> rtn;
            for(size_t i=0; i<DIM; i++){
                rtn(i,0)=(*func_ptr[i])(x);
            }
            return rtn;
        }
#endif // WITHOUT_EIGEN
        virtual FuncPtr<T,DIM> operator[](const size_t idx) const
        {
            return func_ptr[idx];
        }
        virtual std::array<FuncPtr<T,DIM>,DIM> to_array() const
        {
            return func_ptr;
        }

    };

    template<typename T, int DIM>
    VecFuncPtr<T,DIM> jacobian(const FuncPtr<T,DIM>& f)
    {
        std::array<FuncPtr<T,DIM>,DIM> rtn;
        for(size_t i=0; i<DIM; i++){
            rtn[i]=(*f).derivative(i);
        }
        return VecFuncPtr<T,DIM>(rtn);
    }


    ///
    /// utulity-- calculate Hessian Functor
    ///
    template<typename T, int DIM>
    class MatFuncPtr {
    private:
        std::array<std::array<FuncPtr<T,DIM>,DIM>,DIM> func_ptr;
    public:
        MatFuncPtr(const std::array<std::array<FuncPtr<T,DIM>,DIM>,DIM>& _func_ptr) : func_ptr(_func_ptr) {}
#ifdef WITHOUT_EIGEN
        virtual std::array<std::array<T,DIM>,DIM> operator()(const std::array<T,DIM>& x) const
        {
            std::array<std::array<T,DIM>,DIM> rtn;
            for(size_t i=0; i<DIM; i++){
                for(size_t j=0; j<DIM; j++){
                    rtn[i][j]=(*func_ptr[i][j])(x);
                }
            }
            return rtn;
        }
#else
        virtual Eigen::Matrix<T,DIM,DIM> operator()(const Eigen::Matrix<T,DIM,1>& x) const
        {
            Eigen::Matrix<T,DIM,DIM> rtn;
            for(size_t i=0; i<DIM; i++){
                for(size_t j=0; j<DIM; j++){
                    rtn(i,j)=(*func_ptr[i][j])(x);
                }
            }
            return rtn;
        }
#endif // WITHOUT_EIGEN
        virtual std::array<FuncPtr<T,DIM>,DIM> operator[](const size_t idx) const
        {
            return func_ptr[idx];
        }

    };

    template<typename T, int DIM>
    MatFuncPtr<T,DIM> hessian(const FuncPtr<T,DIM>& f)
    {
        std::array<std::array<FuncPtr<T,DIM>,DIM>,DIM> rtn;
        auto jac=jacobian(f);
        for(size_t i=0; i<DIM; i++){
            rtn[i]=jacobian(jac[i]).to_array();
        }
        return MatFuncPtr<T,DIM>(rtn);
    }

#ifdef WITHOUT_EIGEN
    ///
    /// utulity-- convert string for display
    ///
    template<typename T, int DIM>
    std::string to_string(const std::array<T,DIM>& vec)
    {
        std::ostringstream oss;
        for(auto e: vec){
            oss << " " << std::setw(8) << e << std::endl;
        }
        return oss.str();
    }

    template<typename T, int DIM>
    std::string to_string(const std::array<std::array<T,DIM>,DIM>& mat)
    {
        std::ostringstream oss;
        for(auto row: mat){
            for(auto e: row){
                oss << " " << std::setw(8) << e ;
            }
            oss << std::endl;
        }
        return oss.str();
    }
#endif // WITHOUT_EIGEN
}
#endif // AUTOMATICDIFFERENTIATION_HPP_INCLUDED
