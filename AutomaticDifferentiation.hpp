#ifndef AUTOMATICDIFFERENTIATION_HPP_INCLUDED
#define AUTOMATICDIFFERENTIATION_HPP_INCLUDED

#ifndef WITHOUT_EIGEN
#include <eigen3/Eigen/Core>
#endif // WITHOUT_EIGEN
#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include <sstream>
#include <iomanip>

#define OPERATOR_FUNC_DEFINE_VEC(operant,operant_name) \
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

#define MATH_FUNC_DEFINE_VEC(operant,operant_name) \
    template<typename T> \
    FuncPtr<T> operant(const FuncPtr<T>& functor){ \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, functor)); \
    }


namespace AutomaticDifferentiation_Vector {

    template<typename T>
    class Functor {
    public:
        virtual T operator()(const std::vector<T>& x) const =0; /// The same usage as a mathmatical function.
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
        virtual T operator()(const std::vector<T>& x) const
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
        virtual T operator()(const std::vector<T>& x) const
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
        virtual T operator()(const std::vector<T>& x) const
        {
            return x[index];
        }
        FuncPtr<T> derivative(size_t idx) const
        {
            if(idx==index) return FuncPtr<T>(new Constant<T>(1.0));
            else           return FuncPtr<T>(new Constant<T>(0.0));
        };
    };

    OPERATOR_FUNC_DEFINE_VEC(+,SUM);
    OPERATOR_FUNC_DEFINE_VEC(-,DIFFERENCE);
    OPERATOR_FUNC_DEFINE_VEC(*,PRODUCT);
    OPERATOR_FUNC_DEFINE_VEC(/,QUOTIENT);

    MATH_FUNC_DEFINE_VEC(cos,COS);
    MATH_FUNC_DEFINE_VEC(sin,SIN);
    MATH_FUNC_DEFINE_VEC(tan,TAN);
    MATH_FUNC_DEFINE_VEC(acos,ACOS);
    MATH_FUNC_DEFINE_VEC(asin,ASIN);
    MATH_FUNC_DEFINE_VEC(atan,ATAN);
    MATH_FUNC_DEFINE_VEC(exp,EXP);
    MATH_FUNC_DEFINE_VEC(log,LOG);
    MATH_FUNC_DEFINE_VEC(sqrt,SQRT);

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
    /// utulity-- calculate Jacobian Functor
    ///
    template<typename T>
    class VecFuncPtr {
    private:
        std::vector<FuncPtr<T>> func_ptr;
    public:
        VecFuncPtr(const std::vector<FuncPtr<T>>& _func_ptr) : func_ptr(_func_ptr) {}
        virtual std::vector<T> operator()(const std::vector<T>& x) const
        {
            std::vector<T> rtn;
            for(size_t i=0; i<func_ptr.size(); i++){
                rtn.push_back(std::move((*func_ptr[i])(x)));
            }
            return rtn;
        }
        virtual FuncPtr<T> operator[](const size_t idx) const
        {
            return func_ptr[idx];
        }
        virtual std::vector<FuncPtr<T>> to_vector() const
        {
            return func_ptr;
        }
    };

    template<typename T>
    VecFuncPtr<T> jacobian(const FuncPtr<T>& f, const size_t dim)
    {
        std::vector<FuncPtr<T>> rtn;
        for(size_t i=0; i<dim; i++){
            rtn.push_back(move((*f).derivative(i)));
        }
        return VecFuncPtr<T>(rtn);
    }


    ///
    /// utulity-- calculate Hessian Functor
    ///
    template<typename T>
    class MatFuncPtr {
    private:
        std::vector<std::vector<FuncPtr<T>>> func_ptr;
    public:
        MatFuncPtr(const std::vector<std::vector<FuncPtr<T>>>& _func_ptr) : func_ptr(_func_ptr) {}
        virtual std::vector<std::vector<T>> operator()(const std::vector<T>& x) const
        {
            size_t dim=func_ptr.size();
            std::vector<std::vector<T>> rtn(dim,std::vector<T>(dim));
            for(size_t i=0; i<dim; i++){
                for(size_t j=0; j<dim; j++){
                    rtn[i][j]=(*func_ptr[i][j])(x);
                }
            }
            return rtn;
        }
        virtual std::vector<FuncPtr<T>> operator[](const size_t idx) const
        {
            return func_ptr[idx];
        }

    };

    template<typename T>
    MatFuncPtr<T> hessian(const FuncPtr<T>& f, const size_t dim)
    {
        std::vector<std::vector<FuncPtr<T>>> rtn;
        auto jac=jacobian(f,dim);
        for(size_t i=0; i<dim; i++){
            rtn.push_back(std::move(jacobian(jac[i],dim).to_vector()));
        }
        return MatFuncPtr<T>(rtn);
    }

    ///
    /// utulity-- convert string for display
    ///
    template<typename T>
    std::string to_string(const std::vector<T>& vec)
    {
        std::ostringstream oss;
        for(auto e: vec){
            oss << " " << std::setw(8) << e << std::endl;
        }
        return oss.str();
    }

    template<typename T>
    std::string to_string(const std::vector<std::vector<T>>& mat)
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
}



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

#define MATH_FUNC_DEFINE(operant,operant_name) \
    template<typename T, int DIM> \
    FuncPtr<T,DIM> operant(const FuncPtr<T,DIM>& functor){ \
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::operant_name, functor)); \
    }


namespace AutomaticDifferentiation_Array {

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

    enum FuncType { SUM, DIFFERENCE, PRODUCT, QUOTIENT, MINUS, COS, SIN, TAN, ACOS, ASIN, ATAN, EXP, LOG, SQRT };

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

    MATH_FUNC_DEFINE(cos,COS);
    MATH_FUNC_DEFINE(sin,SIN);
    MATH_FUNC_DEFINE(tan,TAN);
    MATH_FUNC_DEFINE(acos,ACOS);
    MATH_FUNC_DEFINE(asin,ASIN);
    MATH_FUNC_DEFINE(atan,ATAN);
    MATH_FUNC_DEFINE(exp,EXP);
    MATH_FUNC_DEFINE(log,LOG);
    MATH_FUNC_DEFINE(sqrt,SQRT);

    template<typename T, int DIM>
    FuncPtr<T,DIM> operator-(const FuncPtr<T,DIM>& functor){
        return FuncPtr<T,DIM>(new Operator<T,DIM>(FuncType::MINUS, functor));
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

/// I recommend "AutomaticDifferentiation_Vector".
/// AutomaticDifferentiation_Array will be removed.
namespace AutomaticDifferentiation=AutomaticDifferentiation_Vector;

#endif // AUTOMATICDIFFERENTIATION_HPP_INCLUDED
