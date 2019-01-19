#ifndef AUTOMATICDIFFERENTIATION_HPP_INCLUDED
#define AUTOMATICDIFFERENTIATION_HPP_INCLUDED

#include <eigen3/Eigen/Core>
#include <vector>
#include <memory>
#include <cmath>

#define OPERATOR_FUNC_DEFINE(operant,operant_name) \
    template<typename T> \
    FuncPtr<T> operator operant (const FuncPtr<T>& left, const FuncPtr<T>& right) \
    { \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, left, right)); \
    } \
    template<typename T> \
    FuncPtr<T> operator operant (const FuncPtr<T>& left, const T& right_val) \
    { \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, left, FuncPtr<T>(new Constant<T>(right_val)))); \
    } \
    template<typename T> \
    FuncPtr<T> operator operant (const T& left_val, const FuncPtr<T>& right) \
    { \
        return FuncPtr<T>(new Operator<T>(FuncType::operant_name, FuncPtr<T>(new Constant<T>(left_val)), right)); \
    }

#define MATH_FUNC_DEFINE(operant,operant_name) \
    template<typename T> \
    FuncPtr<T> operant(const FuncPtr<T>& functor) \
    { \
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
        template<typename TYPE>
        friend std::string toString(const FuncPtr<TYPE>& functor);
        template<typename TYPE>
        friend inline void simplification(FuncPtr<TYPE>& functor);
    public:
        Constant(const T& _val) : val(_val) {}
        virtual T operator()(const T* x) const { return val; }
        virtual FuncPtr<T> derivative(size_t idx) const { return FuncPtr<T>(new Constant(0.0)); }
    };

    enum FuncType { SUM, DIFFERENCE, PRODUCT, QUOTIENT, MINUS, COS, SIN, TAN, ACOS, ASIN, ATAN, EXP, LOG, SQRT, COSH, SINH, TANH };

    template<typename T>
    class Operator : public Functor<T> {
        FuncType func_type;
        FuncPtr<T> left;
        FuncPtr<T> right;
        template<typename TYPE>
        friend std::string toString(const FuncPtr<TYPE>& functor);
        template<typename TYPE>
        friend inline void simplification(FuncPtr<TYPE>& functor);
    public:
        Operator(FuncType _func_type, const FuncPtr<T>& _left, const FuncPtr<T>& _right=nullptr) : func_type(_func_type), left(_left), right(_right) {}
        virtual T operator()(const T* x) const;
        virtual FuncPtr<T> derivative(size_t idx) const;
    };

    template<typename T>
    class Variable : public Functor<T> {
    private:
        const size_t index;
        template<typename TYPE>
        friend std::string toString(const FuncPtr<TYPE>& functor);
        template<typename TYPE>
        friend inline void simplification(FuncPtr<TYPE>& functor);
    public:
        Variable(const size_t _index) : index(_index) {}
        virtual T operator()(const T* x) const { return x[index]; }
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
    MATH_FUNC_DEFINE(cosh,COSH);
    MATH_FUNC_DEFINE(sinh,SINH);
    MATH_FUNC_DEFINE(tanh,TANH);

    template<typename T>
    FuncPtr<T> operator-(const FuncPtr<T>& functor){ return FuncPtr<T>(new Operator<T>(FuncType::MINUS, functor)); }


    ///
    /// utulity-- show Functor
    ///
    template<typename T>
    std::string toString(const FuncPtr<T>& functor);

    ///
    /// Simplification of Functor
    ///
    template<typename T>
    void simplification(FuncPtr<T>& functor);

    ///
    /// utulity-- create variables
    ///
    template<typename T>
    std::vector<FuncPtr<T>> createVariables(size_t dim);

    ///
    /// utulity-- create zero constant value
    ///
    template<typename T>
    T zero(){ return T(0.0); }

    template<>
    inline FuncPtr<double> zero(){ return FuncPtr<double>(new Constant<double>(0.0)); }

    ///
    /// utulity-- create one constant value
    ///
    template<typename T>
    T one(){ return T(1.0); }

    template<>
    inline FuncPtr<double> one(){ return FuncPtr<double>(new Constant<double>(1.0)); }

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
        const FuncPtr<T>& operator()(const size_t row, const size_t col=0) const { return func_ptr[nCols*row+col]; }
        MatFuncPtr<T> operator*(const MatFuncPtr<T>& other) const;
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> operator()(const T* x) const;
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

template<typename T>
inline T AutomaticDifferentiation::Operator<T>::operator()(const T* x) const
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
    case FuncType::COSH:
        return std::cosh( (*left)(x) );
        break;
    case FuncType::SINH:
        return std::sinh( (*left)(x) );
        break;
    case FuncType::TANH:
        return std::tanh( (*left)(x) );
        break;
    default:
        throw std::string("Not defined operator in ")+__func__;
    }
}

template<typename T>
inline AutomaticDifferentiation::FuncPtr<T> AutomaticDifferentiation::Operator<T>::derivative(size_t idx) const
{
    FuncPtr<T> rtn;
    switch(func_type){
    case FuncType::SUM:
        rtn= (*left).derivative(idx) + (*right).derivative(idx);
        break;
    case FuncType::DIFFERENCE:
        rtn= (*left).derivative(idx) - (*right).derivative(idx);
        break;
    case FuncType::PRODUCT:
        rtn= (*left).derivative(idx) * right + left * (*right).derivative(idx);
        break;
    case FuncType::QUOTIENT:
        rtn= (*left).derivative(idx) / right - left * (*right).derivative(idx) / right / right;
        break;
    case FuncType::MINUS:
        rtn= -(*left).derivative(idx);
        break;
    case FuncType::COS:
        rtn= - (*left).derivative(idx) * sin(left);
        break;
    case FuncType::SIN:
        rtn= (*left).derivative(idx) * cos(left);
        break;
    case FuncType::TAN:
        rtn= (*left).derivative(idx) / cos(left) / cos(left);
        break;
    case FuncType::ACOS:
        rtn= - (*left).derivative(idx) / sqrt(1.0 - left*left);
        break;
    case FuncType::ASIN:
        rtn= (*left).derivative(idx)  / sqrt(1.0 - left*left);
        break;
    case FuncType::ATAN:
        rtn= (*left).derivative(idx) / (left*left + 1.0);
        break;
    case FuncType::EXP:
        rtn= (*left).derivative(idx) * exp(left);
        break;
    case FuncType::LOG:
        rtn= (*left).derivative(idx) / (left);
        break;
    case FuncType::SQRT:
        rtn= (*left).derivative(idx) / 2.0 / sqrt(left);
        break;
    case FuncType::COSH:
        rtn= (*left).derivative(idx) * sinh(left);
        break;
    case FuncType::SINH:
        rtn= (*left).derivative(idx) * cosh(left);
        break;
    case FuncType::TANH:
        rtn= (*left).derivative(idx) / cosh(left) / cosh(left);
        break;
    default:
        throw std::string("Not defined operator in ")+__func__;
    }
    simplification(rtn);
    return rtn;
}


///
/// utulity-- show Functor
///
template<typename T>
inline std::string AutomaticDifferentiation::toString(const FuncPtr<T>& functor)
{
    std::shared_ptr<Operator<T>> op=std::dynamic_pointer_cast<Operator<T>>(functor);
    std::shared_ptr<Variable<T>> va=std::dynamic_pointer_cast<Variable<T>>(functor);
    std::shared_ptr<Constant<T>> co=std::dynamic_pointer_cast<Constant<T>>(functor);
    if(op!=nullptr){
        switch(op->func_type){
        case FuncType::SUM:
            return "("+toString(op->left)+")+("+toString(op->right)+")";
            break;
        case FuncType::DIFFERENCE:
            return "("+toString(op->left)+")-("+toString(op->right)+")";
            break;
        case FuncType::PRODUCT:
            return "("+toString(op->left)+")*("+toString(op->right)+")";
            break;
        case FuncType::QUOTIENT:
            return "("+toString(op->left)+")/("+toString(op->right)+")";
            break;
        case FuncType::MINUS:
            return "-("+toString(op->left)+")";
            break;
        case FuncType::COS:
            return "cos("+toString(op->left)+")";
            break;
        case FuncType::SIN:
            return "sin("+toString(op->left)+")";
            break;
        case FuncType::TAN:
            return "tan("+toString(op->left)+")";
            break;
        case FuncType::ACOS:
            return "acos("+toString(op->left)+")";
            break;
        case FuncType::ASIN:
            return "asin("+toString(op->left)+")";
            break;
        case FuncType::ATAN:
            return "atan("+toString(op->left)+")";
            break;
        case FuncType::EXP:
            return "exp("+toString(op->left)+")";
            break;
        case FuncType::LOG:
            return "log("+toString(op->left)+")";
            break;
        case FuncType::SQRT:
            return "sqrt("+toString(op->left)+")";
            break;
        case FuncType::COSH:
            return "cosh("+toString(op->left)+")";
            break;
        case FuncType::SINH:
            return "sinh("+toString(op->left)+")";
            break;
        case FuncType::TANH:
            return "tanh("+toString(op->left)+")";
            break;
        default:
            throw std::string("Not defined operator in ")+__func__;
        }
    }
    if(va!=nullptr){
        return std::string("x_")+std::to_string(va->index);
    }
    if(co!=nullptr){
        return std::to_string(co->val);
    }
    return "ERROR";
}

///
/// Simplification of Functor
///
template<typename T>
inline void AutomaticDifferentiation::simplification(FuncPtr<T>& functor)
{
    std::shared_ptr<Operator<T>> op=std::dynamic_pointer_cast<Operator<T>>(functor);
    if(op!=nullptr){
        /// call recursively
        simplification(op->left);
        if(op->right!=nullptr) simplification(op->right);

        /// conbination of sum and minus -->> difference
        if(op->func_type==FuncType::SUM){
            std::shared_ptr<Operator<T>> right=std::dynamic_pointer_cast<Operator<T>>(op->right);
            if(right!=nullptr && right->func_type==FuncType::MINUS){
                op->func_type=FuncType::DIFFERENCE;
                op->right=right->left;
            }
        }

        /// conbination of difference and minus -->> sum
        if(op->func_type==FuncType::DIFFERENCE){
            std::shared_ptr<Operator<T>> right=std::dynamic_pointer_cast<Operator<T>>(op->right);
            if(right!=nullptr && right->func_type==FuncType::MINUS){
                op->func_type=FuncType::SUM;
                op->right=right->left;
            }
        }

        /// sum and difference: zero -->> other
        if(op->func_type==FuncType::SUM || op->func_type==FuncType::DIFFERENCE){
            std::shared_ptr<Constant<T>> left=std::dynamic_pointer_cast<Constant<T>>(op->left);
            std::shared_ptr<Constant<T>> right=std::dynamic_pointer_cast<Constant<T>>(op->right);
            if(left!=nullptr && left->val==T(0.0)){
                if(op->func_type==FuncType::SUM){
                    functor=op->right;
                }else{ /// DIFFERENCE
                    functor=FuncPtr<T>(new Operator<T>(FuncType::MINUS, op->right, nullptr));
                }
            }else if(right!=nullptr && right->val==T(0.0)){
                functor=op->left;
            }
        }

        /// product
        if(op->func_type==FuncType::PRODUCT){
            std::shared_ptr<Constant<T>> left=std::dynamic_pointer_cast<Constant<T>>(op->left);
            std::shared_ptr<Constant<T>> right=std::dynamic_pointer_cast<Constant<T>>(op->right);
            if(left!=nullptr){
                if(left->val==T(1.0)){
                    functor=op->right;
                }else if(left->val==T(0.0)){
                    functor=op->left; /// set zero
                }
            }else if(right!=nullptr){
                if(right->val==T(1.0)){
                    functor=op->left;
                }else if(right->val==T(0.0)){
                    functor=op->right; /// set zero
                }
            }
        }

        /// quotient
        if(op->func_type==FuncType::QUOTIENT){
            std::shared_ptr<Constant<T>> left=std::dynamic_pointer_cast<Constant<T>>(op->left);
            if(left!=nullptr && left->val==T(0.0)){
                functor=op->left; /// set zero
            }
        }
    }
}

///
/// utulity-- create variables
///
template<typename T>
inline std::vector<AutomaticDifferentiation::FuncPtr<T>> AutomaticDifferentiation::createVariables(size_t dim)
{
    std::vector<FuncPtr<T>> rtn;
    for(size_t i=0; i<dim; i++){
        rtn.emplace_back(new Variable<T>(i));
    }
    return rtn;
}

template<typename T>
inline Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> AutomaticDifferentiation::MatFuncPtr<T>::operator()(const T* x) const
{
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rtn(nRows,nCols);
    for(size_t i=0; i<nRows; i++){
        for(size_t j=0; j<nCols; j++){
            rtn(i,j)=(*(*this)(i,j))(x);
//            rtn(i,j)=(*func_ptr[nCols*i+j])(x);
        }
    }
    return rtn;
}

template<typename T>
inline AutomaticDifferentiation::MatFuncPtr<T> AutomaticDifferentiation::MatFuncPtr<T>::operator*(const AutomaticDifferentiation::MatFuncPtr<T>& other) const{
    assert(nCols==other.nRows);
    MatFuncPtr<T> rtn(nRows,other.nCols);
    for(size_t row=0; row<rtn.nRows; row++){
        for(size_t col=0; col<rtn.nCols; col++){
            rtn(row,col)=FuncPtr<T>(new Operator<T>(FuncType::PRODUCT,(*this)(row,0),other(0,col)));
            for(size_t k=1; k<nCols; k++){
                FuncPtr<T> tmp(new Operator<T>(FuncType::PRODUCT,(*this)(row,k),other(k,col)));
                rtn(row,col)=FuncPtr<T>(new Operator<T>(FuncType::SUM,rtn(row,col),tmp));
            }
            simplification(rtn(row,col));
        }
    }
    return rtn;
}

#endif // AUTOMATICDIFFERENTIATION_HPP_INCLUDED
