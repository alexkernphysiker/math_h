// this file is distributed under
// LGPLv3 license
#ifndef _____EXTENDED_UNCERTAINTIES_ESTIMATION_____________
#define _____EXTENDED_UNCERTAINTIES_ESTIMATION_____________
#include <math.h>
#include "error.h"
#include "sigma.h"
#include "randomfunc.h"
namespace MathTemplates
{
    template<size_t d, class n = double>class value_f;
    template<class n = double>class value_f_chain;
    template<typename numt = double>
    class abstract_value_with_uncertainty_numeric :public abstract_value_with_uncertainty<numt>
    {
    public:
        typedef numt NumberType;
        typedef numt NumberType_ForNumeric;//this is needed to differ template operators declared here from ones declared for value<>
        template<size_t d, class n>friend class value_f;
        template<class n>friend class value_f_chain;
    protected:
        virtual numt mmeasure()const = 0;
    private:
        mutable bool m_cache_valid;
        mutable StandardDeviation<numt> m_cache;
        static size_t mm_number;
        const StandardDeviation<numt>& VAL()const
        {
            if (!m_cache_valid) {
                m_cache = StandardDeviation<numt>();
                for (size_t i = 0;i < mm_number;i++)m_cache << (mmeasure());
                m_cache_valid = true;
            }
            return m_cache;
        }
    public:
        abstract_value_with_uncertainty_numeric() :m_cache_valid(false) {}
        static inline size_t micro_measurements_count() { return abstract_value_with_uncertainty_numeric::mm_number; }
        static inline void set_micro_measurements_count(size_t v) { abstract_value_with_uncertainty_numeric::mm_number = v; }
        virtual ~abstract_value_with_uncertainty_numeric() {}
        virtual const numt& val()const override
        {
            return VAL().val();
        }
        virtual const numt& uncertainty()const override
        {
            return VAL().uncertainty();
        }
    };
    template<typename numt>
    size_t abstract_value_with_uncertainty_numeric<numt>::mm_number = 10000;
    //measured values
    template<class numt = double>
    class value_numeric_const :public abstract_value_with_uncertainty_numeric<numt> {
    private:
        numt m_val;
    protected:
        virtual numt mmeasure()const override { return m_val; }
    public:
        value_numeric_const(const value_numeric_const& source) = default;
        value_numeric_const(const numt& v) :abstract_value_with_uncertainty_numeric<numt>(), m_val(v) {}
        virtual ~value_numeric_const() {}
    };
    template<class numt = double, class DISTR = RandomGauss<numt>>
    class value_numeric_distr :public abstract_value_with_uncertainty_numeric<numt> {
    private:
        DISTR m_distr;
    protected:
        virtual numt mmeasure()const override {
            return m_distr();
        }
    public:
        value_numeric_distr(const value_numeric_distr& S) :abstract_value_with_uncertainty_numeric<numt>(), m_distr(S.m_distr) {}
        value_numeric_distr(const DISTR& D) :abstract_value_with_uncertainty_numeric<numt>(), m_distr(D) {}
        template<typename... Args>
        value_numeric_distr(Args... args) : abstract_value_with_uncertainty_numeric<numt>(), m_distr(args...) {}
        virtual ~value_numeric_distr() {}
    };
    template<class numt>
    inline value_numeric_distr<numt, RandomGauss<numt>> value_numeric(const abstract_value_with_uncertainty<numt>& a) {
        return value_numeric_distr<numt, RandomGauss<numt>>(a.val(), a.uncertainty());
    }
    //calculated values
    template<class numt>
    class value_f_chain :public abstract_value_with_uncertainty_numeric<numt> {
    private:
        std::function<numt(const Chain<numt>&)> m_func;
        Chain<const abstract_value_with_uncertainty_numeric<numt>*> m_chain;
    protected:
        virtual numt mmeasure()const override {
            Chain<numt> P;
            for (const auto A : m_chain)P.push_back(A->mmeasure());
            return m_func(P);
        }
    public:
        value_f_chain(const value_f_chain& source) = default;
        template<class VType>
        value_f_chain(std::function<numt(const Chain<numt>&)> F, const Chain<VType>& chain)
            :abstract_value_with_uncertainty_numeric<numt>(), m_func(F) {
            for (const auto& a : chain)m_chain.push_back(&a);
        }
        virtual ~value_f_chain() {}
    };
    template<class numt, class VType>
    inline value_f_chain<numt> FUNC(std::function<numt(const Chain<numt>&)> F, const Chain<VType>& chain) {
        return value_f_chain<numt>(F, chain);
    }
    template<class numt>
    class value_f<0, numt> :public abstract_value_with_uncertainty_numeric<numt> {
    private:
        std::function<numt()> m_func;
    protected:
        virtual numt mmeasure()const override { return m_func(); }
    public:
        value_f(const value_f& s) = default;
        inline value_f(std::function<numt()> F)
            :abstract_value_with_uncertainty_numeric<numt>(), m_func(F) {}
        virtual ~value_f() {}
    };


    template<class numt>
    class value_f<1, numt> :public value_f<0, numt> {
    public:
        value_f(const value_f& s) = default;
        template<class ArgType1>
        inline value_f(std::function<numt(const numt&)> F, ArgType1 a)
            :value_f<0, numt>([a, F]() {
            const abstract_value_with_uncertainty_numeric<numt>* A = &a;
            return F(A->mmeasure());
                }) {
        }
        virtual ~value_f() {}
    };
    template<class numt>
    class value_f<2, numt> :public value_f<1, numt> {
    public:
        value_f(const value_f& s) = default;
        template<class ArgType1, class ArgType2>
        inline value_f(std::function<numt(const numt&, const numt&)> F, ArgType1 a, ArgType2 b)
            :value_f<1, numt>([a, F](const numt& y) {
            const abstract_value_with_uncertainty_numeric<numt>* A = &a;
            return F(A->mmeasure(), y);
                }, b) {
        }
        virtual ~value_f() {}
    };


    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> SQRT(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return sqrt(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> SQR(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return x * x;}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> EXP(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return exp(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> LOG(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return log(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> SIN(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return sin(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> COS(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return cos(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> TAN(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return tan(x);}, a);
    }
    template<class ArgType1>
    inline value_f<1, typename ArgType1::NumberType_ForNumeric> ATAN(ArgType1 a) {
        return value_f<1, typename ArgType1::NumberType_ForNumeric>([](const typename ArgType1::NumberType_ForNumeric& x) {return atan(x);}, a);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> ATAN2(ArgType1 a, ArgType2 b) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return atan2(x, y);}, a, b);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> POW(ArgType1 a, ArgType2 b
    ) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return pow(x, y);}, a, b);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> operator+(ArgType1 a, ArgType2 b) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return x + y;}, a, b);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> operator-(ArgType1 a, ArgType2 b) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return x - y;}, a, b);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> operator*(ArgType1 a, ArgType2 b) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return x * y;}, a, b);
    }
    template<class ArgType1, class ArgType2>
    inline value_f<2, typename ArgType1::NumberType_ForNumeric> operator/(ArgType1 a, ArgType2 b) {
        return value_f<2, typename ArgType1::NumberType_ForNumeric>([](
            const typename ArgType1::NumberType_ForNumeric& x, const typename ArgType1::NumberType_ForNumeric& y
            ) {return x / y;}, a, b);
    }

    template<size_t dim, class numt>
    class value_f :public value_f<dim - 1, numt> {
    public:
        value_f(const value_f& s) :value_f<dim - 1, numt>(s) {}
        template<class FUNC, class ArgType, typename... Args>
        inline value_f(FUNC F, ArgType a, Args... args)
            : value_f<dim - 1, numt>([a, F](auto...z) {
            const abstract_value_with_uncertainty_numeric<numt>* A = &a;
            return F(A->mmeasure(), z...);
                }, args...) {
        }
        virtual ~value_f() {}
    };
    template<class FF, class ArgType, typename... Args>
    inline value_f<sizeof...(Args) + 1, typename ArgType::NumberType_ForNumeric> FUNC(FF F, ArgType a, Args... args) {
        return value_f<sizeof...(Args) + 1, typename ArgType::NumberType_ForNumeric>(F, a, args...);
    }

};
#endif
