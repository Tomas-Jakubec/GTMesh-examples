#pragma once

#include <lib/NumericHelpers.hpp>


namespace optimisation {

    
    template<unsigned int Dim, typename F>
    Vertex<Dim> steepestGradient(Vertex<Dim> x, F&& f, double treshold, double maxDist = 0.01) {
        auto gradF = gradientOf(f, maxDist * 1e-3);
        auto x_i = x;
        double fx_i = f(x);
        auto g_i = -gradF(x_i);
        int cnt = 0;
        while(++cnt > 2000 && g_i.normEuclid() > treshold){
            Vector<Dim> df = -gradF(x_i);
            double gamma = 1.0;//maxDist / (df).normEuclid();
            double fx_i1 = 0;
            Vertex<Dim> x_i1;
            // line search
            do {
                x_i1 = x_i + gamma * df;
                fx_i1 = f(x_i1);
                gamma *= 0.7;
            } while (fx_i1 > fx_i);
            fx_i = fx_i1;
            x_i = x_i1;
            g_i = -gradF(x_i);
        }
        return x_i;
    }

    template<unsigned int Dim, typename F>
    Vertex<Dim> DFP(Vertex<Dim> x, F&& f, double maxDist = 0.01) {
        auto gradF = gradientOf(f, maxDist * 1e-3);

        auto x_i = x;

        double fx_i = f(x);
        auto g_i = gradF(x_i);
        auto S = unitMatrix<Dim>();
        int cnt = 0;
        while(++cnt > 100 && g_i.normEuclid() > 10e-5){
            auto d = -1.0 * multiply(S,  g_i);
            DBGVAR(S, d, g_i);

            double gamma = 1.0;
            double fx_i1 = 0;
            Vertex<Dim> x_i1;
            // line search
            do {
                x_i1 = x_i + gamma * d;
                fx_i1 = f(x_i1);
                gamma *= 0.7;
            } while (fx_i1 > fx_i);
            fx_i = fx_i1;
            auto g_i1 = gradF(x_i1);
            auto p = x_i1 - x_i;
            auto q = g_i1 - g_i;
            x_i = x_i1;
            g_i = g_i1;
            auto p_q = multiply(p,q);
            if (abs(p_q) < 1e-30) {
                return x_i;
            }
            auto S_q = multiply(S,q);
            S -= tenzorMultiply(S_q, S_q) / multiply(q, S_q) - tenzorMultiply(p,p) / multiply(p, p);
        }
        return x_i;
    }

    template<unsigned int Dim, typename F>
    Vertex<Dim> BFGS(Vertex<Dim> x, F&& f, double maxDist = 0.01, const double treshold = 1e-5) {
        auto gradF = gradientOf(f, maxDist * 1e-3);

        auto x_i = x;

        double fx_i = f(x);
        Vector<Dim> g_i = gradF(x_i);
        auto S = unitMatrix<Dim>();
        int cnt = 0;
        while(++cnt < 100 && g_i.normEuclid() > treshold){
            auto d = -1.0 * multiply(S,  g_i);

            double gamma = 1.0;//maxDist / d.normEuclid();
            double fx_i1 = 0;
            
            Vertex<Dim> x_i1;
            // line search
            do {
                x_i1 = x_i + gamma * d;
                fx_i1 = f(x_i1);
                gamma *= 0.7;
            } while (fx_i1 > fx_i);
            fx_i = fx_i1;
            auto g_i1 = gradF(x_i1);
            auto p = x_i1 - x_i;
            auto q = g_i1 - g_i;
            x_i = x_i1;
            g_i = g_i1;
            auto p_q = multiply(p,q);
            if (abs(p_q) < 1e-30) {
                return x_i;
            }
            auto S_q = multiply(S,q);
            auto deltaS = (tenzorMultiply(p, S_q) + tenzorMultiply(S_q, p) - (1 + (multiply(q, S_q)/p_q)) * tenzorMultiply(p,p)) / p_q;
            S -= deltaS;
        }
        return x_i;
    }

}