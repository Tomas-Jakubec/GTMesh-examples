#ifndef RKMSOLVER_HPP
#define RKMSOLVER_HPP

#include <GTMesh/UnstructuredMesh/UnstructuredMesh.h>

template <typename Functor, typename ...T, unsigned int Dimension>
void performVectorOperation(Functor&& f, MeshDataContainer<T, Dimension>& ...args) {
    const auto& firstVector = std::get<0>(std::forward_as_tuple(args...));
    for (std::size_t i = 0; i < firstVector.template getDataByPos<0>().size(); ++i) {
        f((args.template getDataByPos<0>()[i])...);
    }
}


template <typename Functor, typename ...T, unsigned int Dimension>
double performVectorReductionMax(Functor&& f, MeshDataContainer<T, Dimension>& ...args) {
    const auto& firstVector = std::get<0>(std::forward_as_tuple(args...));
    double res = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < firstVector.template getDataByPos<0>().size(); ++i) {
        double tmp_res = f((args.template getDataByPos<0>()[i])...);
        if (res < tmp_res) {
            res = tmp_res;
        }
    }
    return res;
}

template <typename Problem, std::enable_if_t<HasDefaultArithmeticTraits<typename Problem::ResultType>::value, bool> = true>
void RKMSolver(Problem& problem,
               MeshDataContainer<typename Problem::ResultType, Problem::MeshType::meshDimension()>& compData,//x_ini
               double tau_ini, double startTime, double finalT, double delta)
{
    constexpr unsigned int MeshDimension = Problem::MeshType::meshDimension();
    using container_type = MeshDataContainer<typename Problem::ResultType, MeshDimension>;
    container_type Ktemp(compData);
    container_type K1(compData), K2(compData), K3(compData), K5(compData), K4(compData);

    double tau = tau_ini;
    double time = startTime;
    bool run = true;

    while (time < finalT) {
        if (time + tau > finalT) {
            tau = finalT - time;
        }
        problem.calculateRHS(time, compData, K1);

        performVectorOperation([&tau](auto &Ktemp, auto& compData, auto& K1){
            Ktemp = compData + (tau * (1.0 / 3.0) * K1);},
        Ktemp, compData, K1);

        problem.calculateRHS(time, Ktemp, K2);

        performVectorOperation([&tau](auto &Ktemp, auto& compData, auto& K1, auto& K2){
            Ktemp = compData + (tau * (1.0 / 6.0) * (K1 + K2));},
        Ktemp, compData, K1, K2);

        problem.calculateRHS(time, Ktemp, K3);

        performVectorOperation([&tau](auto &Ktemp, auto& compData, auto& K1, auto& K3){
        Ktemp = compData + (tau * (0.125 * K1 + 0.375 * K3));
        },Ktemp, compData, K1, K3);

        problem.calculateRHS(time, Ktemp, K4);

        performVectorOperation([&tau](auto &Ktemp, auto& compData, auto& K1, auto& K3, auto& K4){
        Ktemp = compData + (tau * ((0.5 * K1) - (1.5 * K3) + (2.0 * K4)));},
        Ktemp, compData, K1, K3, K4);


        problem.calculateRHS(time, Ktemp, K5);
        double error = performVectorReductionMax([](auto& K1, auto& K3, auto& K4, auto& K5)->double{
        return max(abs(0.2 * K1 - 0.9 * K3 + 0.8 * K4 - 0.1 * K5));},
        K1, K3, K4, K5);

        error *= tau * (1.0 / 3.0);

        if (error < delta) {
            performVectorOperation([&tau](auto& compData, auto& K1, auto& K4, auto& K5){
                compData += tau * (1.0 / 6.0) * (((K1 + K5)) + (4.0 * K4));},
                                   compData, K1, K4, K5);
            time += tau;
            if (error == 0.0) continue;
        }
        tau *= std::pow(delta/error, 0.2) * 0.8;
    }
}


#endif // RKMSOLVER_HPP
