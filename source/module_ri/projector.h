//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_H
#define PROJECTOR_H

//get_SR
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"

//diagonalize
#include "diago_lcao_test.h"
//#include "module_hsolver/test/diago_elpa_utils.h"
#include "module_hamilt_general/matrixblock.h"
#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#endif


// TK is SR
template <typename TK, typename TR>
class projector
{
public:
    projector(){};
    ~projector()
    {
        //delete SR;
        //delete Tk;
    };

    void get_SR(LCAO_Hamilt& UHM);
    void diag_S();
    std::vector<double> S2T(const double pm_epl);
private:
    HamiltTEST<double> SR;
    std::vector<double> TR;
    std::vector<double> eigenvalue;
    psi::Psi<double> eigenvector;
};

#include "projector.hpp"
#endif // PROJECTOR_H
