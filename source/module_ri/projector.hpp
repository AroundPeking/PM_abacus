//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_HPP
#define PROJECTOR_HPP

#include "projector.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
//#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/overlap_new.h"
#include "module_hsolver/hsolver_lcao.h"

#include <iostream>

template <typename TK, typename TR>
void projector<TK, TR>::get_Sk(const LCAO_Matrix* LM_in)
{
    //this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(LM_in, kv_in);

    this->Sk = &LM_in->Sloc2;

    int N = 26;

    std::cout << "Overlap Matrix (Sk):" << std::endl;
    std::cout << "dim:" << std::endl;
    std::cout << this->Sk->size() << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << (*this->Sk)[i * N + j] << " ";
        }
        std::cout << std::endl;
    }
    ModuleBase::QUIT();
}
#endif
