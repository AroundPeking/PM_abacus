//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_basis/module_ao/ORB_control.h"

template <typename TK, typename TR>
class projector
{
public:
    projector(){};
    ~projector()
    {
        delete SR;
        //delete Tk;
    };

    void get_SR(LCAO_Hamilt& UHM, const K_Vectors& kv_in);
    //void S2T();
    //void diagonalize_S(ESolver_KS_LCAO<TK, TR>& esolver);
private:
    hamilt::HContainer<TR>* SR = nullptr;
    ORB_control orb_con;    //Basis_LCAO
    Record_adj RA;
    //HContainer<TR>* Tk = nullptr;
};

#include "projector.hpp"
#endif // PROJECTOR_H
