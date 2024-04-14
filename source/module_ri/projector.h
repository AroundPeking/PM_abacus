//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

template <typename TK, typename TR>
class projector
{
public:
    projector(){};
    ~projector()
    {
        delete Sk;
        //delete Tk;
    };

    void get_Sk(const LCAO_Matrix* LM_in);
    //void S2T();
    //void diagonalize_S(ESolver_KS_LCAO<TK, TR>& esolver);
private:
    //const LCAO_Matrix* LM_in;
    const std::vector<std::complex<double>>* Sk;
    //HContainer<TR>* Tk = nullptr;
};

#include "projector.hpp"
#endif // PROJECTOR_H
