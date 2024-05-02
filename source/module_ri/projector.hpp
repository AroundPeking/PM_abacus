//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_HPP
#define PROJECTOR_HPP

#include "projector.h"

#include <iostream>

template <typename TK, typename TR>
void projector<TK, TR>::get_SR(LCAO_Hamilt& UHM, const K_Vectors& kv_in)
{
    //std::cout << "Overlap Matrix (SR):" << std::endl;

    auto &SR_sparse_ptr = UHM.genH.LM->SR_sparse;
    auto &all_R_coor_ptr = UHM.genH.LM->all_R_coor;
    double *line = nullptr;

    for (auto &R_coor : all_R_coor_ptr)
    {
	    if ((R_coor.x == 0) && (R_coor.y == 0) && (R_coor.z == 0))
	    {
		    //std::cout << "pm" << std::endl;
        	//std::cout << R_coor.x << R_coor.y << R_coor.z << std::endl;
        	auto sR = SR_sparse_ptr[R_coor];
        	line = new double[GlobalV::NLOCAL];
        	for(int row = 0; row < GlobalV::NLOCAL; ++row)
    		{
        		ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);
        		auto iter = sR.find(row);
        		if (iter != sR.end())
        		{
            			for (auto &value : iter->second)
            			{
                			line[value.first] = value.second;
            			}
        		}

        		Parallel_Reduce::reduce_all(line, GlobalV::NLOCAL);

        		if(GlobalV::DRANK == 0)
        		{
            			for (int col = 0; col < GlobalV::NLOCAL; ++col)
            			{
                			std::cout << " " << std::fixed << std::scientific << std::setprecision(8) << line[col];
            			}
        		}
			    std::cout << std::endl;

            }

            delete[] line;
            line = nullptr;
        } 
    }
    //ModuleBase::QUIT();
}
#endif
