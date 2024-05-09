//=======================
// AUTHOR : Huanjing Gong
// DATE :   2024-3-26
//=======================

#ifndef PROJECTOR_HPP
#define PROJECTOR_HPP

#include "projector.h"

#include <iostream>
#include <iomanip>

template <typename TK, typename TR>
void projector<TK, TR>::get_SR(LCAO_Hamilt& UHM)
{
    //std::cout << "Overlap Matrix (SR):" << std::endl;

    auto &SR_sparse_ptr = UHM.genH.LM->SR_sparse;
    auto &all_R_coor_ptr = UHM.genH.LM->all_R_coor;
    double *line = nullptr;
    int local_size = GlobalV::NLOCAL * GlobalV::NLOCAL;
    std::vector<double> s;
    std::vector<double> I;
    std::vector<double> s_local;
    std::vector<double>	I_local;
    int ISRC = 0, info;
    int nprows, npcols, myprow, mypcol, icontxt;
    int nlocal, nbands, nb2d;

    s.resize(local_size);
    I.resize(local_size);
    s_local.resize(local_size);	
    I_local.resize(local_size);
    nlocal = GlobalV::NLOCAL;
    nbands = GlobalV::NBANDS;
    nb2d = GlobalV::NB2D;
    MPI_Bcast(&nlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    LCAO_DIAGO_TEST::process_2d(nprows, npcols, myprow, mypcol, icontxt);
    this->SR.nrow = LCAO_DIAGO_TEST::na_rc(nlocal, nb2d, nprows, myprow); // the number of row of the new_matrix in each process
    this->SR.ncol = LCAO_DIAGO_TEST::na_rc(nlocal, nb2d, npcols, mypcol); // the number of column of the new_matrix in each process

    this->eigenvector.resize(1, this->SR.ncol, this->SR.nrow);
    this->eigenvalue.resize(nlocal, 0.0);
    descinit_(this->SR.desc, &nlocal, &nlocal, &nb2d, &nb2d, &ISRC, &ISRC, &icontxt, &(this->SR.nrow), &info);
    if (info != 0)
    {
	    printf("Invalid blacs-distribution. Abort!\n");
            exit(1);
    }

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
                			s[row * GlobalV::NLOCAL + col] = line[col];
							I[row * GlobalV::NLOCAL + col] = 0.0;
							if(col == row)
							{
								I[row * GlobalV::NLOCAL + col] = 1.0;
							}
							//std::cout << std::fixed << std::setprecision(6) << s[row * GlobalV::NLOCAL + col] << ",";
            			}
						//std::cout << std::endl;
        		}

            }
			LCAO_DIAGO_TEST::distribute_data<double>(s.data(),s_local.data(),nlocal,nb2d,this->SR.nrow,this->SR.ncol,icontxt);
			LCAO_DIAGO_TEST::distribute_data<double>(I.data(),I_local.data(),nlocal,nb2d,this->SR.nrow,this->SR.ncol,icontxt);
			MPI_Barrier(MPI_COMM_WORLD);
			this->SR.s_local = I_local;
			this->SR.h_local = s_local;

            delete[] line;
            line = nullptr;
        } 
    }
    //ModuleBase::QUIT();
}

template <typename TK, typename TR>
void projector<TK, TR>::diag_S()
{
	hsolver::DiagH<double>* dh = 0;
    dh = new hsolver::DiagoElpa<double>;

	//std::cout << "Before diag..." << std::endl;
	dh->diag(&(this->SR), this->eigenvector, eigenvalue.data());
	std::cout << "row eigenvalues: " << std::endl;
	for (int i = 0; i < GlobalV::NLOCAL; ++i)
	{
		std::cout << this->eigenvalue[i] << std::endl;
	}
	//std::cout << "eigenvectors: " << std::endl;
	//for (int ir = 0; ir < eigenvector.get_nbasis(); ir++)
	//{
	//	std::cout << eigenvector(0, 25, ir) << std::endl;
	//}
	delete dh;
}

template <typename TK, typename TR>
std::vector<double> projector<TK, TR>::S2T(const double pm_epl)
{
	int nlocal = GlobalV::NLOCAL;
	std::vector<double> Dm, tmp;
	Dm.resize(nlocal * nlocal);
	tmp.resize(nlocal * nlocal, 0.0);
	this->TR.resize(nlocal * nlocal, 0.0);

	for (int i = 0; i < nlocal; ++i)
	{
		if(this->eigenvalue[i] < pm_epl)
		{
			this->eigenvalue[i] = 0.0;
		}
		else
		{
			this->eigenvalue[i] = 1.0;
		}
	}
	std::cout << "modified eigenvalues: " << std::endl;
	for (int i = 0; i < GlobalV::NLOCAL; ++i)
	{
		std::cout << this->eigenvalue[i] << std::endl;
	}
	for(int row = 0; row < nlocal; ++row)
    {	
    	for (int col = 0; col < nlocal; ++col)
    	{
			Dm[row * nlocal + col] = 0.0;
			if(col == row)
			{
				Dm[row * nlocal + col] = this->eigenvalue[row];
			}
			//std::cout<<Dm[row * nlocal + col]<<",";
    	}
	}
	
	for(int row = 0; row < nlocal; ++row)
    {	
    	for(int col = 0; col < nlocal; ++col)
    	{
			for(int ir = 0; ir < nlocal; ++ir)
			{
				for(int ic = 0; ic < nlocal; ++ic)
				{
					tmp[row * nlocal + ir] += this->eigenvector(0, ic, row) * Dm[ic * nlocal + ir];
				}
				this->TR[row * nlocal + col] +=  tmp[row * nlocal + ir] * this->eigenvector(0, ir, col);
				tmp[row * nlocal + ir] = 0.0;
				// for(int ic = 0; ic < nlocal; ++ic)
				// {
				// 	tmp[row * nlocal + ir] += this->eigenvector(0, row, ic) * Dm[ic * nlocal + ir];
				// }
				// this->TR[row * nlocal + col] +=  tmp[row * nlocal + ir] * this->eigenvector(0, col, ir);
				// tmp[row * nlocal + ir] = 0.0;
			}
			std::cout << std::fixed << std::setprecision(6) << this->TR[row * nlocal + col] << ",";
    	}
		std::cout << std::endl;
	}
	return this->TR;
}

#endif
