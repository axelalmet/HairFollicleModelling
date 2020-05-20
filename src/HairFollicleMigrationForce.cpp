/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "HairFollicleMigrationForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
HairFollicleMigrationForce<DIM>::HairFollicleMigrationForce()
    : AbstractForce<DIM>(),
    mMigrationForceStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
HairFollicleMigrationForce<DIM>::~HairFollicleMigrationForce()
{
}

// Method to get the strength of the migration force
template<unsigned DIM>
double HairFollicleMigrationForce<DIM>::GetMigrationForceStrength()
{
    return mMigrationForceStrength;
}

// Method to set the strength of migration force
template<unsigned DIM>
void HairFollicleMigrationForce<DIM>::SetMigrationForceStrength(double migrationForceStrength)
{
    mMigrationForceStrength = migrationForceStrength;
}

template<unsigned DIM>
void HairFollicleMigrationForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // Get the force parameters
    double migration_force_strength = GetMigrationForceStrength(); // Strength of migration force

    // Initialise the migration force vector
    c_vector<double, DIM> migration_force;
    migration_force[0] = 0.0;
    migration_force[1] = -migration_force_strength;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Only apply this migration force to fibroblasts
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType(); // Get the cell type

        if (p_cell_type->IsType<TransitCellProliferativeType>())
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(migration_force); 
        }
        // else if (p_cell_type->IsType<DifferentiatedCellProliferativeType>())
        // {
        //     // Quick hack to kick the diff cells out of the HF base
        //     double y = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

        //     if (y < 4.0)
        //     {
        //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        //         rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(-1.0 * migration_force); 
        //     }
        // }
    }
}

template<unsigned DIM>
void HairFollicleMigrationForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<MigrationForceStrength>"<<  mMigrationForceStrength << "</MigrationForceStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class HairFollicleMigrationForce<1>;
template class HairFollicleMigrationForce<2>;
template class HairFollicleMigrationForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HairFollicleMigrationForce)
