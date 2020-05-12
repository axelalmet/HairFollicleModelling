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

#include "HairFollicleDifferentiationTrackingModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "MovableStemCellProgenyProliferativeType.hpp"
#include "NonMovableStemCellProgenyProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
HairFollicleDifferentiationTrackingModifier<DIM>::HairFollicleDifferentiationTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mNicheBulgeRadius(DOUBLE_UNSET),
      mNicheBulgeCentre(zero_vector<double>(DIM)),
      mMovableProgenyDifferentiationProbability(DOUBLE_UNSET),
      mBaseRadius(DOUBLE_UNSET),
      mTransitDifferentiationHeight(DOUBLE_UNSET)
{
}

template<unsigned DIM>
HairFollicleDifferentiationTrackingModifier<DIM>::~HairFollicleDifferentiationTrackingModifier()
{
}

/*
* Get and set for the size of the stem cell bulge
*/
template<unsigned DIM>
double HairFollicleDifferentiationTrackingModifier<DIM>::GetNicheBulgeRadius()
{
    return mNicheBulgeRadius;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetNicheBulgeRadius(double nicheBulgeRadius)
{
    mNicheBulgeRadius = nicheBulgeRadius;
}

/*
* Get and set for location of the centre of the stem cell bulge
*/
template<unsigned DIM>
c_vector<double, DIM> HairFollicleDifferentiationTrackingModifier<DIM>::GetNicheBulgeCentre()
{
    return mNicheBulgeCentre;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetNicheBulgeCentre(c_vector<double, DIM> nicheBulgeCentre)
{
    mNicheBulgeCentre = nicheBulgeCentre;
}

/*
* Get and set for probability that stem cell differentates into movable progeny.
*/
template<unsigned DIM>
double HairFollicleDifferentiationTrackingModifier<DIM>::GetMovableProgenyDifferentiationProbability()
{
    return mMovableProgenyDifferentiationProbability;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetMovableProgenyDifferentiationProbability(double movableProgenyDifferentiationProbability)
{
    mMovableProgenyDifferentiationProbability = movableProgenyDifferentiationProbability;
}

/*
* Get and set for HF base radius
*/
template<unsigned DIM>
double HairFollicleDifferentiationTrackingModifier<DIM>::GetBaseRadius()
{
    return mBaseRadius;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetBaseRadius(double baseRadius)
{
    mBaseRadius = baseRadius;
}

/**
 * Get the height at which TA cells differentiate
 */
template<unsigned DIM>
double HairFollicleDifferentiationTrackingModifier<DIM>::GetTransitDifferentiationHeight()
{
    return mTransitDifferentiationHeight;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetTransitDifferentiationHeight(double transitDifferentiationHeight)
{
    mTransitDifferentiationHeight = transitDifferentiationHeight;
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // First check the current attachments to the basement membrane, which depends on the cell type.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only consider stem cells
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        // Get the node index
        c_vector<double, DIM> current_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        
        if (p_cell_type->template IsType<StemCellProliferativeType>()) // For stem cells
        {          
            if (current_location[0] >= mNicheBulgeCentre[0]) // Only celsl that are within the hair follicle should differentiate
            {
                if (norm_2(current_location - mNicheBulgeCentre) > mNicheBulgeRadius) // If the cell has fallen outside of the bulge, it will differentiate
                {
                    double progeny_fate = RandomNumberGenerator::Instance()->ranf();

                    boost::shared_ptr<AbstractCellProperty> p_movable_type =
                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<MovableStemCellProgenyProliferativeType>();
                    boost::shared_ptr<AbstractCellProperty> p_non_movable_type =
                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<NonMovableStemCellProgenyProliferativeType>();

                    if (progeny_fate < mMovableProgenyDifferentiationProbability)
                    {
                        cell_iter->SetCellProliferativeType(p_movable_type);
                    }
                    else
                    {
                        cell_iter->SetCellProliferativeType(p_non_movable_type);
                    }
                }
            }
        }
        // Also track the progenies for differentiation
        else if ( p_cell_type->template IsType<MovableStemCellProgenyProliferativeType>())
        {
            // If progeny cells are in the lower half of the base and fall out of the annulus of progeny cells.
            if (current_location[1] < 0.0)
            {
                if (norm_2(current_location) < mBaseRadius - 1.5)
                {
                    boost::shared_ptr<AbstractCellProperty> p_transit_type =
                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<TransitCellProliferativeType>();
                    cell_iter->SetCellProliferativeType(p_transit_type);
                }
            }
        }
        else if (p_cell_type->template IsType<TransitCellProliferativeType>() ) // For TA cells, if they've gone past a certain height, let's differentiate them.
        {   

            if (current_location[1] > mTransitDifferentiationHeight)
            {
                boost::shared_ptr<AbstractCellProperty> p_diff_type =
                cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<DifferentiatedCellProliferativeType>();
                cell_iter->SetCellProliferativeType(p_diff_type);
            }
        }
        
        

    }
}

template<unsigned DIM>
void HairFollicleDifferentiationTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class HairFollicleDifferentiationTrackingModifier<1>;
template class HairFollicleDifferentiationTrackingModifier<2>;
template class HairFollicleDifferentiationTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HairFollicleDifferentiationTrackingModifier)
