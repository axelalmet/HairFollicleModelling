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

#include "HairFollicleUniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "MovableStemCellProgenyProliferativeType.hpp"
#include "NonMovableStemCellProgenyProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "DermalSheathCellProliferativeType.hpp"
#include "Debug.hpp"

HairFollicleUniformG1GenerationalCellCycleModel::HairFollicleUniformG1GenerationalCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mGeneration(0),
      mMaxTransitGenerations(3),
      mQuiescentVolumeFraction(0.0),
      mEquilibriumVolume(0.25*M_PI),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0)
{
}

HairFollicleUniformG1GenerationalCellCycleModel::HairFollicleUniformG1GenerationalCellCycleModel(const HairFollicleUniformG1GenerationalCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel),
    mGeneration(rModel.mGeneration),
     mMaxTransitGenerations(rModel.mMaxTransitGenerations),
     mQuiescentVolumeFraction(rModel.mQuiescentVolumeFraction),
     mEquilibriumVolume(rModel.mEquilibriumVolume),
     mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
     mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* HairFollicleUniformG1GenerationalCellCycleModel::CreateCellCycleModel()
{
    return new HairFollicleUniformG1GenerationalCellCycleModel(*this);
}

void HairFollicleUniformG1GenerationalCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 4*p_gen->ranf(); // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = GetTransitCellG1Duration() + 2*p_gen->ranf(); // U[4,6] for default parameters (mTransitG1CellDuration) according to Meineke
    }
    // Doesn't proliferate for diff cells, stem progenies or fibroblasts
    else 
    {
        mG1Duration = DBL_MAX;
    }

}

void HairFollicleUniformG1GenerationalCellCycleModel::UpdateCellCyclePhase()
{

    if ((mQuiescentVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
    {
        EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
    }

    // Get cell volume
    double cell_volume = mpCell->GetCellData()->GetItem("volume");

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell volume
        double dt = SimulationTime::Instance()->GetTimeStep();
        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

        if (cell_volume < quiescent_volume)
        {
            // Update the duration of the current period of contact inhibition.
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
            mG1Duration += dt;
        }
        else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
        }
    }

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if ((!mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        &&(!mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>()) )
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }

}

void HairFollicleUniformG1GenerationalCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double HairFollicleUniformG1GenerationalCellCycleModel::GetQuiescentVolumeFraction() const
{
    return mQuiescentVolumeFraction;
}

void HairFollicleUniformG1GenerationalCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double HairFollicleUniformG1GenerationalCellCycleModel::GetEquilibriumVolume() const
{
    return mEquilibriumVolume;
}

double HairFollicleUniformG1GenerationalCellCycleModel::GetCurrentQuiescentDuration() const
{
    return mCurrentQuiescentDuration;
}

double HairFollicleUniformG1GenerationalCellCycleModel::GetCurrentQuiescentOnsetTime() const
{
    return mCurrentQuiescentOnsetTime;
}

void HairFollicleUniformG1GenerationalCellCycleModel::ResetForDivision()
{
    mGeneration++;
    if (mGeneration > mMaxTransitGenerations)
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mGeneration = 0;
    }
    AbstractSimplePhaseBasedCellCycleModel::ResetForDivision();
}

void HairFollicleUniformG1GenerationalCellCycleModel::InitialiseDaughterCell()
{
    /*
     * If the parent cell is a stem cell then its generation was reset
     * to zero when ResetForDivision() was called. The daughter cell's
     * generation must therefore be incremented here.
     */
    if (mGeneration == 0)
    {
        mGeneration = 1;
    }
    /*
    * We modify the generation-based cell cycle model so that stem cells divide symmetrically
    * and TA cells divide symmetrically for a fixed number of divisions.
    */
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        boost::shared_ptr<AbstractCellProperty> p_stem_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_stem_type);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);

        if (mGeneration > mMaxTransitGenerations)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
        }
    }
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

void HairFollicleUniformG1GenerationalCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned HairFollicleUniformG1GenerationalCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void HairFollicleUniformG1GenerationalCellCycleModel::SetMaxTransitGenerations(unsigned maxTransitGenerations)
{
    mMaxTransitGenerations = maxTransitGenerations;
}

unsigned HairFollicleUniformG1GenerationalCellCycleModel::GetMaxTransitGenerations() const
{
    return mMaxTransitGenerations;
}

void HairFollicleUniformG1GenerationalCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HairFollicleUniformG1GenerationalCellCycleModel)
