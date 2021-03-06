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

#ifndef HairFollicleUniformG1GenerationalCellCycleModel_HPP_
#define HairFollicleUniformG1GenerationalCellCycleModel_HPP_

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A fixed-generation-based cell cycle model for the HF model, where
 * stem cells divide symmetrically and TA cells divide for a fixed
 * number of generations. We extend UniformG1GenerationalCellCycleModel
 * to include contact inhibition as well.
 * 
 */
class HairFollicleUniformG1GenerationalCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mGeneration;
        archive & mMaxTransitGenerations;
        archive & mQuiescentVolumeFraction;
        archive & mEquilibriumVolume;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;
    }

protected:

    /** The generation of this cell (cells with a StemCellProliferativeType have a generation of 0) */
    unsigned mGeneration;

    /** How many generations a transit cell lives for before becoming fully differentiated. */
    unsigned mMaxTransitGenerations;

    /**
     * The fraction of the cells' equilibrium volume in G1 phase below which these cells are quiescent.
     */
    double mQuiescentVolumeFraction;

    /**
     * The cell equilibrium volume while in G1 phase.
     */
    double mEquilibriumVolume;

    /**
     * The time when the current period of quiescence began.
     */
    double mCurrentQuiescentOnsetTime;

    /**
     * How long the current period of quiescence has lasted.
     * Has units of hours.
     */
    double mCurrentQuiescentDuration;

    /**
     * Set the duration of G1 phase. This method is called on each cell at the
     * start of a simulation, and for both daughter cells immediately following
     * cell division.
     *
     * If the cell associated with this cell-cycle model has stem proliferative
     * type, then the G1 phase duration is drawn from the uniform distribution
     * U[14,18]. If the cell has transit proliferative type (semi-differentiated),
     * then the G1 phase duration is drawn from the uniform distribution U[4,6].
     * These two distributions, proposed by Meineke et al (doi:10.1046/j.0960-7722.2001.00216.x),
     * reflect indirect biological observations that stem cells cycle more
     * slowly than their progeny.
     *
     * If the cell is differentiated, then the G1 phase duration is set to DBL_MAX,
     * so that the cell will never reach the end of G1 phase.
     */
    void SetG1Duration();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     *
     * The only way for external code to create a copy of a cell-cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    HairFollicleUniformG1GenerationalCellCycleModel(const HairFollicleUniformG1GenerationalCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    HairFollicleUniformG1GenerationalCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * @param quiescentVolumeFraction
     */
    void SetQuiescentVolumeFraction(double quiescentVolumeFraction);

    /**
     * @return mQuiescentVolumeFraction
     */
    double GetQuiescentVolumeFraction() const;

    /**
     * @param equilibriumVolume
     */
    void SetEquilibriumVolume(double equilibriumVolume);

    /**
     * @return mEquilibriumVolume
     */
    double GetEquilibriumVolume() const;

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration() const;

    /**
     * @return mCurrentQuiescentOnsetTime
     */
    double GetCurrentQuiescentOnsetTime() const;

    /** Overridden ResetForDivision() method. */
    void ResetForDivision();

    /**
     * Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

    /**
     * Sets the cell's generation.
     *
     * @param generation the cell's generation
     */
    void SetGeneration(unsigned generation);

    /**
     * @return the cell's generation.
     */
    unsigned GetGeneration() const;

    /**
     * Set mMaxTransitGenerations.
     *
     * @param maxTransitGenerations the new value of mMaxTransitGenerations
     */
    void SetMaxTransitGenerations(unsigned maxTransitGenerations);

    /**
     * @return mMaxTransitGenerations
     */
    unsigned GetMaxTransitGenerations() const;

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HairFollicleUniformG1GenerationalCellCycleModel)

#endif /*HairFollicleUniformG1GenerationalCellCycleModel_HPP_*/
