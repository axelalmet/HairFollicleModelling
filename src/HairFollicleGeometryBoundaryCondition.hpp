#ifndef HAIRFOLLICLEGEOMETRYBOUNDARYCONDITION_HPP_
#define HAIRFOLLICLEGEOMETRYBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A boundary condition to impose the hair follicle geometry,
 * imposed by two parameters, the "scale" and "width" of the hair follicle,
 * which really determine the shape of the bottom region.
 */
template<unsigned DIM>
class HairFollicleGeometryBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /**
     * Scale of the hair follicle base.
     */
    double mHairFollicleBaseScale;

    /**
     * Width of the hair follicle base
     */
    double mHairFollicleBaseRadius;

    /**
     * Width of the hair follicle top
     */
    double mHairFollicleTopWidth;

    /**
     * Size of stem cell niche bulge
     */
    double mNicheBulgeRadius;

    /*
     * Centre of stem cell niche
     */
    c_vector<double, DIM> mNicheBulgeCentre;

    /*
     * Max height for the progeny of the HF
     */
    double mMaxHeight;

    /*
     * Height when TA cells differentiate
     */
    double mTransitDifferentiationHeight;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param hairFollicleBaseScale determines how the local minima of the HF base is
     * @param hairFollicleBaseRadius determines the width of the HF bulge
     * @param hairFollicleTopWidth determines the width of the HF top
     * @param nicheBulgeRadius size of stem cell niche
     * @param nicheBulgeCentre centre of stem cell niche (assumed to be circular in shape)
     * @param maxHeight maximal height for the non-differentiated HF cells
     * @param transitDifferentiationHeight height when TA cells differentiate (to separate from diff'd cells)
     */
    HairFollicleGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                           double hairFollicleBaseScale,
                           double hairFollicleBaseRadius,
                           double hairFollicleTopWidth, 
                           double nicheBulgeRadius,
                           c_vector<double, DIM> nicheBulgeCentre,
                           double maxHeight,
                           double transitDifferentiationHeight);

    /**
     * @return #mHairFollicleScale.
     */
    double rGetHairFollicleBaseScale() const;

    /**
     * @return #mHairFollicleWidth.
     */
    double rGetHairFollicleBaseRadius() const;

    /**
     * @return #mHairFollicleWidth.
     */
    double rGetHairFollicleTopWidth() const;

    /**
     * @return #mNicheBulgeRadius
     */
    double rGetNicheBulgeRadius() const;

        /**
     * @return #mNicheBulgeCentre.
     */
    const c_vector<double, DIM>& rGetNicheBulgeCentre() const;

    /**
     * @return #mMaxHeight
     */
    double rGetMaxHeight() const;

    /**
     * @return #mTransitDifferentiationHeight
     */
    double rGetTransitDifferentiationHeight() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HairFollicleGeometryBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const HairFollicleGeometryBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive doubles 
    double base_scale = t->rGetHairFollicleBaseScale();
    ar << base_scale;

    double base_width = t->rGetHairFollicleBaseRadius();
    ar << base_width;

    double top_width = t->rGetHairFollicleTopWidth();
    ar << top_width;

    double niche_bulge_radius = t->rGetNicheBulgeRadius();
    ar << niche_bulge_radius;

    c_vector<double, DIM> niche_bulge_centre = t->rGetNicheBulgeCentre();

    for (unsigned i = 0; i < DIM; i++)
    {
        ar << niche_bulge_centre[i];
    }

    double max_height = t->rGetMaxHeight();
    ar << max_height;

    double transit_height = t->rGetTransitDifferentiationHeight();
    ar << transit_height;
}

/**
 * De-serialize constructor parameters and initialize a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, HairFollicleGeometryBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    double base_scale;
    ar >> base_scale;

    double base_radius;
    ar >> base_radius;

    double top_width;
    ar >> top_width;

    double niche_bulge_radius;
    ar >> niche_bulge_radius;

    c_vector<double, DIM> niche_bulge_centre;
    for (unsigned i = 0; i < DIM; i++)
    {
        ar >> niche_bulge_centre[i];
    }

    double max_height;
    ar >> max_height;

    double transit_height;
    ar >> transit_height;

    // Invoke inplace constructor to initialise instance
    ::new(t)HairFollicleGeometryBoundaryCondition<DIM>(p_cell_population, base_scale, base_radius, top_width, niche_bulge_radius, niche_bulge_centre, max_height, transit_height);
}
}
} // namespace ...

#endif /*HairFollicleGeometryBoundaryCondition_HPP_*/
