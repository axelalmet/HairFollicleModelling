#include "HairFollicleGeometryBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCellProperty.hpp"
#include "StemCellProliferativeType.hpp"
#include "MovableStemCellProgenyProliferativeType.hpp"
#include "NonMovableStemCellProgenyProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "DermalSheathCellProliferativeType.hpp"
#include "Debug.hpp"

/* Boundary condition that imposes HF geometry shape.
 */
template<unsigned DIM>
HairFollicleGeometryBoundaryCondition<DIM>::HairFollicleGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
		double hairFollicleBaseScale,
		double hairFollicleBaseRadius,
        double hairFollicleTopWidth,
        double nicheBulgeRadius,
        c_vector<double, DIM> nicheBulgeCentre,
        double maxHeight)
		: AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
		  mHairFollicleBaseScale(hairFollicleBaseScale),
          mHairFollicleBaseRadius(hairFollicleBaseRadius),
          mHairFollicleTopWidth(hairFollicleTopWidth),
          mNicheBulgeRadius(nicheBulgeRadius),
          mNicheBulgeCentre(nicheBulgeCentre),
          mMaxHeight(maxHeight)
		  {
		  }

template<unsigned DIM>
double HairFollicleGeometryBoundaryCondition<DIM>::rGetHairFollicleBaseScale() const
{
	return mHairFollicleBaseScale;
}

template<unsigned DIM>
double HairFollicleGeometryBoundaryCondition<DIM>::rGetHairFollicleBaseRadius() const
{
	return mHairFollicleBaseRadius;
}

template<unsigned DIM>
double HairFollicleGeometryBoundaryCondition<DIM>::rGetHairFollicleTopWidth() const
{
	return mHairFollicleTopWidth;
}

template<unsigned DIM>
double HairFollicleGeometryBoundaryCondition<DIM>::rGetNicheBulgeRadius() const
{
	return mNicheBulgeRadius;
}

template<unsigned DIM>
const c_vector<double, DIM>& HairFollicleGeometryBoundaryCondition<DIM>::rGetNicheBulgeCentre() const
{
	return mNicheBulgeCentre;
}

template<unsigned DIM>
double HairFollicleGeometryBoundaryCondition<DIM>::rGetMaxHeight() const
{
	return mMaxHeight;
}

template<unsigned DIM>
void HairFollicleGeometryBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
	///\todo Move this to constructor. If this is in the constructor then Exception always throws.
	if (dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(this->mpCellPopulation)==NULL)
	{
		EXCEPTION("HairFollicleGeometryBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
	}

	assert((dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
			|| (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation)) );


	if (DIM != 1)
	{
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                cell_iter != this->mpCellPopulation->End();
                ++cell_iter)
        {

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            c_vector<double, DIM> current_location = p_node->rGetLocation();

            double x = current_location[0];
            double y = current_location[1];

            /* 
             * I think it's easier to do this based on cell type, because the HF population
             * structure is so hierarchical. What that means is:
             * Stem cells: confined to semi-circular bulge defined by its centre and radius that sits
             *              outside of HF, but when in HF, they sit within a horizontal strip defined 
             *              by the radius.
             * Movable/non-movable progeny:
             *              Confined to two vertical strips about width of HF and annulus that 
             *              defines the circular HF base. NB Movable progeny should only be in region x < 0.
             * Dermal papilla fibroblasts: 
             *              Confined to lower region of HF base, that is defined by the circular HF base
             *              and a reflected parabola.
             * Dermal sheath fibroblasts:
             *              Line the entire hair follicle, sit below the niche bulge and line up with the DP.
             * Matrix TA cells:
             *              Confined within HF base but above the parabolic region where the DP cells sit.
             * 
             * Differentiated cells: 
             *              Sit above TA cells and confined within the vertical boundaries imposed by the
             *              SC progeny.
             */

            // If we look at stem cells, we judge them a little bit differently
            if (cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
            {
                c_vector<double, DIM> nearest_point = current_location;

                if (x < mNicheBulgeCentre[0])
                {                    
                    // Confine stem cells to the bulge niche
                    if (norm_2(current_location - mNicheBulgeCentre) > mNicheBulgeRadius)
                    {
                        nearest_point[0] = mNicheBulgeCentre[0] + (x - mNicheBulgeCentre[0]) * mNicheBulgeRadius / norm_2(current_location - mNicheBulgeCentre);
                        nearest_point[1] = mNicheBulgeCentre[1] + (y - mNicheBulgeCentre[1]) * mNicheBulgeRadius / norm_2(current_location - mNicheBulgeCentre);
                    }
                }
                else
                {
                    if (x > -mHairFollicleTopWidth + 1.0)
                    {
                        nearest_point[0] = -mHairFollicleTopWidth + 1.0;
                    }

                    if (y > mNicheBulgeCentre[1] + mNicheBulgeRadius)
                    {
                        nearest_point[1] = mNicheBulgeCentre[1] + mNicheBulgeRadius;
                    }
                    
                }
                
                // Set the new location
                p_node->rGetModifiableLocation() = nearest_point;
            }
            else if  ( (cell_iter->GetCellProliferativeType()->template IsType<MovableStemCellProgenyProliferativeType>())||
                        (cell_iter->GetCellProliferativeType()->template IsType<NonMovableStemCellProgenyProliferativeType>()) )
            {
                c_vector<double, DIM> nearest_point = current_location;

                // If we are above the circle that defines the HF base.
                if (y >=  pow(pow(mHairFollicleBaseRadius, 2.0) - pow(mHairFollicleTopWidth, 2.0), 0.5))
                {
                    // Make sure that the cells are confined within the vertical strip along the HF
                    if (x < 0.0)
                    {
                        if (x < -mHairFollicleTopWidth)
                        {
                            nearest_point[0] = -mHairFollicleTopWidth;
                        }
                        else if (x > -mHairFollicleTopWidth + 1.0)
                        {
                            nearest_point[0] = -mHairFollicleTopWidth + 1.0;
                        }

                        if (y > mNicheBulgeCentre[1]) // We need to make sure progeny stay outside of the bulge.
                        {
                            if ( y < mNicheBulgeCentre[1] + mNicheBulgeRadius)
                            {
                                nearest_point[1] = mNicheBulgeCentre[1] + mNicheBulgeRadius;
                            }
                        }
                        else
                        {
                            if ( y > mNicheBulgeCentre[1] - mNicheBulgeRadius)
                            {
                                nearest_point[1] = mNicheBulgeCentre[1] - mNicheBulgeRadius;
                            }
                        }
                    }
                    else
                    {
                        if (x > mHairFollicleTopWidth)
                        {
                            nearest_point[0] = mHairFollicleTopWidth;
                        }
                        else if (x < mHairFollicleTopWidth - 1.0)
                        {
                            nearest_point[0] = mHairFollicleTopWidth - 1.0;
                        }

                    }
                    if (y > mMaxHeight) // Impose a hard barrier on the progeny (not on diff cells, because we want hair to shoot out)
                    {
                        nearest_point[1] = mMaxHeight;
                    }
                }
                else if (y >=  pow(pow(mHairFollicleBaseRadius - 1.25, 2.0) - pow(mHairFollicleTopWidth - 1, 2.0), 0.5))
                {
                    // Make sure that the cells are confined within the vertical strip along the HF
                    if (x < 0.0)
                    {
                        if (x > -mHairFollicleTopWidth + 1.0)
                        {
                            nearest_point[0] = -mHairFollicleTopWidth + 1.0;
                        }
                        else if (x < -mHairFollicleTopWidth)
                        {
                            if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
                            {
                                nearest_point[0] = x * mHairFollicleBaseRadius / norm_2(current_location);
                                nearest_point[1] = y * mHairFollicleBaseRadius / norm_2(current_location);
                            }
                        }

                        if (y > mNicheBulgeCentre[1]) // We need to make sure progeny stay outside of the bulge.
                        {
                            if ( y < mNicheBulgeCentre[1] + mNicheBulgeRadius)
                            {
                                nearest_point[1] = mNicheBulgeCentre[1] + mNicheBulgeRadius;
                            }
                        }
                        else
                        {
                            if ( y > mNicheBulgeCentre[1] - mNicheBulgeRadius)
                            {
                                nearest_point[1] = mNicheBulgeCentre[1] - mNicheBulgeRadius;
                            }
                        }
                    }
                    else
                    {
                        if (x < mHairFollicleTopWidth - 1.0)
                        {
                            nearest_point[0] = mHairFollicleTopWidth - 1.0;
                        }
                        else if (x > mHairFollicleTopWidth)
                        {
                            if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
                            {
                                nearest_point[0] = x * mHairFollicleBaseRadius / norm_2(current_location);
                                nearest_point[1] = y * mHairFollicleBaseRadius / norm_2(current_location);
                            }
                        }

                    }
                    if (y > mMaxHeight) // Impose a hard barrier on the progeny (not on diff cells, because we want hair to shoot out)
                    {
                        nearest_point[1] = mMaxHeight;
                    }
                }
                else if (y >= 0.0) // If we are in the upper half of the HF base
                {
                    if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
                    {
                        nearest_point[0] = x * mHairFollicleBaseRadius / norm_2(current_location);
                        nearest_point[1] = y * mHairFollicleBaseRadius / norm_2(current_location);
                    }
                    // else if (pow(x, 2.0) + pow(y, 2.0) < pow(mHairFollicleBaseRadius - 1.25, 2.0) )
                    // {
                    //     nearest_point[0] = x * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                    //     nearest_point[1] = y * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                    // }
                }
                else // Should be in te lower half of the base
                {   
                    // if ( (x < -sqrt(fabs(y)/mHairFollicleBaseScale))||(x > sqrt(fabs(y)/mHairFollicleBaseScale) ) )
                    // {
                        if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
                        {
                            nearest_point[0] = x * mHairFollicleBaseRadius / norm_2(current_location);
                            nearest_point[1] = y * mHairFollicleBaseRadius / norm_2(current_location);
                        }
                        // else if (pow(x, 2.0) + pow(y, 2.0) < pow(mHairFollicleBaseRadius - 1.25, 2.0) )
                        // {
                        //     nearest_point[0] = x * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                        //     nearest_point[1] = y * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                        // }

                    // }

                    if ( (x < 0.0)&&(x > -sqrt(fabs(y)/mHairFollicleBaseScale)) )
                    {
                        nearest_point[0] = -sqrt(fabs(y)/mHairFollicleBaseScale);
                    }
                    else if ( (x > 0.0)&&(x < sqrt(fabs(y)/mHairFollicleBaseScale)) )
                    {
                        nearest_point[0] = sqrt(fabs(y)/mHairFollicleBaseScale);
                    }
                    
                    // Shouldn't enter the DP region
                    // if (y < - mHairFollicleBaseScale * x * x)
                    // {
                    //     nearest_point[1] = - mHairFollicleBaseScale * x * x;
                    // }
                }

                // Set the new location
                p_node->rGetModifiableLocation() = nearest_point;
            }
            // Dermal papilla cells (modelled as fibroblasts) should be confined to the lower parabola of the stem cell niche
            // Essentially, DP cells are static in this model.
            else if (cell_iter->GetCellProliferativeType()->template IsType<FibroblastCellProliferativeType>())
            {
                c_vector<double, DIM> nearest_point = current_location;

                if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius + 0.5, 2.0)) // If cell has gone too far below and falls outside of the disc.
                {
                    nearest_point[0] = x * (mHairFollicleBaseRadius + 0.5)/ norm_2(current_location);
                    nearest_point[1] = y * (mHairFollicleBaseRadius + 0.5)/ norm_2(current_location);
                }

                // Cell can't be too high or too wide
                if (y > -mHairFollicleBaseScale * x * x)
                {
                    nearest_point[1] = -mHairFollicleBaseScale * x * x;
                }
                
                // Set the new location
                p_node->rGetModifiableLocation() = nearest_point;
            }
            // Dermal sheaths line the HF 
            else if (cell_iter->GetCellProliferativeType()->template IsType<DermalSheathCellProliferativeType>())
            {
                c_vector<double, DIM> nearest_point = current_location;

                // If we're above the circle that defines the HF base.
                if (y >= pow(pow(mHairFollicleBaseRadius + 0.5 , 2.0) - pow(mHairFollicleTopWidth + 0.5, 2.0), 0.5))
                {
                    if ( (x < 0.0)&&(y > mNicheBulgeCentre[1] - mNicheBulgeRadius) )
                    {
                        nearest_point[1] = mNicheBulgeCentre[1] - mNicheBulgeRadius;
                    }

                    if (x < -mHairFollicleTopWidth - 0.5)
                    {
                        nearest_point[0] = -mHairFollicleTopWidth - 0.5;

                        if (y > mNicheBulgeCentre[1] - mNicheBulgeRadius)
                        {
                            nearest_point[1] = mNicheBulgeCentre[1] - mNicheBulgeRadius;
                        }
                    }
                    else if (x > mHairFollicleTopWidth + 0.5)
                    {
                        nearest_point[0] = mHairFollicleTopWidth + 0.5;
                    }

                    if ( (x > 0.0)&&(y > mMaxHeight) )
                    {
                        nearest_point[1] = mMaxHeight;
                    }
                }
                else
                {

                    if (y > 0.0)
                    {
                        if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius + 0.5, 2.0))
                        {
                            nearest_point[0] = x * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                        }
                        else if (pow(x, 2.0) + pow(y, 2.0) < pow(mHairFollicleBaseRadius + 0.5, 2.0))
                        {
                            nearest_point[0] = x * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                        }
                    }
                    else
                    {

                        if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius + 0.5, 2.0))
                        {
                            nearest_point[0] = x * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                        }
                        else if (pow(x, 2.0) + pow(y, 2.0) < pow(mHairFollicleBaseRadius + 0.5, 2.0))
                        {
                            nearest_point[0] = x * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius + 0.5) / norm_2(current_location);
                        }

                        if ( (x < 0.0)&&(x >= -pow(fabs(y)/mHairFollicleBaseScale, 0.5)) ) 
                        {
                            nearest_point[0] = -pow(fabs(y)/mHairFollicleBaseScale, 0.5);
                        }
                        else if ( (x > 0.0)&&(x <= pow(fabs(y)/mHairFollicleBaseScale, 0.5)) )
                        {
                            nearest_point[0] = pow(fabs(y)/mHairFollicleBaseScale, 0.5);
                        }
                        
                    }
                }

                // Set the new location
                p_node->rGetModifiableLocation() = nearest_point;
            }
            else // Should be a differentiated and/or TA cells
            {
                c_vector<double, DIM> nearest_point = current_location;

                if (y >= pow(pow(mHairFollicleBaseRadius - 1.5, 2.0) - pow(mHairFollicleTopWidth - 1.0, 2.0), 0.5))
                {
    
                    if ( (x < 0.0)&&(x < -mHairFollicleTopWidth + 1.0) )
                    {
                        nearest_point[0] = -mHairFollicleTopWidth + 1.0;
                    }
                    else if ( (x > 0.0)&&(x > mHairFollicleTopWidth - 1.0))
                    {
                        nearest_point[0] = mHairFollicleTopWidth - 1.0;
                    }
                }
                // If we are below the HF outer root sheath, within the base, we should make sure the TA/diff cells
                // don't mix with the stem cell progeny.
                else if (y >= 0.0)
                {
                    if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius - 1.25, 2.0))
                    {
                            nearest_point[0] = x * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                    }
                    else if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
                    {
                            nearest_point[0] = x * (mHairFollicleBaseRadius) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius) / norm_2(current_location);
                    }
                }
                else // Cells should be in lower half of base
                {
                    // Confine cells so that they're not outside of the base.
                    if ( (x < -sqrt(fabs(y)/mHairFollicleBaseScale))||(x > sqrt(fabs(y)/mHairFollicleBaseScale) ) )
                    {
                        if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius - 1.25, 2.0))
                        {
                            nearest_point[0] = x * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                            nearest_point[1] = y * (mHairFollicleBaseRadius - 1.25) / norm_2(current_location);
                        }
                    } 
                    else // Make sure they don't enter DP region as well
                    {
                        if (y < -mHairFollicleBaseScale * x * x)
                        {
                            nearest_point[1] = -mHairFollicleBaseScale * x * x;
                        }
                    }
                }

                // Set the new location
                p_node->rGetModifiableLocation() = nearest_point;  
            }
        }
    }
	else
	{
		// DIM == 1
		NEVER_REACHED;
		//PlaneBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
	}
}

template<unsigned DIM>
bool HairFollicleGeometryBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
	bool condition_satisfied = true;

	// if (DIM == 1)
	// {
	// 	EXCEPTION("PlaneBoundaryCondition is not implemented in 1D");
	// }
	// else
	// {
    //     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
    //             cell_iter != this->mpCellPopulation->End();
    //             ++cell_iter)
    //     {
    //         unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
    //         Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

    //         c_vector<double, DIM> current_location = p_node->rGetLocation();

    //         double x = current_location[0];
    //         double y = current_location[1];
            
    //         if (y < 0.0)
    //         {
    //             // If the cell has gone past the base.
    //             if (y < -mHairFollicleBaseScale * x * x ) //If we have moved past the 
    //             {
    //                 condition_satisfied = false;
    //                 break;
    //             }
    //             if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
    //             {
    //                 condition_satisfied = false;
    //                 break;
    //             }
    //         }
    //         else
    //         {
    //             if (y > mHairFollicleBaseRadius - 1.0)
    //             {
    //                 if ( (x < - mHairFollicleTopWidth)||(x > mHairFollicleTopWidth) )
    //                 {
    //                     condition_satisfied = false;
    //                     break;
    //                 }
    //             }
    //             else
    //             {  
    //                 if (pow(x, 2.0) + pow(y, 2.0) > pow(mHairFollicleBaseRadius, 2.0))
    //                 {
    //                     condition_satisfied = false;
    //                     break;
    //                 }

    //             }
    //         }
	// 	}
	// }

	return condition_satisfied;
}

template<unsigned DIM>
void HairFollicleGeometryBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<HairFollicleBaseScale>" << mHairFollicleBaseScale << "</HairFollicleBaseScale>\n";
	*rParamsFile << "\t\t\t<HairFollicleBaseRadius>" << mHairFollicleBaseRadius << "</HairFollicleBaseRadius>\n";
	*rParamsFile << "\t\t\t<HairFollicleTopWidth>" << mHairFollicleTopWidth << "</HairFollicleTopWidth>\n";
	*rParamsFile << "\t\t\t<NicheBulgeRadius>" << mNicheBulgeRadius << "</NicheBulgeRadius>\n";

    *rParamsFile << "\t\t\t<NicheBulgeCentre>";
    for (unsigned i = 0; i != DIM-1U; i++)
    {
        *rParamsFile << mNicheBulgeCentre[i] << ",";
    }
	*rParamsFile << mNicheBulgeCentre[DIM - 1] << "</NicheBulgeCentre>\n";

    *rParamsFile << "\t\t\t<MaxHeight>" << mMaxHeight << "</MaxHeight>\n";

	// Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class HairFollicleGeometryBoundaryCondition<1>;
template class HairFollicleGeometryBoundaryCondition<2>;
template class HairFollicleGeometryBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HairFollicleGeometryBoundaryCondition)
