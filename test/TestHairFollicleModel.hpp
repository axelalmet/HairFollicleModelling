#ifndef TESTHAIRFOLLICLEMODEL_HPP_
#define TESTHAIRFOLLICLEMODEL_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "NodesOnlyMesh.hpp" // Nodes only mesh
#include "GeneralisedLinearSpringForce.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "RepulsionForce.hpp" // Repulsion-only force (really only applied to OS models)
#include "HairFollicleMigrationForce.hpp" // Migration force to keep TA cells in base.
#include "NoCellCycleModel.hpp" // Cell cycle for fibroblasts that is dependent on exposure to wound-derived growth factors
#include "HairFollicleUniformG1GenerationalCellCycleModel.hpp" // Simple generation-based cell cycle model that draws G1 times from uniform distribution
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "StemCellProliferativeType.hpp" // HF stem cell
#include "NonMovableStemCellProgenyProliferativeType.hpp" // HF non-movable stem cell progeny types
#include "MovableStemCellProgenyProliferativeType.hpp" // HF movable stem cell progeny types
#include "TransitCellProliferativeType.hpp" // HF matrix TA cells
#include "DifferentiatedCellProliferativeType.hpp" // HF differentiated cells
#include "FibroblastCellProliferativeType.hpp" // Fibroblasts to represent the dermal papillae
#include "DermalSheathCellProliferativeType.hpp" // Fibroblasts to represent the dermal papillae
#include "CellLabel.hpp" // Add a cell label to tag certain cells
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "HairFollicleDifferentiationTrackingModifier.hpp" // Enforce differentiation in the right places
#include "VolumeTrackingModifier.hpp" // Track the cell volumes, which is needed for contact inhibition.
#include "HairFollicleGeometryBoundaryCondition.hpp" // Boundary to maintain the geometry of the HF.
#include "CellAncestorWriter.hpp" // Track cell ancestors
// #include "VolumeTrackingModifier.hpp" // Track cell volumes for contact inhibition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "HairFollicleModel";
static const double M_DT = 0.005;
static const double M_END_TIME = 150.0;
static const double M_SAMPLING_TIMESTEP = 25.0/M_DT;

/*
* A test model to test a simple fixed hair follicle model in a fixed geometry
* with four cell types: stem cells, movable stem cell progenies, non-movable
* stem cell progenies, and basal transit-amplifying cells
*/
class TestHairFollicleModel : public AbstractCellBasedTestSuite
{
public:
    void TestHairFollicle()
    {
        //Set the number of cells across and down for the array
        unsigned cells_across = 14;
        unsigned cells_up = 35;

        // Certain parameters needed for the HF geometry
        double hf_base_radius = 5.0;
        double hf_base_scale = 1.0;
        double hf_top_width = 3.0;

        // Parameters that determine where the stem cell bulge is 
        double niche_bulge_radius = 1.75; // Size of the stem cell bulge.
        c_vector<double, 2> niche_bulge_centre;
        niche_bulge_centre[0] = -0.8*hf_top_width - 1.0;
        niche_bulge_centre[1] = 20.0;

        // Set some parameters for node-based cell populations
        double radius_of_interaction = 1.5; // Radius of interaction to determine neighbourhoods
        double division_separation = 0.1; // Initial resting length upon division

        // Some mechanical parameters
        double spring_stiffness = 30.0;
        double migration_force_strength = 2.5;

        // Set the proportion of movable to nonmovable stem cell progeny
        double movable_progeny_probability = 0.5;

        // Set the height at which we impose TA differentiation
        // double transit_differentiation_height = -100.0;

        HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
        MutableMesh<2, 2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

        // Need to translate the mesh across
        p_generating_mesh->Translate(-7.0, -hf_base_radius - 1.0);

        // Remove cells outside of the intiial HF geometry
        for (AbstractMesh<2, 2>::NodeIterator node_iter = p_generating_mesh->GetNodeIteratorBegin();
                node_iter != p_generating_mesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            c_vector<double,2> node_location = node_iter->rGetLocation();
            double x = node_location[0];
            double y = node_location[1];

            if (y >= 0.0)
            {

                if (y <= pow(pow(hf_base_radius - 1.5, 2.0) - pow(hf_top_width - 1.0, 2.0), 0.5) + 1.0)
                {
                    if (pow(x, 2.0) + pow(y, 2.0) > pow(hf_base_radius + 1.0, 2.0))
                    {
                        p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                    }
                    
                }
                else
                {
                    if (x > hf_top_width + 1.0)
                    {
                        p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                    }
                    else
                    {
                        if (y > niche_bulge_centre[1] + niche_bulge_radius)
                        {
                            if (x < -hf_top_width)
                            {
                                p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                            }
                        }
                        if (y < niche_bulge_centre[1] - niche_bulge_radius)
                        {

                            if (x < -hf_top_width - 1.0)
                            {
                                p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                            }
                        }
                        else
                        {
                            if (x < niche_bulge_centre[0] - pow(pow(niche_bulge_radius, 2.0) - pow(y - niche_bulge_centre[1], 2.0), 0.5))
                            {
                                p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                            }
                        }
                    }
                }
                

            }
            else
            {
                if (pow(x, 2.0) + pow(y, 2.0) > pow(hf_base_radius + 1.0, 2.0))
                {
                    p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                }
                // else if (x > -pow(fabs(y)/hf_base_scale, 0.5)&&(x < pow(fabs(y)/hf_base_scale, 0.5)))
                // {
                //     p_generating_mesh->DeleteNodePriorToReMesh(node_index);
                // }
            }
        }
        p_generating_mesh->ReMesh();

        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, radius_of_interaction);

        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_movable_type(CellPropertyRegistry::Instance()->Get<MovableStemCellProgenyProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_nonmovable_type(CellPropertyRegistry::Instance()->Get<NonMovableStemCellProgenyProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_dp_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_ds_type(CellPropertyRegistry::Instance()->Get<DermalSheathCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_cell_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        std::vector<CellPtr> cells; //Create vector of cells

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++) // Iterator for periodic mesh
        {
            double birth_time = -30.0 * RandomNumberGenerator::Instance()->ranf();

            // Set contact-inhibition-based cell cycle model with fixed-generation-based divisions for TA cells.
            HairFollicleUniformG1GenerationalCellCycleModel* p_cycle_model = new HairFollicleUniformG1GenerationalCellCycleModel();
            p_cycle_model->SetEquilibriumVolume(0.25*M_PI);
            p_cycle_model->SetQuiescentVolumeFraction(0.8);
            p_cycle_model->SetMaxTransitGenerations(2);
            p_cycle_model->SetStemCellG1Duration(38.0);
            p_cycle_model->SetTransitCellG1Duration(26.0);
            p_cycle_model->SetDimension(2);
        
            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));

            p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(p_nonmovable_type); // Initialise as non-movable progeny

            // Set the volume
            p_cell->GetCellData()->SetItem("volume", 0.25*M_PI);
            
            cells.push_back(p_cell);
        }

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells); // Used for non-periodic
        cell_population.SetMeinekeDivisionSeparation(division_separation);

        // Mark ancestors
        cell_population.SetCellAncestorsToLocationIndices();

        // Add cell writers to track clones
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Let's specify the cell types and work out the maximum height for the HF barrier.
        double max_height = -hf_base_radius;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

            // Let's change some of the non-movable cells on the left to be movable progeny type
            if ( (x < 0.0)&&(y < niche_bulge_centre[1] - niche_bulge_radius) )
            {
                double movable_progeny_fate = RandomNumberGenerator::Instance()->ranf();
                if (movable_progeny_fate < movable_progeny_probability)
                {
                    cell_iter->SetCellProliferativeType(p_movable_type);
                }
            }
            // We need to set the stem cell bulge, TA population, DP (fibroblast) and differentiated populations.
            if (y > 0.0)
            {
                if (y > pow(pow(hf_base_radius - 1.5, 2.0) - pow(hf_top_width - 1.0, 2.0), 0.5))
                {
                    if (y > max_height)
                    {
                        max_height = y;
                    }
            
                    if (pow(x - niche_bulge_centre[0], 2.0) + pow(y - niche_bulge_centre[1], 2.0) <= pow(niche_bulge_radius, 2.0))
                    {
                        cell_iter->SetCellProliferativeType(p_stem_type);
                    }
                    if ( (x > niche_bulge_centre[0])&&(x < -hf_top_width + 1.6)
                            &&(y >= niche_bulge_centre[1] - niche_bulge_radius)&&(y <= niche_bulge_centre[1] + niche_bulge_radius) )
                    {
                        cell_iter->SetCellProliferativeType(p_stem_type);

                        // Try and tag one cell
                        if ( (y < niche_bulge_centre[1] - niche_bulge_radius + 1.25)&&(x > -hf_top_width + 0.5)&&(x < -hf_top_width + 1.75) )
                        {
                            cell_iter->AddCellProperty(p_cell_label);
                        }

                    }

                    if ((x > -hf_top_width + 1.6)&&(x < hf_top_width - 1.6))
                    {
                        cell_iter->SetCellProliferativeType(p_diff_type);
                    }
                    else if (x > hf_top_width)
                    {
                        cell_iter->SetCellProliferativeType(p_ds_type);
                    }
                    else if ( (x < -hf_top_width)&&(y < niche_bulge_centre[1] - niche_bulge_radius) )
                    {
                        cell_iter->SetCellProliferativeType(p_ds_type);
                    }
                    
                }
                else
                {
                        if ( (y <= 1.0)&&(y > -hf_base_scale * x * x) )
                        {
                            if (pow(x, 2.0) + pow(y, 2.0) < pow(hf_base_radius - 1.5, 2.0) )
                            {
                                    cell_iter->SetCellProliferativeType(p_transit_type);
                            }
                        }
                        else
                        {
                            if (pow(x, 2.0) + pow(y, 2.0) < pow(hf_base_radius - 1.5, 2.0) )
                            {
                                cell_iter->SetCellProliferativeType(p_diff_type);
                            }
                        }

                    if (pow(x, 2.0) + pow(y, 2.0) > pow(hf_base_radius - 0.1, 2.0))
                    {
                        cell_iter->SetCellProliferativeType(p_ds_type);
                    }
                }
            }
            else
            {

                if (pow(x, 2.0) + pow(y, 2.0) < pow(hf_base_radius - 1.5, 2.0))
                {
                    cell_iter->SetCellProliferativeType(p_transit_type);
                }

                if (pow(x, 2.0) + pow(y, 2.0) > pow(hf_base_radius - 0.1, 2.0))
                {
                    cell_iter->SetCellProliferativeType(p_ds_type);
                }
                
                if (y <= -hf_base_scale * x * x)
                {
                    cell_iter->SetCellProliferativeType(p_dp_type);
                }
            }

        }

        // Initialise simulator class
        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory
        std::stringstream out;
        out << "/Movable_" << movable_progeny_probability << "/TestContactInhibition/";
        std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(M_DT);
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add sort-of-linear spring force
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(radius_of_interaction);
        simulator.AddForce(p_spring_force); 

        // Add repulsion-only force
        // MAKE_PTR(RepulsionForce<2>, p_repulsion_force);
        // p_repulsion_force->SetMeinekeSpringStiffness(spring_stiffness);
        // p_repulsion_force->SetCutOffLength(radius_of_interaction);
        // simulator.AddForce(p_repulsion_force); 

        // Add migration force to retain TA cells
        MAKE_PTR(HairFollicleMigrationForce<2>, p_migration_force);
        p_migration_force->SetMigrationForceStrength(migration_force_strength);
        simulator.AddForce(p_migration_force);

        // Create a modifier to track which cells are attached to the basement membrane.
        MAKE_PTR(HairFollicleDifferentiationTrackingModifier<2>, p_differentiation_tracking_modifier);
        p_differentiation_tracking_modifier->SetNicheBulgeCentre(niche_bulge_centre);
        p_differentiation_tracking_modifier->SetNicheBulgeRadius(niche_bulge_radius);
        p_differentiation_tracking_modifier->SetMovableProgenyDifferentiationProbability(movable_progeny_probability);
        p_differentiation_tracking_modifier->SetTransitDifferentiationHeight(4.0);
        p_differentiation_tracking_modifier->SetBaseRadius(0.8*hf_base_radius);
		simulator.AddSimulationModifier(p_differentiation_tracking_modifier);

        // // Create a modifier to track cell volumes
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Add a boundary condition to maintain the hair follicle geometry
        MAKE_PTR_ARGS(HairFollicleGeometryBoundaryCondition<2>, p_bc, (&cell_population, 
                                                                        hf_base_scale, 
                                                                        0.9*hf_base_radius,
                                                                        0.9*hf_top_width,
                                                                        niche_bulge_radius,
                                                                        niche_bulge_centre,
                                                                        max_height));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTSCARFORMATIONINCROSSSECTIONALGEOMETRY */