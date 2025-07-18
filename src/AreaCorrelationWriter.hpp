/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef AREACORRELATIONWRITER_HPP_
#define AREACORRELATIONWRITER_HPP_

#include "AbstractCellPopulationCountWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <map>
#include <vector>
#include "UblasVectorInclude.hpp"
#include "Node.hpp"

/**
 * A class for writing the area correlation between cell neighbours.
 * Only cells that are not on the tissue boundary are considered in this correlation
 * in order to avoid boundary effects.
 * The formula used for the correlation is
 * (<A_i*A_j> - <A>^2)/var(A)
 * for all pairs of cells <i,j> in the tissue. `A' is the area of a single cell
 * in this formula.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AreaCorrelationWriter : public AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
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
        archive & boost::serialization::base_object<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    AreaCorrelationWriter();

    /**
     * Overridden WriteHeader() method.
     *
     * Write the header to file.
     *
     * @param pCellPopulation a pointer to the population to be written.
     */
    virtual void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * This will throw an exception
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * This will throw an exception
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * This will throw an exception
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * This will throw an exception.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Calculate the area correlation for this tissue and write to file
     *
     * @param pCellPopulation  The cell population
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Helper function to calculate the mean and variance of the cell area for cells
     * that are not on the tissue boundary.
     *
     * @param pCellPopulation, the population
     *
     * @returns area_statistics, the first entry is the mean, the second entry is the
     * variance
     */
    c_vector<double, 2> GetMeanInternalAreaAndVariance(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation );

    /**
     * Helper function to find all index pairs for cells that are adjacent to each other
     * and such that neither cell in each pair is on the tissue boundary. Each pair
     * will appear exactly once in the vector.
     *
     * @param pCellPopulation, the population
     *
     * @return pairs_vector, vector of pairs of indices.
     */
    std::vector< c_vector<unsigned,2> > GetAllInternalCellNeighbourIndexPairs(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation );
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AreaCorrelationWriter)

#endif /*AREACORRELATIONWRITER_HPP_*/
