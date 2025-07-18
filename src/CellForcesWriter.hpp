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

#ifndef CELLFORCESWRITER_HPP_
#define CELLFORCESWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <vector>
#include "UblasVectorInclude.hpp"
#include "Node.hpp"

/**
 * A class written using the visitor pattern for writing forces around cells
 * in a vertex model to a file.
 *
 * The output file is called CellForces.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which in this cell is a dummy.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellForcesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CellForcesWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    double GetAreaForceContribution(CellPtr pCell, VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
    double GetLineTensionForceContribution(CellPtr pCell, VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
    double GetPerimeterForceContribution(CellPtr pCell, VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    double GetLineTensionParameter(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexBasedCellPopulation<SPACE_DIM>& rVertexCellPopulation);
    double GetAreaElasticityParameter();
    double GetPerimeterContractilityParameter();
    double GetLineTensionParameter();
    double GetBoundaryLineTensionParameter();

    virtual void WriteTimeStamp();
    virtual void WriteNewline();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellForcesWriter)

#endif /* CELLFORCESWRITER_HPP_ */
