/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef VERTEXMODELDATAWRITER_HPP_
#define VERTEXMODELDATAWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing the polygon class
 * (i.e. number of edges) of each cell to file. This class should only be
 * used when simulating a vertex-based model.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexModelDataWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
     *
     * @param directory the path to the directory in to which this class should write
     * @param cellDatafilename the filename
     */
    VertexModelDataWriter();

    /**
     * Method for VTK output - doesn't do anything yet
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden method to write the current time stamp to the file.
     * This class does NOT do this by default.
     */
    void WriteTimeStamp();

    /**
     * Overridden method to add a newline character to the stream.
     * This class does NOT do this by default.
     */
    virtual void WriteNewline();

    /**
     * Helper function to determine the average area of cell neighbours for this cell
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    double GetAverageCellAreaOfNeighbours(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Helper function to determine the average number of neighbours for this cell
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    double GetAverageNeighbourNumberOfNeighbours(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Helper function to determine if any neighbours of the cell are on the boundary of the population.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    bool IsCellOnInnerBoundary(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexModelDataWriter)

#endif /*VERTEXMODELDATAWRITER_HPP_*/
