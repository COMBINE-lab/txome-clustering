/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _GATB_TOOLS_FRONTLINE_HPP_
#define _GATB_TOOLS_FRONTLINE_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Terminator.hpp>
#include <set>
#include <queue>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

// types using in advanced traversal functions
struct NodeNt
{
    Node             node;
    kmer::Nucleotide nt;

    NodeNt (const Node& node, const kmer::Nucleotide& nt)  : node(node), nt(nt)  {}

    bool operator< (const NodeNt &other) const
    {
        // need to define a strict weak ordering
        if (node.kmer   != other.node.kmer)    {  return (node.kmer   < other.node.kmer);    }
        if (node.strand != other.node.strand)  {  return (node.strand < other.node.strand);  }
        return (nt < other.nt);
    }
};

/********************************************************************************/

// auxiliary class that is used by MonumentTraversal and deblooming
class Frontline
{
public:

    /** Constructor. */
    Frontline (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode,
        const Node&       previousNode,
        std::set<Node>*   all_involved_extensions = 0
    );

    /** Constructor. */
    Frontline (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode
    );

    /** */
    virtual ~Frontline() {}

    /** */
    bool go_next_depth();

    size_t size  () const  {  return _frontline.size();  }
    size_t depth () const  {  return _depth;             }

    NodeNt front () { return _frontline.front(); }

    enum reason
    {
        NONE,
        ALREADY_FRONTLINED,
        IN_BRANCHING_DEPTH,
        IN_BRANCHING_BREADTH,
        IN_BRANCHING_OTHER,
        MARKED
    };
    reason stopped_reason;

protected:

    virtual bool check (const Node& node)  { return true; }

    Direction _direction;

    const Graph& _graph;

    Terminator&  _terminator;

    typedef std::queue<NodeNt> queue_nodes;
    queue_nodes _frontline;

    int  _depth;

    std::set<Node>* _all_involved_extensions;

    std::set<Node::Value> _already_frontlined; // making it simpler now
};

/********************************************************************************/

// auxiliary class that is used by MonumentTraversal and deblooming
class FrontlineBranching : public Frontline
{
public:

    /** Constructor. */
    FrontlineBranching (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode,
        const Node&       previousNode,
        std::set<Node>*   all_involved_extensions
    );

    /** Constructor. */
    FrontlineBranching (
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode
    );

private:

    bool check (const Node& node);
};

// a middle ground between Frontline and FrontlineBranching:
// check whether all nodes inside the frontline must be reachable from startingNode
class FrontlineReachable : public Frontline
{
public:

    /** Constructor. */
    FrontlineReachable(
        Direction         direction,
        const Graph&      graph,
        Terminator&       terminator,
        const Node&       startingNode,
        const Node&       previousNode,
        std::set<Node>*   all_involved_extensions
    );

private:

    bool check (const Node& node);
};


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_FRONTLINE_HPP_ */

