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

/** \file Pool.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for pooling objects
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/* Cette class dÈfinit une pool memoire pour allocation rapide de la table de hachage
 * utilisee quand seed >14 */
template <typename graine_type, typename value_type=int>  class Pool
{
public:

    typedef u_int32_t  cell_ptr_t;

    struct cell
    {
        graine_type graine;
        cell_ptr_t  suiv;
        value_type val;
    };

    /** Default constructor.
     * \param[in] tai :  2^22  16 M cells *16 o    blocs de 256 Mo
     * \param[in] N : 2^10  soit 4 G cells max
     * */
    Pool (size_t tai=4194304, size_t N=1024) : TAI_POOL(tai), N_POOL(N)
    {
        n_pools = 0; n_cells=0;
        //allocation table de pool :
        tab_pool = (cell**)  MALLOC (N_POOL*sizeof(cell*) );

        tab_pool[0]=0; n_pools++; // la premiere pool est NULL, pour conversion null_internal -> null

        //allocation de la premiere pool :
        pool_courante =(cell*)  MALLOC (TAI_POOL*sizeof(cell) );
        tab_pool[n_pools] = pool_courante;
        n_pools++;
    }

    /**  Destructeur  */
    ~Pool()
    {
        // la pool 0 est NULL
        for(size_t i=1;i<n_pools;i++)  {  FREE ( tab_pool[i] );  }

        FREE (tab_pool);
    }

    /**  allocate cell, return internal pointer type ( 32bits) */
    cell_ptr_t  allocate_cell()
    {
        cell_ptr_t internal_adress = 0;
        // ncells = nb de cells deja utilisees
        if (n_cells <TAI_POOL)
        {
            internal_adress  = n_pools -1;    // low 10 bits  : pool number
            internal_adress |= n_cells << 10; // 22 high bits : cell number in pool
            n_cells ++;

            return internal_adress;
        }
        else // la piscine est pleine, on en alloue une nouvelle
        {
            if (n_pools>= N_POOL)
            {
                // will happen when  4G cells are allocated, representing 64 Go
                throw system::Exception ("Internal memory allocator is full!");
            }
            pool_courante =(cell*)  MALLOC (TAI_POOL*sizeof(cell) );
            tab_pool[n_pools] = pool_courante;
            n_pools++;
            n_cells = 1;

            internal_adress = n_pools -1; // low 8 bits  : pool number
            // 22 high bits are 0

            return internal_adress;
        }
    }

    /** */
    cell*  internal_ptr_to_cell_pointer(cell_ptr_t internal_ptr)
    {
        unsigned int numpool =  internal_ptr & 1023;
        unsigned int numcell =  internal_ptr >> 10;

        return (tab_pool[numpool] + numcell);
    }

    /** vide toutes piscines (garde juste une pool vide)  */
    void  clear ()
    {
        for(size_t i=2;i<n_pools;i++) // garde la premiere pool pour usage futur
        {
            FREE ( tab_pool[i] );
        }

        //on repasse sur premiere pool
        pool_courante = tab_pool[1];
        n_cells=0;
        n_pools=2;
    }

private:

    /** table de cell, pour usage courant */
    cell* pool_courante;

    /** stockage de tous les pointeurs pool */
    cell** tab_pool;

    /** nombre de piscines remplies */
    unsigned int n_pools;

    /**  niveau de remplissage de la piscine courante */
    unsigned int n_cells;

    size_t TAI_POOL;
    size_t N_POOL;
};

/********************************************************************************/

//make it an allocator usable by std vector ?
class MemAllocator
{
public:

    //clear all previous allocs, and alloc pool capacity
    void reserve(u_int64_t size)
    {
        if(size ==0 && mainbuffer !=NULL)
        {
            FREE (mainbuffer);
            capacity = used_space = 0;
            mainbuffer = NULL ;
        }

        /** We add a little bit of memory in case "align" method is called often. */
        size_t extraMem = 1024;

        mainbuffer = (char*) MALLOC(size+extraMem);
        capacity   = size;
        used_space = 0;
    }

    //should be thread safe
    char* pool_malloc(u_int64_t requested_size, const char* message="")
    {
        u_int64_t synced_used_space = __sync_fetch_and_add(&used_space, requested_size);

        if (requested_size> (capacity - synced_used_space))
        {
            __sync_fetch_and_add(&used_space, -requested_size);

            throw system::Exception ("Pool allocation failed for %lld bytes (%s). Current usage is %lld and capacity is %lld",
                requested_size, message, used_space, capacity
            );

            return NULL;
        }

        return mainbuffer + synced_used_space;
    }

    /** Force alignment. */
    void align (u_int8_t alignBytes)
    {
        size_t offset = alignBytes-1 + sizeof(char*);
        char* current = mainbuffer + used_space;
        char* buffer  =  (char*)(((size_t)(current)+offset)&~(alignBytes-1));

        /** We update the used space. */
        used_space = (u_int64_t)(buffer-mainbuffer);
    }

    u_int64_t getCapacity ()  {  return capacity;   }

    u_int64_t getUsedSpace()  {  return used_space; }


    void free_all()
    {
        used_space = 0;
    }

    MemAllocator() : mainbuffer(NULL),capacity(0),used_space(0),_synchro(0)
    {
        setSynchro (system::impl::System::thread().newSynchronizer());
    }

    ~MemAllocator()
    {
        if (mainbuffer != NULL)  {  FREE (mainbuffer);  }
        setSynchro (0);
    }

    system::ISynchronizer* getSynchro()  { return _synchro; }

private :
    char*     mainbuffer;
    u_int64_t capacity; //in bytes
    u_int64_t used_space;

    system::ISynchronizer* _synchro;
    void setSynchro (system::ISynchronizer* synchro) { SP_SETATTR(synchro); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_ */
