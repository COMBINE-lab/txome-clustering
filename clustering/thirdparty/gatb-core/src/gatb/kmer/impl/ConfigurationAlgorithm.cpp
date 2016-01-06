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

#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/kmer/impl/LinearCounter.hpp>

#include <cmath>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/

#define DEBUG(a)  //printf a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

// estimated the number of distinct kmers in a dataset
// wrapper around a Linear Counter. Adapted from Kmergenie code.
// why not a hyperloglog? it seems that the transition from the 32bit-hash initial implementation to 64bits, and supporting billions of elements, is nontrivial, so i didn't bother
// probably deserves to be in its own file
template<size_t span>
class EstimateNbDistinctKmers
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::ModelCanonical  ModelCanonical;
    typedef typename Kmer<span>::ModelDirect     ModelDirect;
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
    typedef typename Model::Kmer                        KmerType;

    /** */
    void estimate()
    {
        nb_distinct_kmers =(unsigned long)( (float)(linearCounter->count( )) * ((float)nbKmersTotal / (float)nbProcessedKmers)); // dubious linear extrapolation, that's all I got

        abs_error = abs((long)(nb_distinct_kmers-previous_nb_distinct_kmers));

        previous_nb_distinct_kmers = nb_distinct_kmers;
    }

    /** */
    void operator() (Sequence& sequence)
    {
        /** We build the kmers from the current sequence. */
        if (model.build (sequence.getData(), kmers) == false)  {  throw "reached EOF"; return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            linearCounter->add((kmers[i].value()));


            // heuristics to stop early, i found that it's inaccurate with low coverage (e.g. on dsk/test/FiftyK.fastq)
            /*
            if (nbProcessedReads % eval_every_N_reads == 0 )
            {

                // let's see if the estimation converges..
                // the following stopping condition will grossly over-estimate the number of distinct kmers
                // but I expect the correct result to be in the same order of magnitude
                // and better to overestimate than underestimate (for both dsk and kmergenie)
                   estimate();
                   bool debug = true;
                   if (debug)
                       printf("linear estimator at %ld kmers, number of distinct kmers estimated now: %ld, abs error: %ld\n",nbProcessedKmers, nb_distinct_kmers, abs_error);
                   if (abs_error < previous_nb_distinct_kmers/20) // 5% error
                   {
                       throw "LinearCounter converged"; // well, "converged" is a big word
                       return;
                   }
                   if (!linearCounter->is_accurate())
                   {
                   printf("LinearCounter is inaccurate";
                   return;
                   }

            }*/

        }
        nbProcessedKmers += kmers.size();
        nbProcessedReads++;
        //if (nbProcessedReads % 100000 == 0) printf("nb: %ld\n",nbProcessedReads);

        // disabled progress
        //if (nbCurProgressKmers > 500000)   {  _progress.inc (nbCurProgressKmers);  nbCurProgressKmers = 0;  }
    }

    EstimateNbDistinctKmers (Model& model, u_int32_t max_memory, unsigned long nb_kmers_total, tools::dp::IteratorListener* progress)
        : model(model),  eval_every_N_reads(10000000),   nbKmersTotal(nb_kmers_total),
        nbProcessedKmers(0), nbCurProgressKmers(0), previous_nb_distinct_kmers(0), nbProcessedReads(0), abs_error(0)
        //, _progress  (progress,System::thread().newSynchronizer())
    {
        unsigned long size_linearCounter; // (in bits)
        /* let's set it to just use half of all memory available at most, ok? this isn't very robust for huge dataset, so to be tested*/
        /* if it's a tiny dataset, let's set it to total number of kmers */
        size_linearCounter = std::min( nb_kmers_total, (unsigned long) (max_memory*8*1024*1024/2) );
        linearCounter =  new LinearCounter<span>(size_linearCounter);
    }

    unsigned long getEstimation()
    {
        estimate();
        // soo.. if it's not accurate, let's assume we have a hugeload of kmers, and let's be safe, we return the total number of kmers
        if (!linearCounter->is_accurate())
        {
            cout << "Warning: linear counter was not accurate, returning worst-case estimation of number of distinct kmers";
            return nbKmersTotal;
        }
        return nb_distinct_kmers;
    }

private:

    /** Local resources. */
    Model&    model;
    unsigned long nbProcessedReads, nbProcessedKmers;
    unsigned long nbCurProgressKmers;
    unsigned long nbKmersTotal;
    unsigned long abs_error;
    vector<KmerType> kmers;
    LinearCounter<span> *linearCounter;
    int eval_every_N_reads;
    unsigned long previous_nb_distinct_kmers, nb_distinct_kmers;
};


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
ConfigurationAlgorithm<span>::ConfigurationAlgorithm (bank::IBank* bank, IProperties* input)
    : Algorithm("configuration", -1, input), _bank(0), _input (0)
{
    setBank  (bank);
    setInput (input);

    _config._kmerSize           = input->getInt (STR_KMER_SIZE);
    _config._minim_size         = input->getInt (STR_MINIMIZER_SIZE);
    _config._repartitionType    = input->getInt (STR_REPARTITION_TYPE);
    _config._minimizerType      = input->getInt (STR_MINIMIZER_TYPE);

    parse (input->getStr (STR_SOLIDITY_KIND), _config._solidityKind);

    _config._max_disk_space     = input->getInt (STR_MAX_DISK);
    _config._max_memory         = input->getInt (STR_MAX_MEMORY);
    _config._nbCores            = input->get(STR_NB_CORES) ? input->getInt(STR_NB_CORES) : 0;

    _config._abundance = CountRange (input->getInt (STR_KMER_ABUNDANCE_MIN), input->getInt (STR_KMER_ABUNDANCE_MAX));

    if (_config._nbCores == 0)  { _config._nbCores = system::impl::System::info().getNbCores(); }

    _config._nb_partitions_in_parallel = _config._nbCores;

    _config._partitionType = 0;

    _config._nb_bits_per_kmer = Type::getSize();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
ConfigurationAlgorithm<span>::~ConfigurationAlgorithm ()
{
    setBank  (0);
    setInput (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void ConfigurationAlgorithm<span>::execute ()
{
    float load_factor = 0.7;

   /** We get some information about the bank. */

    /** By default, we want to have mmers of size 8. However (for unit tests for instance),
     * we may need to have kmer sizes less than 8; in such a case, we set by convention m=k-1. */
    if (_config._minim_size == 0)  {   _config._minim_size = 8;  }

    _config._minim_size = std::min ((int)_config._kmerSize-1, (int)_config._minim_size);

    // optimism == 0 mean that we guarantee worst case the memory usage,
    // any value above assumes that, on average, any distinct k-mer will be seen 'optimism+1' times
    int optimism = 0; // 0: guarantees to always work; above 0: risky

    /** We get some information about the bank. */
    _bank->estimate (_config._estimateSeqNb, _config._estimateSeqTotalSize, _config._estimateSeqMaxSize);

    /** Some checks. */
    if (_config._estimateSeqNb==0)  { throw Exception ("Empty bank"); }

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space_min = 2000;
    _config._available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    size_t meanSeqLen = (size_t) ( (double) _config._estimateSeqTotalSize / (double) _config._estimateSeqNb);
    size_t usedSeqLen = meanSeqLen > _config._kmerSize ? meanSeqLen : _config._estimateSeqMaxSize;

    int64_t kmersNb = (usedSeqLen - _config._kmerSize + 1) * _config._estimateSeqNb;

    /** We have to be sure that the kmers number is ok. */
    if (kmersNb <= 0)  {  throw Exception ("Configuration failed : estimated that longest sequence is %ld nt but kmer size is %ld", _config._estimateSeqMaxSize, _config._kmerSize);     }

    /** The estimated kmers number is ok. */
    _config._kmersNb  = kmersNb;

    _config._volume =  _config._kmersNb * sizeof(Type) / MBYTE;  // in MBytes

    if (_config._volume == 0)   { _config._volume = 1; }    // tiny files fix

    u_int64_t volume_minim = _config._volume * 0.5 *1.2  ; //0.5 for using kxmers   1.2 if bad repartition of minimizers ( todo sampling to assert ram usage)

    if (volume_minim == 0)   { volume_minim = 1; }    // tiny files fix

    /** We get max(75%, 100% - X GB) */
    if (_config._max_disk_space == 0)  { _config._max_disk_space = std::max ((75*_config._available_space)/100, _config._available_space-available_space_min);  }
    if (_config._max_disk_space == 0)  { _config._max_disk_space = 10000; }

    if (_config._max_memory == 0)  {  _config._max_memory = System::info().getMemoryProject(); }
    if (_config._max_memory == 0)  {  _config._max_memory = 1000; }

    /* make sure to not use more mem than system, when max_memory has default value (useful for docker images) */
    if (_config._max_memory == 2000)  {
        unsigned long system_mem = System::info().getMemoryPhysicalTotal() / MBYTE;
        if (_config._max_memory > (system_mem * 2) / 3)
        {
            _config._max_memory = (system_mem * 2) / 3;
            cout << "Warning: default memory usage (2000 MB) is close or above system max, setting memory to: " << _config._max_memory << " MB" << endl;
        }
    }

    assert (_max_disk_space > 0);

    _config._nb_passes = ( (_config._volume/3) / _config._max_disk_space ) + 1; //minim, approx volume /3
    //_nb_passes = 1; //do not constrain nb passes on disk space anymore (anyway with minim, not very big)
    //increase it only if ram issue

    //printf("_volume  %lli volume_minim %lli _max_disk_space %lli  _nb_passes init %i  \n", _volume,volume_minim,_max_disk_space,_nb_passes);
    size_t max_open_files = System::file().getMaxFilesNumber() / 2;
    u_int64_t volume_per_pass;
    float est_volume_distinct_ratio;

#if 0
    /* disabled by default; this was an experiment */
    if (_flagEstimateNbDistinctKmers)
    {
        /* we estimate the volume of distinct kmers vs total number of kmers.
         * we store it in the variable "est_volume_distinct_ratio"
         * to compute it, we need a linear counter, let's call it now */

        TIME_INFO (getTimeInfo(), "estimate_distinct_kmers");
        Iterator<Sequence>* itSeq = _bank->iterator();
        LOCAL (itSeq);

        //_progress->setMessage (progressFormat0); // not touching progress here anymore
        Model model (_kmerSize, _minim_size);
        EstimateNbDistinctKmers<span> estimate_nb_distinct_kmers_function(model, _max_memory, kmersNb, _progress);

        /** We launch the iteration of the sequences iterator with the created functors. */
        try {
            itSeq->iterate (estimate_nb_distinct_kmers_function);
        }
        catch (const char* except)
        {

        }
        _estimatedDistinctKmerNb = estimate_nb_distinct_kmers_function.getEstimation();
        est_volume_distinct_ratio = (float) _estimatedDistinctKmerNb / (float)kmersNb;
        //est_volume_distinct_ratio = 1; // for debug
        /* est_volume_distinct_ratio == 1 mean that we guarantee worst case the memory usage,
           the value mean that, on average, a k-mer will be seen 'est_volume_distinct_ratio' times */
        // if wrongly estimated, the error 'OAHash: max rehashes..' can happen
        printf ("LinearCounter done, estimated %ld number of distinct kmers, ratio to total number of kmers: %.2f\n", (long)_estimatedDistinctKmerNb, est_volume_distinct_ratio);
    }
#endif

    do  {

        assert (_nb_passes > 0);
        volume_per_pass = volume_minim / _config._nb_passes;

        assert (_max_memory > 0);
        //printf("volume_per_pass %lli  _nbCores %zu _max_memory %i \n",volume_per_pass, _nbCores,_max_memory);

        if (_config._partitionType == 1) // adjust partition size for hash table
        {
            _config._nb_partitions = (u_int32_t) ceil((float) _config._nb_partitions / load_factor);
            _config._nb_partitions = ((_config._nb_partitions * OAHash<Type>::size_entry()) + sizeof(Type)-1) / sizeof(Type); // also adjust for hash overhead
#if 0
            if (_flagEstimateNbDistinctKmers)
            {
                // use our estimation of number of distinct kmers to refine number of partitions
                // it's essentially a way to set optimism optimally
                // i'm not enabling it because computing it is slow, and reward was too small
                _nb_partitions = std::max ((u_int32_t) ceil( (float) _nb_partitions *  est_volume_distinct_ratio  * 1.3 ), (u_int32_t)1);  // 1.3 is for security
            }
            else
#endif
            {
                _config._nb_partitions = std::max ((_config._nb_partitions/(optimism+1)), (u_int32_t)1);
            }
        }
        else
        {
            // _nb_partitions  = ( (volume_per_pass*_nbCores) / _max_memory ) + 1;
            _config._nb_partitions  = ( ( volume_per_pass* _config._nb_partitions_in_parallel) / _config._max_memory ) + 1;

            //printf("nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
            //_nb_partitions = max_open_files; break;
        }
             if (_config._nb_partitions >= max_open_files && _config._nb_partitions_in_parallel >1)     { _config._nb_partitions_in_parallel  = _config._nb_partitions_in_parallel /2;  }
        else if (_config._nb_partitions >= max_open_files && _config._nb_partitions_in_parallel == 1)   { _config._nb_passes++;  }
        else                                                                            { break;         }

        //printf("update nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
    } while (1);

    if (_config._nb_partitions < 50 &&  (max_open_files - _config._nb_partitions  > 30) ) _config._nb_partitions += 30; //some more does not hurt

    //round nb parti to upper multiple of _nb_partitions_in_parallel if possible
    int  incpart = _config._nb_partitions_in_parallel - _config._nb_partitions % _config._nb_partitions_in_parallel;
    incpart = incpart % _config._nb_partitions_in_parallel;
    if(((int)max_open_files - (int)_config._nb_partitions  > incpart)) _config._nb_partitions+= incpart ;

    //_nb_partitions_in_parallel = 1 ;

    //then put _nbCores_per_partition

    _config._nbCores_per_partition =  _config._nbCores / _config._nb_partitions_in_parallel ; //how to best distrib available cores ?

    // with this formula we'll sometimes use less than _nbCores (maybe enforce _nb_partitions_in_parallel is power of two ?)
    DEBUG (("ConfigurationAlgorithm<span>::execute  _nbCores %zu  _nb_partitions_in_parallel %zu   _nbCores_per_partition %zu  nb part %u nb passes %i\n",
        _config._nbCores, _config._nb_partitions_in_parallel, _config._nbCores_per_partition, _config._nb_partitions, _config._nb_passes
    ));

    assert(_nbCores_per_partition > 0);

    /** We set the config as set. */
    _config._isComputed = true;

    /** We collect some statistics. */
    getInfo()->add (1, _config.getProperties());
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class ConfigurationAlgorithm <KSIZE_1>;
template class ConfigurationAlgorithm <KSIZE_2>;
template class ConfigurationAlgorithm <KSIZE_3>;
template class ConfigurationAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

