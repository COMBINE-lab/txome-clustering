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

#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/kmer/impl/Sequence2SuperKmer.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/debruijn/impl/ContainerNode.hpp>

#include <gatb/tools/storage/impl/StorageTools.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* progressFormat0 = "DSK: Collecting stats on %s ";


/********************************************************************************/

template<size_t span>
class MmersFrequency
{
public:
    /** Shortcut. */
    typedef typename RepartitorAlgorithm<span>::ModelDirect  ModelDirect;
    typedef typename ModelDirect::Kmer                       KmerTypeDirect;

    void operator() (Sequence& sequence)
    {
        /** We first check whether we got mmers from the sequence or not. */
        if (_minimodel.build (sequence.getData(), _mmers) == false)  { return; }

        /** We loop over the mmers of the sequence. */
        for (size_t i=0; i<_mmers.size(); i++)
        {
            if (_mmers[i].isValid() == false)
                continue;

            /** increment m-mer count */
            _m_mer_counts[_mmers[i].value().getVal()] ++;
        }

        if (_nbProcessedMmers > 500000)   {  _progress.inc (_nbProcessedMmers);  _nbProcessedMmers = 0;  }
    }

    /** Constructor. */
    MmersFrequency (int mmerSize,  IteratorListener* progress,  uint32_t* m_mer_counts)
    :
        _minimodel(mmerSize), _progress (progress,System::thread().newSynchronizer()),
      _m_mer_counts(m_mer_counts),  _nbProcessedMmers(0)
    {
        u_int64_t nbminim = (uint64_t)pow(4.0,mmerSize);

        for (u_int64_t i = 0; i < nbminim; i++)  {   _m_mer_counts[i] = 0;  }
    }

protected:

    ModelDirect             _minimodel;
    vector<KmerTypeDirect>  _mmers;
    ProgressSynchro         _progress;
    uint32_t*               _m_mer_counts;
    size_t                  _nbProcessedMmers;
};

/********************************************************************************/
/* This functor class takes a Sequence as input, splits it into super kmers and
 * get information about the distribution of minimizers.
 */
template<size_t span>
class SampleRepart  : public Sequence2SuperKmer<span>
{
    //ie ce sera posible d avoir plus dinfo , estim ram max par exemple ?
public:

    /** Shortcut. */
    typedef typename Sequence2SuperKmer<span>::Type             Type;
    typedef typename Sequence2SuperKmer<span>::ModelCanonical   ModelCanonical;
    typedef typename Sequence2SuperKmer<span>::Model            Model;
    typedef typename Model::Kmer                           KmerType;
    typedef typename Kmer<span>::SuperKmer                 SuperKmer;

    /** */
    void processSuperkmer (SuperKmer& superKmer)
    {
        DEBUG (("SampleRepart: should count superk %i \n", superKmer.size()));

        if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid() ) //check if falls into pass
        {
            bool prev_which = superKmer[0].which();
            size_t kx_size = 0;

            /** Shortcut. */
            size_t superKmerLen = superKmer.size();

            /** We increase superkmer counter the current minimizer. */
            _local_pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmerLen);

            /** We loop over the kmer of the superkmer (except the first one).
             *  We update the pInfo each time we find a kxmer in the superkmer. */
            for (size_t ii=1 ; ii < superKmerLen; ii++)
            {
                /** A kxmer is defined by having successive canonical kmers. Here, we just care that
                 * successive kmer values are on the same strand. */
                if (superKmer[ii].which() != prev_which || kx_size >= _kx) // kxmer_size = 1 //cost should diminish with larger kxmer
                {
                    /** We increase the number of kxmer found for the current minimizer. */
                    _local_pInfo.incKxmer_per_minimBin (superKmer.minimizer);
                    kx_size = 0;
                }
                else
                {
                    kx_size++;
                }

                prev_which = superKmer[ii].which() ;
            }

            /** We add the pending kxmer to the bin. */
            _local_pInfo.incKxmer_per_minimBin (superKmer.minimizer);
        }
    }

    /** Constructor. */
    SampleRepart (
        Model&            model,
        Configuration&    config,
        size_t            nbPasses,
        size_t            currentPass,
        size_t            nbPartitions,
        IteratorListener* progress,
        BankStats&        bankStats,
        PartiInfo<5>&     pInfo
    )
    :   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats)
        ,_kx(4), _extern_pInfo(pInfo), _local_pInfo(config._nb_partitions, model.getMmersModel().getKmerSize())
    {
    }

    /** Destructor. */
    ~SampleRepart ()
    {
        //add to global parti_info
        _extern_pInfo += _local_pInfo;
    }


private:
    size_t        _kx;
    PartiInfo<5>& _extern_pInfo;
    PartiInfo<5>  _local_pInfo;
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
RepartitorAlgorithm<span>::RepartitorAlgorithm (
    IBank* bank,
    Group& group,
    const Configuration& config
)
    :  Algorithm("repartition"), _config(config), _bank(bank), _group(group), _freq_order(0)
{
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
RepartitorAlgorithm<span>::~RepartitorAlgorithm ()
{
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
void RepartitorAlgorithm<span>::execute ()
{
    /** We compute the distribution of the minimizers. As a result, we will have a hash function
     * that gives a hash code for a minimizer value.
     * IMPORTANT ! we have to give the passes number because it has impact on the computation. */
    Repartitor repartitor (_config._nb_partitions, _config._minim_size, _config._nb_passes);

    /* now is a good time to switch to frequency-based minimizers if required:
      because right after we'll start using minimizers to compute the distribution
      of superkmers in bins */
    if (_config._minimizerType == 1)  {  computeFrequencies (repartitor);  }

    computeRepartition (repartitor);
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
void RepartitorAlgorithm<span>::computeFrequencies (Repartitor& repartitor)
{
    DEBUG (("RepartitorAlgorithm<span>::computeFrequencies\n"));

    u_int64_t estimateSeqNb;
    u_int64_t estimateSeqTotalSize;
    u_int64_t estimateSeqMaxSize;

    _bank->estimate (estimateSeqNb, estimateSeqTotalSize, estimateSeqMaxSize);

    u_int64_t nbseq_sample = std::max ( u_int64_t (estimateSeqNb * 0.05) ,u_int64_t( 1000000ULL) ) ;

    u_int64_t rg = ((u_int64_t)1 << (2*_config._minim_size));
    //cout << "\nAllocating " << ((rg*sizeof(uint32_t))/1024) << " KB for " << _minim_size <<"-mers frequency counting (" << rg << " elements total)" << endl;
    uint32_t *m_mer_counts = new uint32_t[rg];

    Model model (_config._kmerSize, _config._minim_size);

    /** We create a sequence iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    // can we reuse the it_sample variable above?
    Iterator<Sequence>* it_sample = createIterator (
        new TruncateIterator<Sequence> (*itSeq, nbseq_sample),
        nbseq_sample,
        "Approximating frequencies of minimizers"
    );
    LOCAL (it_sample);

    /** We compute an estimation of minimizers frequencies from a part of the bank. */
    // actually.. let's try with the whole thing (itSeq instead of it_sample)
    getDispatcher()->iterate (it_sample,  MmersFrequency<span> (
        _config._minim_size, 0 /*_progress*/, m_mer_counts)
    );

    // single threaded, for debugging
    /*MmersFrequency<span> mmersfrequency(model, _progress, bstatsDummy, m_mer_counts);
    it_sample->iterate(mmersfrequency);*/

    /* sort frequencies */
    for (u_int64_t i(0); i < rg; i++)
    {
        if (m_mer_counts[i] > 0)
            _counts.push_back(make_pair(m_mer_counts[i],i));
    }
    delete[] m_mer_counts;

    sort (_counts.begin(), _counts.end());

    /* assign frequency to minimizers */
    _freq_order = new uint32_t[rg];

    for (u_int64_t i = 0; i < rg ; i++)
        _freq_order[i] = rg; // set everything not seen to highest value (not a minimizer)

    for (unsigned int i = 0; i < _counts.size(); i++)
    {
        _freq_order[_counts[i].second] = i;
    }

    // small but necessary trick: the largest minimizer has to have largest rank, as it's used as the default "largest" value
    _freq_order[rg-1] = rg-1;

    repartitor.setMinimizerFrequencies (_freq_order);
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
void RepartitorAlgorithm<span>::computeRepartition (Repartitor& repartitor)
{
    DEBUG (("RepartitorAlgorithm<span>::computeRepartition\n"));

    /** We create a kmer model; using the frequency order if we're in that mode */
    Model model (_config._kmerSize, _config._minim_size, typename Kmer<span>::ComparatorMinimizerFrequency(), _freq_order);

    int mmsize = model.getMmersModel().getKmerSize();

    PartiInfo<5> sample_info (_config._nb_partitions, mmsize);

    /** We create a sequence iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    u_int64_t nbseq_sample = std::max ( u_int64_t (_config._estimateSeqNb * 0.05) ,u_int64_t( 1000000ULL) ) ;

    string bankShortName = System::file().getBaseName(_bank->getId());

    /** We create an iterator over a truncated part of the input bank. */
    Iterator<Sequence>* it_sample = createIterator (
        new TruncateIterator<Sequence> (*itSeq, nbseq_sample),
        nbseq_sample,
        Stringify::format (progressFormat0, bankShortName.c_str()).c_str()
    );
    LOCAL (it_sample);

    BankStats bstatsDummy;

    size_t  currentPass = 0;

    /** We compute a distribution of Superkmers from a part of the bank. */
    getDispatcher()->iterate (it_sample,  SampleRepart<span> (
        model,
        _config,
        1, // we don't care about the actual number of passes, we just use 1
        0, // we don't care about the actual number of passes, the current one is 0
        _config._nb_partitions,
        0,
        bstatsDummy,
        sample_info
    ));

    if (_config._minimizerType == 1)
    {
        repartitor.justGroup (sample_info, _counts);
    }
    else
    {
        repartitor.computeDistrib (sample_info);
        if (_config._repartitionType == 1)
        {
            repartitor.justGroupLexi (sample_info); // For bcalm, i need the minimizers to remain in order. so using this suboptimal but okay repartition
        }
    }

    /** We save the distribution (may be useful for debloom for instance). */
    repartitor.save (getGroup());
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class RepartitorAlgorithm <KSIZE_1>;
template class RepartitorAlgorithm <KSIZE_2>;
template class RepartitorAlgorithm <KSIZE_3>;
template class RepartitorAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
