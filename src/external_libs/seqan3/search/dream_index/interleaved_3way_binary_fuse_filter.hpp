/*!\file
 * \author Jens-Uwe Ulrich <jens-uwe.ulrich AT hpi.de>
 * \brief Provides seqan3::interleaved_3way_binary_fuse_filter.
 */

#pragma once
#define __STDC_FORMAT_MACROS

#include <algorithm>
#include <bit>

#include <sdsl/bit_vectors.hpp>
#include <seqan3/core/detail/strong_type.hpp>

#include <bitset>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace seqan3
{

/*!\brief The IXF binning directory. A data structure that efficiently answers set-membership queries for multiple bins.
 * \ingroup search_dream_index
 * \implements seqan3::cerealisable
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) that can answer these queries efficiently.
 *
 * ### Interleaved XOR Filter (IXF)
 *
 * TODO: add description
 * The implementation of the XOR filter construction and query is based on the fastfilter_cpp github repository
 * https://github.com/FastFilter/fastfilter_cpp
 *
 * The Interleaved XOR Filter now applies the concept of a XOR Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in `b` data sets/bins.
 * Conceptually, a XOR Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting `b` XOR Filters are then interleaved such that the `i`'th l-bits if each XOR Filter are
 * adjacent to each other, with l being 8 or 16 . For l=2: 
 * ```
 *  XOR Filter 0         XOR Filter 1        XOR Filter 2        XOR Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where `x.y` denotes the `y`'th bit of the `x`'th XOR Filter.
 * ```
 * Interleaved XOR Filter
 * |0.0|0.1|1.0|1.1|2.0|2.1|3.0|3.1|0.2|0.3|1.2|1.3|2.2|2.3|3.2|3.3|
 * ```
 * A query can now be searched in all `b` bins by computing the `h` hash functions, retrieving the `h` sub-bitvectors of
 * length `b * l` starting at the positions indicated by the hash functions. The bitwise XOR of these sub-bitvectors yields
 * the binningvector, a bitvector of length `b * l` where the `i*l`'th bits indicate set membership in the `i`'th bin iff all 
 * `l` bits of `i` are 0.
 *
 * ### Querying
 * To query the Interleaved XOR Filter for a value, call seqan3::interleaved_3way_binary_fuse_filter::membership_agent() and use
 * the returned seqan3::interleaved_3way_binary_fuse_filter::membership_agent.
 *
 * To count the occurrences of a range of values in the Interleaved XOR Filter, call
 * seqan3::interleaved_3way_binary_fuse_filter::counting_agent() and use
 * the returned seqan3::interleaved_3way_binary_fuse_filter::counting_agent_type.
 *
 *
 * ### Thread safety
 *
 * The Interleaved XOR Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time). Furthermore, XOR Filters are immutable data structures
 * by itself.
 *
 * 
 */
template <typename FingerprintType = uint8_t>
//!\cond
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
//!\endcond
class interleaved_3way_binary_fuse_filter
{
private:
    struct t2val {
        uint64_t t2 = 0;
        uint64_t t2count = 0;
    };

    using t2val_t = t2val;
    static constexpr int blockShift = 18;

    using data_type = sdsl::int_vector<>;

    //!\brief The number of bins specified by the user.
    size_t bins{};
    //!\brief The number of bins stored in the IXF (next multiple of 64 of `bins`).
    //size_t technical_bins{};
    //!\brief The size of each bin in bits.
    size_t bin_size_{};
    //!\brief number of elements that can be stored in each bin.
    size_t max_bin_elements{};
    //!\brief The number of 64-bit integers needed to store `bins` many bits (e.g. `bins = 50` -> `bin_words = 1`).
    size_t bin_words{};
    //!\brief Segment length for binary fuse filter sizing
    size_t segment_length{};
    //!\brief Segment count for binary fuse filter sizing
    size_t segment_count{};
    //!\brief Segment count length for binary fuse filter sizing
    size_t segment_count_length{};
    //!\brief The int vector of fingerprints
    data_type data{};
    //!\brief number of bits used to store one hashed item in the XOR filter
    size_t ftype{};
    //!\brief seed for  hashing
    size_t seed{13572355802537770549ULL};
    //!\brief number of bins to query in parallel
    size_t bins_per_batch{};
    //!\brief Stash for each bin to handle difficult-to-place elements (binary fuse filter feature)
    std::vector<std::vector<uint64_t>> bin_stashes{};
    //!\brief Maximum number of elements allowed in stash per bin
    uint32_t max_stash{1};

    /*!\brief Utilizes hashing of a 64-bit key.

     * \param   64-bit key to hash
     * \returns a 64-bit hash value for the given key
     * 
     */
    inline constexpr uint64_t murmur64(uint64_t h) const
    {
        h += seed;
        h ^= h >> 33;
        h *= UINT64_C(0xff51afd7ed558ccd);
        h ^= h >> 33;
        h *= UINT64_C(0xc4ceb9fe1a85ec53);
        h ^= h >> 33;
        return h;
    }

#ifdef __SIZEOF_INT128__
    static inline constexpr uint64_t binary_fuse_mulhi(uint64_t a, uint64_t b) {
        return (uint64_t)(((__uint128_t)a * b) >> 64U);
    }
#else
    static inline constexpr uint64_t binary_fuse_mulhi(uint64_t a, uint64_t b) {
        const uint64_t a0 = (uint32_t) a;
        const uint64_t a1 = a >> 32;
        const uint64_t b0 = (uint32_t) b;
        const uint64_t b1 = b >> 32;
        const uint64_t p11 = a1 * b1;
        const uint64_t p01 = a0 * b1;
        const uint64_t p10 = a1 * b0;
        const uint64_t p00 = a0 * b0;

        const uint64_t middle = p10 + (p00 >> 32) + (uint32_t) p01;
        return p11 + (middle >> 32) + (uint32_t)((p01 & 0xFFFFFFFF) + (uint32_t) middle) < (uint32_t) p01;
    }
#endif

    /*!\brief   Fingerprint funtion of the underlying filter
     * \param   64-bit hash value
     *  
     * \returns an unsigned 8-bit or 16-bit unsigned integer fingerprint for 
     * storing hash values in the filter
     * 
     */
    inline constexpr FingerprintType fingerprint(const uint64_t hash) const 
    {
        //TODO this might have a role in preventing FP spike
        /*assert(hash != 0);
        FingerprintType h = (FingerprintType) hash ^ (hash >> 32);
        if (h == 0)
            h = (FingerprintType) hash ^ (hash >> 16);
        if (h == 0)
            h = (FingerprintType) hash ^ (hash >> 8);
        
        assert(h != 0);

        return h;*/
        return (FingerprintType) hash ^ (hash >> 32);
    }

    /*!\brief Calculate hash indices for a given 64-bit hash value
     * \param hash 64-bit hash value
     * \returns array of 3 indices mapped to segments
     */
    inline constexpr std::array<size_t, 3> get_hashes(uint64_t hash) const
    {
        std::array<size_t, 3> ans;
        uint64_t hi = binary_fuse_mulhi(hash, segment_count_length);
        ans[0] = (uint32_t)hi;
        ans[1] = ans[0] + segment_length;
        ans[2] = ans[1] + segment_length;
        ans[1] ^= (uint32_t)(hash >> 18) & (segment_length - 1);
        ans[2] ^= (uint32_t)(hash) & (segment_length - 1);
        return ans;
    }

    /*!\brief Calculate segment length for binary fuse filter sizing (3-way)
     * \param size Number of elements
     * \returns Segment length as power of 2
     */
    inline constexpr uint32_t calculate_segment_length(uint32_t size) const
    {
        return ((uint32_t)1) << (unsigned)(std::floor(std::log((double)(size)) / std::log(3.33) + 2.25));
    }

    /*!\brief Calculate size factor for binary fuse filter (3-way)
     * \param size Number of elements
     * \returns Size factor for capacity calculation
     */
    inline constexpr double calculate_size_factor(uint32_t size) const
    {
        return std::max(1.125, 0.875 + 0.25 * std::log(1000000.0) / std::log((double)size));
    }

    public:
    /*!\brief Check if a key exists in a bin's stash using binary search
     * \param bin Bin index
     * \param key Key to search for
     * \returns true if key is in stash, false otherwise
     */
    inline bool stash_contains(size_t bin, uint64_t key) const
    {
        if (bin >= bin_stashes.size() || bin_stashes[bin].empty())
            return false;
        
        const auto& stash = bin_stashes[bin];
        return std::binary_search(stash.begin(), stash.end(), key);
    }
    private:

    
    /*!\brief Count number of occurences of each index in a bin
     * \param tmp pointer to an sdsl::int_vector<> data structure
     * \param b Block number in which we count
     * \param len corresponds to the number of keys that match to the index
     * \param t2vals data structure to store counts for each index
     *  
     */
    void applyBlock(uint64_t* tmp, int b, int len, t2val_t* t2vals)
    {
        for (int i = 0; i < len; i += 2) {
            uint64_t x = tmp[(b << blockShift) + i];
            int index = (int) tmp[(b << blockShift) + i + 1];
            t2vals[index].t2count++;
            t2vals[index].t2 ^= x;
        }
    }

    /*!\brief Remove index from index occurence array and add index to single entry set
     * \param tmp pointer to an sdsl::int_vector<> data structure
     * \param b Block number in which we count
     * \param len corresponds to the number of keys that match to the index
     * \param t2vals data structure to store counts for each index
     * \param alone reference to an array of indexes that occure only once in t2vals
     * \param alonePos number of indexes that occur only once in t2vals
     *  
     * remove all indexes in tmp array from index occurence array and remove corresponding hashes by XOR-ing with the hash value
     *  add index to single entry set if an index now occurs only once in the remaining set of indexes
     */
    int applyBlock2(uint64_t* tmp, int b, int len, t2val_t*  t2vals, sdsl::int_vector<>& alone, int alonePos) 
    {
        for (int i = 0; i < len; i += 2) {
            uint64_t hash = tmp[(b << blockShift) + i];
            int index = (int) tmp[(b << blockShift) + i + 1];
            int oldCount = t2vals[index].t2count;
            if (oldCount >= 1) {
                int newCount = oldCount - 1;
                t2vals[index].t2count = newCount;
                if (newCount == 1) {
                    alone[alonePos++] = index;
                }
                t2vals[index].t2 ^= hash;
            }
        }
        return alonePos;
    }

    /*!\brief Find all index positions to which only one hashed key maps
     * \param elements list of keys to store in the XOR filter
     * \param t2vals data structure to store counts for each index
     * \param alone_positions reference to an array of indexes that occure only once
     *  
     */
    int find_alone_positions(std::vector<size_t>& elements, std::vector<t2val_t>& t2vals, sdsl::int_vector<>& alone_positions)
    {
        std::fill(t2vals.begin(), t2vals.end(), t2val_t{ 0,0 });
        // number of elements / 2^18 => if more than 2^18 elements, we need 2 blocks
        int blocks = 1 + (bin_size_ >> blockShift);
        sdsl::int_vector<> tmp = sdsl::int_vector(blocks << blockShift, 0, 64);       
        sdsl::int_vector<> tmpc = sdsl::int_vector(blocks,0);
        
        for(size_t k : elements) {
            uint64_t hash = murmur64(k);
            std::array<size_t, 3> h_arr = get_hashes(hash);
            for (int hi = 0; hi < 3; hi++) {
                int index = h_arr[hi];
                int b = index >> blockShift;
                int i2 = tmpc[b];
                tmp[(b << blockShift) + i2] = hash;
                tmp[(b << blockShift) + i2 + 1] = index;
                tmpc[b] += 2;
                if (i2 + 2 == (1 << blockShift)) {
                    applyBlock(tmp.data(), b, i2 + 2, t2vals.data());
                    tmpc[b] = 0;
                }
            }
        }
       
        // count occurences of index positions for all computed hash values
        for (int b = 0; b < blocks; b++) {
            applyBlock(tmp.data(), b, tmpc[b], t2vals.data());
        }

        // pick only index positions where only one unique hash value points to => those are our start positions
        int alonePos = 0;
        for (size_t i = 0; i < bin_size_; i++) 
        {
            if (t2vals[i].t2count == 1) {
                alone_positions[alonePos++] = i;
            }
        }
        return alonePos;
    }

    /*!\brief Determine the order in which keys are added to the filter
     * \param reverse_order resulting stack of keys to insert in the filter
     * \param reverse_h stack of hash function numbers
     * \param t2vals data structure to store counts for each index
     * \param alone_positions reference to an array of indexes that occure only once
     * \param size number of keys to insert
     * \param alone_pos number of indexes that occur only once in t2vals
     * 
     * \returns the number of keys in the stack
     *  
     */
    size_t fill_stack(sdsl::int_vector<>& reverse_order, 
                      sdsl::int_vector<>& reverse_h, 
                      std::vector<t2val_t>& t2vals, 
                      sdsl::int_vector<>& alone_positions,
                      size_t size, 
                      int alone_pos)
    {
        int blocks = 1 + (bin_size_ >> blockShift);
        sdsl::int_vector<> tmp = sdsl::int_vector(blocks << blockShift, 0, 64);
        sdsl::int_vector<> tmpc = sdsl::int_vector(blocks, 0);
        size_t reverse_order_pos = 0;
        int best_block = 0;
       
        while (reverse_order_pos < size)
        {
            if (alone_pos == 0) 
            {
                // we need to apply blocks until we have an entry that is alone
                // (that is, until alonePos > 0)
                // so, find a large block (the larger the better)
                // but don't need to search very long
                // start searching where we stopped the last time
                // (to make it more even)
                
                for (int i = 0, b = best_block, best = -1; i < blocks; i++)
                {
                    if (b >= blocks)
                    {
                        b = 0;
                    }
                    if (tmpc[b] > best) 
                    {
                        best = tmpc[b];
                        best_block = b;
                        if (best > (1 << (blockShift - 1))) 
                        {
                            // sufficiently large: stop
                            break;
                        }
                    }
                }
                
                if (tmpc[best_block] > 0) {
                    alone_pos = applyBlock2(tmp.data(), best_block, tmpc[best_block], t2vals.data(), alone_positions, alone_pos);
                    tmpc[best_block] = 0;
                }
                // applying a block may not actually result in a new entry that is alone
                if (alone_pos == 0) {
                    for (int b = 0; b < blocks && alone_pos == 0; b++) {
                        if (tmpc[b] > 0) {
                            alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                            tmpc[b] = 0;
                        }
                    }
                }
            }
            if (alone_pos == 0) {
                break;
            }
            
            int i = alone_positions[--alone_pos];
            int b = i >> blockShift;
            if (tmpc[b] > 0) {
                alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                tmpc[b] = 0;
            }
            
            uint8_t found = -1;
            if (t2vals[i].t2count == 0 || t2vals[i].t2count > 1) {
                continue;
            }
            long hash = t2vals[i].t2;
            
            std::array<size_t, 3> h_arr = get_hashes(hash);
            for (int hi = 0; hi < 3; hi++) {
                int h = h_arr[hi];
                if (h == i) {
                    found = (uint8_t) hi;
                    t2vals[i].t2count = 0;
                } else {
                    int b = h >> blockShift;
                    int i2 = tmpc[b];
                    tmp[(b << blockShift) + i2] = hash;
                    tmp[(b << blockShift) + i2 + 1] = h;
                    tmpc[b] +=  2;
                    if (tmpc[b] >= 1 << blockShift) {
                        alone_pos = applyBlock2(tmp.data(), b, tmpc[b], t2vals.data(), alone_positions, alone_pos);
                        tmpc[b] = 0;
                    }
                }
            }
            
            reverse_order[reverse_order_pos] = hash;
            reverse_h[reverse_order_pos] = found;
            reverse_order_pos++;
        }
        return reverse_order_pos;
    }

    /*!\brief Adds a 8 or 16 bit value to the xor filter for each hashed key
     * \param reverse_order stack of keys to insert in the filter
     * \param reverse_h stack of hash function numbers
     * \param bin number of the interleaved XOR filter
     * 
     *  
     */
    void fill_filter(sdsl::int_vector<>& reverse_order, sdsl::int_vector<>& reverse_h, uint bin)
    {
        
        for (int i = reverse_order.size() - 1; i >= 0; i--) 
        {
            // the hash of the key we insert next
            uint64_t hash = reverse_order[i];
            int found = reverse_h[i];
            // which entry in the table we can change
            int change = -1;
            // we set table[change] to the fingerprint of the key,
            // unless the other two entries are already occupied
            FingerprintType xor2 = fingerprint(hash);
            std::array<size_t, 3> h_arr = get_hashes(hash);
            for (int hi = 0; hi < 3; hi++) 
            {
                size_t h = h_arr[hi];
                if (found == hi) 
                {
                    change = h;
                } 
                else 
                {
                    uint64_t idx = bins * h;
                    idx += bin;
                    xor2 ^= data[idx];

                }
            }
            uint64_t idx = bins * change;
            idx += bin;
           
            data[idx] = xor2;
        }
    }
    
    /*!\brief Adds all keys of all bins in one go
     * \param elements two-dimensional list of bins and their keys
     *  
     */
    void add_elements(std::vector<std::vector<size_t>>& elements)
    {
        //std::cerr << "add_elements called with " << elements.size() << " bins" << std::flush;
        // stack sigma
        // order in which elements will be inserted into their corresponding bins
        std::vector<sdsl::int_vector<>> reverse_orders;
        // order in which hash seeds are used for element insertion
        std::vector<sdsl::int_vector<>> reverse_hs;
        size_t reverse_order_pos;
        // repeat until all xor filters can be build with same seed
        uint32_t iteration = 0;
        while (true)
        {
            iteration++;
            if (iteration % 10 == 0) {
                std::cerr << " [retry " << iteration << "]" << std::flush;
            }
            // sets the same seed for all xor filters
            int i = 0;
            bool success = true;
            reverse_hs.clear();
            reverse_order_pos = 0;
            reverse_orders.clear();
            for (std::vector<size_t> vec : elements)
            {
                //if (vec.size() == 0)
                //    continue;

                std::vector<t2val_t> t2vals_vec(bin_size_);
                sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
                int alone_position_nr = find_alone_positions(vec, t2vals_vec, alone_positions);
               
                sdsl::int_vector<> rev_order_i = sdsl::int_vector(vec.size(),0,64);
                sdsl::int_vector<> reverse_hi = sdsl::int_vector(vec.size(),0,8);
                reverse_order_pos = fill_stack(std::ref(rev_order_i), std::ref(reverse_hi), t2vals_vec, alone_positions, vec.size(), alone_position_nr);
                reverse_orders.emplace_back(std::move(rev_order_i));
                reverse_hs.emplace_back(std::move(reverse_hi));
                if (reverse_order_pos != vec.size())
                {
                    success = false;
                    set_seed();
                    break;
                }

                i++;
            }
            if (success)
                break;
        }
        
        for (size_t i = 0; i < elements.size(); ++i)
        {
            //if (elements[i].size() == 0)
            //    continue;
            fill_filter(reverse_orders[i], reverse_hs[i], i);
        }
        
    }

public:

    class membership_agent; // documented upon definition below
    template <std::integral value_t>
    class counting_agent_type; // documented upon definition below

    /*!\name Constructors, destructor and assignment
     */
    interleaved_3way_binary_fuse_filter() = default; //!< Defaulted.
    interleaved_3way_binary_fuse_filter(interleaved_3way_binary_fuse_filter const &) = default; //!< Defaulted.
    interleaved_3way_binary_fuse_filter & operator=(interleaved_3way_binary_fuse_filter const &) = default; //!< Defaulted.
    interleaved_3way_binary_fuse_filter(interleaved_3way_binary_fuse_filter &&) = default; //!< Defaulted.
    interleaved_3way_binary_fuse_filter & operator=(interleaved_3way_binary_fuse_filter &&) = default; //!< Defaulted.
    ~interleaved_3way_binary_fuse_filter() = default; //!< Defaulted.

    /*!\brief Construct an Interleaved XOR Filter from given sets of keys
     * \param elements two-dimensional list of bins and keys to store in the IXF
     *
     * \attention This constructor can only be used to construct Interleaved XOR Filters from 
     * the given sets of keys. Further adding of bins/keys will corrupt the data structure.
     *
     * \details
     *
     * ### Example
     * 
     */
    interleaved_3way_binary_fuse_filter(std::vector<std::vector<size_t>>& elements, uint32_t max_stash_ = 1)
    {
        max_stash = max_stash_;
        bins = elements.size();
        int idx = 0;
        for (std::vector<size_t> v : elements)
        {
            if (v.size() > max_bin_elements)
                max_bin_elements = v.size();
        }
        
        // Binary fuse filter sizing
        uint32_t arity = 3;
        segment_length = max_bin_elements == 0 ? 4 : calculate_segment_length(max_bin_elements);
        if (segment_length > 262144) segment_length = 262144;
        double size_factor = max_bin_elements <= 1 ? 0 : calculate_size_factor(max_bin_elements);
        uint32_t capacity = max_bin_elements <= 1 ? 0 : (uint32_t)(std::round((double)max_bin_elements * size_factor));
        uint32_t init_segment_count = (capacity + segment_length - 1) / segment_length - (arity - 1);
        uint32_t array_length = (init_segment_count + arity - 1) * segment_length;
        this->segment_count = (array_length + segment_length - 1) / segment_length;
        if (this->segment_count <= arity - 1) this->segment_count = 1;
        else this->segment_count = this->segment_count - (arity - 1);
        array_length = (this->segment_count + arity - 1) * segment_length;
        this->segment_count_length = this->segment_count * segment_length;
        bin_size_ = array_length;
        
        ftype = CHAR_BIT * sizeof(FingerprintType);
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (bin_size_ == 0)
            throw std::logic_error{"The size of a bin must be > 0."};


        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
        }
        else
        {
            bin_words = (bins + 15) >> 2; // = ceil(bins/4)
        }

        
        data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
        bin_stashes.resize(bins);

        bool success = false;
        bool has_failed = false;
        
        while (!success) {
            success = true;
            
            // Use add_bin_elements() which supports stash for difficult-to-place hashes
            for (size_t bin_idx = 0; bin_idx < elements.size(); ++bin_idx)
            {
                success = add_bin_elements(bin_idx, elements[bin_idx]);

                if (!success) {
                    has_failed = true;
                    std::cerr << " [Failed on IBFF3, bin id: (" << bin_idx << ")] ... " << std::flush;
                    clear();
                    for (auto & stash : bin_stashes) {
                        stash.clear();
                    }
                    set_seed();
                    break;
                }
            }
        }
        if (has_failed) {
            std::cerr << "Success!\n";
        }        
    }

    /*!\brief Construct an Interleaved XOR Filter object for later adding of bins
     * \param bins_ The number of bins.
     * \param size  maximum number of elements to store in a single filter
     *
     * \attention This constructor should be used if large numbers of bins and keys shall be
     * stored in an interleaved XOR filter. It creates the necessary data structure based
     * on the given number of bins and elements, which can be added later on by utilizing
     * add_bin_elements() function
     *
     * \details
     *
     * ### Example
     * 
     */
    interleaved_3way_binary_fuse_filter(size_t bins_, size_t max_bin_elements_)
    {
        this->bins = bins_;
        this->max_bin_elements = max_bin_elements_;
        
        // Binary fuse filter sizing
        uint32_t arity = 3;
        segment_length = max_bin_elements == 0 ? 4 : calculate_segment_length(max_bin_elements);
        if (segment_length > 262144) segment_length = 262144;
        double size_factor = max_bin_elements <= 1 ? 0 : calculate_size_factor(max_bin_elements);
        uint32_t capacity = max_bin_elements <= 1 ? 0 : (uint32_t)(std::round((double)max_bin_elements * size_factor));
        uint32_t init_segment_count = (capacity + segment_length - 1) / segment_length - (arity - 1);
        uint32_t array_length = (init_segment_count + arity - 1) * segment_length;
        this->segment_count = (array_length + segment_length - 1) / segment_length;
        if (this->segment_count <= arity - 1) this->segment_count = 1;
        else this->segment_count = this->segment_count - (arity - 1);
        array_length = (this->segment_count + arity - 1) * segment_length;
        this->segment_count_length = this->segment_count * segment_length;
        bin_size_ = array_length;
        
        ftype = CHAR_BIT * sizeof(FingerprintType);
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (max_bin_elements == 0)
            throw std::logic_error{"The number of elements to store must be > 0."};


        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
        }
        else
        {
            bin_words = (bins + 15) >> 2; // = ceil(bins/4)
        }

        data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
        bin_stashes.resize(bins);
    }

    bool operator==(interleaved_3way_binary_fuse_filter const & rhs) const noexcept
    {
        return bins == rhs.bins && 
               bin_size_ == rhs.bin_size_ && 
               seed == rhs.seed && 
               data == rhs.data &&
               bin_stashes == rhs.bin_stashes;
    }

    bool operator!=(interleaved_3way_binary_fuse_filter const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    bool operator<(interleaved_3way_binary_fuse_filter const & rhs) const noexcept
    {
        if (bins != rhs.bins) return bins < rhs.bins;
        if (bin_size_ != rhs.bin_size_) return bin_size_ < rhs.bin_size_;
        if (seed != rhs.seed) return seed < rhs.seed;
        if (data.size() != rhs.data.size()) return data.size() < rhs.data.size();
        if (data != rhs.data) return std::lexicographical_compare(data.begin(), data.end(), rhs.data.begin(), rhs.data.end());
        return bin_stashes < rhs.bin_stashes;
    }

    /*!\brief Returns the total number of elements currently in the stash across all bins.
     * \returns Total stash count
     */
    size_t stash_usage() const
    {
        size_t count = 0;
        for (const auto& stash : bin_stashes)
            count += stash.size();
        return count;
    }

    /*!\\brief Returns a const reference to all bin stashes for extraction.
     * \\returns Const reference to the vector of stashes (one per bin)
     */
    std::vector<std::vector<uint64_t>> const & get_bin_stashes() const
    {
        return bin_stashes;
    }

    /*!\brief add all keys of a vector/list to the given bin number
     * \param bin       bins number of the interleaved XOR filter
     * \param elements  elements to store in the given bin
     * \returns `true` if elements could be added successfully, `false` otherwise.
     * 
     * \attention If this method returns `false`, you have to clear the filter (@see clear()),
     * reset the seed (@see set_seed()) and repeat adding all elements of all bins. 
     *
     * \details
     *
     * ### Example
     * 
     *  seqan3::interleaved_3way_binary_fuse_filter<> ixf(100, 2000);
	 *  std::vector<uint64_t> elems{};
	 *  while (true)
	 *  {
     *      bool success = true;
	 *	    for (int e = 0; e < 100 ; ++e)
	 *	    {
	 *		    std::vector<uint64_t> tmp{};
	 *		    for (uint64_t i = 0; i < 2000; ++i)
	 *		    {
	 * 			    uint64_t key = (e*2000) + i;
	 *			    tmp.emplace_back(key);
	 *		    }
	 *		    success = ixf.add_bin_elements(e, tmp);
	 *		    if (!success)
	 *		    {
	 *			    ixf.clear();
	 *			    ixf.set_seed();
	 *			    break;
	 *		    }
	 *		    if (e == 2)
	 *			    elems=std::move(tmp);
	 *	    }
     *
	 *	    if (success)
	 *		    break;
	 *  }
     * 
     */
    bool add_bin_elements(size_t bin, std::vector<size_t>& elements)
    {
        // Binary fuse filter construction with victim selection and stash support
        
        // Key-level deduplication: sort and remove duplicate keys before graph construction.
        // Duplicate keys produce identical hashes that XOR to 0 at all positions,
        // making them impossible to peel and forcing unnecessary stashing.
        std::sort(elements.begin(), elements.end());
        auto last_unique = std::unique(elements.begin(), elements.end());
        elements.erase(last_unique, elements.end());
        
        size_t size = elements.size();
        //std::cerr << "Starting add_bin_elements for bin " << bin << " with " << size << " elements" << std::flush;
        if (size == 0) return true;
        
        int64_t max_seed_retries = 1; // -1 for infinite retries
        
        // Retry indefinitely if max_Seed_retries = -1, else retry max_seed_retries times
        for (int64_t seed_retry = 0; max_seed_retries < 0 || seed_retry < max_seed_retries; ++seed_retry) {
            if (seed_retry > 0) {
                set_seed();
                std::cerr << " [seed retry " << seed_retry << "]" << std::flush;
            }
            
            // Clear stash for this bin
            bin_stashes[bin].clear();
            
            // Allocate auxiliary arrays
            std::vector<t2val_t> t2vals_vec(bin_size_);
            sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
            sdsl::int_vector<> rev_order = sdsl::int_vector(size, 0, 64);
            sdsl::int_vector<> reverse_h = sdsl::int_vector(size, 0, 8);
            std::vector<uint8_t> key_status(size, 0); // 0 = not processed, 1 = processed/stashed
            
            // Find initial alone positions
            int alone_position_nr = find_alone_positions(elements, t2vals_vec, alone_positions);
            
            // Graph-level duplicate detection using XOR-cancellation.
            // After graph construction, if two identical hashes mapped to the same positions,
            // they XOR to 0. We detect this by checking: if at any of the 3 positions the
            // t2 hash is 0 and the count is exactly 2, those two entries are duplicates.
            // We remove the duplicate from the graph and track it separately.
            uint32_t duplicates = 0;
            for (size_t ki = 0; ki < size; ki++) {
                uint64_t hash = murmur64(elements[ki]);
                std::array<size_t, 3> h_arr = get_hashes(hash);
                uint32_t h0 = h_arr[0], h1 = h_arr[1], h2 = h_arr[2];
                
                // Check if all 3 positions have t2 XOR'd to 0 (possible duplicate pair)
                if ((t2vals_vec[h0].t2 & t2vals_vec[h1].t2 & t2vals_vec[h2].t2) == 0) {
                    if ((t2vals_vec[h0].t2 == 0 && t2vals_vec[h0].t2count == 2)
                     || (t2vals_vec[h1].t2 == 0 && t2vals_vec[h1].t2count == 2)
                     || (t2vals_vec[h2].t2 == 0 && t2vals_vec[h2].t2count == 2)) {
                        // This is a hash-level duplicate — remove from graph
                        key_status[ki] = 1;
                        duplicates++;
                        t2vals_vec[h0].t2count--;
                        t2vals_vec[h0].t2 ^= hash;
                        if (t2vals_vec[h0].t2count == 1) alone_positions[alone_position_nr++] = h0;
                        t2vals_vec[h1].t2count--;
                        t2vals_vec[h1].t2 ^= hash;
                        if (t2vals_vec[h1].t2count == 1) alone_positions[alone_position_nr++] = h1;
                        t2vals_vec[h2].t2count--;
                        t2vals_vec[h2].t2 ^= hash;
                        if (t2vals_vec[h2].t2count == 1) alone_positions[alone_position_nr++] = h2;
                    }
                }
            }
            
            // Try to fill stack with peeling + victim selection
            size_t stacksize = 0;
            size_t key_cursor = 0;
            uint32_t stash_count = 0;
            uint32_t iteration = 0;
            bool construction_failed = false;
        
            while (true) {
                iteration++;
                if (iteration % 1000 == 0) {
                std::cerr << "Bin " << bin << " iteration " << iteration 
                          << ": stacksize=" << stacksize << ", stash=" << stash_count 
                          << ", key_cursor=" << key_cursor << ", size=" << size 
                          << ", alone_pos=" << alone_position_nr << std::flush;
            }
            
            // Phase 1: Peeling - process alone positions
            while (alone_position_nr > 0) {
                int i = alone_positions[--alone_position_nr];
                
                if (t2vals_vec[i].t2count == 0 || t2vals_vec[i].t2count > 1) {
                    continue;
                }
                
                uint64_t hash = t2vals_vec[i].t2;
                uint8_t found = -1;
                
                // Find which hash function points to this position
                std::array<size_t, 3> h_arr = get_hashes(hash);
                for (int hi = 0; hi < 3; hi++) {
                    int h = h_arr[hi];
                    if (h == i) {
                        found = (uint8_t) hi;
                        t2vals_vec[i].t2count = 0;
                    } else {
                        // Update other positions
                        if (t2vals_vec[h].t2count > 0) {
                            t2vals_vec[h].t2count--;
                            t2vals_vec[h].t2 ^= hash;
                            if (t2vals_vec[h].t2count == 1) {
                                alone_positions[alone_position_nr++] = h;
                            }
                        }
                    }
                }
                
                rev_order[stacksize] = hash;
                reverse_h[stacksize] = found;
                stacksize++;
            }
            
            // Check if done
            if (stacksize + duplicates + stash_count == size) {
                break; // Success!
            }
            
            // Phase 2: Victim selection - find a key to stash
            bool found_victim = false;
            while (key_cursor < size) {
                if (key_status[key_cursor] == 0) {
                    uint64_t key = elements[key_cursor];
                    uint64_t hash = murmur64(key);
                    std::array<size_t, 3> h_arr = get_hashes(hash);
                    uint32_t h0 = h_arr[0];
                    uint32_t h1 = h_arr[1];
                    uint32_t h2 = h_arr[2];
                    
                    // Check if this is a valid victim (all positions still have count >= 1)
                    if (t2vals_vec[h0].t2count >= 1 && t2vals_vec[h1].t2count >= 1 && t2vals_vec[h2].t2count >= 1) {
                        // Stash this key
                        bin_stashes[bin].push_back(key);
                        stash_count++;
                        
                        // Check stash limit
                        if (stash_count > max_stash) {
                            //std::cerr << " [Bin " << bin << " stash overflow: " << stash_count << " > " << max_stash << "]" << std::flush;
                            construction_failed = true;
                            break; // Break inner loop to retry with new seed
                        }
                        
                        key_status[key_cursor] = 1;
                        
                        // Remove from graph
                        t2vals_vec[h0].t2count--;
                        t2vals_vec[h0].t2 ^= hash;
                        if (t2vals_vec[h0].t2count == 1) alone_positions[alone_position_nr++] = h0;
                        
                        t2vals_vec[h1].t2count--;
                        t2vals_vec[h1].t2 ^= hash;
                        if (t2vals_vec[h1].t2count == 1) alone_positions[alone_position_nr++] = h1;
                        
                        t2vals_vec[h2].t2count--;
                        t2vals_vec[h2].t2 ^= hash;
                        if (t2vals_vec[h2].t2count == 1) alone_positions[alone_position_nr++] = h2;
                        
                        found_victim = true;
                        break; // Go back to peeling
                    } else {
                        // Already peeled, mark as processed
                        key_status[key_cursor] = 1;
                    }
                }
                key_cursor++;
            }
            
            if (construction_failed) {
                break; // Exit while(true) to retry with new seed
            }
            
            if (!found_victim && (stacksize + duplicates + stash_count != size)) {
                std::cerr << " [construction incomplete, retrying]" << std::flush;
                construction_failed = true;
                break; // Exit while(true) to retry with new seed
            }
            
            // If we found a victim and we are not done yet, continue the loop back to peeling
            if (found_victim && (stacksize + duplicates + stash_count < size)) {
                continue;
            }
            
            // Construction succeeded - exit while(true)
            break;
        } // end while(true)
        
        if (construction_failed) {
            //std::cerr << " Retry " << (seed_retry + 2);
            //if (max_seed_retries > 0) std::cerr << "/" << max_seed_retries;
            //std::cerr << " with new seed..." << std::endl;
            continue; // Try next seed in for loop
        }
        
        // Success! Sort and deduplicate stash
        if (!bin_stashes[bin].empty()) {
            std::sort(bin_stashes[bin].begin(), bin_stashes[bin].end());
            auto last = std::unique(bin_stashes[bin].begin(), bin_stashes[bin].end());
            bin_stashes[bin].erase(last, bin_stashes[bin].end());
        }
        
        // Resize to actual stacksize so fill_filter doesn't process trailing zeros
        rev_order.resize(stacksize);
        reverse_h.resize(stacksize);

        // Fill filter with successfully peeled elements
        fill_filter(rev_order, reverse_h, bin);
        
        // Print statistics for builds that use stash
        if (!bin_stashes[bin].empty()) {
            std::cerr << "Bin " << bin << " stash build stats: "
                      << "Filter built for: " << max_bin_elements << " elements, "
                      << "Hashes inserted: " << size << ", "
                      << "Hashes in filter: " << stacksize << ", "
                      << "Hashes in stash: " << bin_stashes[bin].size() << std::endl;
        }
        
        return true;
    } // end seed retry loop
    
    // All retries exhausted
    /*
    std::cerr << "BUILD FAILURE - Bin " << bin << ": All ";
    if (max_seed_retries > 0) std::cerr << max_seed_retries << " ";
    std::cerr << "seed retries exhausted\n"
              << "  Filter built for: " << max_bin_elements << " elements\n"
              << "  Hashes attempted to insert: " << size << std::endl; */
    return false;
}
   
    /*!\brief Clears the entire IXF.
     * \
     * \attention This function only sets the entries in the filter to 0. It will not
     * resize the filter or reset the number of bins.
     *
     * \details
     *
     * ### Example
     *
     */
    void clear()
    {
        for (size_t idx = 0; idx < data.size(); ++idx)
            data[idx] = 0;
    }

    /*!\brief Sets a new seed for the hashing function
     * \
     * \attention Setting a new seed requires clearing the entire filter and repeat
     * adding all elements of all bins.
     *
     * \details
     *
     * ### Example
     *
     */
    void set_seed()
    {
        ::std::random_device random;
        seed = random();
        seed <<= 32;
        seed |= random();
    }

    size_t get_seed() const
    {
        return seed;
    }

    void set_max_stash(uint32_t value)
    {
        max_stash = value;
    }

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns a seqan3::interleaved_3way_binary_fuse_filter::membership_agent to be used for lookup.
     * `seqan3::interleaved_3way_binary_fuse_filter::membership_agent`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     */
    class membership_agent membership_agent() const
    {
        using AgentType = class membership_agent;
        return AgentType{*this};
    }

    /*!\brief Returns a seqan3::interleaved_3way_binary_fuse_filter::counting_agent_type to be used for counting.
     * `seqan3::interleaved_3way_binary_fuse_filter::counting_agent_type`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     */
    template <typename value_t = uint32_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
    

    /*!\brief Returns the number of bins that the Interleaved XOR Filter manages.
     * \returns The number of bins.
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }

    /*!\brief Returns the size of a single bin that the Interleaved XOR Filter manages.
     * \returns The size in bits of a single bin.
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return data.bit_size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Test for equality.
     * \param[in] lhs A `seqan3::interleaved_3way_binary_fuse_filter`.
     * \param[in] rhs `seqan3::interleaved_3way_binary_fuse_filter` to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    friend bool operator==(interleaved_3way_binary_fuse_filter const & lhs, interleaved_3way_binary_fuse_filter const & rhs) noexcept
    {
        return std::tie(lhs.bins, lhs.bins_per_batch, lhs.bin_size_, lhs.ftype, lhs.bin_words, lhs.segment_length, lhs.segment_count, lhs.segment_count_length, lhs.max_bin_elements, lhs.max_stash, lhs.bin_stashes,
                        lhs.seed, lhs.data) ==
               std::tie(rhs.bins, rhs.bins_per_batch, rhs.bin_size_, rhs.ftype, rhs.bin_words, rhs.segment_length, rhs.segment_count, rhs.segment_count_length, rhs.max_bin_elements, rhs.max_stash, rhs.bin_stashes,
                        rhs.seed, rhs.data);
    }

    /*!\brief Test for inequality.
     * \param[in] lhs A `seqan3::interleaved_3way_binary_fuse_filter`.
     * \param[in] rhs `seqan3::interleaved_3way_binary_fuse_filter` to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(interleaved_3way_binary_fuse_filter const & lhs, interleaved_3way_binary_fuse_filter const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(bins);
        archive(bin_size_);
        archive(bins_per_batch);
        archive(bin_words);
        archive(ftype);
        if (ftype != (CHAR_BIT * sizeof(FingerprintType)))
        {
            throw std::logic_error{"The interleaved XOR filter was built with a fingerprint size of " + std::to_string(ftype) +
                                   " but it is being read into an interleaved XOR filter with fingerprint of size " +
                                   std::to_string(CHAR_BIT * sizeof(FingerprintType)) + "."};
        }
        archive(segment_length);
        archive(segment_count);
        archive(segment_count_length);
        archive(max_bin_elements);
        archive(max_stash);
        archive(seed);
        archive(data);
        archive(bin_stashes);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the seqan3::interleaved_3way_binary_fuse_filter.
 * 
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/membership_agent_construction.cpp
 */
template <typename FingerprintType>
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
class interleaved_3way_binary_fuse_filter<FingerprintType>::membership_agent
{
private:
    //!\brief The type of the augmented seqan3::interleaved_3way_binary_fuse_filter.
    using ixf_t = interleaved_3way_binary_fuse_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_3way_binary_fuse_filter.
    ixf_t const * ixf_ptr{nullptr};

public:
    class binning_bitvector;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default; //!< Defaulted.
    membership_agent(membership_agent const &) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default; //!< Defaulted.
    ~membership_agent() = default; //!< Defaulted.

    /*!\brief Construct a membership_agent from a seqan3::interleaved_3way_binary_fuse_filter.
     * \private
     * \param ixf The seqan3::interleaved_3way_binary_fuse_filter.
     */
    explicit membership_agent(ixf_t const & ixf) :
        ixf_ptr(std::addressof(ixf)), result_buffer(ixf.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    binning_bitvector result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of a given value.
     * \param[in] value The raw value to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_3way_binary_fuse_filter::membership_agent for each thread.
     */

     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) & noexcept
    {
        assert(ixf_ptr != nullptr);
        assert(result_buffer.size() == ixf_ptr->bin_count());

        size_t bins = ixf_ptr->bin_count();
        uint8_t ftype = ixf_ptr->ftype;
        uint8_t bins_per_batch = ixf_ptr->bins_per_batch;
        uint64_t hash = ixf_ptr->murmur64(value);

        FingerprintType f = ixf_ptr->fingerprint(hash);
        std::array<size_t, 3> h_arr = ixf_ptr->get_hashes(hash);
        uint64_t h0 = h_arr[0];
        uint64_t h1 = h_arr[1];
        uint64_t h2 = h_arr[2];


        h0 = (h0*bins) * ftype;
        h1 = (h1*bins) * ftype;
        h2 = (h2*bins) * ftype;


        // concatenate (64/ftype) the fingerprint
        uint64_t fc64 = f;
        for (uint8_t b = 1; b < bins_per_batch ; ++b)
        {
            fc64 = (fc64 << ftype) | f;
        }

        for (size_t batch = 0; batch < ixf_ptr->bin_words; ++batch)
        {
            size_t batch_start = batch * 64;
            
            uint64_t v = fc64 ^ ixf_ptr->data.get_int(h0 + batch_start, 64) 
                              ^ ixf_ptr->data.get_int(h1 + batch_start, 64) 
                              ^ ixf_ptr->data.get_int(h2 + batch_start, 64);

            size_t used_bins = batch * bins_per_batch;
            uint8_t bits = 0;
            for (size_t bin = 0; bin < bins_per_batch; ++bin)
            {
                if (used_bins + bin == result_buffer.size())
                    break;

                uint64_t tmp = v << ((bins_per_batch - (bin+1)) * ftype ); 
                uint8_t tmpb = std::bitset<8>(tmp >> (64-ftype)).none() << bin;
                bits |= tmpb;
                   
            }
            result_buffer.data.set_int(used_bins, bits, ftype);
        }

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    //!\}
     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) && noexcept = delete;
};

//!\brief A bitvector representing the result of a call to `bulk_contains` of the seqan3::interleaved_bloom_filter.
template <typename FingerprintType>
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
class interleaved_3way_binary_fuse_filter<FingerprintType>::membership_agent::binning_bitvector
{
private:
    //!\brief The underlying datatype to use.
    using data_type = sdsl::bit_vector;
    //!\brief The bitvector.
    data_type data{};

    friend class membership_agent;

/*    template <std::integral value_t>
    friend class counting_agent_type;

    template <std::integral value_t>
    friend class counting_vector;
*/
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_bitvector() = default; //!< Defaulted.
    binning_bitvector(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector(binning_bitvector &&) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector &&) = default; //!< Defaulted.
    ~binning_bitvector() = default; //!< Defaulted.

    //!\brief Construct with given size.
    explicit binning_bitvector(size_t const size) :
        data(size)
    {}
    //!\}

    //!\brief Returns the number of elements.
    size_t size() const noexcept
    {
        return data.size();
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the first element of the container.
    auto begin() noexcept
    {
        return data.begin();
    }

    //!\copydoc begin()
    auto begin() const noexcept
    {
        return data.begin();
    }

    //!\brief Returns an iterator to the element following the last element of the container.
    auto end() noexcept
    {
        return data.end();
    }

    //!\copydoc end()
    auto end() const noexcept
    {
        return data.end();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Test for equality.
    friend bool operator==(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return lhs.data == rhs.data;
    }

    //!\brief Test for inequality.
    friend bool operator!=(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
     //!\brief Return the i-th element.
    auto operator[](size_t const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    //!\copydoc operator[]()
    auto operator[](size_t const i) const noexcept
    {
        assert(i < size());
        return data[i];
    }

    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}
};


/*!\brief Manages counting ranges of values for the seqan3::interleaved_3way_binary_fuse_filter.
 * \attention Calling seqan3::interleaved_3way_binary_fuse_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <typename FingerprintType>
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
template <std::integral value_t>
class interleaved_3way_binary_fuse_filter<FingerprintType>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::interleaved_3way_binary_fuse_filter.
    using ixf_t = interleaved_3way_binary_fuse_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_3way_binary_fuse_filter.
    ixf_t const * ixf_ptr{nullptr};

    //!\brief Store a seqan3::interleaved_3way_binary_fuse_filter::membership_agent to call `bulk_contains`.
    class ixf_t::membership_agent membership_agent;

public:
    class counting_vector;
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    explicit counting_agent_type(ixf_t const & ixf) :
        ixf_ptr(std::addressof(ixf)), membership_agent(ixf), result_buffer(ixf.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    counting_vector result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector const & bulk_count(value_range_t && values) & noexcept
    {
        assert(ixf_ptr != nullptr);
        assert(result_buffer.size() == ixf_ptr->bin_count());

        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, 0);

        for (auto && value : values)
            result_buffer += membership_agent.bulk_contains(value);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}

};

/*!\brief A data structure that behaves like a std::vector and can be used to consolidate the results of multiple calls
 *        to seqan3::interleaved_3way_binary_fuse_filter::membership_agent::bulk_contains.
 * \ingroup search_dream_index
 * \tparam value_t The type of the count. Must model std::integral.
 *
 * \details
 *
 * When using the seqan3::interleaved_3way_binary_fuse_filter::membership_agent::bulk_contains operation, a common use case is to
 * add up, for example, the results for all k-mers in a query. This yields, for each bin, the number of k-mers of a
 * query that are in the respective bin. Such information can be used to apply further filtering or abundance estimation
 * based on the k-mer counts.
 *
 * The seqan3::counting_vector offers an easy way to add up the individual
 * seqan3::interleaved_3way_binary_fuse_filter::membership_agent::binning_bitvector by offering an `+=` operator.
 *
 * The `value_t` template parameter should be chosen in a way that no overflow occurs if all calls to `bulk_contains`
 * return a hit for a specific bin. For example, `uint8_t` will suffice when processing short Illumina reads, whereas
 * long reads will require at least `uint32_t`.
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_vector.cpp
 */
template <typename FingerprintType>
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
template<std::integral value_t>
class interleaved_3way_binary_fuse_filter<FingerprintType>::counting_agent_type<value_t>::counting_vector : public std::vector<value_t>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t>;
    using ixf_t = interleaved_3way_binary_fuse_filter<FingerprintType>;

    friend class counting_agent_type<value_t>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_vector() = default; //!< Defaulted.
    counting_vector(counting_vector const &) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector const &) = default; //!< Defaulted.
    counting_vector(counting_vector &&) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector &&) = default; //!< Defaulted.
    ~counting_vector() = default; //!< Defaulted.

    using base_t::base_t;
    class ixf_t::membership_agent membership_agent;
    //!\}

    /*!\brief Bin-wise adds the bits of a seqan3::interleaved_3way_binary_fuse_filter::membership_agent::binning_bitvector.
     * \tparam rhs_t The type of the right-hand side.
     *         Must be seqan3::interleaved_3way_binary_fuse_filter::membership_agent::binning_bitvector.
     * \param rhs The seqan3::interleaved_3way_binary_fuse_filter::membership_agent::binning_bitvector.
     * \attention The counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    template <typename rhs_t>
/*    //!\cond
        requires std::same_as<rhs_t, membership_agent::binning_bitvector>
    //!\endcond
*/
    counting_vector & operator+=(rhs_t const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        // Each iteration can handle 64 bits, so we need to iterate `((rhs.size() + 63) >> 6` many times
        for (size_t batch = 0, bin = 0; batch < ((rhs.size() + 63) >> 6); bin = 64 * ++batch)
        {
            size_t tmp = rhs.raw_data().get_int(batch * 64); // get 64 bits starting at position `batch * 64`
            if (tmp ^ (1ULL<<63)) // This is a special case, because we would shift by 64 (UB) in the while loop.
            {
                while (tmp > 0)
                {
                    // Jump to the next 1 and increment the corresponding vector entry.
                    uint8_t step = std::countr_zero(tmp);
                    bin += step++;
                    tmp >>= step;
                    ++(*this)[bin++];
                }
            }
            else
            {
                ++(*this)[bin + 63];
            }
        }
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan3::counting_vector`s.
     * \param rhs The other seqan3::counting_vector.
     * \attention The seqan3::counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }
};

} // namespace seqan3
