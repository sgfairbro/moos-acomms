/************************************************************/
/*    NAME: Toby Schneider                                  */
/*    ORGN: MIT                                             */
/*    FILE: arithmetic.cpp                                  */
/*    DATE: 2012-03-12                                      */
/************************************************************/

#include "arithmetic.h"
#include "goby/util/as.h" // for goby::util::as<>()

#define ARITHMETIC_DEBUG

goby::acomms::Bitset ArithmeticCodec::encode_repeated(const std::vector<goby::int32>& wire_values)
{
    std::pair<double, double> encoded_range(0, 1);
    for(int i = 0, n = max_repeat(); i < n; ++i)
    {
        // set out of range or nonexistent symbols to 20. Inelegant, but we have to do something
        // (real system should use NaN or other symbol to indicate these failures)
        int symbol = (i < wire_values.size()) ? wire_values[i] : 20;
        if(!symbol_probabilities_.count(symbol))
            symbol = 20;

#ifdef ARITHMETIC_DEBUG
        std::cout << "Symbol is: " << symbol << std::endl;
        std::cout << "Old range is: [" << encoded_range.first << "," << encoded_range.second << ")" << std::endl;
#endif
        
        // find symbol range
        double spread = (encoded_range.second - encoded_range.first);

        // narrow range based on symbol given
        std::map<int, double>::const_iterator symbol_it = symbol_probabilities_.find(symbol);
        for(std::map<int, double>::const_iterator it = symbol_probabilities_.begin();
            it != symbol_it;
            ++it)
        {
            encoded_range.first += spread*it->second;
        }
        encoded_range.second = encoded_range.first + spread*symbol_it->second;
        
#ifdef ARITHMETIC_DEBUG
        std::cout << "New range is: [" << std::setprecision(15) << encoded_range.first << "," << encoded_range.second << ")" << std::endl;
#endif
    }    

    // convert range to bitset
    goby::acomms::Bitset data_bits = range_to_bits(encoded_range);

#ifdef ARITHMETIC_DEBUG
    std::cout << "Data bits: " << data_bits << std::endl;
    double cumulative_probability = 1;
    for(int i = 0, n = wire_values.size(); i < n; ++i)
    {
        goby::int32 wire_value = wire_values[i];
        if(!symbol_probabilities_.count(wire_value))
            throw(goby::acomms::DCCLException("Symbol " + goby::util::as<std::string>(wire_value) + " not in probability map"));            
            
        double probability = symbol_probabilities_[wire_value];
        cumulative_probability *= probability;
    }
    
    std::cout << "Actual size: " << data_bits.size() << ", entropy (lower bound) size: " << goby::util::ceil_log2(1/cumulative_probability) << std::endl;
#endif
    
    return data_bits;
}
    
std::vector<goby::int32> ArithmeticCodec::decode_repeated(goby::acomms::Bitset* bits)
{
    // `bits` contains min_size_repeated number of bits (here, 3)
    
    std::vector<goby::int32> out;
    int last_symbol = -1; // some symbol that isn't in our map

    std::pair<double, double> decoded_range(0, 1);

    // loop until we have (in this case) 5 symbols
    while(out.size() < max_repeat())
    {
        double spread = (decoded_range.second - decoded_range.first);

        // see what our minimum bitset gives us for a range. If we're lucky we're done
        std::pair<double, double> current_range = bits_to_range(*bits);            

        // place a marker at the lower bound of the known decoded range
        std::pair<double, double> search_range(decoded_range);

        // loop over all the symbols moving our search range
        // until the computed range `current_range` lies *completely* within the search range

        //  lX     search_range.first
        //  uX     search_range.second
        //  l*     current_range.first
        //  u*     current_range.second
        //  ||     decoded_range
        //
        // e.g.
        //  |       l*  lX      u*  uX                |        
        for(std::map<int, double>::const_iterator it = symbol_probabilities_.begin(),
                end = symbol_probabilities_.end();
            it != end;
            ++it)
        {
            double search_spread = spread*it->second;
            search_range.second = search_range.first + search_spread;
            
            // is the lower bound of the current range inside our symbol's range?
            if(current_range.first >= search_range.first &&
               current_range.first < search_range.second)
            {
                // |    lX   l*    uX                  |
                // is the upper bound *also* inside our symbol's range?
                if (current_range.second <= search_range.second)
                {
                    // |    lX   l*  u*  uX                 |
                    // we've found a symbol!
                    last_symbol = it->first;
                    out.push_back(last_symbol);
                    decoded_range = search_range;
#ifdef ARITHMETIC_DEBUG
                    std::cout << "found symbol: " << last_symbol << std::endl;
                    std::cout << "range: [" << decoded_range.first << "," << decoded_range.second << ")" << std::endl;
#endif
                    break;
                }
                else
                {
                    // |    lX   l*   uX   u*               |
                    // no go, not enough bits
                    bits->get_more_bits(1);
                    break;
                }
            }
            else
            {
                // |    lX      uX     l*   u*         |
                // no overlap at all yet, keep moving the marker
                search_range.first += search_spread;
            }
        }
    }

    return out;
}
    
unsigned ArithmeticCodec::size_repeated(const std::vector<goby::int32>& field_values)
{
    return encode_repeated(field_values).size();
}        

std::pair<double, double> ArithmeticCodec::bits_to_range(goby::acomms::Bitset bits)
{
    // lower bound is the given bitset followed by 000000...
    double lower = 0;
    for(int i = 0, n = bits.size(); i < n; ++i)
        lower += bits[i]/(pow(2.0, i+1));

    // upper bound is the given bitset followed by 111111...
    // remember 0.111111... = 1 in binary, like 0.999999... = 1 in decimal
    double upper = lower + 1/(pow(2.0, bits.size()));
        
    return std::make_pair(lower,upper);
}

goby::acomms::Bitset ArithmeticCodec::range_to_bits(std::pair<double, double> range)
{
    double lower = range.first, upper = range.second;
    goby::acomms::Bitset lower_bits, upper_bits;

    // keep adding bits until the binary representation of lower and upper
    // differ OR we have found an exact representation for both decimals
    while(lower_bits == upper_bits && !(lower == 0 && upper == 0))
    {
        double lower_int_part, upper_int_part;
        lower *= 2;
        lower = modf(lower, &lower_int_part);
        lower_bits.push_back(lower_int_part >= 1 ? 1 : 0);
            
        upper *= 2;
        upper = modf(upper, &upper_int_part);

        if(upper_int_part == 2) // deal with 1 case
            upper = 1;
        
        upper_bits.push_back(upper_int_part >= 1 ? 1 : 0);
    }

    // if lower is an exact representation return it
    if(lower == 0)
    {
        return lower_bits;
    }
    else // keep adding bits until we uniquely identify the range
    {
        while(bits_to_range(upper_bits).second > range.second)
            upper_bits.push_back(0);

        while(bits_to_range(lower_bits).first < range.first)
            lower_bits.push_back(1);

        // return the smaller of the two representations of the range
        return (lower_bits.size() < upper_bits.size()) ? lower_bits : upper_bits;
    }
        
}

void ArithmeticCodec::validate()
{
    DCCLFieldCodecBase::require(DCCLFieldCodecBase::dccl_field_options().max_repeat() <= LARGEST_MAX_REPEAT_SIZE,
                                "(goby.field).dccl.max_repeat must be less than or equal to " +
                                goby::util::as<std::string>(static_cast<int>(LARGEST_MAX_REPEAT_SIZE)));
}
