/************************************************************/
/*    NAME: Toby Schneider                                  */
/*    ORGN: MIT                                             */
/*    FILE: arithmetic.cpp                                  */
/*    DATE: 2012-03-12                                      */
/************************************************************/

#include "arithmetic.h"
#include "goby/util/as.h" // for goby::util::as<>()

#define ARITHMETIC_DEBUG

goby::acomms::Bitset ArithmeticCodec::encode_repeated(const vector<goby::int32>& wire_values)
{

#ifdef ARITHMETIC_DEBUG
    cout << "Wire values: ";
    for (int i = 0; i < wire_values.size(); i++){
        cout << wire_values[i] << " ";
    }
    cout << endl; 
#endif

    double upper_bound = 1, lower_bound = 0, prev_lower = 0, prev_upper = 0;
    double prob_fx_n_minus_one = 0;
    double prob_fx_n = symbol_probabilities_[wire_values[0]];

    //Use the sequence Sayood gives to calculate the upper and lower bounds
    for (int i = 0; i < wire_values.size(); i++){

        prob_fx_n = symbol_probabilities_[wire_values[i]];

        prev_lower = lower_bound; 
        prev_upper = upper_bound;

        lower_bound = prev_lower + (prev_upper - prev_lower)*prob_fx_n_minus_one; 
        upper_bound = prev_lower + (prev_upper - prev_lower)*prob_fx_n;

        prob_fx_n_minus_one = prob_fx_n; 
        
    }
#ifdef ARITHMETIC_DEBUG

    cout << lower_bound << endl;
    cout << upper_bound << endl;

#endif

    pair<double, double> encoded_range(lower_bound, upper_bound); 

    return range_to_bits(encoded_range);
    
}
    
std::vector<goby::int32> ArithmeticCodec::decode_repeated(goby::acomms::Bitset* bits)
{

    pair<double, double> range = bits_to_range(*bits); 
    double tag = (range.first + range.second) / 2;

#ifdef ARITHMETIC_DEBUG
    cout << "Decoded lower: " << range.first << endl; 
    cout << "Decoded upper: " << range.second << endl; 
    cout << "Tag: " << tag << endl;
#endif

    vector<goby::int32> decoded; 
    double upper_bound = 1, lower_bound = 0, prev_lower = 0, prev_upper = 0, diff = 0;
    double prob_fx_n, prob_fx_n_minus_one; 

    while (decoded.size() < max_repeat()){

        prev_lower = lower_bound; 
        prev_upper = upper_bound; 

        diff = prev_upper - prev_lower;

        prob_fx_n_minus_one = 0; 

        for(map<int, double>::const_iterator it = symbol_probabilities_.begin(),
                end = symbol_probabilities_.end();
            it != end;
            ++it)
        {

            prob_fx_n = it->second;  
            lower_bound = prev_lower + (prev_upper - prev_lower)*prob_fx_n_minus_one; 
            upper_bound = prev_lower + (prev_upper - prev_lower)*prob_fx_n;

            if (tag > lower_bound && tag < upper_bound){
#ifdef ARITHMETIC_DEBUG
                cout << "Lower bound: " << lower_bound << endl;
                cout << "Symbol: " << it->first << endl;
                cout << "Probability: " << it->second << endl;
                cout << "Upper bound: " << upper_bound << endl;
#endif
                decoded.push_back(it->first); 
                break; 
            }
            else
                prob_fx_n_minus_one = prob_fx_n; 


        }

    }

    return decoded; 

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
