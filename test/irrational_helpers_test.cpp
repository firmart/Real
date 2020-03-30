#include <catch2/catch.hpp>
#include <real/real.hpp>
#include <real/irrational_helpers.hpp>
#include <string>

TEST_CASE( "Compute irrational numbers", "[irrational]" ) {

    SECTION("2000 bits of pi computed with BBP algorithm") {
        std::string pi("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000000011011100000111001101000100101001000000100100111000001000100010100110011111001100011101000000001000001011101111101010011000111011000100111001101100100010010100010100101000001000011110011000111000110100000001001101110111101111100101010001100110110011110011010011101001000011000110110011000000101011000010100110110111110010010111110001010000110111010011111110000100110101011011010110110101010001110000100100010111100100100001011011010101110110011000100101111001111110110001101111010001001100010000101110100110100110001101111110110101101011000010111111111101011100101101101111010000000110101101111110110111101110001110000110101111111011010110101000100110011111101001011010111010011111001001000001000101111100010010110001111111100110010010010010100001100110010100011110110011100100010110110011110111000010000000000111110010111000101000010110001110111111000001011001100011011010010010000011011000011100010101011101001110011010011010010001011000111111101010001111110100100100110011110101111110000011011001010101110100100011110111001010001110101101100101100001110001100010111100110101011000100000100001010101001010111011100111101101010100101001000001110111000010010110100101100110110101100111000011000011010101001110010010101011110010011000000001001111000101110100011011000000100011001010000110000010000101111100001100101001000001011110010001100010111000110110110011100011101111100011100111100111011100101100000110000000111010000110000000111001101100100111100000111010001011101100000001111010001010001111101101011100010101011101111100000110111101001100010100101100100111011110001010111100101111110110100101010101100000010111000110000011100110010101010010010111110011101010100101010110101011100101000101011101001000100110000110001001100011111010000001010001000000010101011100101000111001011010100010101010101011000100001011011010110100110011000101110000110100000100010100000111101000110011101010000101010100");
        for(int i = 1; i <= pi.length(); i++){
            CHECK(boost::real::irrational::pi_binary_get_nth_digit(i) == (pi[i-1] - '0'));
        }
    }
}