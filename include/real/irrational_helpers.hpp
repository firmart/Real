#ifndef BOOST_REAL_IRRATIONAL_HELPERS_HPP
#define BOOST_REAL_IRRATIONAL_HELPERS_HPP

#include <vector>
#include <algorithm>
#include <real/real.hpp>
#include <math.h>

namespace boost {
    namespace real {
        namespace irrational {
            // used for extra precision. should be replaced with something more definitive in the future.
            inline const int PLACEHOLDER = 10; 

            /**
             * @brief The function returns the n-th digit of the champernowne number in the
             * binary integer version.
             *
             * @param n - The number digit index.
             * @return The value of the champernowne number n-th digit (either 0 or 1)
             */
            int champernowne_binary_get_nth_digit(unsigned int n) {
                std::vector<int> binary = {1};
                unsigned int index = 0;

                while (index < n) {

                    if (std::all_of(binary.begin(), binary.end(), [](int d) -> bool { return d == 1; })) {

                        for (int i = (int)binary.size() - 1; i >= 0; i--) {
                            binary[i] = 0;
                        }
                        binary.insert(binary.begin(), 1);

                    } else {

                        for (int i = (int)binary.size() - 1; i >= 0; i--) {
                            if (binary[i] == 0) {
                                for (int j = (int)binary.size() - 1; j > i; j--) {
                                    binary[j] = 0;
                                }
                                binary[i] = 1;
                                break;
                            }
                        }
                    }

                    index += binary.size();
                }

                return binary[binary.size() - 1 - (index - n)];
            }

            /// @TODO: figure out how to avoid unnecessary recalculation by saving
            // some information in real_algorithm

            // e^x = sum_{k=0}^\infty x^k / k! = 1 + x + x^2/2! + x^3/3! + x^4/4!
            // calculates e^x, where x = m/l
            template<int m, int l = 1>
            int exponential_get_nth_digit(unsigned int n) {
                std::vector<int> const one {1};
                std::vector<int> const zero {0};

                exact_number one_n {one, 1};
                exact_number last_term {one, 1}; // x^k_n / k_n!
                exact_number current_value; // summation up to k_n

                exact_number k;
                exact_number k_fac;
                exact_number const mn (m);
                exact_number const ln (l);
                std::cout << mn.as_string() << ", " << ln.as_string() << '\n';
                exact_number x = mn;

                /// @TODO look into possible precision problems
                x.divide_vector(ln, n+PLACEHOLDER);

                // prepare to calculate k=2 term
                k = one_n + one_n;
                last_term = x;
                current_value = one_n + x;

                // keep getting terms from the taylor series until the terms go below our precision bound
                /// @TODO: ensure cast doesn't overflow
                while(last_term.exponent >= 1-(int)(n+PLACEHOLDER)) {
                    last_term *= x;
                    last_term.divide_vector(k, n+PLACEHOLDER);
                    current_value += last_term;

                    k = k + one_n;
                }

                return current_value.digits[n];
            }

            // TODO optimize: fast binary exponentiation
            static long mul_exp2_mod(long p, long exp, long m) {
                long result = p;
                for(long i = 0; i < exp; i++) {
                    result = (result * 2) % m;
                }
                return result % m;
            }

            /**
             * @brief The function returns the n-th digit of binary expansion of pi.
             * @param n - The number digit index.
             * @return The n-th digit of pi's binary expansion. 
             */

            int pi_binary_get_nth_digit(unsigned int n) {
                if (n == 1) return 1;
                else if (n == 2) return 1;
                n -= 3; // BBP algorithm computes the (n+1)-th *fractional bit*
                double result = 0;
                for(int k = 0; k <= (long)n/4; k++){
                    long factor = (n - 4*k);
                    result += (double)mul_exp2_mod(4, factor, 8*k + 1)/(8*k+1);
                    result -= (double)mul_exp2_mod(2, factor, 8*k + 4)/(8*k+4);
                    result -= (double)mul_exp2_mod(1, factor, 8*k + 5)/(8*k+5);
                    result -= (double)mul_exp2_mod(1, factor, 8*k + 6)/(8*k+6);
                }
                if (result < 0) result -= floor(result);
                return (int)(result * 2);
            }

            // e^x = sum_{k=0}^\inf = x^0/0! + x^1/1! + x^2/2! + x^3/3! + x^3/6! + ...

            //class exponential {
            //    private:
            //    // would be nice to interoperate between long, int, and boost::real::real,
            //    // and have ctors from the integral types
            //    boost::real::real k_prev = boost::real::real_explicit("0");
            //    boost::real::real const * const x_ptr;
            //    boost::real::real last_term; // x^kn / kn!
            //    boost::real::real current_value; // summation from k0 to k_n, with precision digits

            //    public:
            //    exponential(boost::real::real &x) : x_ptr(&x) {
            //        last_term = boost::real::real ("1");
            //        current_value = boost::real::real ("1");
            //    };

            //    int get_nth_digit(unsigned int n) {
            //        boost::real::real one = boost::real::real_explicit("1");
            //        // if n < k_prev, reset

            //        boost::real::real min_bound;
            //        std::get<boost::real::real_explicit>(min_bound.get_real_number()).digits = {1};
            //        std::get<boost::real::real_explicit>(min_bound.get_real_number()).exponent = 1-n;

            //        // keep getting terms from the taylor series until the terms go below our precision bound
            //        while((last_term > min_bound) || (last_term == min_bound)) {
            //            last_term = last_term * (*x_ptr) / (k_prev + one);
            //            current_value = current_value + last_term;
            //        }

            //        return 0;
            //    }

            //};
        }
    }
}

#endif //BOOST_REAL_IRRATIONAL_HELPERS_HPP
