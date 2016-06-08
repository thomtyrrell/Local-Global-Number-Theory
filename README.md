# Local-Global-Number-Theory
Code for studying the Brauer-Manin obstruction, Hasse principle, and other local/global problems in arithmetic geometry and number theory.  All files with a .sage or .sws extension are written in Python.

## monogenic_number_fields.sage
The method general_discriminant(K) takes a number field $K$ as input and returns an integer polynomial in degree(K) variables that represents the discriminant of the general element (divided by the discriminant of K) in the ring of integers of $K$.  If this polynomial vanishes modulo a prime, it implies that $K$ is not monogenic.  Ex.  K = cubic field defined by root of x^3 - x^2 - 2*x - 8.  

## mordell_weil_sieve.sage
This code can be used to study the Brauer-Manin Obstruction.  It is an implementation of a "Mordell-Weil Sieve" for some special curves of genus 2.  For more information, see <a href="http://jmilne.org/math/Students/b.pdf">"The Brauer-Manin Obstruction for Curves" by Victor Scharaschkin</a>.

## bremner.sage
This code implements a construction due to <a href="http://www.sciencedirect.com/science/article/pii/S0022314X97921892">Andrew Bremner in "Some Interesting Curves of genus 2 to 7"</a>.  Assumes MordellWeilSieve.sage is stored in the same directory.

## Modular Forms and Twists of Rank 0.sws
This sage worksheet breaks down one example in the theory of 3/2-weight modular forms and their relationship to ranks of elliptic curves.  For more information, see <a href=http://www.mathcs.emory.edu/~ono/publications-cv/pdfs/014.pdf>"Rank zero quadratic twists of modular elliptic curves" by Ken Ono</a>.  