
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *    M O D E R N   I N T R O D U C T O R Y   A N A L Y S I S
 *
 *    Computer Programs with E D U C A T I O N A L   V E R S I O N S
 *
 *    irrational roots:

 *    19 July 2019:  irrational roots Version "1"
 *    Irrational Roots 1:  irroots1_tenths.cpp
 *                   Using location principle via synthetic division/substitution
 *                   to "bracket" the root between succussive "tenths"
 *
 *     We will use global variable FIGURES as in irroots3.cpp (irroots),
 *     likewise, setting with FIGURES = setSigFig(getSigFig);
 *
  double stepsize(INTERVAL I)  {
    double diff = I.second - I.first;

    if ( diff > 1) return 1;        // changed from diff >= 1
    else return 0.1*diff;            // changed from return diff
  }

 *     std::setprecision(LP - 1); // 0 |---> integer
 *
 *     std::round(std::pow(10, LP - 1)*tmp)/std::pow(10, LP);
 *
 *** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <cmath>    // for std::abs(double)
#include <string>  // for dynamic tab
#include <utility>  // for std::pair
#include <vector>   // for vector of pairs to display Intervals
#include <iomanip>
#include <algorithm>    // std::copy
#include <cmath>
#include <istream>  // to use cin.get()

typedef std::pair<double, double>  INTERVAL;

int FIGURES = 7;

int LP = 1;   // Initialize LP to 1
               // It is incremented AFTER each succussive intervals
               // is added to std::vector<INTERVAL> intervals
               // and subsequent calls to setInterval.
               // Hence, LP = intervals.size() or just increment LP in
               // pair<INTERVAL, double> info = setIntervals(x2, h, intervals);
               // which corresponds to setIntervals(X, dx, V) {}
               // which returns pair<INTERVAL, double> data.

int t = 0;
char choice = 'N';

std::string Ntab(int n);
double stepsize(INTERVAL I);

void printIntervals(std::vector<std::pair<double, double>> V);
std::pair<INTERVAL, double> getInterval();

std::vector<double> SyntheticSubstitution_0(int n, std::vector<double> P);

// For displaying each row of coefficients, R = SyntheticSubstitution_0, for each
// interval set by setInterval(x2, h, intervals) inside DO-WHILE Loop
double SyntheticSubstitution_X(double X, std::vector<double> P_vec, std::vector<double> f_vec, int n);

std::pair<INTERVAL, double> setInterval( double X,
                                         double dx,
                                         std::vector<INTERVAL> V);

 double bracket_Root_Between_Tenths(
                       std::pair<INTERVAL, double> data,
                       std::vector<double> P_vec,
                       std::vector<double> f_vec, int n);

// working with FIGURES and isCloseToEqual() in DO-WHILE Loop
int setSigFig(int F);
int getSigFig();

/* Dubiously wrong way

inline bool isEqual(double x, double y)     {
       return x == y; }
 */

inline bool isEqual(double x, double y)
 {
   const double epsilon = 1e-5;  /* 1e-5 =  0.00001 */;
   return std::abs(x - y) <= epsilon * std::abs(x);
   // see Knuth section 4.2.2 pages 217-218
 }

inline bool isCloseToEqual(double x, double y)     {

        if ( ( std::abs(x - y) <= std::pow(10, -FIGURES) * std::abs(x) )
          &&
             ( std::abs(y - x) <= std::pow(10, -FIGURES) * std::abs(y))
          )
          // these might as well be equal
          return true;
          else return false;
  }


int main()  {
  int degree;
  double coefficient;
  double x1, x2, x; // endpoints of interval to be searched, plus main x
                   // initialize both x and global nextX
                  // to set the inner loop in motion
  double h;      // subdivision size
  double highest, y, y1 = 0;
  std::vector<INTERVAL> intervals{};
  std::pair<INTERVAL, double> info;

  std::cout << "\nWhat is the degree of your polynomial? ";
  std::cin >> degree;

  std::vector<double> P(degree+1);
  std::vector<double> R(degree+1);
  // create row2 on the fly ...


// 1. Build P = f(x)
  for (int j = degree; j >= 1; j--) {
    std::cout << "coefficient of x^" << j << ": ";
    std::cin >> coefficient;
    P.at(j) = coefficient;
  }
  std::cout << "\nconstant term: ";
  std::cin >> coefficient;
  P.at(0) = coefficient;


// 2. Get Interval (a, b), a < b and compute size of subinterval dx = h
//    Using [(a,b), dx] in query with application of Location principle
//    using endpoints a to b subintervals dx.
  std::pair<INTERVAL, double> query;
  query = getInterval();

  // 3. "query" least number of significant figures required."
  FIGURES = getSigFig();

  // 4. Apply location principle (to build table displaying f(x) ---> f(0))
  // The first time is to locate tow consecutibe integers where f(x1)*f(x0) < 0

      R = SyntheticSubstitution_0(degree, P);

// 5. Finish making table entries

/***********************************************************************/


x1 = query.first.first;
x2 = query.first.second;
h = query.second;


for (x = x1; x < x2+h/2; x += h) {            // Make Initial Table

   highest = y = P.at(degree);
   R.at(degree) = y;
   std::cout << "\t" << x << "\t" << highest << "\t";
   for (int j = degree; j >= 1; j--)  {
       y = y*x + P.at(j - 1);  // synthetic division
       R.at(j) = y;
       std::cout << R.at(j) << '\t';
    }
    std::cout << '\n';


     if ( !isCloseToEqual(y, 0.0) )  {   // if y != 0

      if ( ( isCloseToEqual(y*y1, 0.0) )  || (y*y1 > 0.0)  )  {
                y1 = y;
      }
      else if (y*y1 < 0)  {                                   // SIGN CHANGE !


     std::cout << "\nThere is a zero between " << x - h
               << " and " << x << ".\n\n";

     y1 = y;
     t += 1;

      while (LP <= FIGURES)
      {
       // STEP 6: set interval:
       info = setInterval(x, h, intervals);

       // Now we will be able to access x1 = info.first.first, x2 = info.first.second,
       // and info.second (), which is dx and becomes the next h.
       x1 = info.first.first;
       x2 = info.first.second;
       h = info.second;

       x = bracket_Root_Between_Tenths(info, P, R, degree);
       x = std::round(std::pow(10, LP - 1)*x)/std::pow(10, LP - 1);

     }
     break;  // Break out of initial table and bgin DO_WHILE with fresh R.

   }  // (end) SIGN CHANGE //******// y*y1 < 0 |---> sign change

 }     // if y != 0
 else {
       std::cout << '\n' << std::fixed  << std::setprecision(LP-1)
                 << x << " is a zero.\n";
       t += 1;
       y1 = y;
    }

 } // for (make table)

 if (t == 0)  {
       std::cout << "\nNo zeros found.\n\n";
     }
     else std::cout << "\n\n";

  // FINAL DISPLAY:
  std::cout << "\nThere were " << LP << " applications of the Location "
            << " Principle via synthetic divisions.\n";
  printIntervals(intervals);
  std::cout << "\nOur efforts yield: x = ";
 std::cout << std::fixed << std::setprecision(LP - 2)
          << x << "\n\n";

}  // end main()

std::string Ntab(int n)  {
  std::string s = "";
  for (int i = 0; i < n; ++i)   s += '\t';
  return s;
}


void printIntervals(std::vector<std::pair<double, double>> V)   {
  std::cout << "\n\tZ E R O S:\n";
    for (auto& v : V)
      std::cout << "\nThere is a zero between " << std::fixed
                << std::setprecision(0) << v.first
                << " and " << v.second << ".\n";

}

double stepsize(INTERVAL I)  {
  double diff = I.second - I.first;

  if ( diff > 1) return 1;        // changed from diff >= 1
  else return 0.1*diff;            // changed from return diff
}


std::pair<INTERVAL, double> getInterval()  {
  INTERVAL i;
  std::cout << "\nendpoints of interval to be searched: ";
  std::cin >> i.first >> i.second;
  //std::cout << "\nsubdivision size: ";
  double step = stepsize(i);
  //std::cin >> step;
  std::pair<INTERVAL, double> results;
  results.first = i;
  results.second = step;
  return results;
}

std::vector<double>
SyntheticSubstitution_0(int n, std::vector<double> P)
{
 std::vector<double> R_local(n+1);
 std::cout << "\n\n\tx   |" << Ntab(n) << "     | f(x)\n";
 std::cout << "---------------------------------------------------";
 double y = P.at(n);
 std::cout << "\n\t0   |" << "\t" << y << "\t";
 double x = 0.0;
     for (int j = n; j >= 1; j--)  {
       y = y*x + P.at(j - 1);  // synthetic division
       R_local.at(j) = y;
       std::cout << R_local.at(j) << '\t';
     }
     std::cout << '\n';
     return R_local;
 }

double SyntheticSubstitution_X(double X, std::vector<double> P, std::vector<double> f_vec, int n)
  {
   int ROUND = LP - 1;
   double y = (std::round(std::pow(10, ROUND)*P.at(n)) / std::pow(10, ROUND));
   double x = (std::round(std::pow(10, ROUND)*X) / std::pow(10, ROUND));

   std::cout << "\n\t" << std::fixed  << std::setprecision(ROUND)
             << x << "   |   " << y << "\t"; // still print old highest first

       for (int j = n; j >= 1; j--)  {
         y = y*x + P.at(j - 1);  // synthetic division
         f_vec.at(j) = y;
         double tmp = f_vec.at(j);
         tmp = std::round(std::pow(10, ROUND)*tmp)/std::pow(10, ROUND);
         std::cout << std::fixed  << std::setprecision(ROUND)  << tmp << '\t';
       }
       std::cout << '\n';
       return y;
   }

std::pair<INTERVAL, double> setInterval(double X,
                                         double dx,
                                         std::vector<INTERVAL> V)  {
    std::pair<INTERVAL, double> data;
    INTERVAL I;
    I.first = X - dx;
    I.second = X;
    V.push_back(I);
    data.first = I;
    data.second = stepsize(I);
    // UPDATE LP for version using "tenths": irroots1 := irroots1_tenths.cpp
    LP++;
    return data;
 }

 double bracket_Root_Between_Tenths(
                       std::pair<INTERVAL, double> data,
                       std::vector<double> P_vec,
                       std::vector<double> f_vec, int n)
 {
   double root;
   double x1 = data.first.first;
   double x2 = data.first.second;
   double y, y1 = 0;
   double h = data.second;

   for (double x = x1; x < x2+h/2; x += h) {    // Make a Table

         y = SyntheticSubstitution_X(x, P_vec, f_vec, n);

         if ( !isCloseToEqual(y, 0.0) )  {   // if y != 0

          if ( ( isCloseToEqual(y*y1, 0.0) )  || (y*y1 > 0.0)  )  {
                    y1 = y;
          }
          else if (y*y1 < 0)  {          // SIGN CHANGE !

         std::cout << "\nThere is a zero between " << std::fixed
                   << std::setprecision(LP - 1) << x - h
                   << " and " << x << ".\n";
         y1 = y;
         root = x;
         break;
       }
     }     // END if y != 0
     else {
           std::cout << "\n" << x << " is a zero.\n";
           y1 = y;
           root = x;
        }
    }
   return root;
 }

int getSigFig()  {
  int S;
     std::cout << "\nEnter the least number of significant figures required: ";
     std::cin >> S;
     return S;
  }

// unnecessary?
int setSigFig(int F)  {
   int G = 32;
   if (F <= 16) G = 16;
   if (F <= 8)  G = 8;
   if (F <= 4)  G = 4;
   return G;
 }
