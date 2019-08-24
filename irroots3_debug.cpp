
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *    M O D E R N   I N T R O D U C T O R Y   A N A L Y S I S
 *
 *    Computer Programs with E D U C A T I O N A L   V E R S I O N S
 *
 *    irrational roots using: repeated synthetic division
 *    17 July 2019:  irrational roots Version 3.5
 *    Irrational Roots irroots3.cpp: Brute Force Version Using
 *                 Repeated Synthetic Division Count [global RSD]
 *                 Local variables PREC_X and PREC_Y in function
 *                 for setting precision.
 *
 *     There cannot be any less than 1 significant figure(s).
 *     Therefore, set fig = 1; precision(0) implies whole number
 *     n digits of precision to the right of radix is "figures - 1"
 *     That is, 8 figures implies need for 8 -1 = 7 digits of precision.
 *
  double stepsize(INTERVAL I)  {
    double diff = I.second - I.first;

    if ( diff > 1) return 1;        // changed from diff >= 1
    else return 0.1*diff;            // changed from return diff
  }

 *     std::setprecision(digits of accuracy after radix); // 0 |---> integer
 *
 *     nextX = std::round(std::pow(10, number of significant figures)*thisX)/ \
 *                          std::pow(10, fig);
 *
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
#include <map>
#include <iterator>


typedef std::pair<double, double>  INTERVAL;

int FIGURES = 7;         // 6 works nicely
                   // Setting to 1 or 2 gives 4 figures of accuracy
                   // Setting to 3 or 4 gives 8 figures of accuracy
                   // Setting to 5, 6, 7, 8 gives 16 figures
                   // Setting to 9..15 gives 32 figures
                   // Set higher for a cheap thrill.


int RSD = 0;   // Initialize RSD (Repeated Synthetic Division) to 0.
               // It is incremented AFTER each "repeated synthetic division"
               // and subsequent "Next Approximation of the Zero."
int t = 0;
double nextX = -6321;
double preX;
char choice = 'N';

std::string Ntab(int n);
double stepsize(INTERVAL I);

void printIntervals(std::vector<std::pair<double, double>> V);
std::pair<INTERVAL, double> getInterval();
std::vector<double> applyLocationPrinciple(int n, std::vector<double> P);
double repeated_Synthetic_Division(int n,
                                   std::vector<double> P,
                                   std::vector<double> R1,
                                   double highest, double X,
                                   double f);
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
  double h;  // subdivision size
  double highest, y, y1 = 0;
  std::vector<INTERVAL> intervals{};

  std::cout << "\nWhat is the degree of your polynomial? ";
  std::cin >> degree;

  //double P[degree+1];
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

// 3. Apply location principle (to build table displaying f(x) ---> f(0)):

    R = applyLocationPrinciple(degree, P);


// 4. Finish making table entries

/***********************************************************************/

x1 = query.first.first;
x2 = query.first.second;
h = query.second;

for (x = x1; x < x2; x += h) {

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
                      << " and " << x << ".\n";
            y1 = y;
            t += 1;

            INTERVAL I;
            I.first = x - h;
            I.second = x;
            intervals.push_back(I);
            std::cout << '\n';


// build second row showing f'(x)
std::vector<double> Q(R);
nextX = repeated_Synthetic_Division(degree, Q, R, highest, x, R.at(1));

std::cout << std::fixed;
std::cout << "\n\nCHOOSE X: ";  // RSD now 1

std::cout << std::fixed  << std::setprecision(std::pow(2, RSD))
          << nextX << "\t\n"; // 1.83

//Approximate X: ";
nextX = (std::round((double)std::pow(10, 1)*nextX)) / std::pow(10, 1);
preX = x;
x = nextX;   // IMPORTANT to update main() value of x

std::cout << "\nWith precision 1, x = nextX = " << std::fixed
          << std::setprecision(1) << x << " = " << nextX << '\n';

// display nextX

std::cout << "\nThe zero is approximately " << nextX << "\n\n";
std::cout << "\nIs this approximation satisfactory <Y/N> ? :  ";
std::cin >> choice;

 if ( (choice == 'y') || (choice == 'Y') )  {

     std::cout << "\n\nClose enough for government work!\n";
     std::cout << "\nThat approximation is good to " << RSD+1
               << " figures.\n\n";  // RSD = 1, so 2 figures
     return 0;
  }
  else
  {
    std::cout << "\nOur two consecutive integers have been found.\n";
    std::cout << "\nPrevious value of x = " << preX << "\tnextX = " << nextX;
    //x = nextX;
    std::cout << "\nWe've updated x = nextX = x - f(x)/f'(x).\n\n";
    std::cout << "\ncurrent value of x = " << x;

    do  {
        std::cout << "\nSince dx = " << h << " and x2 = " << x2
                  << " we're done\napplying the Location Principle, and now proceed"
                  << " to \napply Newton's method using repeated synthetic division.\n\n";
        std::cout << std::endl << std::endl;
        std::cout << "\t Press ENTER  [ <-----||| ]  to continue.\n";
        std::cin.clear();  // not necessary
        // The following is sufficient:
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
       } while (std::cin.get() != '\n');

       /*
       You place continue; at whatever point in the flow you
       want it to move on to the next iteration of the loop.
       */

    //continue;  // In this case, it is not what we want.
    break;
  }

}  // (end) SIGN CHANGE //******// y*y1 < 0 |---> sign change

}     // if y != 0
else {
      std::cout << "\n" << x << " is a zero.\n";
      t += 1;
      y1 = y;
   }

} // for (make table)
/*************************************************************************/
//std::cout << "\nOutside \"make table\"\n";

while ( (choice == 'n') || (choice == 'N') )  {
  std::vector<double> f_vec(degree+1);  // like R1
  // IMPORTANT: We must round these results correctly or else the error
  // will carry itself through the rest of the program
  y = P.at(degree);

  for (int j = degree; j >= 1; j--)  {
     y = y*x + P.at(j - 1);  // synthetic division
     f_vec.at(j) = y;
   }
   std::cout << '\n';

   std::vector<double> f_as_P(f_vec);

   double nX = repeated_Synthetic_Division(degree, f_as_P, f_vec, highest, x, f_vec.at(1));

   std::cout << std::fixed;  // RSD now 2
   std::cout << "\nCHOOSE X: " << std::setprecision((int)std::pow(2, RSD)) << nX << '\n';

   std::cout << "\nApproximate X: ";

  if (RSD == 2) {
    nextX = std::round(std::pow(10, (int)std::pow(2, RSD))*nX)/ \
            std::pow(10, (int)std::pow(2, RSD));

    std::cout << "\tWith precision (2^RSD) = (" << (int)std::pow(2, RSD)
              << "), nextX = " << std::fixed
              << std::setprecision((int)std::pow(2, RSD))
              << nextX << '\t';
  }
  else {
    nextX = std::round(std::pow(10, (int)std::pow(2, RSD))*nX)/ \
            std::pow(10, (int)std::pow(2, RSD));

    std::cout << "\tWith precision (2^RSD-1) = (" << (int)std::pow(2, RSD) - 1
              << "), nextX = " << std::fixed
              << std::setprecision((int)std::pow(2, RSD)-1)
              << nextX << '\t';
  }

   preX = x;
   x = nextX;

   std::cout << "\nThe zero is approximately " << x << "\n\n";
   std::cout << "\nIs this approximation satisfactory <Y/N> ? :  ";
   std::cin >> choice;
   if ( isCloseToEqual(preX, nextX) ) break;
   }


/*************************************************************************/

 if ( ( isCloseToEqual(x, nextX) ) && ( (choice == 'n') || (choice == 'N')) )
 {
      std::cout << "\n\nThis approximation is already good to "
                << std::fixed  << std::setprecision(std::pow(2, RSD))
                << (int)std::pow(2, RSD) << " figures.\n";
                 // RSD = 2, 3, 4 ... so 4, 8, 16 ... figures
  }
 if (t == 0)  {
       std::cout << "\nNo zeros found.\n\n";
     }
     else std::cout << "\n\n";

  // FINAL DISPLAY:
  std::cout << "\n\nThere were " << RSD << " repeated synthetic divisions.\n";
  printIntervals(intervals);
  std::cout << "\nOur efforts yield: x = ";
  switch (RSD) {
    case 0: std::cout << std::fixed
                      << std::setprecision(0)
                      << nextX << "\n\n";
            break;

    case 1: std::cout << std::fixed
                      << std::setprecision(1)
                      << nextX << "\n\n";
            break;

    case 2: std::cout << std::fixed
                      << std::setprecision(4)
                      << nextX << "\n\n";
            break;

    case 3: std::cout << std::fixed
                      << std::setprecision(7)
                      << nextX << "\n\n";
            break;

     case 4: std::cout << std::fixed
                       << std::setprecision(15)
                       << nextX << "\n\n";
             break;

     case 5:  std::cout << std::fixed
                       << std::setprecision(31)
                       << nextX << "\n\n";
              break;

    default:   std::cout << std::fixed
                      << std::setprecision(std::pow(2, RSD) - 1)
                      << nextX << "\n\n";
                break;
  }
}


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

std::vector<double> applyLocationPrinciple(int n, std::vector<double> P)   {

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

double repeated_Synthetic_Division(int n,
                std::vector<double> P, // y = y * x + P[j-1]
                std::vector<double> R1,  // stores row 1 f(x)
                double highest,        // coefficient of highest degree term
                double X,             // x
                double f)          // y = f(x)
{

 int PREC_X, PREC_Y, ROUND;  // where Y = f(x)

 switch (RSD) {
   case 0: PREC_X = PREC_Y = RSD;
           ROUND = 1;
           break;

   case 1: PREC_X = 1;
           PREC_Y = 3;
           ROUND = 3;
           break;

   case 2: PREC_X = 4;
           PREC_Y = 8;
           ROUND = 7;
           break;

   case 3: PREC_X = 7;
           PREC_Y = 15;
           ROUND = 15;
           break;

    case 4: PREC_X = 15;
            PREC_Y = 31;
            ROUND = 31;
            break;

    case 5: PREC_X = 31;
            PREC_Y = 63;
            ROUND = 63;

   default:  break;

 }

 std::vector<double> row2(n);
 double x = X;
 double y = highest;
 double tmp;
 std::cout << "\n\n x = " << x << "   |" << Ntab(n)
           << "     | f(" << x << ")\n";
  std::cout << "---------------------------------------------------\n\t\t";
  tmp = highest;
  tmp = (std::round(std::pow(10, ROUND)*tmp) / \
                    std::pow(10, ROUND));    // precision 0 for whole numbers

  std::cout << std::fixed  << std::setprecision(PREC_Y)  << tmp << '\t';

 for (int j = n; j >= 1; j--)  {
   tmp = R1.at(j);
   tmp = (std::round(std::pow(10, ROUND)*tmp) / \
                     std::pow(10, ROUND));

   std::cout << std::fixed  << std::setprecision(PREC_Y)  << tmp << '\t';
  }
 std::cout << "\n\n x = " << std::setprecision(PREC_X) << x << Ntab(n) << "     | f'("
           << x << ")\n";    // derivative of f
 std::cout << "---------------------------------------------------";

 double dy;
 if (isEqual(x, 0.0))
   std::cout << "\n\t0   |\t" << y << "\t";
 else std::cout << "\n\t   |\t" << std::fixed
                << std::setprecision(PREC_Y)  <<highest << "\t";

  for (int j = n+1; j > 1; j--)  {  // "n+1" down to 2: first elements just place-holders
                                    // we have highest degree coefficient already
                                    // so Q.at("n+1" - 1) = Q.at(n) = old P.at(n-1)
       y = y*x + P.at(j - 1);  // synthetic division
       row2.at(j-2) = y;       // IMPORTANT:  changed to j-2 from j-1
       }
      dy = row2.at(1);         // IMPORTANT: chnaged from (2) to (1): dy = row2[1]

      for (int i = n-1 ; i > 1; --i)  {  // IMPORTANT: changed from n to n-1

        tmp = row2.at(i);
        tmp = (std::round((double)std::pow(10, ROUND)*tmp) / std::pow(10, ROUND));
        std::cout << std::fixed  << std::setprecision(PREC_Y)  << tmp << '\t';
      }
      tmp = row2.at(1);
      tmp = (std::round((double)std::pow(10, ROUND)*tmp) / \
            std::pow(10, ROUND));
      std::cout << std::fixed  << std::setprecision(PREC_Y)  << tmp << '\t';
     std::cout << "\n\nrow2.at(1) = " << tmp << " = dy = " << dy << std::endl;
     std::cout << std::endl;
                                   // CHECK to see if f'(x) == 0
     if ( !(isCloseToEqual(dy, 0.0)) ) {
       std::cout << std::fixed  << std::setprecision(PREC_X)
                 << "\n\nNext X = \n  x - f(x)/f'(x) = " << x << " - [f("
                 << x << ")/f'(" << x << ")]\n    = " << x
                 << std::fixed  << std::setprecision(PREC_Y)
                 << " - (" << f << "/" << dy << ") = ";
     double X = x - f/dy;

     std::cout << std::fixed  << std::setprecision(PREC_Y+2)  << X << '\t';
     RSD++;
     return X;
     }
     else {
       std::cout << "\nf'(" << x << ") = " << dy << " == 0 ? ";
       do  {
           std::cout << "\nThis method requires that neither the first nor second\n"
                     << " derivatives be zero.\n";
           std::cout << std::endl << std::endl;
           std::cout << "\t Press ENTER  [ <-----||| ]  to continue.\n";
           std::cin.clear();  // not necessary
           // The following is sufficient:
           std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
          } while (std::cin.get() != '\n');

       return x;
     }
  }
