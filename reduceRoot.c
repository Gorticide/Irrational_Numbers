/************************************************************************
 *** Factorize the number under the root, pick out the factors that   ***
 *** come out in pairs and leave the rest under the root.             ***
 ***                                                                  ***
 *** √800 = √(2 x 2 x 2 x 2 x 5 x 2 x 5) = √(22 x 22 x 52 x 2)        ***
 ***      = (2 x 2 x 5)√2 = 20√2                                      ***
 ***                                                                  ***
 ************************************************************************/

#include<stdio.h>
/* #include<conio.h>  // Microsoft DOS */ 
#include<curses.h>    // Linux: compile with -lncurses  for getch()
#include <stdlib.h> 

int main(int argc, char *argv[]) {

    int h = 2;
    int inside_root;
    int outside_root = 1;

    /*
     printf("Enter integer under the radicand: ");
     scanf("%d", & n);
     */

    inside_root = atoi(argv[1]);;

    while ( (h * h) <= inside_root)  {
    if (inside_root % (h * h) == 0)   { // # inside_root evenly divisible by h^2
        inside_root /= (h * h);   // extract square factor
        outside_root *= h;        // term outside radicand multiplied by factor    
     }
    else   h += 1;
   }
    printf("\n\n( %d )*SQRT[ %d ]\n\n", outside_root, inside_root);
    getch();
    return 0;
}