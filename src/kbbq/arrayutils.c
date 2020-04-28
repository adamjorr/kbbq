#include "include/arrayutils.h"

static void increment_at_i(int a[], unsigned int i[], unsigned int n){
	for (unsigned int j = 0; j < n; j++){
		a[i[j]] += 1;
	}
}

// static void increment_at_b(int a[], bool i[], unsigned int n){
// 	for (unsigned int j = 0; j < n; j++){
// 		if(i[j]){
// 			a[j] += 1;
// 		}
// 	}
// }
