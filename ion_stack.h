#ifndef _MCERD_ION_STACK_H_
#define _MCERD_ION_STACK_H_

#include "general.h"

Ion *next_ion(Global *global, Ion *cur_ion, Ion *ions_moving);
Ion *prev_ion(Global *global, Ion *cur_ion, Ion *ions_moving);
void copy_ions(Ion *ion, Target *target, int dest, int src, int copy_stopping);
void cascades_create_additional_ions(Global *global, Detector *detector, Target *target, Ion **ion);

#endif /* _MCERD_ION_STACK_H_ */
