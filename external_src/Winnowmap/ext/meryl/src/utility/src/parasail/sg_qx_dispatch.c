/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/cpuid.h"

/* forward declare the dispatcher functions */
parasail_function_t parasail_sg_qx_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_table_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_table_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_table_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_table_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_table_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_table_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_table_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_table_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_table_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_table_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_table_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_table_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_table_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_rowcol_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_stats_rowcol_diag_8_dispatcher;
parasail_function_t parasail_sg_qx_trace_scan_64_dispatcher;
parasail_function_t parasail_sg_qx_trace_scan_32_dispatcher;
parasail_function_t parasail_sg_qx_trace_scan_16_dispatcher;
parasail_function_t parasail_sg_qx_trace_scan_8_dispatcher;
parasail_function_t parasail_sg_qx_trace_striped_64_dispatcher;
parasail_function_t parasail_sg_qx_trace_striped_32_dispatcher;
parasail_function_t parasail_sg_qx_trace_striped_16_dispatcher;
parasail_function_t parasail_sg_qx_trace_striped_8_dispatcher;
parasail_function_t parasail_sg_qx_trace_diag_64_dispatcher;
parasail_function_t parasail_sg_qx_trace_diag_32_dispatcher;
parasail_function_t parasail_sg_qx_trace_diag_16_dispatcher;
parasail_function_t parasail_sg_qx_trace_diag_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_table_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_table_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_rowcol_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_stats_rowcol_striped_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_scan_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_scan_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_scan_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_scan_profile_8_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_striped_profile_64_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_striped_profile_32_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_striped_profile_16_dispatcher;
parasail_pfunction_t parasail_sg_qx_trace_striped_profile_8_dispatcher;

/* declare and initialize the pointer to the dispatcher function */
#if 0
parasail_function_t * parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_8_dispatcher;
parasail_function_t * parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_64_dispatcher;
parasail_function_t * parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_32_dispatcher;
parasail_function_t * parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_16_dispatcher;
parasail_function_t * parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_8_dispatcher;
parasail_function_t * parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_64_dispatcher;
#endif
parasail_function_t * parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_32_dispatcher;
parasail_function_t * parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_16_dispatcher;
#if 0
parasail_function_t * parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_8_dispatcher;
parasail_function_t * parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_64_dispatcher;
parasail_function_t * parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_32_dispatcher;
parasail_function_t * parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_16_dispatcher;
parasail_function_t * parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_8_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_64_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_32_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_16_dispatcher;
parasail_pfunction_t * parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_8_dispatcher;
#endif

/* dispatcher function implementations */

#if 0
parasail_result_t* parasail_sg_qx_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_stats_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_stats_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_stats_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_stats_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx_stats_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx_stats_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx_stats_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx_stats_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx_stats_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx_stats_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx_stats_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx_stats_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_table_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_table_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_table_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_table_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_table_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_table_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_table_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_table_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx_table_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx_table_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx_table_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx_table_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx_table_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_table_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx_table_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_table_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx_table_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_table_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx_table_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_table_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_table_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_stats_table_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_table_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_stats_table_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_table_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_stats_table_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_table_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_stats_table_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_table_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx_stats_table_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx_stats_table_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx_stats_table_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx_stats_table_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx_stats_table_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx_stats_table_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx_stats_table_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx_stats_table_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_table_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_rowcol_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_rowcol_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_rowcol_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_rowcol_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_rowcol_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_rowcol_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_rowcol_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_rowcol_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx_rowcol_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx_rowcol_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx_rowcol_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx_rowcol_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx_rowcol_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx_rowcol_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx_rowcol_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx_rowcol_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_rowcol_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_stats_rowcol_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_rowcol_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_stats_rowcol_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_rowcol_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_stats_rowcol_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_rowcol_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_stats_rowcol_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_stats_rowcol_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx_stats_rowcol_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx_stats_rowcol_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx_stats_rowcol_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx_stats_rowcol_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx_stats_rowcol_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx_stats_rowcol_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx_stats_rowcol_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx_stats_rowcol_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_stats_rowcol_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_trace_scan_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_64_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_trace_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_trace_scan_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_32_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_trace_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_trace_scan_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_16_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_trace_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_trace_scan_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_8_pointer = parasail_sg_qx_scan;
    }
    return parasail_sg_qx_trace_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_striped_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx_trace_striped_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
#endif

parasail_result_t* parasail_sg_qx_trace_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx_trace_striped_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx_trace_striped_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

#if 0
parasail_result_t* parasail_sg_qx_trace_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx_trace_striped_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx_trace_diag_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_trace_diag_64_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx_trace_diag_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_trace_diag_32_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx_trace_diag_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_trace_diag_16_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx_trace_diag_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_trace_diag_8_pointer = parasail_sg_qx;
    }
    return parasail_sg_qx_trace_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_profile_64_pointer = parasail_sg_qx_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_profile_32_pointer = parasail_sg_qx_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_profile_16_pointer = parasail_sg_qx_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_scan_profile_8_pointer = parasail_sg_qx_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_profile_64_pointer = parasail_sg_qx_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_profile_32_pointer = parasail_sg_qx_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_profile_16_pointer = parasail_sg_qx_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_striped_profile_8_pointer = parasail_sg_qx_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_profile_64_pointer = parasail_sg_qx_stats_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_profile_32_pointer = parasail_sg_qx_stats_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_profile_16_pointer = parasail_sg_qx_stats_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_scan_profile_8_pointer = parasail_sg_qx_stats_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_profile_64_pointer = parasail_sg_qx_stats_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_profile_32_pointer = parasail_sg_qx_stats_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_profile_16_pointer = parasail_sg_qx_stats_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_striped_profile_8_pointer = parasail_sg_qx_stats_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_profile_64_pointer = parasail_sg_qx_table_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_table_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_profile_32_pointer = parasail_sg_qx_table_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_table_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_profile_16_pointer = parasail_sg_qx_table_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_table_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_scan_profile_8_pointer = parasail_sg_qx_table_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_table_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_table_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_profile_64_pointer = parasail_sg_qx_table_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_table_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_profile_32_pointer = parasail_sg_qx_table_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_table_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_profile_16_pointer = parasail_sg_qx_table_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_table_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_table_striped_profile_8_pointer = parasail_sg_qx_table_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_table_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_table_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = parasail_sg_qx_stats_table_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = parasail_sg_qx_stats_table_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = parasail_sg_qx_stats_table_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = parasail_sg_qx_stats_table_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = parasail_sg_qx_stats_table_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = parasail_sg_qx_stats_table_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = parasail_sg_qx_stats_table_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = parasail_sg_qx_stats_table_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_table_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_table_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = parasail_sg_qx_rowcol_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = parasail_sg_qx_rowcol_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = parasail_sg_qx_rowcol_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = parasail_sg_qx_rowcol_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = parasail_sg_qx_rowcol_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = parasail_sg_qx_rowcol_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = parasail_sg_qx_rowcol_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = parasail_sg_qx_rowcol_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_rowcol_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_rowcol_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = parasail_sg_qx_stats_rowcol_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = parasail_sg_qx_stats_rowcol_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = parasail_sg_qx_stats_rowcol_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = parasail_sg_qx_stats_rowcol_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = parasail_sg_qx_stats_rowcol_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = parasail_sg_qx_stats_rowcol_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = parasail_sg_qx_stats_rowcol_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = parasail_sg_qx_stats_rowcol_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_stats_rowcol_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_stats_rowcol_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_profile_64_pointer = parasail_sg_qx_trace_scan_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_trace_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_profile_32_pointer = parasail_sg_qx_trace_scan_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_trace_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_profile_16_pointer = parasail_sg_qx_trace_scan_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_trace_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_scan_profile_8_pointer = parasail_sg_qx_trace_scan_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_trace_scan_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_trace_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_64_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_avx2_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_sse41_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_sse2_128_64;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_altivec_128_64;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_profile_64_pointer = parasail_sg_qx_trace_striped_profile_neon_128_64;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_profile_64_pointer = NULL;
    }
    return parasail_sg_qx_trace_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_32_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_sse2_128_32;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_altivec_128_32;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_profile_32_pointer = parasail_sg_qx_trace_striped_profile_neon_128_32;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_profile_32_pointer = NULL;
    }
    return parasail_sg_qx_trace_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_16_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_sse2_128_16;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_altivec_128_16;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_profile_16_pointer = parasail_sg_qx_trace_striped_profile_neon_128_16;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_profile_16_pointer = NULL;
    }
    return parasail_sg_qx_trace_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_8_dispatcher(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_sse2_128_8;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_altivec_128_8;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        parasail_sg_qx_trace_striped_profile_8_pointer = parasail_sg_qx_trace_striped_profile_neon_128_8;
    }
    else
#endif
    {
        parasail_sg_qx_trace_striped_profile_8_pointer = NULL;
    }
    return parasail_sg_qx_trace_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}
#endif

/* implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct impl */

#if 0
parasail_result_t* parasail_sg_qx_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_table_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_table_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_rowcol_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_stats_rowcol_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_scan_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_striped_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
#endif

parasail_result_t* parasail_sg_qx_trace_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

#if 0
parasail_result_t* parasail_sg_qx_trace_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_diag_64_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_trace_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return parasail_sg_qx_trace_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

parasail_result_t* parasail_sg_qx_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_table_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_table_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_rowcol_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_stats_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_stats_rowcol_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_scan_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_scan_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_scan_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_scan_profile_8_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_striped_profile_64_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_striped_profile_32_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_striped_profile_16_pointer(profile, s2, s2Len, open, gap);
}

parasail_result_t* parasail_sg_qx_trace_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    return parasail_sg_qx_trace_striped_profile_8_pointer(profile, s2, s2Len, open, gap);
}
#endif
