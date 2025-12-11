// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef BENT_CHANNELLING_DEFAULT_CONFIG_H
#define BENT_CHANNELLING_DEFAULT_CONFIG_H



// --- DEFAULT METHOD -------------------------------------------------
// 2 → M2 
// 3 → M3
// 4 → M4 (recommended)
#define BENTCH_DEFAULT_METHOD   4

// --- DEFAULT ORDER --------------------------------------------------
// 12 → Yoshida 12th order (your standard)
#define BENTCH_DEFAULT_ORDER    4

// --- DEFAULT VARIANT ------------------------------------------------
// 1 → v1 (your main implementation)
#define BENTCH_DEFAULT_VARIANT  2



// --- DEFAULT BENDING RADIUS -----------------------------------------
// R = large number → effectively straight crystal
#define BENTCH_DEFAULT_R        10

// --- DEFAULT NUMBER OF STEPS ----------------------------------------
// 1 → no substepping
#define BENTCH_DEFAULT_NSTEPS   10


// --- DEFAULT PHYSICS PARAMETERS (Silicon (110)) --------------------
#define BENTCH_DEFAULT_U0       21.7681     // eV
#define BENTCH_DEFAULT_UMAX     23.9037     // eV
#define BENTCH_DEFAULT_ATF      1.94e-11    // m
#define BENTCH_DEFAULT_UT       7.5e-12     // m
#define BENTCH_DEFAULT_DP       1.92e-10    // m
#define BENTCH_DEFAULT_ALPHA 0.722452
#define BENTCH_DEFAULT_BETA  0.573481
// --- CHANNEL GEOMETRY ----------------------------------------------
#define BENTCH_DEFAULT_XC       (0.5 * BENTCH_DEFAULT_DP-0.6565*BENTCH_DEFAULT_UT)


// --- SMALL NONZERO BENDING PARAMETER -------------------------------
#define BENTCH_DEFAULT_ETA      (1.0 / BENTCH_DEFAULT_R)

#endif
