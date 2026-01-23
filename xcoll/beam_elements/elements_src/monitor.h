// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_MONITOR_H
#define XCOLL_MONITOR_H

#ifndef C_LIGHT
#define C_LIGHT 299792458.0
#endif

#include <xtrack/headers/track.h>
#include <xtrack/headers/atomicadd.h>


/*gpufun*/
void ParticleStatsMonitor_track_local_particle(ParticleStatsMonitorData el, LocalParticle* part0){
    ParticleStatsMonitorData_set__cached(el, 0);
    int64_t const start_at_turn = ParticleStatsMonitorData_get_start_at_turn(el);
    int64_t const part_id_start = ParticleStatsMonitorData_get_part_id_start(el);
    int64_t const part_id_end   = ParticleStatsMonitorData_get_part_id_end(el);
    double const frev = ParticleStatsMonitorData_get_frev(el);
    double const sampling_frequency = ParticleStatsMonitorData_get_sampling_frequency(el);

    ParticleStatsMonitorRecord record = ParticleStatsMonitorData_getp_data(el);

    int64_t max_slot = ParticleStatsMonitorRecord_len_count(record);

    // Decode _selector
    int16_t const _selector         = ParticleStatsMonitorData_get__selector(el);
    int16_t const monitor_x         =  _selector % 2;
    int16_t const monitor_px        = (_selector >> 1) % 2;
    int16_t const monitor_y         = (_selector >> 2) % 2;
    int16_t const monitor_py        = (_selector >> 3) % 2;
    int16_t const monitor_zeta      = (_selector >> 4) % 2;
    int16_t const monitor_pzeta     = (_selector >> 5) % 2;
    int16_t const monitor_delta     = (_selector >> 6) % 2;
    int16_t const monitor_mean      = (_selector >> 7) % 2;
    int16_t const monitor_variance  = (_selector >> 8) % 2;

    //start_per_particle_block (part0->part);
        int64_t particle_id = LocalParticle_get_particle_id(part);
        if (part_id_end < 0 || (part_id_start <= particle_id && particle_id < part_id_end)){

            // zeta is the absolute path length deviation from the reference particle: zeta = (s - beta0*c*t)
            // but without limits, i.e. it can exceed the circumference (for coasting beams)
            // as the particle falls behind or overtakes the reference particle
            double const zeta = LocalParticle_get_zeta(part);
            double const beta0 = LocalParticle_get_beta0(part);
            double const at_turn = LocalParticle_get_at_turn(part);

            double x = 0;
            double px = 0;
            double y = 0;
            double py = 0;
            double pzeta = 0;
            double delta = 0;

            // Compute sample index
            int64_t slot = round(sampling_frequency * ( (at_turn-start_at_turn)/frev - zeta/beta0/C_LIGHT ));

            if (slot >= 0 && slot < max_slot){
                GPUGLMEM int64_t *count = ParticleStatsMonitorRecord_getp1_count(record, slot); atomicAdd(count, 1);

                // Read coordinates only if needed
                if (monitor_x){ x = LocalParticle_get_x(part); }
                if (monitor_px){ px = LocalParticle_get_px(part); }
                if (monitor_y){ y = LocalParticle_get_y(part); }
                if (monitor_py){ py = LocalParticle_get_py(part); }
                if (monitor_pzeta){
                    double const ptau = LocalParticle_get_ptau(part);
                    pzeta = ptau/beta0;
                }
                if (monitor_delta){ delta = LocalParticle_get_delta(part); }

                // Store sums of the coordinates only if needed
                if (monitor_mean){
                    if (monitor_x){ GPUGLMEM double *x_sum1 = ParticleStatsMonitorRecord_getp1_x_sum1(record, slot); atomicAdd(x_sum1, x); }
                    if (monitor_px){ GPUGLMEM double *px_sum1 = ParticleStatsMonitorRecord_getp1_px_sum1(record, slot); atomicAdd(px_sum1, px); }
                    if (monitor_y){ GPUGLMEM double *y_sum1 = ParticleStatsMonitorRecord_getp1_y_sum1(record, slot); atomicAdd(y_sum1, y); }
                    if (monitor_py){ GPUGLMEM double *py_sum1 = ParticleStatsMonitorRecord_getp1_py_sum1(record, slot); atomicAdd(py_sum1, py); }
                    if (monitor_zeta){ GPUGLMEM double *zeta_sum1 = ParticleStatsMonitorRecord_getp1_zeta_sum1(record, slot); atomicAdd(zeta_sum1, zeta); }
                    if (monitor_pzeta){ GPUGLMEM double *pzeta_sum1 = ParticleStatsMonitorRecord_getp1_pzeta_sum1(record, slot); atomicAdd(pzeta_sum1, pzeta); }
                    if (monitor_delta){ GPUGLMEM double *delta_sum1 = ParticleStatsMonitorRecord_getp1_delta_sum1(record, slot); atomicAdd(delta_sum1, delta); }
                }

                // Store squared sums of the coordinates only if needed
                if (monitor_variance){
                    if (monitor_x){
                        GPUGLMEM double *x_x_sum2 = ParticleStatsMonitorRecord_getp1_x_x_sum2(record, slot); atomicAdd(x_x_sum2, x*x);
                        if (monitor_px){ GPUGLMEM double *x_px_sum2 = ParticleStatsMonitorRecord_getp1_x_px_sum2(record, slot); atomicAdd(x_px_sum2, x*px);}
                        if (monitor_y){ GPUGLMEM double *x_y_sum2 = ParticleStatsMonitorRecord_getp1_x_y_sum2(record, slot); atomicAdd(x_y_sum2, x*y);}
                        if (monitor_py){ GPUGLMEM double *x_py_sum2 = ParticleStatsMonitorRecord_getp1_x_py_sum2(record, slot); atomicAdd(x_py_sum2, x*py);}
                        if (monitor_zeta){ GPUGLMEM double *x_zeta_sum2 = ParticleStatsMonitorRecord_getp1_x_zeta_sum2(record, slot); atomicAdd(x_zeta_sum2, x*zeta);}
                        if (monitor_pzeta){ GPUGLMEM double *x_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_x_pzeta_sum2(record, slot); atomicAdd(x_pzeta_sum2, x*pzeta);}
                        if (monitor_delta){ GPUGLMEM double *x_delta_sum2 = ParticleStatsMonitorRecord_getp1_x_delta_sum2(record, slot); atomicAdd(x_delta_sum2, x*delta);}
                    }
                    if (monitor_px){
                        GPUGLMEM double *px_px_sum2 = ParticleStatsMonitorRecord_getp1_px_px_sum2(record, slot); atomicAdd(px_px_sum2, px*px);
                        if (monitor_y){ GPUGLMEM double *px_y_sum2 = ParticleStatsMonitorRecord_getp1_px_y_sum2(record, slot); atomicAdd(px_y_sum2, px*y); }
                        if (monitor_py){ GPUGLMEM double *px_py_sum2 = ParticleStatsMonitorRecord_getp1_px_py_sum2(record, slot); atomicAdd(px_py_sum2, px*py); }
                        if (monitor_zeta){ GPUGLMEM double *px_zeta_sum2 = ParticleStatsMonitorRecord_getp1_px_zeta_sum2(record, slot); atomicAdd(px_zeta_sum2, px*zeta); }
                        if (monitor_pzeta){ GPUGLMEM double *px_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_px_pzeta_sum2(record, slot); atomicAdd(px_pzeta_sum2, px*pzeta); }
                        if (monitor_delta){ GPUGLMEM double *px_delta_sum2 = ParticleStatsMonitorRecord_getp1_px_delta_sum2(record, slot); atomicAdd(px_delta_sum2, px*delta); }
                    }
                    if (monitor_y){
                        GPUGLMEM double *y_y_sum2 = ParticleStatsMonitorRecord_getp1_y_y_sum2(record, slot); atomicAdd(y_y_sum2, y*y);
                        if (monitor_py){ GPUGLMEM double *y_py_sum2 = ParticleStatsMonitorRecord_getp1_y_py_sum2(record, slot); atomicAdd(y_py_sum2, y*py); }
                        if (monitor_zeta){ GPUGLMEM double *y_zeta_sum2 = ParticleStatsMonitorRecord_getp1_y_zeta_sum2(record, slot); atomicAdd(y_zeta_sum2, y*zeta); }
                        if (monitor_pzeta){ GPUGLMEM double *y_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_y_pzeta_sum2(record, slot); atomicAdd(y_pzeta_sum2, y*pzeta); }
                        if (monitor_delta){ GPUGLMEM double *y_delta_sum2 = ParticleStatsMonitorRecord_getp1_y_delta_sum2(record, slot); atomicAdd(y_delta_sum2, y*delta); }
                    }
                    if (monitor_py){
                        GPUGLMEM double *py_py_sum2 = ParticleStatsMonitorRecord_getp1_py_py_sum2(record, slot); atomicAdd(py_py_sum2, py*py);
                        if (monitor_zeta){ GPUGLMEM double *py_zeta_sum2 = ParticleStatsMonitorRecord_getp1_py_zeta_sum2(record, slot); atomicAdd(py_zeta_sum2, py*zeta); }
                        if (monitor_pzeta){ GPUGLMEM double *py_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_py_pzeta_sum2(record, slot); atomicAdd(py_pzeta_sum2, py*pzeta); }
                        if (monitor_delta){ GPUGLMEM double *py_delta_sum2 = ParticleStatsMonitorRecord_getp1_py_delta_sum2(record, slot); atomicAdd(py_delta_sum2, py*delta); }
                    }
                    if (monitor_zeta){
                        GPUGLMEM double *zeta_zeta_sum2 = ParticleStatsMonitorRecord_getp1_zeta_zeta_sum2(record, slot); atomicAdd(zeta_zeta_sum2, zeta*zeta);
                        if (monitor_pzeta){ GPUGLMEM double *zeta_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_zeta_pzeta_sum2(record, slot); atomicAdd(zeta_pzeta_sum2, zeta*pzeta); }
                        if (monitor_delta){ GPUGLMEM double *zeta_delta_sum2 = ParticleStatsMonitorRecord_getp1_zeta_delta_sum2(record, slot); atomicAdd(zeta_delta_sum2, zeta*delta); }
                    }
                    if (monitor_pzeta){
                        GPUGLMEM double *pzeta_pzeta_sum2 = ParticleStatsMonitorRecord_getp1_pzeta_pzeta_sum2(record, slot); atomicAdd(pzeta_pzeta_sum2, pzeta*pzeta);
                        if (monitor_delta){ GPUGLMEM double *pzeta_delta_sum2 = ParticleStatsMonitorRecord_getp1_pzeta_delta_sum2(record, slot); atomicAdd(pzeta_delta_sum2, pzeta*delta); }
                    }
                    if (monitor_delta){
                        GPUGLMEM double *delta_delta_sum2 = ParticleStatsMonitorRecord_getp1_delta_delta_sum2(record, slot); atomicAdd(delta_delta_sum2, delta*delta);
                    }
                }
            }
        }
	//end_per_particle_block
}

#endif /* XCOLL_MONITOR_H */