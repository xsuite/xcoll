// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EMITTANCE_MONITOR_H
#define XCOLL_EMITTANCE_MONITOR_H

#ifndef C_LIGHT
#define C_LIGHT 299792458.0
#endif

/*gpufun*/
void EmittanceMonitor_track_local_particle(EmittanceMonitorData el, LocalParticle* part0){
    EmittanceMonitorData_set__cached(el, 0);
    int64_t const start_at_turn = EmittanceMonitorData_get_start_at_turn(el);
    int64_t const part_id_start = EmittanceMonitorData_get_part_id_start(el);
    int64_t const part_id_end   = EmittanceMonitorData_get_part_id_end(el);
    double const frev = EmittanceMonitorData_get_frev(el);
    double const sampling_frequency = EmittanceMonitorData_get_sampling_frequency(el);

    EmittanceMonitorRecord record = EmittanceMonitorData_getp_data(el);

    int64_t max_slot = EmittanceMonitorRecord_len_count(record);

    //start_per_particle_block(part0->part)
        int64_t particle_id = LocalParticle_get_particle_id(part);
        if (part_id_end < 0 || (part_id_start <= particle_id && particle_id < part_id_end)){

            // zeta is the absolute path length deviation from the reference particle: zeta = (s - beta0*c*t)
            // but without limits, i.e. it can exceed the circumference (for coasting beams)
            // as the particle falls behind or overtakes the reference particle
            double const zeta = LocalParticle_get_zeta(part);
            double const at_turn = LocalParticle_get_at_turn(part);
            double const beta0 = LocalParticle_get_beta0(part);

            // compute sample index
            int64_t slot = round(sampling_frequency * ( (at_turn-start_at_turn)/frev - zeta/beta0/C_LIGHT ));

            if (slot >= 0 && slot < max_slot){
                double const x  = LocalParticle_get_x(part);
                double const px = LocalParticle_get_px(part);
                double const y  = LocalParticle_get_y(part);
                double const py = LocalParticle_get_py(part);

                /*gpuglmem*/ double *count      = EmittanceMonitorRecord_getp1_count(record, slot);      atomicAdd(count, 1);
                /*gpuglmem*/ double *x_sum1     = EmittanceMonitorRecord_getp1_x_sum1(record, slot);     atomicAdd(x_sum1, x);
                /*gpuglmem*/ double *px_sum1    = EmittanceMonitorRecord_getp1_px_sum1(record, slot);    atomicAdd(px_sum1, px);
                /*gpuglmem*/ double *y_sum1     = EmittanceMonitorRecord_getp1_y_sum1(record, slot);     atomicAdd(y_sum1, y);
                /*gpuglmem*/ double *py_sum1    = EmittanceMonitorRecord_getp1_py_sum1(record, slot);    atomicAdd(py_sum1, py);
                /*gpuglmem*/ double *x_x_sum2   = EmittanceMonitorRecord_getp1_x_x_sum2(record, slot);   atomicAdd(x_x_sum2, x*x);
                /*gpuglmem*/ double *x_px_sum2  = EmittanceMonitorRecord_getp1_x_px_sum2(record, slot);  atomicAdd(x_px_sum2, x*px);
                /*gpuglmem*/ double *x_y_sum2   = EmittanceMonitorRecord_getp1_x_y_sum2(record, slot);   atomicAdd(x_y_sum2, x*y);
                /*gpuglmem*/ double *x_py_sum2  = EmittanceMonitorRecord_getp1_x_py_sum2(record, slot);  atomicAdd(x_py_sum2, x*py);
                /*gpuglmem*/ double *px_px_sum2 = EmittanceMonitorRecord_getp1_px_px_sum2(record, slot); atomicAdd(px_px_sum2, px*px);
                /*gpuglmem*/ double *px_y_sum2  = EmittanceMonitorRecord_getp1_px_y_sum2(record, slot);  atomicAdd(px_y_sum2, px*y);
                /*gpuglmem*/ double *px_py_sum2 = EmittanceMonitorRecord_getp1_px_py_sum2(record, slot); atomicAdd(px_py_sum2, px*py);
                /*gpuglmem*/ double *y_y_sum2   = EmittanceMonitorRecord_getp1_y_y_sum2(record, slot);   atomicAdd(y_y_sum2, y*y);
                /*gpuglmem*/ double *y_py_sum2  = EmittanceMonitorRecord_getp1_y_py_sum2(record, slot);  atomicAdd(y_py_sum2, y*py);
                /*gpuglmem*/ double *py_py_sum2 = EmittanceMonitorRecord_getp1_py_py_sum2(record, slot); atomicAdd(py_py_sum2, py*py);            }
        }
	//end_per_particle_block
}

#endif /* XCOLL_EMITTANCE_MONITOR_H */