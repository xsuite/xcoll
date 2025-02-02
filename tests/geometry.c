#include <stdint.h>
#ifndef XOBJ_TYPEDEF_ArrNFloat64
#define XOBJ_TYPEDEF_ArrNFloat64
typedef   struct ArrNFloat64_s * ArrNFloat64;
 static inline ArrNFloat64 ArrNFloat64_getp(ArrNFloat64 restrict  obj){
  int64_t offset=0;
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ArrNFloat64_len(ArrNFloat64 restrict  obj){
  int64_t offset=0;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ArrNFloat64_get(const ArrNFloat64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ArrNFloat64_set(ArrNFloat64 restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ArrNFloat64_getp1(ArrNFloat64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_ArrNInt64
#define XOBJ_TYPEDEF_ArrNInt64
typedef   struct ArrNInt64_s * ArrNInt64;
 static inline ArrNInt64 ArrNInt64_getp(ArrNInt64 restrict  obj){
  int64_t offset=0;
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ArrNInt64_len(ArrNInt64 restrict  obj){
  int64_t offset=0;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ArrNInt64_get(const ArrNInt64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ArrNInt64_set(ArrNInt64 restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ArrNInt64_getp1(ArrNInt64 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_ArrNUint32
#define XOBJ_TYPEDEF_ArrNUint32
typedef   struct ArrNUint32_s * ArrNUint32;
 static inline ArrNUint32 ArrNUint32_getp(ArrNUint32 restrict  obj){
  int64_t offset=0;
  return (ArrNUint32)(( char*) obj+offset);
}
 static inline int64_t ArrNUint32_len(ArrNUint32 restrict  obj){
  int64_t offset=0;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline uint32_t ArrNUint32_get(const ArrNUint32 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*4;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void ArrNUint32_set(ArrNUint32 restrict  obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=16+i0*4;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* ArrNUint32_getp1(ArrNUint32 restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=16+i0*4;
  return ( uint32_t*)(( char*) obj+offset);
}
#endif
#ifndef XOBJ_TYPEDEF_RecordIndex
#define XOBJ_TYPEDEF_RecordIndex
typedef   struct RecordIndex_s * RecordIndex;
 static inline RecordIndex RecordIndex_getp(RecordIndex restrict  obj){
  int64_t offset=0;
  return (RecordIndex)(( char*) obj+offset);
}
 static inline int64_t RecordIndex_get_capacity(const RecordIndex restrict  obj){
  int64_t offset=0;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void RecordIndex_set_capacity(RecordIndex restrict  obj, int64_t value){
  int64_t offset=0;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* RecordIndex_getp_capacity(RecordIndex restrict  obj){
  int64_t offset=0;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline uint32_t RecordIndex_get_num_recorded(const RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void RecordIndex_set_num_recorded(RecordIndex restrict  obj, uint32_t value){
  int64_t offset=0;
  offset+=8;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* RecordIndex_getp_num_recorded(RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline uint32_t RecordIndex_get__dummy(const RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void RecordIndex_set__dummy(RecordIndex restrict  obj, uint32_t value){
  int64_t offset=0;
  offset+=16;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* RecordIndex_getp__dummy(RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline int64_t RecordIndex_get_buffer_id(const RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void RecordIndex_set_buffer_id(RecordIndex restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=24;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* RecordIndex_getp_buffer_id(RecordIndex restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( int64_t*)(( char*) obj+offset);
}
#endif


 static inline
int64_t RecordIndex_get_slot(RecordIndex record_index){

    if (record_index == NULL){
        return -2;
    }
    int64_t capacity = RecordIndex_get_capacity(record_index);
      uint32_t* num_recorded = RecordIndex_getp_num_recorded(record_index);

    if(*num_recorded >= capacity){
        return -1;
    }

    uint32_t slot;

//    slot = atomic_add(num_recorded, 1);   //only_for_context opencl
//    slot = atomicAdd(num_recorded, 1);    //only_for_context cuda

//    #pragma omp atomic capture            //only_for_context cpu_openmp
    slot = (*num_recorded)++;             //only_for_context cpu_serial cpu_openmp

    if (slot >= capacity){
        *num_recorded = capacity;
        return -1;
    }
    return (int64_t)slot;
}

#ifndef XOBJ_TYPEDEF_ParticlesData
#define XOBJ_TYPEDEF_ParticlesData
typedef   struct ParticlesData_s * ParticlesData;
 static inline ParticlesData ParticlesData_getp(ParticlesData restrict  obj){
  int64_t offset=0;
  return (ParticlesData)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_get__capacity(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__capacity(ParticlesData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp__capacity(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_get__num_active_particles(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__num_active_particles(ParticlesData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=16;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp__num_active_particles(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_get__num_lost_particles(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__num_lost_particles(ParticlesData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=24;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp__num_lost_particles(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_get_start_tracking_at_element(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_start_tracking_at_element(ParticlesData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=32;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp_start_tracking_at_element(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline double ParticlesData_get_q0(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=40;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_q0(ParticlesData restrict  obj, double value){
  int64_t offset=0;
  offset+=40;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp_q0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=40;
  return ( double*)(( char*) obj+offset);
}
 static inline double ParticlesData_get_mass0(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=48;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_mass0(ParticlesData restrict  obj, double value){
  int64_t offset=0;
  offset+=48;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp_mass0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=48;
  return ( double*)(( char*) obj+offset);
}
 static inline double ParticlesData_get_t_sim(const ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=56;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_t_sim(ParticlesData restrict  obj, double value){
  int64_t offset=0;
  offset+=56;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp_t_sim(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=56;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_p0c(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=280;
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_p0c(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=280;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_p0c(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=280;
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_p0c(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=280;
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_p0c(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=280;
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_gamma0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+64);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_gamma0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+64);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_gamma0(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+64);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_gamma0(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+64);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_gamma0(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+64);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_beta0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+72);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_beta0(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+72);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_beta0(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+72);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_beta0(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+72);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_beta0(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+72);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_s(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_s(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_s(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_s(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_s(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_zeta(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_zeta(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_zeta(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_zeta(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_zeta(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_x(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_x(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_x(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_x(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_x(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_y(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_y(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_y(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_y(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_y(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_px(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_px(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_px(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_px(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_px(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_py(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_py(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_py(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_py(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_py(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_ptau(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_ptau(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_ptau(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_ptau(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_ptau(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_delta(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_delta(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_delta(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_delta(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_delta(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_rpp(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_rpp(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_rpp(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_rpp(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_rpp(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_rvv(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_rvv(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_rvv(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_rvv(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_rvv(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_chi(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_chi(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_chi(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_chi(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_chi(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_charge_ratio(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_charge_ratio(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_charge_ratio(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_charge_ratio(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_charge_ratio(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_weight(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_weight(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_weight(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_weight(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_weight(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_ax(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_ax(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_ax(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_ax(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_ax(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 ParticlesData_getp_ay(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_ay(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double ParticlesData_get_ay(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_ay(ParticlesData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* ParticlesData_getp1_ay(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_pdg_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_pdg_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_pdg_id(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_pdg_id(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_pdg_id(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_particle_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_particle_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_particle_id(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_particle_id(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_particle_id(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_at_element(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_at_element(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_at_element(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_at_element(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_at_element(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_at_turn(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_at_turn(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_at_turn(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_at_turn(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_at_turn(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_state(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_state(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_state(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_state(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_state(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 ParticlesData_getp_parent_particle_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len_parent_particle_id(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t ParticlesData_get_parent_particle_id(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set_parent_particle_id(ParticlesData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* ParticlesData_getp1_parent_particle_id(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNUint32 ParticlesData_getp__rng_s1(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  return (ArrNUint32)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len__rng_s1(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline uint32_t ParticlesData_get__rng_s1(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*4;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__rng_s1(ParticlesData restrict  obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*4;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* ParticlesData_getp1__rng_s1(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*4;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline ArrNUint32 ParticlesData_getp__rng_s2(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  return (ArrNUint32)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len__rng_s2(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline uint32_t ParticlesData_get__rng_s2(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*4;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__rng_s2(ParticlesData restrict  obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*4;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* ParticlesData_getp1__rng_s2(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*4;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline ArrNUint32 ParticlesData_getp__rng_s3(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  return (ArrNUint32)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len__rng_s3(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline uint32_t ParticlesData_get__rng_s3(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*4;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__rng_s3(ParticlesData restrict  obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*4;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* ParticlesData_getp1__rng_s3(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*4;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline ArrNUint32 ParticlesData_getp__rng_s4(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  return (ArrNUint32)(( char*) obj+offset);
}
 static inline int64_t ParticlesData_len__rng_s4(ParticlesData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline uint32_t ParticlesData_get__rng_s4(const ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*4;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void ParticlesData_set__rng_s4(ParticlesData restrict  obj, int64_t i0, uint32_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*4;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* ParticlesData_getp1__rng_s4(ParticlesData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*4;
  return ( uint32_t*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2024.                 //
// ######################################### //

#ifndef XTRACK_BASE_RNG_H
#define XTRACK_BASE_RNG_H

//#include <stdint.h> //only_for_context none

// Combined LCG-Thausworthe generator from (example 37-4):
// https://developer.nvidia.com/gpugems/gpugems3/part-vi-gpu-computing/chapter-37-efficient-random-number-generation-and-application
#define MASK 0xffffffffUL
#define TAUSWORTHE(s,a,b,c,d) ((((s) &c) <<d) &MASK) ^ (((((s) <<a) &MASK)^(s)) >>b)
#define LCG(s,A,C) ((((A*(s)) &MASK) + C) &MASK)

 static inline
uint32_t rng_get_int32 (uint32_t *s1, uint32_t *s2, uint32_t *s3, uint32_t *s4 ){
  *s1 = TAUSWORTHE (*s1, 13, 19, 4294967294UL, 12);  // p1=2^31-1
  *s2 = TAUSWORTHE (*s2, 2, 25, 4294967288UL, 4);    // p2=2^30-1
  *s3 = TAUSWORTHE (*s3, 3, 11, 4294967280UL, 17);   // p3=2^28-1
  *s4 = LCG(*s4, 1664525, 1013904223UL);             // p4=2^32

  // Combined period is lcm(p1,p2,p3,p4) ~ 2^121
  return ((*s1) ^ (*s2) ^ (*s3) ^ (*s4));
}

#ifndef TWO_TO_32
#define TWO_TO_32 4294967296.0
#endif

 static inline
double rng_get (uint32_t *s1, uint32_t *s2, uint32_t *s3, uint32_t *s4 ){

  return rng_get_int32(s1, s2, s3, s4) / TWO_TO_32; // uniform in [0, 1) 1e10 resolution

}

 static inline
void rng_set (uint32_t *s1, uint32_t *s2, uint32_t *s3, uint32_t *s4, uint32_t s ){
  *s1 = LCG (s, 69069, 0);
  if (*s1 < 2) *s1 += 2UL;
  *s2 = LCG (*s1, 69069, 0);
  if (*s2 < 8) *s2 += 8UL;
  *s3 = LCG (*s2, 69069, 0);
  if (*s3 < 16) *s3 += 16UL;
  *s4 = LCG (*s3, 69069, 0);

  /* "warm it up" */
  rng_get (s1, s2, s3, s4);
  rng_get (s1, s2, s3, s4);
  rng_get (s1, s2, s3, s4);
  rng_get (s1, s2, s3, s4);
  rng_get (s1, s2, s3, s4);
  rng_get (s1, s2, s3, s4);
  return;
}

#endif /* XTRACK_BASE_RNG_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2024.                 //
// ######################################### //

#ifndef XTRACK_PARTICLES_RNG_H
#define XTRACK_PARTICLES_RNG_H

 
void Particles_initialize_rand_gen(ParticlesData particles,
      uint32_t* seeds, int n_init){

for (int ii=0; ii<n_init; ii++){ //autovectorized


    uint32_t s1, s2, s3, s4, s;
    s = seeds[ii];

    rng_set(&s1, &s2, &s3, &s4, s);

    ParticlesData_set__rng_s1(particles, ii, s1);
    ParticlesData_set__rng_s2(particles, ii, s2);
    ParticlesData_set__rng_s3(particles, ii, s3);
    ParticlesData_set__rng_s4(particles, ii, s4);

}//end autovectorized


}

#endif /* XTRACK_PARTICLES_RNG_H */


 typedef struct {
                 int64_t  _capacity;
                 int64_t  _num_active_particles;
                 int64_t  _num_lost_particles;
                 int64_t  start_tracking_at_element;
                 double  q0;
                 double  mass0;
                 double  t_sim;
      double* p0c;
      double* gamma0;
      double* beta0;
      double* s;
      double* zeta;
      double* x;
      double* y;
      double* px;
      double* py;
      double* ptau;
      double* delta;
      double* rpp;
      double* rvv;
      double* chi;
      double* charge_ratio;
      double* weight;
      double* ax;
      double* ay;
      int64_t* pdg_id;
      int64_t* particle_id;
      int64_t* at_element;
      int64_t* at_turn;
      int64_t* state;
      int64_t* parent_particle_id;
      uint32_t* _rng_s1;
      uint32_t* _rng_s2;
      uint32_t* _rng_s3;
      uint32_t* _rng_s4;
                 int64_t ipart;
                 int64_t endpart;
      int8_t* io_buffer;
} LocalParticle;


             static inline
            void LocalParticle_add_to_p0c(LocalParticle* part, double value){
#ifndef FREEZE_VAR_p0c
  part->p0c[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_gamma0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_gamma0
  part->gamma0[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_beta0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_beta0
  part->beta0[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_s(LocalParticle* part, double value){
#ifndef FREEZE_VAR_s
  part->s[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_zeta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_zeta
  part->zeta[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_x(LocalParticle* part, double value){
#ifndef FREEZE_VAR_x
  part->x[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_y(LocalParticle* part, double value){
#ifndef FREEZE_VAR_y
  part->y[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_px(LocalParticle* part, double value){
#ifndef FREEZE_VAR_px
  part->px[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_py(LocalParticle* part, double value){
#ifndef FREEZE_VAR_py
  part->py[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_ptau(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ptau
  part->ptau[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_delta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_delta
  part->delta[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_rpp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rpp
  part->rpp[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_rvv(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rvv
  part->rvv[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_chi(LocalParticle* part, double value){
#ifndef FREEZE_VAR_chi
  part->chi[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_charge_ratio(LocalParticle* part, double value){
#ifndef FREEZE_VAR_charge_ratio
  part->charge_ratio[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_weight(LocalParticle* part, double value){
#ifndef FREEZE_VAR_weight
  part->weight[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_ax(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ax
  part->ax[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_ay(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ay
  part->ay[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_pdg_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_pdg_id
  part->pdg_id[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_particle_id
  part->particle_id[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_at_element(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_element
  part->at_element[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_at_turn(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_turn
  part->at_turn[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_state(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_state
  part->state[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to_parent_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_parent_particle_id
  part->parent_particle_id[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to__rng_s1(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s1
  part->_rng_s1[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to__rng_s2(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s2
  part->_rng_s2[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to__rng_s3(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s3
  part->_rng_s3[part->ipart] += value;
#endif
}


             static inline
            void LocalParticle_add_to__rng_s4(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s4
  part->_rng_s4[part->ipart] += value;
#endif
}


 static inline
int64_t LocalParticle_get__capacity(LocalParticle* part){
  return part->_capacity;
}
 static inline
int64_t LocalParticle_get__num_active_particles(LocalParticle* part){
  return part->_num_active_particles;
}
 static inline
int64_t LocalParticle_get__num_lost_particles(LocalParticle* part){
  return part->_num_lost_particles;
}
 static inline
int64_t LocalParticle_get_start_tracking_at_element(LocalParticle* part){
  return part->start_tracking_at_element;
}
 static inline
double LocalParticle_get_q0(LocalParticle* part){
  return part->q0;
}
 static inline
double LocalParticle_get_mass0(LocalParticle* part){
  return part->mass0;
}
 static inline
double LocalParticle_get_t_sim(LocalParticle* part){
  return part->t_sim;
}
 static inline
double LocalParticle_get_p0c(LocalParticle* part){
  return part->p0c[part->ipart];
}
 static inline
double LocalParticle_get_gamma0(LocalParticle* part){
  return part->gamma0[part->ipart];
}
 static inline
double LocalParticle_get_beta0(LocalParticle* part){
  return part->beta0[part->ipart];
}
 static inline
double LocalParticle_get_s(LocalParticle* part){
  return part->s[part->ipart];
}
 static inline
double LocalParticle_get_zeta(LocalParticle* part){
  return part->zeta[part->ipart];
}
 static inline
double LocalParticle_get_x(LocalParticle* part){
  return part->x[part->ipart];
}
 static inline
double LocalParticle_get_y(LocalParticle* part){
  return part->y[part->ipart];
}
 static inline
double LocalParticle_get_px(LocalParticle* part){
  return part->px[part->ipart];
}
 static inline
double LocalParticle_get_py(LocalParticle* part){
  return part->py[part->ipart];
}
 static inline
double LocalParticle_get_ptau(LocalParticle* part){
  return part->ptau[part->ipart];
}
 static inline
double LocalParticle_get_delta(LocalParticle* part){
  return part->delta[part->ipart];
}
 static inline
double LocalParticle_get_rpp(LocalParticle* part){
  return part->rpp[part->ipart];
}
 static inline
double LocalParticle_get_rvv(LocalParticle* part){
  return part->rvv[part->ipart];
}
 static inline
double LocalParticle_get_chi(LocalParticle* part){
  return part->chi[part->ipart];
}
 static inline
double LocalParticle_get_charge_ratio(LocalParticle* part){
  return part->charge_ratio[part->ipart];
}
 static inline
double LocalParticle_get_weight(LocalParticle* part){
  return part->weight[part->ipart];
}
 static inline
double LocalParticle_get_ax(LocalParticle* part){
  return part->ax[part->ipart];
}
 static inline
double LocalParticle_get_ay(LocalParticle* part){
  return part->ay[part->ipart];
}
 static inline
int64_t LocalParticle_get_pdg_id(LocalParticle* part){
  return part->pdg_id[part->ipart];
}
 static inline
int64_t LocalParticle_get_particle_id(LocalParticle* part){
  return part->particle_id[part->ipart];
}
 static inline
int64_t LocalParticle_get_at_element(LocalParticle* part){
  return part->at_element[part->ipart];
}
 static inline
int64_t LocalParticle_get_at_turn(LocalParticle* part){
  return part->at_turn[part->ipart];
}
 static inline
int64_t LocalParticle_get_state(LocalParticle* part){
  return part->state[part->ipart];
}
 static inline
int64_t LocalParticle_get_parent_particle_id(LocalParticle* part){
  return part->parent_particle_id[part->ipart];
}
 static inline
uint32_t LocalParticle_get__rng_s1(LocalParticle* part){
  return part->_rng_s1[part->ipart];
}
 static inline
uint32_t LocalParticle_get__rng_s2(LocalParticle* part){
  return part->_rng_s2[part->ipart];
}
 static inline
uint32_t LocalParticle_get__rng_s3(LocalParticle* part){
  return part->_rng_s3[part->ipart];
}
 static inline
uint32_t LocalParticle_get__rng_s4(LocalParticle* part){
  return part->_rng_s4[part->ipart];
}


             static inline
            void LocalParticle_set_p0c(LocalParticle* part, double value){
#ifndef FREEZE_VAR_p0c
  part->p0c[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_gamma0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_gamma0
  part->gamma0[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_beta0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_beta0
  part->beta0[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_s(LocalParticle* part, double value){
#ifndef FREEZE_VAR_s
  part->s[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_zeta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_zeta
  part->zeta[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_x(LocalParticle* part, double value){
#ifndef FREEZE_VAR_x
  part->x[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_y(LocalParticle* part, double value){
#ifndef FREEZE_VAR_y
  part->y[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_px(LocalParticle* part, double value){
#ifndef FREEZE_VAR_px
  part->px[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_py(LocalParticle* part, double value){
#ifndef FREEZE_VAR_py
  part->py[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_ptau(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ptau
  part->ptau[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_delta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_delta
  part->delta[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_rpp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rpp
  part->rpp[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_rvv(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rvv
  part->rvv[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_chi(LocalParticle* part, double value){
#ifndef FREEZE_VAR_chi
  part->chi[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_charge_ratio(LocalParticle* part, double value){
#ifndef FREEZE_VAR_charge_ratio
  part->charge_ratio[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_weight(LocalParticle* part, double value){
#ifndef FREEZE_VAR_weight
  part->weight[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_ax(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ax
  part->ax[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_ay(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ay
  part->ay[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_pdg_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_pdg_id
  part->pdg_id[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_particle_id
  part->particle_id[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_at_element(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_element
  part->at_element[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_at_turn(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_turn
  part->at_turn[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_state(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_state
  part->state[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set_parent_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_parent_particle_id
  part->parent_particle_id[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set__rng_s1(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s1
  part->_rng_s1[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set__rng_s2(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s2
  part->_rng_s2[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set__rng_s3(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s3
  part->_rng_s3[part->ipart] = value;
#endif
}

             static inline
            void LocalParticle_set__rng_s4(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s4
  part->_rng_s4[part->ipart] = value;
#endif
}


             static inline
            void LocalParticle_scale_p0c(LocalParticle* part, double value){
#ifndef FREEZE_VAR_p0c
  part->p0c[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_gamma0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_gamma0
  part->gamma0[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_beta0(LocalParticle* part, double value){
#ifndef FREEZE_VAR_beta0
  part->beta0[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_s(LocalParticle* part, double value){
#ifndef FREEZE_VAR_s
  part->s[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_zeta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_zeta
  part->zeta[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_x(LocalParticle* part, double value){
#ifndef FREEZE_VAR_x
  part->x[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_y(LocalParticle* part, double value){
#ifndef FREEZE_VAR_y
  part->y[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_px(LocalParticle* part, double value){
#ifndef FREEZE_VAR_px
  part->px[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_py(LocalParticle* part, double value){
#ifndef FREEZE_VAR_py
  part->py[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_ptau(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ptau
  part->ptau[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_delta(LocalParticle* part, double value){
#ifndef FREEZE_VAR_delta
  part->delta[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_rpp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rpp
  part->rpp[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_rvv(LocalParticle* part, double value){
#ifndef FREEZE_VAR_rvv
  part->rvv[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_chi(LocalParticle* part, double value){
#ifndef FREEZE_VAR_chi
  part->chi[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_charge_ratio(LocalParticle* part, double value){
#ifndef FREEZE_VAR_charge_ratio
  part->charge_ratio[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_weight(LocalParticle* part, double value){
#ifndef FREEZE_VAR_weight
  part->weight[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_ax(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ax
  part->ax[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_ay(LocalParticle* part, double value){
#ifndef FREEZE_VAR_ay
  part->ay[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_pdg_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_pdg_id
  part->pdg_id[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_particle_id
  part->particle_id[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_at_element(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_element
  part->at_element[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_at_turn(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_at_turn
  part->at_turn[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_state(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_state
  part->state[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale_parent_particle_id(LocalParticle* part, int64_t value){
#ifndef FREEZE_VAR_parent_particle_id
  part->parent_particle_id[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale__rng_s1(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s1
  part->_rng_s1[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale__rng_s2(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s2
  part->_rng_s2[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale__rng_s3(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s3
  part->_rng_s3[part->ipart] *= value;
#endif
}


             static inline
            void LocalParticle_scale__rng_s4(LocalParticle* part, uint32_t value){
#ifndef FREEZE_VAR__rng_s4
  part->_rng_s4[part->ipart] *= value;
#endif
}



         static inline
        void LocalParticle_exchange(LocalParticle* part, int64_t i1, int64_t i2){
        
    {
    double temp = part->p0c[i2];
    part->p0c[i2] = part->p0c[i1];
    part->p0c[i1] = temp;
     }
    {
    double temp = part->gamma0[i2];
    part->gamma0[i2] = part->gamma0[i1];
    part->gamma0[i1] = temp;
     }
    {
    double temp = part->beta0[i2];
    part->beta0[i2] = part->beta0[i1];
    part->beta0[i1] = temp;
     }
    {
    double temp = part->s[i2];
    part->s[i2] = part->s[i1];
    part->s[i1] = temp;
     }
    {
    double temp = part->zeta[i2];
    part->zeta[i2] = part->zeta[i1];
    part->zeta[i1] = temp;
     }
    {
    double temp = part->x[i2];
    part->x[i2] = part->x[i1];
    part->x[i1] = temp;
     }
    {
    double temp = part->y[i2];
    part->y[i2] = part->y[i1];
    part->y[i1] = temp;
     }
    {
    double temp = part->px[i2];
    part->px[i2] = part->px[i1];
    part->px[i1] = temp;
     }
    {
    double temp = part->py[i2];
    part->py[i2] = part->py[i1];
    part->py[i1] = temp;
     }
    {
    double temp = part->ptau[i2];
    part->ptau[i2] = part->ptau[i1];
    part->ptau[i1] = temp;
     }
    {
    double temp = part->delta[i2];
    part->delta[i2] = part->delta[i1];
    part->delta[i1] = temp;
     }
    {
    double temp = part->rpp[i2];
    part->rpp[i2] = part->rpp[i1];
    part->rpp[i1] = temp;
     }
    {
    double temp = part->rvv[i2];
    part->rvv[i2] = part->rvv[i1];
    part->rvv[i1] = temp;
     }
    {
    double temp = part->chi[i2];
    part->chi[i2] = part->chi[i1];
    part->chi[i1] = temp;
     }
    {
    double temp = part->charge_ratio[i2];
    part->charge_ratio[i2] = part->charge_ratio[i1];
    part->charge_ratio[i1] = temp;
     }
    {
    double temp = part->weight[i2];
    part->weight[i2] = part->weight[i1];
    part->weight[i1] = temp;
     }
    {
    double temp = part->ax[i2];
    part->ax[i2] = part->ax[i1];
    part->ax[i1] = temp;
     }
    {
    double temp = part->ay[i2];
    part->ay[i2] = part->ay[i1];
    part->ay[i1] = temp;
     }
    {
    int64_t temp = part->pdg_id[i2];
    part->pdg_id[i2] = part->pdg_id[i1];
    part->pdg_id[i1] = temp;
     }
    {
    int64_t temp = part->particle_id[i2];
    part->particle_id[i2] = part->particle_id[i1];
    part->particle_id[i1] = temp;
     }
    {
    int64_t temp = part->at_element[i2];
    part->at_element[i2] = part->at_element[i1];
    part->at_element[i1] = temp;
     }
    {
    int64_t temp = part->at_turn[i2];
    part->at_turn[i2] = part->at_turn[i1];
    part->at_turn[i1] = temp;
     }
    {
    int64_t temp = part->state[i2];
    part->state[i2] = part->state[i1];
    part->state[i1] = temp;
     }
    {
    int64_t temp = part->parent_particle_id[i2];
    part->parent_particle_id[i2] = part->parent_particle_id[i1];
    part->parent_particle_id[i1] = temp;
     }
    {
    uint32_t temp = part->_rng_s1[i2];
    part->_rng_s1[i2] = part->_rng_s1[i1];
    part->_rng_s1[i1] = temp;
     }
    {
    uint32_t temp = part->_rng_s2[i2];
    part->_rng_s2[i2] = part->_rng_s2[i1];
    part->_rng_s2[i1] = temp;
     }
    {
    uint32_t temp = part->_rng_s3[i2];
    part->_rng_s3[i2] = part->_rng_s3[i1];
    part->_rng_s3[i1] = temp;
     }
    {
    uint32_t temp = part->_rng_s4[i2];
    part->_rng_s4[i2] = part->_rng_s4[i1];
    part->_rng_s4[i1] = temp;
     }}



             static inline
              int8_t* LocalParticle_get_io_buffer(LocalParticle* part){
                return part->io_buffer;
            }

            

             static inline
            void Particles_to_LocalParticle(ParticlesData source,
                                            LocalParticle* dest,
                                            int64_t id,
                                            int64_t eid){
  dest->_capacity = ParticlesData_get__capacity(source);
  dest->_num_active_particles = ParticlesData_get__num_active_particles(source);
  dest->_num_lost_particles = ParticlesData_get__num_lost_particles(source);
  dest->start_tracking_at_element = ParticlesData_get_start_tracking_at_element(source);
  dest->q0 = ParticlesData_get_q0(source);
  dest->mass0 = ParticlesData_get_mass0(source);
  dest->t_sim = ParticlesData_get_t_sim(source);
  dest->p0c = ParticlesData_getp1_p0c(source, 0);
  dest->gamma0 = ParticlesData_getp1_gamma0(source, 0);
  dest->beta0 = ParticlesData_getp1_beta0(source, 0);
  dest->s = ParticlesData_getp1_s(source, 0);
  dest->zeta = ParticlesData_getp1_zeta(source, 0);
  dest->x = ParticlesData_getp1_x(source, 0);
  dest->y = ParticlesData_getp1_y(source, 0);
  dest->px = ParticlesData_getp1_px(source, 0);
  dest->py = ParticlesData_getp1_py(source, 0);
  dest->ptau = ParticlesData_getp1_ptau(source, 0);
  dest->delta = ParticlesData_getp1_delta(source, 0);
  dest->rpp = ParticlesData_getp1_rpp(source, 0);
  dest->rvv = ParticlesData_getp1_rvv(source, 0);
  dest->chi = ParticlesData_getp1_chi(source, 0);
  dest->charge_ratio = ParticlesData_getp1_charge_ratio(source, 0);
  dest->weight = ParticlesData_getp1_weight(source, 0);
  dest->ax = ParticlesData_getp1_ax(source, 0);
  dest->ay = ParticlesData_getp1_ay(source, 0);
  dest->pdg_id = ParticlesData_getp1_pdg_id(source, 0);
  dest->particle_id = ParticlesData_getp1_particle_id(source, 0);
  dest->at_element = ParticlesData_getp1_at_element(source, 0);
  dest->at_turn = ParticlesData_getp1_at_turn(source, 0);
  dest->state = ParticlesData_getp1_state(source, 0);
  dest->parent_particle_id = ParticlesData_getp1_parent_particle_id(source, 0);
  dest->_rng_s1 = ParticlesData_getp1__rng_s1(source, 0);
  dest->_rng_s2 = ParticlesData_getp1__rng_s2(source, 0);
  dest->_rng_s3 = ParticlesData_getp1__rng_s3(source, 0);
  dest->_rng_s4 = ParticlesData_getp1__rng_s4(source, 0);
  dest->ipart = id;
  dest->endpart = eid;
}


             static inline
            void LocalParticle_to_Particles(
                                            LocalParticle* source,
                                            ParticlesData dest,
                                            int64_t id,
                                            int64_t set_scalar){
if (set_scalar){
  ParticlesData_set__capacity(dest,      LocalParticle_get__capacity(source));
  ParticlesData_set__num_active_particles(dest,      LocalParticle_get__num_active_particles(source));
  ParticlesData_set__num_lost_particles(dest,      LocalParticle_get__num_lost_particles(source));
  ParticlesData_set_start_tracking_at_element(dest,      LocalParticle_get_start_tracking_at_element(source));
  ParticlesData_set_q0(dest,      LocalParticle_get_q0(source));
  ParticlesData_set_mass0(dest,      LocalParticle_get_mass0(source));
  ParticlesData_set_t_sim(dest,      LocalParticle_get_t_sim(source));
}
  ParticlesData_set_p0c(dest, id,       LocalParticle_get_p0c(source));
  ParticlesData_set_gamma0(dest, id,       LocalParticle_get_gamma0(source));
  ParticlesData_set_beta0(dest, id,       LocalParticle_get_beta0(source));
  ParticlesData_set_s(dest, id,       LocalParticle_get_s(source));
  ParticlesData_set_zeta(dest, id,       LocalParticle_get_zeta(source));
  ParticlesData_set_x(dest, id,       LocalParticle_get_x(source));
  ParticlesData_set_y(dest, id,       LocalParticle_get_y(source));
  ParticlesData_set_px(dest, id,       LocalParticle_get_px(source));
  ParticlesData_set_py(dest, id,       LocalParticle_get_py(source));
  ParticlesData_set_ptau(dest, id,       LocalParticle_get_ptau(source));
  ParticlesData_set_delta(dest, id,       LocalParticle_get_delta(source));
  ParticlesData_set_rpp(dest, id,       LocalParticle_get_rpp(source));
  ParticlesData_set_rvv(dest, id,       LocalParticle_get_rvv(source));
  ParticlesData_set_chi(dest, id,       LocalParticle_get_chi(source));
  ParticlesData_set_charge_ratio(dest, id,       LocalParticle_get_charge_ratio(source));
  ParticlesData_set_weight(dest, id,       LocalParticle_get_weight(source));
  ParticlesData_set_ax(dest, id,       LocalParticle_get_ax(source));
  ParticlesData_set_ay(dest, id,       LocalParticle_get_ay(source));
  ParticlesData_set_pdg_id(dest, id,       LocalParticle_get_pdg_id(source));
  ParticlesData_set_particle_id(dest, id,       LocalParticle_get_particle_id(source));
  ParticlesData_set_at_element(dest, id,       LocalParticle_get_at_element(source));
  ParticlesData_set_at_turn(dest, id,       LocalParticle_get_at_turn(source));
  ParticlesData_set_state(dest, id,       LocalParticle_get_state(source));
  ParticlesData_set_parent_particle_id(dest, id,       LocalParticle_get_parent_particle_id(source));
  ParticlesData_set__rng_s1(dest, id,       LocalParticle_get__rng_s1(source));
  ParticlesData_set__rng_s2(dest, id,       LocalParticle_get__rng_s2(source));
  ParticlesData_set__rng_s3(dest, id,       LocalParticle_get__rng_s3(source));
  ParticlesData_set__rng_s4(dest, id,       LocalParticle_get__rng_s4(source));
}

 static inline
double LocalParticle_get_xp(LocalParticle* part){
    double const px = LocalParticle_get_px(part);
    double const rpp = LocalParticle_get_rpp(part);
    // INFO: this is not the angle, but sin(angle)
    return px*rpp;
}

 static inline
double LocalParticle_get_yp(LocalParticle* part){
    double const py = LocalParticle_get_py(part);
    double const rpp = LocalParticle_get_rpp(part);
    // INFO: this is not the angle, but sin(angle)
    return py*rpp;
}

 static inline
void LocalParticle_set_xp(LocalParticle* part, double xp){
#ifndef FREEZE_VAR_px
    double rpp = LocalParticle_get_rpp(part);
    // INFO: xp is not the angle, but sin(angle)
    LocalParticle_set_px(part, xp/rpp);
#endif
}

 static inline
void LocalParticle_set_yp(LocalParticle* part, double yp){
#ifndef FREEZE_VAR_py
    double rpp = LocalParticle_get_rpp(part);
    // INFO: yp is not the angle, but sin(angle)
    LocalParticle_set_py(part, yp/rpp);
#endif
}

 static inline
void LocalParticle_add_to_xp(LocalParticle* part, double xp){
#ifndef FREEZE_VAR_px
    LocalParticle_set_xp(part, LocalParticle_get_xp(part) + xp);
#endif
}

 static inline
void LocalParticle_scale_xp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_px
    LocalParticle_set_xp(part, LocalParticle_get_xp(part) * value);
#endif
}

 static inline
void LocalParticle_add_to_yp(LocalParticle* part, double yp){
#ifndef FREEZE_VAR_py
    LocalParticle_set_yp(part, LocalParticle_get_yp(part) + yp);
#endif
}

 static inline
void LocalParticle_scale_yp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_py
    LocalParticle_set_yp(part, LocalParticle_get_yp(part) * value);
#endif
}

 static inline
void LocalParticle_set_xp_yp(LocalParticle* part, double xp, double yp){
    double rpp = LocalParticle_get_rpp(part);
#ifndef FREEZE_VAR_px
    LocalParticle_set_px(part, xp/rpp);
#endif
#ifndef FREEZE_VAR_py
    LocalParticle_set_py(part, yp/rpp);
#endif
}

 static inline
void LocalParticle_add_to_xp_yp(LocalParticle* part, double xp, double yp){
    LocalParticle_set_xp_yp(part, LocalParticle_get_xp(part) + xp, LocalParticle_get_yp(part) + yp);
}

 static inline
void LocalParticle_scale_xp_yp(LocalParticle* part, double value_x, double value_y){
    LocalParticle_set_xp_yp(part, LocalParticle_get_xp(part) * value_x, LocalParticle_get_yp(part) * value_y);
}
 static inline
double LocalParticle_get_exact_xp(LocalParticle* part){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);
    double const rpp = 1./sqrt(one_plus_delta*one_plus_delta - px*px - py*py);
    // INFO: this is not the angle, but sin(angle)
    return px*rpp;
}

 static inline
double LocalParticle_get_exact_yp(LocalParticle* part){
    double const py = LocalParticle_get_py(part);
    double const px = LocalParticle_get_px(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);
    double const rpp = 1./sqrt(one_plus_delta*one_plus_delta - px*px - py*py);
    // INFO: this is not the angle, but sin(angle)
    return py*rpp;
}

 static inline
void LocalParticle_set_exact_xp(LocalParticle* part, double xp){
#ifndef FREEZE_VAR_px
    double rpp = LocalParticle_get_rpp(part);
    // Careful! If yp also changes, use LocalParticle_set_exact_xp_yp!
    double const yp = LocalParticle_get_exact_yp(part);
    rpp *= sqrt(1 + xp*xp + yp*yp);
    // INFO: xp is not the angle, but sin(angle)
    LocalParticle_set_px(part, xp/rpp);
#endif
}

 static inline
void LocalParticle_set_exact_yp(LocalParticle* part, double yp){
#ifndef FREEZE_VAR_py
    double rpp = LocalParticle_get_rpp(part);
    // Careful! If xp also changes, use LocalParticle_set_exact_xp_yp!
    double const xp = LocalParticle_get_exact_xp(part);
    rpp *= sqrt(1 + xp*xp + yp*yp);
    // INFO: yp is not the angle, but sin(angle)
    LocalParticle_set_py(part, yp/rpp);
#endif
}

 static inline
void LocalParticle_add_to_exact_xp(LocalParticle* part, double xp){
#ifndef FREEZE_VAR_px
    LocalParticle_set_exact_xp(part, LocalParticle_get_exact_xp(part) + xp);
#endif
}

 static inline
void LocalParticle_scale_exact_xp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_px
    LocalParticle_set_exact_xp(part, LocalParticle_get_exact_xp(part) * value);
#endif
}

 static inline
void LocalParticle_add_to_exact_yp(LocalParticle* part, double yp){
#ifndef FREEZE_VAR_py
    LocalParticle_set_exact_yp(part, LocalParticle_get_exact_yp(part) + yp);
#endif
}

 static inline
void LocalParticle_scale_exact_yp(LocalParticle* part, double value){
#ifndef FREEZE_VAR_py
    LocalParticle_set_exact_yp(part, LocalParticle_get_exact_yp(part) * value);
#endif
}

 static inline
void LocalParticle_set_exact_xp_yp(LocalParticle* part, double xp, double yp){
    double rpp = LocalParticle_get_rpp(part);
    rpp *= sqrt(1 + xp*xp + yp*yp);
#ifndef FREEZE_VAR_px
    LocalParticle_set_px(part, xp/rpp);
#endif
#ifndef FREEZE_VAR_py
    LocalParticle_set_py(part, yp/rpp);
#endif
}

 static inline
void LocalParticle_add_to_exact_xp_yp(LocalParticle* part, double xp, double yp){
    LocalParticle_set_exact_xp_yp(part, LocalParticle_get_exact_xp(part) + xp, LocalParticle_get_exact_yp(part) + yp);
}

 static inline
void LocalParticle_scale_exact_xp_yp(LocalParticle* part, double value_x, double value_y){
    LocalParticle_set_exact_xp_yp(part, LocalParticle_get_exact_xp(part) * value_x, LocalParticle_get_exact_yp(part) * value_y);
}


         static inline
        double LocalParticle_get_energy0(LocalParticle* part){

            double const p0c = LocalParticle_get_p0c(part);
            double const m0  = LocalParticle_get_mass0(part);

            return sqrt( p0c * p0c + m0 * m0 );
        }

         static inline
        void LocalParticle_update_ptau(LocalParticle* part, double new_ptau_value){

            double const beta0 = LocalParticle_get_beta0(part);

            double const ptau = new_ptau_value;

            double const irpp = sqrt(ptau*ptau + 2*ptau/beta0 +1);

            double const new_rpp = 1./irpp;
            LocalParticle_set_delta(part, irpp - 1.);

            double const new_rvv = irpp/(1 + beta0*ptau);
            LocalParticle_set_rvv(part, new_rvv);
            LocalParticle_set_ptau(part, ptau);

            LocalParticle_set_rpp(part, new_rpp );
        }

         static inline
        void LocalParticle_update_delta(LocalParticle* part, double new_delta_value){
            double const beta0 = LocalParticle_get_beta0(part);
            double const delta_beta0 = new_delta_value * beta0;
            double const ptau_beta0  = sqrt( delta_beta0 * delta_beta0 +
                                        2. * delta_beta0 * beta0 + 1. ) - 1.;

            double const one_plus_delta = 1. + new_delta_value;
            double const rvv    = ( one_plus_delta ) / ( 1. + ptau_beta0 );
            double const rpp    = 1. / one_plus_delta;
            double const ptau = ptau_beta0 / beta0;

            LocalParticle_set_delta(part, new_delta_value);

            LocalParticle_set_rvv(part, rvv );
            LocalParticle_set_rpp(part, rpp );
            LocalParticle_set_ptau(part, ptau );

        }

         static inline
        double LocalParticle_get_pzeta(LocalParticle* part){

            double const ptau = LocalParticle_get_ptau(part);
            double const beta0 = LocalParticle_get_beta0(part);

            return ptau/beta0;

        }

         static inline
        void LocalParticle_update_pzeta(LocalParticle* part, double new_pzeta_value){

            double const beta0 = LocalParticle_get_beta0(part);
            LocalParticle_update_ptau(part, beta0*new_pzeta_value);

        }

         static inline
        void increment_at_element(LocalParticle* part0, int64_t const increment){


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

                LocalParticle_add_to_at_element(part, increment);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }


        }

         static inline
        void increment_at_turn(LocalParticle* part0, int flag_reset_s){


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

            LocalParticle_add_to_at_turn(part, 1);
            LocalParticle_set_at_element(part, 0);
            if (flag_reset_s>0){
                LocalParticle_set_s(part, 0.);
            }

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

        }

         static inline
        void increment_at_turn_backtrack(LocalParticle* part0, int flag_reset_s,
                                         double const line_length,
                                         int64_t const num_elements){


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

            LocalParticle_add_to_at_turn(part, -1);
            LocalParticle_set_at_element(part, num_elements);
            if (flag_reset_s>0){
                LocalParticle_set_s(part, line_length);
            }

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

        }

        // check_is_active has different implementation on CPU and GPU

        #define CPU_SERIAL_IMPLEM //only_for_context cpu_serial
//        #define CPU_OMP_IMPLEM //only_for_context cpu_openmp

        #ifdef CPU_SERIAL_IMPLEM

         static inline
        int64_t check_is_active(LocalParticle* part) {
            int64_t ipart=0;
            while (ipart < part->_num_active_particles){
                if (part->state[ipart]<1){
                    LocalParticle_exchange(
                        part, ipart, part->_num_active_particles-1);
                    part->_num_active_particles--;
                    part->_num_lost_particles++;
                }
                else{
                    ipart++;
                }
            }

            if (part->_num_active_particles==0){
                return 0;//All particles lost
            } else {
                return 1; //Some stable particles are still present
            }
        }

        #else // not CPU_SERIAL_IMPLEM
        #ifdef CPU_OMP_IMPLEM

         static inline
        int64_t check_is_active(LocalParticle* part) {
        #ifndef SKIP_SWAPS
            int64_t ipart = part->ipart;
            int64_t endpart = part->endpart;

            int64_t left = ipart;
            int64_t right = endpart - 1;
            int64_t swap_made = 0;
            int64_t has_alive = 0;

            if (left == right) return part->state[left] > 0;

            while (left < right) {
                if (part->state[left] > 0) {
                    left++;
                    has_alive = 1;
                }
                else if (part->state[right] <= 0) right--;
                else {
                    LocalParticle_exchange(part, left, right);
                    left++;
                    right--;
                    swap_made = 1;
                }
            }

            return swap_made || has_alive;
        #else
            return 1;
        #endif
        }

         static inline
        void count_reorganized_particles(LocalParticle* part) {
            int64_t num_active = 0;
            int64_t num_lost = 0;

            for (int64_t i = part->ipart; i < part->endpart; i++) {
                if (part->state[i] <= -999999999) break;
                else if (part->state[i] > 0) num_active++;
                else num_lost++;
            }

            part->_num_active_particles = num_active;
            part->_num_lost_particles = num_lost;
        }

        #else // not CPU_SERIAL_IMPLEM and not CPU_OMP_IMPLEM

         static inline
        int64_t check_is_active(LocalParticle* part) {
            return LocalParticle_get_state(part)>0;
        };

        #endif // CPU_OMP_IMPLEM
        #endif // CPU_SERIAL_IMPLEM

        #undef CPU_SERIAL_IMPLEM //only_for_context cpu_serial
//        #undef CPU_OMP_IMPLEM //only_for_context cpu_openmp


        
                    #ifdef XTRACK_GLOBAL_XY_LIMIT

                     static inline
                    void global_aperture_check(LocalParticle* part0) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

                            double const x = LocalParticle_get_x(part);
                            double const y = LocalParticle_get_y(part);

                        int64_t const is_alive = (int64_t)(
                                  (x >= -XTRACK_GLOBAL_XY_LIMIT) &&
                                  (x <=  XTRACK_GLOBAL_XY_LIMIT) &&
                                  (y >= -XTRACK_GLOBAL_XY_LIMIT) &&
                                  (y <=  XTRACK_GLOBAL_XY_LIMIT) );

                        // I assume that if I am in the function is because
                            if (!is_alive){
                               LocalParticle_set_state(part, -1);
                        }

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

                    }

                    #endif

                     static inline
                    void LocalParticle_add_to_energy(LocalParticle* part, double delta_energy, int pz_only ){
                        double ptau = LocalParticle_get_ptau(part);
                        double const p0c = LocalParticle_get_p0c(part);
                        double const charge_ratio = LocalParticle_get_charge_ratio(part);
                        double const chi = LocalParticle_get_chi(part);
                        double const mass_ratio = chi / charge_ratio;
                        
                        ptau += delta_energy/p0c * mass_ratio;

                        double const old_rpp = LocalParticle_get_rpp(part);

                        LocalParticle_update_ptau(part, ptau);

                        if (!pz_only) {
                            double const new_rpp = LocalParticle_get_rpp(part);
                            double const f = old_rpp / new_rpp;
                            LocalParticle_scale_px(part, f);
                            LocalParticle_scale_py(part, f);
                        }
                    }


                     static inline
                    void LocalParticle_update_p0c(LocalParticle* part, double new_p0c_value){

                        double const mass0 = LocalParticle_get_mass0(part);
                        double const old_p0c = LocalParticle_get_p0c(part);
                        double const old_delta = LocalParticle_get_delta(part);
                        double const old_beta0 = LocalParticle_get_beta0(part);

                        double const ppc = old_p0c * old_delta + old_p0c;
                        double const new_delta = (ppc - new_p0c_value)/new_p0c_value;

                        double const new_energy0 = sqrt(new_p0c_value*new_p0c_value + mass0 * mass0);
                        double const new_beta0 = new_p0c_value / new_energy0;
                        double const new_gamma0 = new_energy0 / mass0;

                        LocalParticle_set_p0c(part, new_p0c_value);
                        LocalParticle_set_gamma0(part, new_gamma0);
                        LocalParticle_set_beta0(part, new_beta0);

                        LocalParticle_update_delta(part, new_delta);

                        LocalParticle_scale_px(part, old_p0c/new_p0c_value);
                        LocalParticle_scale_py(part, old_p0c/new_p0c_value);

                        LocalParticle_scale_zeta(part, new_beta0/old_beta0);

                    }

                     static inline
                    void LocalParticle_kill_particle(LocalParticle* part, int64_t kill_state) {
                        LocalParticle_set_x(part, 1e30);
                        LocalParticle_set_px(part, 1e30);
                        LocalParticle_set_y(part, 1e30);
                        LocalParticle_set_py(part, 1e30);
                        LocalParticle_set_zeta(part, 1e30);
                        LocalParticle_update_delta(part, -1);  // zero energy
                        LocalParticle_set_state(part, kill_state);
                    }
                 

#ifndef XOBJ_TYPEDEF_DriftData
#define XOBJ_TYPEDEF_DriftData
typedef   struct DriftData_s * DriftData;
 static inline DriftData DriftData_getp(DriftData restrict  obj){
  int64_t offset=0;
  return (DriftData)(( char*) obj+offset);
}
 static inline double DriftData_get_length(const DriftData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void DriftData_set_length(DriftData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* DriftData_getp_length(DriftData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_ELEM_H
#define XTRACK_DRIFT_ELEM_H

 static inline
void Drift_track_local_particle(DriftData el, LocalParticle* part0){

    double length = DriftData_get_length(el);
    #ifdef XSUITE_BACKTRACK
        length = -length;
    #endif


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        Drift_single_particle(part, length);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }


}


#endif /* XTRACK_DRIFT_ELEM_H */


             static inline
            void Drift_track_local_particle_with_transformations(DriftData el, LocalParticle* part0){
    Drift_track_local_particle(el, part0);
}


             
            void Drift_track_particles(
               DriftData el,

                             ParticlesData particles,

                             int64_t flag_increment_at_element,
                  int8_t* io_buffer){

//            #define CONTEXT_OPENMP  //only_for_context cpu_openmp
            #ifdef CONTEXT_OPENMP
                const int64_t capacity = ParticlesData_get__capacity(particles);
                const int num_threads = omp_get_max_threads();

                #ifndef XT_OMP_SKIP_REORGANIZE
                    const int64_t num_particles_to_track = ParticlesData_get__num_active_particles(particles);

                    {
                        LocalParticle lpart;
                        lpart.io_buffer = io_buffer;
                        Particles_to_LocalParticle(particles, &lpart, 0, capacity);
                        check_is_active(&lpart);
                        count_reorganized_particles(&lpart);
                        LocalParticle_to_Particles(&lpart, particles, 0, capacity);
                    }
                #else // When we skip reorganize, we cannot just batch active particles
                    const int64_t num_particles_to_track = capacity;
                #endif

                const int64_t chunk_size = (num_particles_to_track + num_threads - 1)/num_threads; // ceil division
            #endif // CONTEXT_OPENMP

//            #pragma omp parallel for                                                           //only_for_context cpu_openmp
//            for (int64_t batch_id = 0; batch_id < num_threads; batch_id++) {                   //only_for_context cpu_openmp
                LocalParticle lpart;
                lpart.io_buffer = io_buffer;
//                int64_t part_id = batch_id * chunk_size;                                       //only_for_context cpu_openmp
//                int64_t end_id = (batch_id + 1) * chunk_size;                                  //only_for_context cpu_openmp
//                if (end_id > num_particles_to_track) end_id = num_particles_to_track;          //only_for_context cpu_openmp

                int64_t part_id = 0;                    //only_for_context cpu_serial
//                int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
//                int64_t part_id = get_global_id(0);                    //only_for_context opencl
                int64_t end_id = 0; // unused outside of openmp  //only_for_context cpu_serial cuda opencl

                int64_t part_capacity = ParticlesData_get__capacity(particles);
                if (part_id<part_capacity){
                    Particles_to_LocalParticle(particles, &lpart, part_id, end_id);
                    if (check_is_active(&lpart)>0){
              Drift_track_local_particle_with_transformations(el, &lpart);

                    }
                    if (check_is_active(&lpart)>0 && flag_increment_at_element){
                            increment_at_element(&lpart, 1);
                    }
                }
//            } //only_for_context cpu_openmp

            // On OpenMP we want to additionally by default reorganize all
            // the particles.
//            #ifndef XT_OMP_SKIP_REORGANIZE                             //only_for_context cpu_openmp
//            LocalParticle lpart;                                       //only_for_context cpu_openmp
//            lpart.io_buffer = io_buffer;                               //only_for_context cpu_openmp
//            Particles_to_LocalParticle(particles, &lpart, 0, capacity);//only_for_context cpu_openmp
//            check_is_active(&lpart);                                   //only_for_context cpu_openmp
//            #endif                                                     //only_for_context cpu_openmp
        }

#ifndef XOBJ_TYPEDEF_XYShiftData
#define XOBJ_TYPEDEF_XYShiftData
typedef   struct XYShiftData_s * XYShiftData;
 static inline XYShiftData XYShiftData_getp(XYShiftData restrict  obj){
  int64_t offset=0;
  return (XYShiftData)(( char*) obj+offset);
}
 static inline double XYShiftData_get_dx(const XYShiftData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void XYShiftData_set_dx(XYShiftData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XYShiftData_getp_dx(XYShiftData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double XYShiftData_get_dy(const XYShiftData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void XYShiftData_set_dy(XYShiftData restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XYShiftData_getp_dy(XYShiftData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_XYSHIFT_H
#define XTRACK_XYSHIFT_H


 static inline
void XYShift_single_particle(LocalParticle* part, double dx, double dy){

    LocalParticle_add_to_x(part, -dx );
    LocalParticle_add_to_y(part, -dy );
 
}


 static inline
void XYShift_track_local_particle(XYShiftData el, LocalParticle* part0){

    double dx = XYShiftData_get_dx(el);
    double dy = XYShiftData_get_dy(el);

    #ifdef XSUITE_BACKTRACK
        dx = -dx;
        dy = -dy;
    #endif


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        XYShift_single_particle(part, dx, dy);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }


}


#endif /* XTRACK_XYSHIFT_H */


             static inline
            void XYShift_track_local_particle_with_transformations(XYShiftData el, LocalParticle* part0){
    XYShift_track_local_particle(el, part0);
}


             
            void XYShift_track_particles(
               XYShiftData el,

                             ParticlesData particles,

                             int64_t flag_increment_at_element,
                  int8_t* io_buffer){

//            #define CONTEXT_OPENMP  //only_for_context cpu_openmp
            #ifdef CONTEXT_OPENMP
                const int64_t capacity = ParticlesData_get__capacity(particles);
                const int num_threads = omp_get_max_threads();

                #ifndef XT_OMP_SKIP_REORGANIZE
                    const int64_t num_particles_to_track = ParticlesData_get__num_active_particles(particles);

                    {
                        LocalParticle lpart;
                        lpart.io_buffer = io_buffer;
                        Particles_to_LocalParticle(particles, &lpart, 0, capacity);
                        check_is_active(&lpart);
                        count_reorganized_particles(&lpart);
                        LocalParticle_to_Particles(&lpart, particles, 0, capacity);
                    }
                #else // When we skip reorganize, we cannot just batch active particles
                    const int64_t num_particles_to_track = capacity;
                #endif

                const int64_t chunk_size = (num_particles_to_track + num_threads - 1)/num_threads; // ceil division
            #endif // CONTEXT_OPENMP

//            #pragma omp parallel for                                                           //only_for_context cpu_openmp
//            for (int64_t batch_id = 0; batch_id < num_threads; batch_id++) {                   //only_for_context cpu_openmp
                LocalParticle lpart;
                lpart.io_buffer = io_buffer;
//                int64_t part_id = batch_id * chunk_size;                                       //only_for_context cpu_openmp
//                int64_t end_id = (batch_id + 1) * chunk_size;                                  //only_for_context cpu_openmp
//                if (end_id > num_particles_to_track) end_id = num_particles_to_track;          //only_for_context cpu_openmp

                int64_t part_id = 0;                    //only_for_context cpu_serial
//                int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
//                int64_t part_id = get_global_id(0);                    //only_for_context opencl
                int64_t end_id = 0; // unused outside of openmp  //only_for_context cpu_serial cuda opencl

                int64_t part_capacity = ParticlesData_get__capacity(particles);
                if (part_id<part_capacity){
                    Particles_to_LocalParticle(particles, &lpart, part_id, end_id);
                    if (check_is_active(&lpart)>0){
              XYShift_track_local_particle_with_transformations(el, &lpart);

                    }
                    if (check_is_active(&lpart)>0 && flag_increment_at_element){
                            increment_at_element(&lpart, 1);
                    }
                }
//            } //only_for_context cpu_openmp

            // On OpenMP we want to additionally by default reorganize all
            // the particles.
//            #ifndef XT_OMP_SKIP_REORGANIZE                             //only_for_context cpu_openmp
//            LocalParticle lpart;                                       //only_for_context cpu_openmp
//            lpart.io_buffer = io_buffer;                               //only_for_context cpu_openmp
//            Particles_to_LocalParticle(particles, &lpart, 0, capacity);//only_for_context cpu_openmp
//            check_is_active(&lpart);                                   //only_for_context cpu_openmp
//            #endif                                                     //only_for_context cpu_openmp
        }

#ifndef XOBJ_TYPEDEF_SRotationData
#define XOBJ_TYPEDEF_SRotationData
typedef   struct SRotationData_s * SRotationData;
 static inline SRotationData SRotationData_getp(SRotationData restrict  obj){
  int64_t offset=0;
  return (SRotationData)(( char*) obj+offset);
}
 static inline double SRotationData_get_cos_z(const SRotationData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void SRotationData_set_cos_z(SRotationData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* SRotationData_getp_cos_z(SRotationData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double SRotationData_get_sin_z(const SRotationData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void SRotationData_set_sin_z(SRotationData restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* SRotationData_getp_sin_z(SRotationData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_SROTATION_H
#define XTRACK_SROTATION_H

 static inline
void SRotation_track_local_particle(SRotationData el, LocalParticle* part0){

    double sin_z = SRotationData_get_sin_z(el);
    double cos_z = SRotationData_get_cos_z(el);

    #ifdef XSUITE_BACKTRACK
        sin_z = -sin_z;
    #endif


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        SRotation_single_particle(part, sin_z, cos_z);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }


}

#endif /* XTRACK_SROTATION_H */


             static inline
            void SRotation_track_local_particle_with_transformations(SRotationData el, LocalParticle* part0){
    SRotation_track_local_particle(el, part0);
}


             
            void SRotation_track_particles(
               SRotationData el,

                             ParticlesData particles,

                             int64_t flag_increment_at_element,
                  int8_t* io_buffer){

//            #define CONTEXT_OPENMP  //only_for_context cpu_openmp
            #ifdef CONTEXT_OPENMP
                const int64_t capacity = ParticlesData_get__capacity(particles);
                const int num_threads = omp_get_max_threads();

                #ifndef XT_OMP_SKIP_REORGANIZE
                    const int64_t num_particles_to_track = ParticlesData_get__num_active_particles(particles);

                    {
                        LocalParticle lpart;
                        lpart.io_buffer = io_buffer;
                        Particles_to_LocalParticle(particles, &lpart, 0, capacity);
                        check_is_active(&lpart);
                        count_reorganized_particles(&lpart);
                        LocalParticle_to_Particles(&lpart, particles, 0, capacity);
                    }
                #else // When we skip reorganize, we cannot just batch active particles
                    const int64_t num_particles_to_track = capacity;
                #endif

                const int64_t chunk_size = (num_particles_to_track + num_threads - 1)/num_threads; // ceil division
            #endif // CONTEXT_OPENMP

//            #pragma omp parallel for                                                           //only_for_context cpu_openmp
//            for (int64_t batch_id = 0; batch_id < num_threads; batch_id++) {                   //only_for_context cpu_openmp
                LocalParticle lpart;
                lpart.io_buffer = io_buffer;
//                int64_t part_id = batch_id * chunk_size;                                       //only_for_context cpu_openmp
//                int64_t end_id = (batch_id + 1) * chunk_size;                                  //only_for_context cpu_openmp
//                if (end_id > num_particles_to_track) end_id = num_particles_to_track;          //only_for_context cpu_openmp

                int64_t part_id = 0;                    //only_for_context cpu_serial
//                int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
//                int64_t part_id = get_global_id(0);                    //only_for_context opencl
                int64_t end_id = 0; // unused outside of openmp  //only_for_context cpu_serial cuda opencl

                int64_t part_capacity = ParticlesData_get__capacity(particles);
                if (part_id<part_capacity){
                    Particles_to_LocalParticle(particles, &lpart, part_id, end_id);
                    if (check_is_active(&lpart)>0){
              SRotation_track_local_particle_with_transformations(el, &lpart);

                    }
                    if (check_is_active(&lpart)>0 && flag_increment_at_element){
                            increment_at_element(&lpart, 1);
                    }
                }
//            } //only_for_context cpu_openmp

            // On OpenMP we want to additionally by default reorganize all
            // the particles.
//            #ifndef XT_OMP_SKIP_REORGANIZE                             //only_for_context cpu_openmp
//            LocalParticle lpart;                                       //only_for_context cpu_openmp
//            lpart.io_buffer = io_buffer;                               //only_for_context cpu_openmp
//            Particles_to_LocalParticle(particles, &lpart, 0, capacity);//only_for_context cpu_openmp
//            check_is_active(&lpart);                                   //only_for_context cpu_openmp
//            #endif                                                     //only_for_context cpu_openmp
        }

#ifndef XOBJ_TYPEDEF_YRotationData
#define XOBJ_TYPEDEF_YRotationData
typedef   struct YRotationData_s * YRotationData;
 static inline YRotationData YRotationData_getp(YRotationData restrict  obj){
  int64_t offset=0;
  return (YRotationData)(( char*) obj+offset);
}
 static inline double YRotationData_get_sin_angle(const YRotationData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void YRotationData_set_sin_angle(YRotationData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* YRotationData_getp_sin_angle(YRotationData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double YRotationData_get_cos_angle(const YRotationData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void YRotationData_set_cos_angle(YRotationData restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* YRotationData_getp_cos_angle(YRotationData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double YRotationData_get_tan_angle(const YRotationData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void YRotationData_set_tan_angle(YRotationData restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* YRotationData_getp_tan_angle(YRotationData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_TRACK_YROTATION_H
#define XTRACK_TRACK_YROTATION_H

 static inline
void YRotation_single_particle(LocalParticle* part, double sin_angle, double cos_angle, double tan_angle){

    double const beta0 = LocalParticle_get_beta0(part);
    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const t = LocalParticle_get_zeta(part)/beta0;
    double const pt = LocalParticle_get_pzeta(part)*beta0;

    double pz = sqrt(1.0 + 2.0*pt/beta0 + pt*pt - px*px - py*py);
    double ptt = 1.0 - tan_angle*px/pz;
    double x_hat = x/(cos_angle*ptt);
    double px_hat = cos_angle*px + sin_angle*pz;
    double y_hat = y + tan_angle*x*py/(pz*ptt);
    double t_hat = t - tan_angle*x*(1.0/beta0+pt)/(pz*ptt);

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_zeta(part,t_hat*beta0);

}

#endif /* XTRACK_TRACK_YROTATION_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_YROTATION_H
#define XTRACK_YROTATION_H

 static inline
void YRotation_track_local_particle(YRotationData el, LocalParticle* part0){

    double sin_angle = YRotationData_get_sin_angle(el);
    double cos_angle = YRotationData_get_cos_angle(el);
    double tan_angle = YRotationData_get_tan_angle(el);

    #ifdef XSUITE_BACKTRACK
        sin_angle = -sin_angle;
        tan_angle = -tan_angle;
    #endif


    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        YRotation_single_particle(part, sin_angle, cos_angle, tan_angle);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }


}

#endif /* XTRACK_YROTATION_H */


             static inline
            void YRotation_track_local_particle_with_transformations(YRotationData el, LocalParticle* part0){
    YRotation_track_local_particle(el, part0);
}


             
            void YRotation_track_particles(
               YRotationData el,

                             ParticlesData particles,

                             int64_t flag_increment_at_element,
                  int8_t* io_buffer){

//            #define CONTEXT_OPENMP  //only_for_context cpu_openmp
            #ifdef CONTEXT_OPENMP
                const int64_t capacity = ParticlesData_get__capacity(particles);
                const int num_threads = omp_get_max_threads();

                #ifndef XT_OMP_SKIP_REORGANIZE
                    const int64_t num_particles_to_track = ParticlesData_get__num_active_particles(particles);

                    {
                        LocalParticle lpart;
                        lpart.io_buffer = io_buffer;
                        Particles_to_LocalParticle(particles, &lpart, 0, capacity);
                        check_is_active(&lpart);
                        count_reorganized_particles(&lpart);
                        LocalParticle_to_Particles(&lpart, particles, 0, capacity);
                    }
                #else // When we skip reorganize, we cannot just batch active particles
                    const int64_t num_particles_to_track = capacity;
                #endif

                const int64_t chunk_size = (num_particles_to_track + num_threads - 1)/num_threads; // ceil division
            #endif // CONTEXT_OPENMP

//            #pragma omp parallel for                                                           //only_for_context cpu_openmp
//            for (int64_t batch_id = 0; batch_id < num_threads; batch_id++) {                   //only_for_context cpu_openmp
                LocalParticle lpart;
                lpart.io_buffer = io_buffer;
//                int64_t part_id = batch_id * chunk_size;                                       //only_for_context cpu_openmp
//                int64_t end_id = (batch_id + 1) * chunk_size;                                  //only_for_context cpu_openmp
//                if (end_id > num_particles_to_track) end_id = num_particles_to_track;          //only_for_context cpu_openmp

                int64_t part_id = 0;                    //only_for_context cpu_serial
//                int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
//                int64_t part_id = get_global_id(0);                    //only_for_context opencl
                int64_t end_id = 0; // unused outside of openmp  //only_for_context cpu_serial cuda opencl

                int64_t part_capacity = ParticlesData_get__capacity(particles);
                if (part_id<part_capacity){
                    Particles_to_LocalParticle(particles, &lpart, part_id, end_id);
                    if (check_is_active(&lpart)>0){
              YRotation_track_local_particle_with_transformations(el, &lpart);

                    }
                    if (check_is_active(&lpart)>0 && flag_increment_at_element){
                            increment_at_element(&lpart, 1);
                    }
                }
//            } //only_for_context cpu_openmp

            // On OpenMP we want to additionally by default reorganize all
            // the particles.
//            #ifndef XT_OMP_SKIP_REORGANIZE                             //only_for_context cpu_openmp
//            LocalParticle lpart;                                       //only_for_context cpu_openmp
//            lpart.io_buffer = io_buffer;                               //only_for_context cpu_openmp
//            Particles_to_LocalParticle(particles, &lpart, 0, capacity);//only_for_context cpu_openmp
//            check_is_active(&lpart);                                   //only_for_context cpu_openmp
//            #endif                                                     //only_for_context cpu_openmp
        }

#ifndef XOBJ_TYPEDEF_InteractionRecordData
#define XOBJ_TYPEDEF_InteractionRecordData
typedef   struct InteractionRecordData_s * InteractionRecordData;
 static inline InteractionRecordData InteractionRecordData_getp(InteractionRecordData restrict  obj){
  int64_t offset=0;
  return (InteractionRecordData)(( char*) obj+offset);
}
 static inline RecordIndex InteractionRecordData_getp__index(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return (RecordIndex)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_get__index_capacity(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__index_capacity(InteractionRecordData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp__index_capacity(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline uint32_t InteractionRecordData_get__index_num_recorded(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__index_num_recorded(InteractionRecordData restrict  obj, uint32_t value){
  int64_t offset=0;
  offset+=16;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* InteractionRecordData_getp__index_num_recorded(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline uint32_t InteractionRecordData_get__index__dummy(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( uint32_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__index__dummy(InteractionRecordData restrict  obj, uint32_t value){
  int64_t offset=0;
  offset+=24;
  *( uint32_t*)(( char*) obj+offset)=value;
}
 static inline  uint32_t* InteractionRecordData_getp__index__dummy(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( uint32_t*)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_get__index_buffer_id(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__index_buffer_id(InteractionRecordData restrict  obj, int64_t value){
  int64_t offset=0;
  offset+=32;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp__index_buffer_id(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_at_element(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=320;
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_at_element(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=320;
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_at_element(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=320;
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_at_element(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=320;
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_at_element(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=320;
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_at_turn(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_at_turn(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_at_turn(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_at_turn(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_at_turn(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+80);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp__inter(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len__inter(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get__inter(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__inter(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1__inter(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+88);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_id_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_id_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_id_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_id_before(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_id_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+96);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_s_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_s_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_s_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_s_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_s_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+104);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_x_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_x_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_x_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_x_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_x_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+112);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_px_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_px_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_px_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_px_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_px_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+120);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_y_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_y_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_y_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_y_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_y_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+128);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_py_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_py_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_py_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_py_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_py_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+136);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_zeta_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_zeta_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_zeta_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_zeta_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_zeta_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+144);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_delta_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_delta_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_delta_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_delta_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_delta_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+152);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_energy_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_energy_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_energy_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_energy_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_energy_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+160);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_mass_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_mass_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_mass_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_mass_before(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_mass_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+168);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_charge_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_charge_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_charge_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_charge_before(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_charge_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+176);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_z_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_z_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_z_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_z_before(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_z_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+184);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_a_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_a_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_a_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_a_before(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_a_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+192);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_pdgid_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_pdgid_before(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_pdgid_before(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_pdgid_before(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_pdgid_before(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+200);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_id_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_id_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_id_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_id_after(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_id_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+208);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_s_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_s_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_s_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_s_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_s_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+216);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_x_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_x_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_x_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_x_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_x_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+224);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_px_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_px_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_px_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_px_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_px_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+232);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_y_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_y_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_y_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_y_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_y_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+240);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_py_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_py_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_py_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_py_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_py_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+248);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_zeta_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_zeta_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_zeta_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_zeta_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_zeta_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+256);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_delta_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_delta_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_delta_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_delta_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_delta_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+264);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_energy_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_energy_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_energy_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_energy_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_energy_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+272);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNFloat64 InteractionRecordData_getp_mass_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+280);
  return (ArrNFloat64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_mass_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+280);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline double InteractionRecordData_get_mass_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+280);
  offset+=16+i0*8;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_mass_after(InteractionRecordData restrict  obj, int64_t i0, double value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+280);
  offset+=16+i0*8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp1_mass_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+280);
  offset+=16+i0*8;
  return ( double*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_charge_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+288);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_charge_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+288);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_charge_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+288);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_charge_after(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+288);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_charge_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+288);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_z_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+296);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_z_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+296);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_z_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+296);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_z_after(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+296);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_z_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+296);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_a_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+304);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_a_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+304);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_a_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+304);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_a_after(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+304);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_a_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+304);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline ArrNInt64 InteractionRecordData_getp_pdgid_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+312);
  return (ArrNInt64)(( char*) obj+offset);
}
 static inline int64_t InteractionRecordData_len_pdgid_after(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+312);
   int64_t* arr = ( int64_t*)(( char*) obj+offset);
  return arr[1];
}
 static inline int64_t InteractionRecordData_get_pdgid_after(const InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+312);
  offset+=16+i0*8;
  return *( int64_t*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set_pdgid_after(InteractionRecordData restrict  obj, int64_t i0, int64_t value){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+312);
  offset+=16+i0*8;
  *( int64_t*)(( char*) obj+offset)=value;
}
 static inline  int64_t* InteractionRecordData_getp1_pdgid_after(InteractionRecordData restrict  obj, int64_t i0){
  int64_t offset=0;
  offset+=*( int64_t*)(( char*) obj+offset+312);
  offset+=16+i0*8;
  return ( int64_t*)(( char*) obj+offset);
}
 static inline double InteractionRecordData_get__sin_rot_s(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=40;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__sin_rot_s(InteractionRecordData restrict  obj, double value){
  int64_t offset=0;
  offset+=40;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp__sin_rot_s(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=40;
  return ( double*)(( char*) obj+offset);
}
 static inline double InteractionRecordData_get__cos_rot_s(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=48;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__cos_rot_s(InteractionRecordData restrict  obj, double value){
  int64_t offset=0;
  offset+=48;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp__cos_rot_s(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=48;
  return ( double*)(( char*) obj+offset);
}
 static inline double InteractionRecordData_get__shift_x(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=56;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__shift_x(InteractionRecordData restrict  obj, double value){
  int64_t offset=0;
  offset+=56;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp__shift_x(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=56;
  return ( double*)(( char*) obj+offset);
}
 static inline double InteractionRecordData_get__shift_y(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=64;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__shift_y(InteractionRecordData restrict  obj, double value){
  int64_t offset=0;
  offset+=64;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp__shift_y(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=64;
  return ( double*)(( char*) obj+offset);
}
 static inline double InteractionRecordData_get__shift_s(const InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=72;
  return *( double*)(( char*) obj+offset);
}
 static inline void InteractionRecordData_set__shift_s(InteractionRecordData restrict  obj, double value){
  int64_t offset=0;
  offset+=72;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* InteractionRecordData_getp__shift_s(InteractionRecordData restrict  obj){
  int64_t offset=0;
  offset+=72;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */


#ifndef XCOLL_INTERACTIONS_H
#define XCOLL_INTERACTIONS_H

#define  XC_UNITIALISED                     0        // NAN  // Do not use

#define  XC_ENTER_JAW_L                    -1        // JI   // point (no children)    Set ds > 0 if entering later
#define  XC_ENTER_JAW_R                    -2        // JI   // point (no children)    Set ds > 0 if entering later
#define  XC_EXIT_JAW                       -3        // JO   // point (no children)
#define  XC_ENTER_JAW                      -4        // JI   // point (no children)    Set ds > 0 if entering later    still here for compatibility
#define  XC_ABSORBED                        1        // A    // point (no children)    Don't use 0 (is default for unitialised)
#define  XC_MULTIPLE_COULOMB_SCATTERING    13        // MCS  // continuous
#define  XC_PN_ELASTIC                     14        // PN   // point (no children)
#define  XC_PP_ELASTIC                     15        // PP   // point (no children)
#define  XC_SINGLE_DIFFRACTIVE             16        // SD   // point (no children)
#define  XC_COULOMB                        17        // C    // point (no children)
#define  XC_CHANNELING                    100        // CH   // continuous
#define  XC_DECHANNELING                  101        // DCH  // point (no children)
#define  XC_VOLUME_REFLECTION_TRANS_CH    102        // VRCH // point (no children)    Transition region around +-xpcrit
#define  XC_VOLUME_REFLECTION             103        // VR   // point (no children)
#define  XC_VOLUME_REFLECTION_TRANS_MCS   104        // VRAM // point (no children)    Transition region around t_P
#define  XC_MULTIPLE_COULOMB_TRANS_VR     105        // AMVR // continuous             Transition region around t_P + 2 xpcrit
#define  XC_VOLUME_CAPTURE                106        // VC   // point (no children)

#endif /* XCOLL_INTERACTIONS_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_STATES_H
#define XCOLL_STATES_H

#define  XC_LOST_ON_EVEREST_BLOCK   -330
#define  XC_LOST_ON_EVEREST_COLL    -331
#define  XC_LOST_ON_EVEREST_CRYSTAL -332
#define  XC_LOST_ON_FLUKA_BLOCK     -333
#define  XC_LOST_ON_FLUKA           -334
#define  XC_LOST_ON_FLUKA_CRYSTAL   -335
#define  XC_LOST_ON_GEANT4_BLOCK    -336
#define  XC_LOST_ON_GEANT4          -337
#define  XC_LOST_ON_GEANT4_CRYSTAL  -338
#define  XC_LOST_ON_ABSORBER        -340

#define  XC_ERR_INVALID_TRACK       -390
#define  XC_ERR_NOT_IMPLEMENTED     -391
#define  XC_ERR_INVALID_XOFIELD     -392
#define  XC_ERR                     -399

#endif /* XCOLL_STATES_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #


#ifndef XCOLL_IMPACTS_H
#define XCOLL_IMPACTS_H

// TODO: do we need to pass RecordIndex?
// probably can do RecordIndex record_index = InteractionRecordData_getp__index(record);  ?
 static inline
int64_t InteractionRecordData_log(InteractionRecordData record, RecordIndex record_index, LocalParticle* parent,
                                  int64_t interaction){
    // This can be used for a point-like interaction where there is no child (or because it's equal to the parent)
    // or to log the parent first, to be followed up with InteractionRecordData_log_child on the same slot

    int64_t i_slot = -1;
    if (record){
        // Get a slot in the record (this is thread safe)
        i_slot = RecordIndex_get_slot(record_index);
        // The returned slot id is negative if record is NULL or if record is full

        if (i_slot>=0){
            InteractionRecordData_set_at_element(record, i_slot, LocalParticle_get_at_element(parent));
            InteractionRecordData_set_at_turn(record, i_slot, LocalParticle_get_at_turn(parent));
            InteractionRecordData_set__inter(record, i_slot, interaction);

            double charge_ratio = LocalParticle_get_charge_ratio(parent);
            double mass_ratio = charge_ratio / LocalParticle_get_chi(parent);
            double energy = ( LocalParticle_get_ptau(parent) + 1 / LocalParticle_get_beta0(parent)
                             ) * mass_ratio * LocalParticle_get_p0c(parent);
            // All fields have to be written, or the arrays will not have the same length
            // TODO: maybe this is not true, as we are setting by slot index? Don't the arrays come pre-initialised?
            InteractionRecordData_set_id_before(record, i_slot, LocalParticle_get_particle_id(parent));
            InteractionRecordData_set_s_before(record,  i_slot, LocalParticle_get_s(parent));
            InteractionRecordData_set_x_before(record,  i_slot, LocalParticle_get_x(parent));
            InteractionRecordData_set_px_before(record, i_slot, LocalParticle_get_px(parent));
            InteractionRecordData_set_y_before(record,  i_slot, LocalParticle_get_y(parent));
            InteractionRecordData_set_py_before(record, i_slot, LocalParticle_get_py(parent));
            InteractionRecordData_set_zeta_before(record,   i_slot, LocalParticle_get_zeta(parent));
            InteractionRecordData_set_delta_before(record,  i_slot, LocalParticle_get_delta(parent));
            InteractionRecordData_set_energy_before(record, i_slot, energy);
            InteractionRecordData_set_mass_before(record,   i_slot, mass_ratio*LocalParticle_get_mass0(parent));
            InteractionRecordData_set_charge_before(record, i_slot, charge_ratio*LocalParticle_get_q0(parent));
            // TODO: particle info
            InteractionRecordData_set_z_before(record, i_slot, -1);
            InteractionRecordData_set_a_before(record, i_slot, -1);
            InteractionRecordData_set_pdgid_before(record, i_slot, -1);

            // TODO: maybe this is not needed
            InteractionRecordData_set_id_after(record, i_slot, -1);
            InteractionRecordData_set_s_after(record, i_slot, -1);
            InteractionRecordData_set_x_after(record, i_slot, -1);
            InteractionRecordData_set_px_after(record, i_slot, -1);
            InteractionRecordData_set_y_after(record, i_slot, -1);
            InteractionRecordData_set_py_after(record, i_slot, -1);
            InteractionRecordData_set_zeta_after(record, i_slot, -1);
            InteractionRecordData_set_delta_after(record, i_slot, -1);
            InteractionRecordData_set_energy_after(record, i_slot, -1);
            InteractionRecordData_set_mass_after(record, i_slot, -1);
            InteractionRecordData_set_charge_after(record, i_slot, -1);
            InteractionRecordData_set_z_after(record, i_slot, -1);
            InteractionRecordData_set_a_after(record, i_slot, -1);
            InteractionRecordData_set_pdgid_after(record, i_slot, -1);
        }
    }
//     printf("Logging %i in slot %i\n", interaction, i_slot);
    return i_slot;
}

 static inline
void InteractionRecordData_log_child(InteractionRecordData record, int64_t i_slot, LocalParticle* child){
    if (record && i_slot>=0){
        double charge_ratio = LocalParticle_get_charge_ratio(child);
        double mass_ratio = charge_ratio / LocalParticle_get_chi(child);
        double energy = ( LocalParticle_get_ptau(child) + 1 / LocalParticle_get_beta0(child)
                         ) * mass_ratio * LocalParticle_get_p0c(child);
        InteractionRecordData_set_id_after(record, i_slot, LocalParticle_get_particle_id(child));
        InteractionRecordData_set_s_after(record,  i_slot, LocalParticle_get_s(child));
        InteractionRecordData_set_x_after(record,  i_slot, LocalParticle_get_x(child));
        InteractionRecordData_set_px_after(record, i_slot, LocalParticle_get_px(child));
        InteractionRecordData_set_y_after(record,  i_slot, LocalParticle_get_y(child));
        InteractionRecordData_set_py_after(record, i_slot, LocalParticle_get_py(child));
        InteractionRecordData_set_zeta_after(record,  i_slot, LocalParticle_get_zeta(child));
        InteractionRecordData_set_delta_after(record, i_slot, LocalParticle_get_delta(child));
        InteractionRecordData_set_energy_after(record, i_slot, energy);
        InteractionRecordData_set_mass_after(record,   i_slot, mass_ratio*LocalParticle_get_mass0(child));
        InteractionRecordData_set_charge_after(record, i_slot, charge_ratio*LocalParticle_get_q0(child));
        // TODO: particle info
        InteractionRecordData_set_z_after(record, i_slot, -1);
        InteractionRecordData_set_a_after(record, i_slot, -1);
        InteractionRecordData_set_pdgid_after(record, i_slot, -1);
//     printf("Slot %i: length %f\n", i_slot, ds);
    }
}

#endif /* XCOLL_IMPACTS_H */

#ifndef XOBJ_TYPEDEF_XcollGeometryData
#define XOBJ_TYPEDEF_XcollGeometryData
typedef   struct XcollGeometryData_s * XcollGeometryData;
 static inline XcollGeometryData XcollGeometryData_getp(XcollGeometryData restrict  obj){
  int64_t offset=0;
  return (XcollGeometryData)(( char*) obj+offset);
}
 static inline double XcollGeometryData_get__sin_rot_s(const XcollGeometryData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryData_set__sin_rot_s(XcollGeometryData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryData_getp__sin_rot_s(XcollGeometryData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryData_get__cos_rot_s(const XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryData_set__cos_rot_s(XcollGeometryData restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryData_getp__cos_rot_s(XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryData_get__shift_x(const XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryData_set__shift_x(XcollGeometryData restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryData_getp__shift_x(XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryData_get__shift_y(const XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryData_set__shift_y(XcollGeometryData restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryData_getp__shift_y(XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryData_get__shift_s(const XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryData_set__shift_s(XcollGeometryData restrict  obj, double value){
  int64_t offset=0;
  offset+=32;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryData_getp__shift_s(XcollGeometryData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SORT_H
#define XCOLL_GEOM_SORT_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"

#ifdef MAX
#undef MAX
#pragma message ("Xcoll geometry: Compiler macro MAX redefined")
#endif
#define MAX(x, y) ({const __typeof__ (x) _x = (x); \
                    const __typeof__ (y) _y = (y); \
                    _x > _y ? _x : _y; })
#ifdef MIN
#undef MIN
#pragma message ("Xcoll geometry: Compiler macro MIN redefined")
#endif
#define MIN(x, y) ({const __typeof__ (x) _x = (x); \
                    const __typeof__ (y) _y = (y); \
                    _x < _y ? _x : _y; })
#ifdef SWAP
#error "Xcoll geometry: Compiler macro SWAP already defined!"
#endif
#define SWAP(d,x,y) ({const __typeof__(*d) _x = MIN(d[x], d[y]); \
                      const __typeof__(*d) _y = MAX(d[x], d[y]); \
                      d[x] = _x; d[y] = _y; })


// Fast methods
// ------------

static inline void sort_array_of_2_double(double* d){
    SWAP(d, 0, 1);
}

static inline void sort_array_of_3_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 1, 2); SWAP(d, 0, 1);
}

static inline void sort_array_of_4_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_5_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 3, 4); SWAP(d, 2, 4); SWAP(d, 2, 3); SWAP(d, 1, 4);
    SWAP(d, 0, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_6_double(double* d){
    SWAP(d, 1, 2); SWAP(d, 4, 5); SWAP(d, 0, 2); SWAP(d, 3, 5); SWAP(d, 0, 1);
    SWAP(d, 3, 4); SWAP(d, 1, 4); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3);
    SWAP(d, 2, 4); SWAP(d, 2, 3);
}

static inline void sort_array_of_7_double(double* d){
    SWAP(d, 1, 2); SWAP(d, 3, 4); SWAP(d, 5, 6); SWAP(d, 0, 2); SWAP(d, 3, 5);
    SWAP(d, 4, 6); SWAP(d, 0, 1); SWAP(d, 4, 5); SWAP(d, 2, 6); SWAP(d, 0, 4);
    SWAP(d, 1, 5); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3); SWAP(d, 2, 4);
    SWAP(d, 2, 3);
}

static inline void sort_array_of_8_double(double* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 4, 5); SWAP(d, 6, 7); SWAP(d, 0, 2);
    SWAP(d, 1, 3); SWAP(d, 4, 6); SWAP(d, 5, 7); SWAP(d, 1, 2); SWAP(d, 5, 6);
    SWAP(d, 0, 4); SWAP(d, 3, 7); SWAP(d, 1, 5); SWAP(d, 2, 6); SWAP(d, 1, 4);
    SWAP(d, 3, 6); SWAP(d, 2, 4); SWAP(d, 3, 5); SWAP(d, 3, 4);
}

static inline void sort_array_of_2_int64(int64_t* d){
    SWAP(d, 0, 1);
}

static inline void sort_array_of_3_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 1, 2); SWAP(d, 0, 1);
}

static inline void sort_array_of_4_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_5_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 3, 4); SWAP(d, 2, 4); SWAP(d, 2, 3); SWAP(d, 1, 4);
    SWAP(d, 0, 3); SWAP(d, 0, 2); SWAP(d, 1, 3); SWAP(d, 1, 2);
}

static inline void sort_array_of_6_int64(int64_t* d){
    SWAP(d, 1, 2); SWAP(d, 4, 5); SWAP(d, 0, 2); SWAP(d, 3, 5); SWAP(d, 0, 1);
    SWAP(d, 3, 4); SWAP(d, 1, 4); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3);
    SWAP(d, 2, 4); SWAP(d, 2, 3);
}

static inline void sort_array_of_7_int64(int64_t* d){
    SWAP(d, 1, 2); SWAP(d, 3, 4); SWAP(d, 5, 6); SWAP(d, 0, 2); SWAP(d, 3, 5);
    SWAP(d, 4, 6); SWAP(d, 0, 1); SWAP(d, 4, 5); SWAP(d, 2, 6); SWAP(d, 0, 4);
    SWAP(d, 1, 5); SWAP(d, 0, 3); SWAP(d, 2, 5); SWAP(d, 1, 3); SWAP(d, 2, 4);
    SWAP(d, 2, 3);
}

static inline void sort_array_of_8_int64(int64_t* d){
    SWAP(d, 0, 1); SWAP(d, 2, 3); SWAP(d, 4, 5); SWAP(d, 6, 7); SWAP(d, 0, 2);
    SWAP(d, 1, 3); SWAP(d, 4, 6); SWAP(d, 5, 7); SWAP(d, 1, 2); SWAP(d, 5, 6);
    SWAP(d, 0, 4); SWAP(d, 3, 7); SWAP(d, 1, 5); SWAP(d, 2, 6); SWAP(d, 1, 4);
    SWAP(d, 3, 6); SWAP(d, 2, 4); SWAP(d, 3, 5); SWAP(d, 3, 4);
}


// Generic methods
// ---------------

int cmpfunc_double(const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
}

static inline void sort_array_of_double(double* arr, int64_t length){
    switch(length){
        case 2:
            sort_array_of_2_double(arr);
            break;
        case 3:
            sort_array_of_3_double(arr);
            break;
        case 4:
            sort_array_of_4_double(arr);
            break;
        case 5:
            sort_array_of_5_double(arr);
            break;
        case 6:
            sort_array_of_6_double(arr);
            break;
        case 7:
            sort_array_of_7_double(arr);
            break;
        case 8:
            sort_array_of_8_double(arr);
            break;
        default:
            qsort(arr, length, sizeof(double), cmpfunc_double);
    }
}

int cmpfunc_int64(const void * a, const void * b) {
   return ( *(int64_t*)a - *(int64_t*)b );
}

static inline void sort_array_of_int64(int64_t* arr, int64_t length){
    switch(length){
        case 2:
            sort_array_of_2_int64(arr);
            break;
        case 3:
            sort_array_of_3_int64(arr);
            break;
        case 4:
            sort_array_of_4_int64(arr);
            break;
        case 5:
            sort_array_of_5_int64(arr);
            break;
        case 6:
            sort_array_of_6_int64(arr);
            break;
        case 7:
            sort_array_of_7_int64(arr);
            break;
        case 8:
            sort_array_of_8_int64(arr);
            break;
        default:
            qsort(arr, length, sizeof(int64_t), cmpfunc_int64);
    }
}

#pragma GCC diagnostic pop
#endif /* XCOLL_GEOM_SORT_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_METHODS_H
#define XCOLL_GEOM_METHODS_H

#define XC_EPSILON 1.e-12
#define XC_S_MAX 1.e21


// This function calculates the overlap between an array and a given interval.
// The array comes in pairs of points, e.g. in-out-in-out... or out-in-out-in...
// IMPORTANT:
// The array and interval are assumed to be sorted!
// Furthermore, the array should have one extra slot allocated at the end, in case it needs to be expanded..
// This is always true for the arrays created by get_s, as we create them with 2*n_segments slots.
 static inline
void calculate_overlap_array_interval(double* arr, int8_t* length, double* interval){
    if (arr[0] > interval[1]){
        // No overlap
        *length = 0;
    }
    if ((*length)%2 == 1){
        // Special case: last interval of array is open until infinity,
        // so we add an extra point at the end to represent infinity.
        arr[*length] = XC_S_MAX*0.1;
        (*length)++;
    } else if (arr[*length-1] < interval[0]){
        // No overlap
        *length = 0;
    }
    int8_t i_start = 0;
    // Find the start of overlap
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[0]){
            if (i%2 == 0){
                // This is the first point of overlap
                i_start = i;
            } else {
                // The vertical restriction is the first point of overlap
                i_start = i-1;
                arr[i_start] = interval[0];
            }
            break;
        }
    }
    // Find the end of overlap
    int8_t i_stop = *length - 1;
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[1]){
            if (i%2 == 0){
                // The previous point is the last point of overlap
                i_stop = i-1;
            } else {
                // The vertical restriction is the first point of overlap
                i_stop = i;
                arr[i_stop] = interval[1];
            }
            break;
        }
    }
    *length = i_stop - i_start + 1;
    if (i_start > 0){
        for (int8_t i=0; i<*length; i++){
            arr[i] = arr[i+i_start];
        }
    }
}


#endif /* XCOLL_GEOM_METHODS_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEGMENTS_H
#define XCOLL_GEOM_SEGMENTS_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


// Parent type for all segments
// ----------------------------
typedef struct Segment_{
    int id;
} Segment_;
typedef Segment_* Segment;


// Line segment between two points (s1, x1) -- (s2, x2)
// -------------------------------------------------

#define XC_LINESEGMENT_ID 0
typedef struct LineSegment_ {
    int id;
    double point1_s;
    double point1_x;
    double point2_s;
    double point2_x;
} LineSegment_;
typedef LineSegment_* LineSegment;

 static inline
LineSegment create_line_segment(double point1_s, double point1_x, double point2_s, double point2_x){
    LineSegment seg = (LineSegment) malloc(sizeof(LineSegment_));
    seg->id = XC_LINESEGMENT_ID;
    seg->point1_s = point1_s;
    seg->point2_s = point2_s;
    seg->point1_x = point1_x;
    seg->point2_x = point2_x;
    return seg;
}


// Half-open line segment from a point (s, x) to +/-infinity along a slope
// -----------------------------------------------------------------------

#define XC_HALFOPENLINESEGMENT_ID 1
// For practical reasons, we store the inverse of the slope (inv_slope = -1/slope)
// of the segment. So inv_slope = 0 means the segment is vertical, inv_slope = 1
// implies a segment at an angle of 135 deg (tilt of +45 degrees), and a horizontal
// half-open segment is not allowed (and also not needed as it would go to infinity
// along the beam).
typedef struct HalfOpenLineSegment_ {
    int id;
    double point_s;
    double point_x;
    double inv_slope;
    int8_t sign;  // Does the segment go to +inf or -inf?
} HalfOpenLineSegment_;
typedef HalfOpenLineSegment_* HalfOpenLineSegment;

 static inline
HalfOpenLineSegment create_halfopen_line_segment(double point_s, double point_x, double inv_slope, int8_t sign){
    HalfOpenLineSegment seg = (HalfOpenLineSegment) malloc(sizeof(HalfOpenLineSegment_));
    seg->id = XC_HALFOPENLINESEGMENT_ID;
    seg->point_s = point_s;
    seg->point_x = point_x;
    seg->inv_slope = inv_slope;
    seg->sign = sign;
    return seg;
}


// Circular arc segment, defined by a centre and radius, and the starting/end angles
// ---------------------------------------------------------------------------------

#define XC_CIRCULARSEGMENT_ID 2
typedef struct CircularSegment_ {
    int id;
    double R;
    double centre_s;
    double centre_x;
    double point1_angle;
    double point2_angle;
} CircularSegment_;
typedef CircularSegment_* CircularSegment;

 static inline
CircularSegment create_circular_segment(double R, double centre_s, double centre_x, double point1_angle, double point2_angle){
    CircularSegment seg = (CircularSegment) malloc(sizeof(CircularSegment_));
    seg->id = XC_CIRCULARSEGMENT_ID;
    seg->R = R;
    seg->centre_s = centre_s;
    seg->centre_x = centre_x;
    seg->point1_angle = point1_angle;
    seg->point2_angle = point2_angle;
    return seg;
}


// Bézier segment, defined by a start and end point P1 and P2, and two control points that define the curve
// --------------------------------------------------------------------------------------------------------

#define XC_BEZIERSEGMENT_ID 3
typedef struct BezierSegment_ {
    int id;
    double point1_s;
    double point1_x;
    double control_point1_s;
    double control_point1_x;
    double point2_s;
    double point2_x;
    double control_point2_s;
    double control_point2_x;
} BezierSegment_;
typedef BezierSegment_* BezierSegment;

 static inline
BezierSegment create_bezier_segment(double point1_s, double point1_x, double control_point1_s, \
                                    double control_point1_x, double point2_s, double point2_x, \
                                    double control_point2_s, double control_point2_x){
    BezierSegment seg = (BezierSegment) malloc(sizeof(BezierSegment_));
    seg->id = XC_BEZIERSEGMENT_ID;
    seg->point1_s = point1_s;
    seg->point1_x = point1_x;
    seg->control_point1_s = control_point1_s;
    seg->control_point1_x = control_point1_x;
    seg->point2_s = point2_s;
    seg->point2_x = point2_x;
    seg->control_point2_s = control_point2_s;
    seg->control_point2_x = control_point2_x;
    return seg;
}


#endif /* XCOLL_GEOM_SEGMENTS_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CROSSING_DRIFT_H
#define XCOLL_GEOM_CROSSING_DRIFT_H


// All segments have a function that calculates the s-coordinate of the crossing(s) with a
// particle trajectory in a drift (given by a line through (part_s, part_x) and a slope part_tan).
// The results are always stored in an array s, and n_hit keeps track of the number of hits over
// multiple segments.


// Line Segments
// -------------

 static inline
void _crossing_drift_line(int8_t* n_hit, double* s, double s1, double s2, double x1, double x2,
                          double part_s, double part_x, double part_tan){
    // We fill in the segment points in the particle trajectory equation; if the results have opposite sign,
    //the two points lie on different sides of the trajectory and hence the segment is crossed.
    double trajectory1 = x1 - part_x - (s1 - part_s)*part_tan;
    double trajectory2 = x2 - part_x - (s2 - part_s)*part_tan;
    if (trajectory1*trajectory2 <= 0){
        // It's a crossing
        if (fabs(s2 - s1) < XC_EPSILON){
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            double poly_tan = (x2 - x1)/(s2 - s1);
            if (fabs(poly_tan - part_tan) < XC_EPSILON){
                // The trajectory is parallel to the segment.
                // TODO: this is situational; we should return s1 if get_s_first and current_s if after current_s
                s[*n_hit] = s1;
                (*n_hit)++;
            } else {
                // Normal crossing of two lines
                s[*n_hit] = (part_x - x1 + s1*poly_tan - part_s*part_tan)/(poly_tan - part_tan);
                (*n_hit)++;
            }
        }
    }
}

 static inline
void crossing_drift_line(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    LineSegment seg = (LineSegment) segment;
    double s1 = seg->point1_s;
    double s2 = seg->point2_s;
    double x1 = seg->point1_x;
    double x2 = seg->point2_x;
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}

 static inline
void crossing_drift_halfopenline(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    HalfOpenLineSegment seg = (HalfOpenLineSegment) segment;
    double s1 = seg->point_s;
    double x1 = seg->point_x;
    // A half-open segment implies one of its points lies at +-inf.
    // In practice we just add a polygon point at the wall overflow (at 1000km for the x-coordinate).
    double x2 = 1.e6*seg->sign;
    double s2;
    if (fabs(seg->inv_slope) < XC_EPSILON){
        s2 = s1;
    } else {
        s2 = s1 - (x2 - x1)*seg->inv_slope;
    }
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}


// Circular Segment
// ----------------

 static inline
void crossing_drift_circular(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    CircularSegment seg = (CircularSegment) segment;
    double R = seg->R;
    double sC = seg->centre_s;
    double xC = seg->centre_x;
    double t1 = seg->point1_angle;
    double t2 = seg->point2_angle;
    // Calculate crossings
    int8_t reversed = 0;
    if (t2 < t1){
        reversed = 1;
    }
    double a = 1 + part_tan*part_tan;
    double bb = sC - part_tan*(part_x - xC - part_tan*part_s); // This is -b/2 with b from the quadratic formula
    double c = sC*sC + (part_x - xC - part_tan*part_s)*(part_x - xC - part_tan*part_s) - R*R;
    double disc = bb*bb - a*c; // This is  2*discriminant**2
    if (disc < 0){
        // No crossing
        return;
    }
    for (int8_t i = 0; i < 2; i++) {
        double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
        double new_s = (bb + sgnD*sqrt(fabs(disc)))/a;
        double new_x = part_x + (new_s - part_s)*part_tan;
        double t = atan2(new_x - xC, new_s - sC);
        if (reversed){
            // t2 < t1, so we are looking at the inverted region of angles
            if (t1 >= t || t >= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        } else {
            if (t1 <= t && t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        }
    }
}


// Bézier Segment
// --------------

 static inline
void _hit_s_bezier(BezierSegment seg, double t, double multiplicity, int8_t* n_hit, double* s){
    double s1 = seg->point1_s;
    double cs1 = seg->control_point1_s;
    double s2 = seg->point2_s;
    double cs2 = seg->control_point2_s;
    double new_s = (1-t)*(1-t)*(1-t)*s1 + 3*(1-t)*(1-t)*t*cs1 + 3*(1-t)*t*t*cs2 + t*t*t*s2;
    for (int8_t i = 0; i < multiplicity; i++) {
        s[*n_hit] = new_s;
        (*n_hit)++;
    }
}

 static inline
void crossing_drift_bezier(void* segment, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    BezierSegment seg = (BezierSegment) segment;
    double s1 = seg->point1_s;
    double x1 = seg->point1_x;
    double cs1 = seg->control_point1_s;
    double cx1 = seg->control_point1_x;
    double s2 = seg->point2_s;
    double x2 = seg->point2_x;
    double cs2 = seg->control_point2_s;
    double cx2 = seg->control_point2_x;
    // The Bézier curve is defined by the parametric equations (with t in [0, 1]):
    // s(t) = (1-t)^3*s1 + 3(1-t)^2*t*cs1 + 3(1-t)*t^2*cs2 + t^3*s2
    // x(t) = (1-t)^3*x1 + 3(1-t)^2*t*cx1 + 3(1-t)*t^2*cx2 + t^3*x2
    double s0 = part_s;
    double x0 = part_x;
    double m = part_tan;
    // Plug the parametric eqs into the drift trajectory x(t) = m*(s(t) - s0) + x0 and solve for t
    // The solutions for t (which we get by Cardano's method) are valid if in [0, 1]
    double a = (m*s1 - x1) - (m*s2 - x2) - 3*(m*cs1 - cx1) + 3*(m*cs2 - cx2);
    double b = 6*(m*cs1 - cx1) - 3*(m*cs2 - cx2) - 3*(m*s1 - x1);
    double c = 3*(m*s1 - x1) - 3*(m*cs1 - cx1);
    double d = (m*s0 - x0) - (m*s1 - x1);
    double t;
    // Edge cases
    if (fabs(a) < XC_EPSILON){
        if (fabs(b) < XC_EPSILON){
            if (fabs(c) < XC_EPSILON){
                if (fabs(d) < XC_EPSILON){
                    // The trajectory is on the Bézier curve
                    // TODO: This cannot happen because we don't have a cubic trajectory.
                    //       Adapt if these ever would be implemented.
                    return;
                } else {
                    // No solutions
                    return;
                }
            } else {
                // This is a linear equation
                t = -d/c;
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        } else {
            // This is a quadratic equation
            double disc = c*c - 4*b*d;
            if (disc < 0){
                // No solutions
                return;
            }
            for (int8_t i = 0; i < 2; i++) {
                double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
                t = (-c + sgnD*sqrt(fabs(disc)))/(2*b);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
            }
        }
    } else {
        // Full cubic equation. Coefficients for the depressed cubic t^3 + p*t + q = 0:
        double p = (3*a*c - b*b)/(3*a*a);
        double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
        double disc = -p*p*p/27 - q*q/4;  // This is the discriminant of the depressed cubic but divided by (4*27)
        if (fabs(disc) < XC_EPSILON){
            if (fabs(p) < XC_EPSILON){
                // One real root with multiplicity 3
                t = -b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 3, n_hit, s);
                }
            } else {
                // Two real roots (one simple and one with multiplicity 2)
                t = 3*q/p - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 1, n_hit, s);
                }
                t = -3*q/(2*p) - b/(3*a);
                if (0 <= t && t <= 1){
                    _hit_s_bezier(seg, t, 2, n_hit, s);
                }
            }
        } else if (disc < 0){
            // One real root
            t = cbrt(-q/2 + sqrt(fabs(disc))) + cbrt(-q/2 - sqrt(fabs(disc))) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        } else {
            // Three real roots
            double phi = acos(3*q/(2*p)*sqrt(fabs(3/p)));
            t = 2*sqrt(fabs(p/3))*cos(phi/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 2*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
            t = 2*sqrt(fabs(p/3))*cos((phi + 4*M_PI)/3) - b/(3*a);
            if (0 <= t && t <= 1){
                _hit_s_bezier(seg, t, 1, n_hit, s);
            }
        }
    }
}


// Array of segments
// -----------------

 static inline
void crossing_drift(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                    double part_s, double part_x, double part_tan){
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            crossing_drift_line(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            crossing_drift_halfopenline(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_CIRCULARSEGMENT_ID){
            crossing_drift_circular(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_BEZIERSEGMENT_ID){
            crossing_drift_bezier(segments[i], n_hit, s, part_s, part_x, part_tan);
        } // TODO: else throw fatal error
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}

 static inline
void crossing_drift_vlimit(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                           double part_s, double part_x, double part_tan_x, \
                           double part_y, double part_tan_y, \
                           double y_min, double y_max){
    if (fabs(part_tan_y) < XC_EPSILON){
        // Trajectory parallel to s axis
        if (part_y < y_min || part_y > y_max){
            // No crossing
            return;
        } else {
            // The particle is completely inside the vertical limits, so only check horizontal
            crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
            return;
        }
    } else {
        crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
        // restrict_s is the region [s0, s1] where the particle is inside the vertical limits
        double* restrict_s = (double*) malloc(2*sizeof(double));
        restrict_s[0] = (y_min - part_y)/part_tan_y + part_s;
        restrict_s[1] = (y_max - part_y)/part_tan_y + part_s;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        calculate_overlap_array_interval(s, n_hit, restrict_s);
        free(restrict_s);
    }
}

 static inline
int max_array_size_drift(Segment* segments, int8_t n_segments){
    int size = 0;
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            size += 1;
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            size += 1;
        } else if (id == XC_CIRCULARSEGMENT_ID){
            size += 2;
        } else if (id == XC_BEZIERSEGMENT_ID){
            size += 3;
        } // TODO: else throw fatal error
    }
    return size;
}

#endif /* XCOLL_GEOM_CROSSING_DRIFT_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_OBJECTS_H
#define XCOLL_GEOM_OBJECTS_H


// Assumption for all objects: the particle at -inf is outside the object (otherwise some comparisons might give wrong results)


// Collimator jaw
// --------------

 static inline
Segment* create_jaw(double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side){
    Segment* segments= (Segment*) malloc(3*sizeof(Segment));
    segments[0] = (Segment) create_halfopen_line_segment(s_U, x_U, tilt_tan, side);
    segments[1] = (Segment) create_line_segment(s_U, x_U, s_D, x_D);
    segments[2] = (Segment) create_halfopen_line_segment(s_D, x_D, tilt_tan, side);
    return segments;
}

 static inline
void destroy_jaw(Segment* segments){
    free((HalfOpenLineSegment) segments[0]);
    free((LineSegment)         segments[1]);
    free((HalfOpenLineSegment) segments[2]);
    free(segments);
}


// Polygon
// -------

 static inline
Segment* create_polygon(double* s_poly, double* x_poly, int8_t num_polys){
    Segment* segments= (Segment*) malloc((unsigned int) num_polys*sizeof(Segment));
    for (int8_t i=0; i<num_polys-1; i++){
        segments[i] = (Segment) create_line_segment(s_poly[i], x_poly[i], s_poly[i+1], x_poly[i+1]);
    }
    segments[num_polys-1] = (Segment) create_line_segment(s_poly[num_polys-1], x_poly[num_polys-1], \
                                                          s_poly[0], x_poly[0]);
    return segments;
}

 static inline
void destroy_polygon(Segment* segments, int8_t num_polys){
    for (int8_t i=0; i<num_polys; i++) {
        free((LineSegment) segments[i]);
    }
    free(segments);
}


// Open polygon
// ------------

 static inline
Segment* create_open_polygon(double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side){
    Segment* segments= (Segment*) malloc((num_polys+1)*sizeof(Segment));
    segments[0] = (Segment) create_halfopen_line_segment(s_poly[0], x_poly[0], tilt_tan, side);
    for (int8_t i=1; i<num_polys; i++){
        segments[i] = (Segment) create_line_segment(s_poly[i-1], x_poly[i-1], s_poly[i], x_poly[i]);
    }
    segments[num_polys] = (Segment) create_halfopen_line_segment(s_poly[num_polys-1], x_poly[num_polys-1], \
                                                                 tilt_tan, side);
    return segments;
}

 static inline
void destroy_open_polygon(Segment* segments, int8_t num_polys){
    free((HalfOpenLineSegment) segments[0]);
    for (int8_t i=1; i<num_polys; i++) {
        free((LineSegment)     segments[i]);
    }
    free((HalfOpenLineSegment) segments[num_polys]);
    free(segments);
}


// Crystal
// -------

// The four corners A, B, C, D are such that AB is the front face, BC the curve furthest from the beam,
// CD the back face, and DA the curve closest to the beam.
 static inline
Segment* create_crystal(double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos){
    Segment* segments= (Segment*) malloc(4*sizeof(Segment));

    // First corner is what defines the crystal position
    double A_s = 0;
    double A_x = jaw_U;

    // Manipulate R in function of sign
    double sgnR = (R > 0) - (R < 0);
    double R_short  = sgnR*(fabs(R) - width);
    double sin_a = length/fabs(R);
    double cos_a = sqrt(1 - length*length/R/R);
    if (fabs(R) < XC_EPSILON){
        // straight crystal - not yet implemented 
        printf("Straight crystal not yet implemented!"); //only_for_context cpu_serial
        fflush(stdout);                 //only_for_context cpu_serial
        return NULL;

    } else if (R < 0){
        // This distinction is needed to keep the crystal at the same location when changing the bend direction
        double R_temp = R_short;
        R_short = R;
        R = R_temp;
    }

    // Bending centre is defined w.r.t. A
    double R_s = A_s - R*tilt_sin;
    double R_x = A_x + R*tilt_cos;

    // Three remaining corner points of crystal
    double B_s = R_s + R_short*tilt_sin;
    double C_s = R_s + fabs(R_short)*sin_a*tilt_cos + R_short*cos_a*tilt_sin;
    double D_s = R_s + fabs(R)*sin_a*tilt_cos + R*cos_a*tilt_sin;
    double B_x = R_x - R_short*tilt_cos;
    double C_x = R_x - cos_a*tilt_cos*R_short + sin_a*tilt_sin*fabs(R_short);
    double D_x = R_x - cos_a*tilt_cos*R + sin_a*tilt_sin*fabs(R);
    double A_t = atan2(A_x - R_x, A_s - R_s);
    double D_t = atan2(D_x - R_x, D_s - R_s);
    double t1 = MIN(A_t, D_t);
    double t2 = MAX(A_t, D_t);

    // Fill segments
    segments[0] = (Segment) create_line_segment(A_s, A_x, B_s, B_x);
    segments[1] = (Segment) create_circular_segment(R, R_s, R_x, t1, t2);
    segments[2] = (Segment) create_line_segment(C_s, C_x, D_s, D_x);
    segments[3] = (Segment) create_circular_segment(R_short, R_s, R_x, t1, t2);

    // printf("R: (%f, %f)   A: (%f, %f)   B: (%f, %f)   C: (%f, %f)   D: (%f, %f)   t1: %f   t2: %f\n", R_s,R_x,A_s,A_x,B_s,B_x,C_s,C_x,D_s,D_x,t1*180/3.141592653589793,t2*180/3.141592653589793); fflush(stdout);
    return segments;
}

 static inline
void destroy_crystal(Segment* segments){
    free((LineSegment)     segments[0]);
    free((CircularSegment) segments[1]);
    free((LineSegment)     segments[2]);
    free((CircularSegment) segments[3]);
    free(segments);
}


#endif /* XCOLL_GEOM_OBJECTS_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H


// Find the s-coordinate of the first crossing of a drift with a set of segments
// Here, first means the smallest s-coordinate (needs to be adapted for back-scattering)
 static inline
double crossing_drift_first(Segment* segments, int8_t n_segments, \
                            double part_s, double part_x, double part_tan){
    // printf("crossing_drift_first\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    if (n_hit==0){
        // No crossing
        free(s);
        // printf("crossing_drift_first done (no crossing)\n");fflush(stdout);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    // printf("crossing_drift_first done\n");fflush(stdout);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments
 static inline
double crossing_drift_after_s(Segment* segments, int8_t n_segments, \
                              double part_s, double part_x, double part_tan, \
                              double after_s){
    // printf("crossing_drift_after_s\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            double result = s[i];
            free(s);
            // printf("crossing_drift_after_s done\n");fflush(stdout);
            return result;
        }
    }
    // No crossing
    free(s);
    // printf("crossing_drift_after_s done (no crossing)\n");fflush(stdout);
    return XC_S_MAX;
}

// Find the s-coordinate of the first crossing of a drift with a set of segments including a vertical restriction
 static inline
double crossing_drift_vlimit_first(Segment* segments, int8_t n_segments, \
                                   double part_s, double part_x, double part_tan_x, \
                                   double part_y, double part_tan_y, \
                                   double y_min, double y_max){
    // printf("crossing_drift_vlimit_first\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    if (n_hit==0){
        // No crossing
        free(s);
        // printf("crossing_drift_vlimit_first done (no crossing)\n");fflush(stdout);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    // printf("crossing_drift_vlimit_first done\n");fflush(stdout);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments including a vertical restriction
 static inline
double crossing_drift_vlimit_after_s(Segment* segments, int8_t n_segments, \
                                     double part_s, double part_x, double part_tan_x, \
                                     double part_y, double part_tan_y, \
                                     double y_min, double y_max, double after_s){
    // printf("crossing_drift_vlimit_after_s\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            double result = s[i];
            free(s);
            // printf("crossing_drift_vlimit_after_s done\n");fflush(stdout);
            return result;
        }
    }
    // No crossing
    free(s);
    // printf("crossing_drift_vlimit_after_s done (no crossing)\n");fflush(stdout);
    return XC_S_MAX;
}


#endif /* XCOLL_GEOM_GET_S_H */
// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_ROT_H
#define XCOLL_GEOM_ROT_H


 static inline
double YRotation_single_particle_rotate_only(LocalParticle* part, double s, double angle){
    double x   = LocalParticle_get_x(part);
    double rpp = LocalParticle_get_rpp(part);
    double sin_y = sin(angle);
    double cos_y = cos(angle);
    LocalParticle_set_x(part, x*cos_y - s*sin_y);
    LocalParticle_add_to_px(part,-angle/rpp);
    return x*sin_y + s*cos_y;  // new s
}

#endif /* XCOLL_GEOM_ROT_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_H
#define XCOLL_COLL_GEOM_H

typedef struct CollimatorGeometry_ {
    // Collimator jaws (with tilts)
    double jaw_LU;
    double jaw_RU;
    double length;
    int8_t side;
    // Jaw angles
    double sin_zL;
    double cos_zL;
    double sin_zR;
    double cos_zR;
    double sin_zDiff;
    double cos_zDiff;
    int8_t jaws_parallel;
    // Tilts
    double sin_yL;
    double cos_yL;
    double sin_yR;
    double cos_yR;
    // Segments
    Segment* segments_L;
    Segment* segments_R;
    // Impact table
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_impacts;
    int8_t record_exits;
} CollimatorGeometry_;
typedef CollimatorGeometry_* CollimatorGeometry;


// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
// Furthermore, the particle is moved to the location where it hits the jaw (drifted to the end if no hit),
//              and transformed to the reference frame of that jaw.
 static inline
int8_t hit_jaws_check_and_transform(LocalParticle* part, CollimatorGeometry restrict cg){
    double part_s = 0, part_x, part_tan;
    int8_t is_hit = 0;
    double s_L = 1.e21;
    double s_R = 1.e21;

    // Find the first hit on the left jaw (if not a right-sided collimator)
    if (cg->side != -1){
        SRotation_single_particle(part, cg->sin_zL, cg->cos_zL);
        part_x = LocalParticle_get_x(part);
#ifdef XCOLL_USE_EXACT
        part_tan = LocalParticle_get_exact_xp(part);
#else
        part_tan = LocalParticle_get_xp(part);
#endif
        s_L = crossing_drift_first(cg->segments_L, 3, part_s, part_x, part_tan);
        if (s_L < XC_S_MAX){
            is_hit = 1;
        } else if (cg->side == 1){
            // If left-sided and no hit, rotate back to lab frame
            SRotation_single_particle(part, -cg->sin_zL, cg->cos_zL);
        }
    }

    // if rightsided:            lab  frame
    // if leftsided  and no hit: lab  frame
    // if leftsided  and hit:    left frame
    // if bothsided  and no hit: left frame
    // if bothsided  and hit:    left frame

    // Find the first hit on the right jaw (if not a left-sided collimator)
    if (cg->side != 1){
        if (cg->side == -1){
            // We didn't rotate to the left frame earlier, so full rotation to the right frame now
            SRotation_single_particle(part, cg->sin_zR, cg->cos_zR);
        } else if (!cg->jaws_parallel){
            // We rotated to the left frame before, so now rotate the difference
            SRotation_single_particle(part, cg->sin_zDiff, cg->cos_zDiff);
        }
        part_x = LocalParticle_get_x(part);
#ifdef XCOLL_USE_EXACT
        part_tan = LocalParticle_get_exact_xp(part);
#else
        part_tan = LocalParticle_get_xp(part);
#endif
        s_R = crossing_drift_first(cg->segments_R, 3, part_s, part_x, part_tan);
        if (s_R < XC_S_MAX && s_R < s_L){
            is_hit = -1;
        } else if (is_hit == 1){
            if (!cg->jaws_parallel){
                // Rotate back to left frame
                SRotation_single_particle(part, -cg->sin_zDiff, cg->cos_zDiff);
            }
        } else {
            // No hit, rotate back to lab frame
            SRotation_single_particle(part, -cg->sin_zR, cg->cos_zR);
        }
    }

    // if rightsided and no hit: lab   frame
    // if rightsided and hit:    right frame
    // if leftsided  and no hit: lab   frame
    // if leftsided  and hit:    left  frame
    // if bothsided  and no hit: lab  frame
    // if bothsided  and hit:    hit   frame

    // Drift to the impact position or end, and move to jaw frame if relevant
    if (is_hit == 1){
        // Move to the impact position
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, s_L);
#else
        Drift_single_particle_expanded(part, s_L);
#endif
        // Shift the reference frame to the jaw corner LU
        XYShift_single_particle(part, cg->jaw_LU, 0);
        LocalParticle_add_to_s(part, -cg->length/2*(1 - cg->cos_yL));
        // Rotate the reference frame to tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), asin(cg->sin_yL));
        LocalParticle_set_s(part, new_s);
        if (cg->record_impacts){
            InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
        }

    } else if (is_hit == -1){
        // Move to the impact position
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, s_R);
#else
        Drift_single_particle_expanded(part, s_R);
#endif
        // Shift the reference frame to the jaw corner RU
        XYShift_single_particle(part, cg->jaw_RU, 0);
        LocalParticle_add_to_s(part, -cg->length/2*(1 - cg->cos_yR));
        // Rotate the reference frame to tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), asin(cg->sin_yR));
        LocalParticle_set_s(part, new_s);
        // Mirror x
        LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
        LocalParticle_scale_exact_xp(part, -1);
#else
        LocalParticle_scale_xp(part, -1);
#endif
        if (cg->record_impacts){
            InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
        }

    } else {
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, cg->length);
#else
        Drift_single_particle_expanded(part, cg->length);
#endif
    }

    return is_hit;
}


 static inline
void hit_jaws_transform_back(int8_t is_hit, LocalParticle* part, CollimatorGeometry restrict cg){
    if (is_hit != 0 && LocalParticle_get_state(part) > 0){
        if (cg->record_exits){
            InteractionRecordData_log(cg->record, cg->record_index, part, XC_EXIT_JAW);
        }
    }
    if (is_hit == 1){
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_yL));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner LU
        XYShift_single_particle(part, -cg->jaw_LU, 0);
        LocalParticle_add_to_s(part, cg->length/2*(1 - cg->cos_yL));
        // If particle survived, drift to end of element
        if (LocalParticle_get_state(part) > 0){
#ifdef XCOLL_USE_EXACT
            Drift_single_particle_exact(part, cg->length - LocalParticle_get_s(part));
#else
            Drift_single_particle_expanded(part, cg->length - LocalParticle_get_s(part));
#endif
        }
        SRotation_single_particle(part, -cg->sin_zL, cg->cos_zL);

    } else if (is_hit == -1){
        // Mirror back
        LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
        LocalParticle_scale_exact_xp(part, -1);
#else
        LocalParticle_scale_xp(part, -1);
#endif
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_yR));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner RU
        XYShift_single_particle(part, -cg->jaw_RU, 0);
        LocalParticle_add_to_s(part, cg->length/2*(1 - cg->cos_yR));
        // If particle survived, drift to end of element
        if (LocalParticle_get_state(part) > 0){
#ifdef XCOLL_USE_EXACT
            Drift_single_particle_exact(part, cg->length - LocalParticle_get_s(part));
#else
            Drift_single_particle_expanded(part, cg->length - LocalParticle_get_s(part));
#endif
        }
        SRotation_single_particle(part, -cg->sin_zR, cg->cos_zR);
    }
}


#endif /* XCOLL_COLL_GEOM_H */

// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_CRY_GEOM_H
#define XCOLL_CRY_GEOM_H

typedef struct CrystalGeometry_ {
    // Crystal inner upstream corner (with tilt)
    double jaw_U;
    double length;
    int8_t side;
    // Jaw angle
    double sin_z;
    double cos_z;
    // Tilt
    double sin_y;
    double cos_y;
    // Crystal geometry
    double bending_radius;
    double bending_angle;
    double width;
    double height;
    double miscut_angle;
    double s_B;    // Bend centre
    double x_B;
    double s_P;    // Miscut centre
    double x_P;
    double t_VImax;
    // Segments
    Segment* segments;
    Segment* segments_jf;  // Segments in jaw frame
    // Impact table
    InteractionRecordData record;
    RecordIndex record_index;
    int8_t record_impacts;
    int8_t record_exits;
} CrystalGeometry_;
typedef CrystalGeometry_* CrystalGeometry;


// This function checks if a particle hits a jaw (and which).
// Return value: 0 (no hit), 1 (hit on left jaw), -1 (hit on right jaw).
// Furthermore, the particle is moved to the location where it hits the jaw (drifted to the end if no hit),
//              and transformed to the reference frame of that jaw.
 static inline
int8_t hit_crystal_check_and_transform(LocalParticle* part, CrystalGeometry restrict cg){
    double part_s = 0, part_x, part_tan_x, part_y, part_tan_y;
    double s = 1.e21;

    // Crystal should be single-sided
    if (cg->side == 0){
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
    }

    SRotation_single_particle(part, cg->sin_z, cg->cos_z);
    part_x = LocalParticle_get_x(part);
#ifdef XCOLL_USE_EXACT
    part_tan_x = LocalParticle_get_exact_xp(part);
#else
    part_tan_x = LocalParticle_get_xp(part);
#endif
    part_y = LocalParticle_get_y(part);
#ifdef XCOLL_USE_EXACT
    part_tan_y = LocalParticle_get_exact_yp(part);
#else
    part_tan_y = LocalParticle_get_yp(part);
#endif
    s = crossing_drift_vlimit_first(cg->segments, 4, part_s, part_x, part_tan_x, part_y, part_tan_y, -cg->height/2, cg->height/2);

    if (s < XC_S_MAX){
        // Hit: Drift to the impact position, and move to jaw frame if relevant
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, s);
#else
        Drift_single_particle_expanded(part, s);
#endif
        // Shift the reference frame to the upstream jaw corner (for a crystal, this is always at s=0)
        XYShift_single_particle(part, cg->jaw_U, 0);
        // Rotate the reference frame to tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        if (cg->side == 1){
            if (cg->record_impacts){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_L);
            }

        } else {
            // Mirror x
            LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
            LocalParticle_scale_exact_xp(part, -1);
#else
            LocalParticle_scale_xp(part, -1);
#endif
            if (cg->record_impacts){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_ENTER_JAW_R);
            }
        }
        return cg->side;

    } else {
        // No hit, rotate back to lab frame
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
        // Drift to end
#ifdef XCOLL_USE_EXACT
        Drift_single_particle_exact(part, cg->length);
#else
        Drift_single_particle_expanded(part, cg->length);
#endif
        return 0;
    }
}


 static inline
void hit_crystal_transform_back(int8_t is_hit, LocalParticle* part, CrystalGeometry restrict cg){
    if (is_hit != 0){
        if (LocalParticle_get_state(part) > 0){
            if (cg->record_exits){
                InteractionRecordData_log(cg->record, cg->record_index, part, XC_EXIT_JAW);
            }
        }
        if (cg->side == -1){
            // Mirror back
            LocalParticle_scale_x(part, -1);
#ifdef XCOLL_USE_EXACT
            LocalParticle_scale_exact_xp(part, -1);
#else
            LocalParticle_scale_xp(part, -1);
#endif
        }
        // Rotate back from tilt
        double new_s = YRotation_single_particle_rotate_only(part, LocalParticle_get_s(part), -asin(cg->sin_y));
        LocalParticle_set_s(part, new_s);
        // Shift the reference frame back from jaw corner U
        XYShift_single_particle(part, -cg->jaw_U, 0);
        // If particle survived, drift to end of element
        if (LocalParticle_get_state(part) > 0){
#ifdef XCOLL_USE_EXACT
            Drift_single_particle_exact(part, cg->length - LocalParticle_get_s(part));
#else
            Drift_single_particle_expanded(part, cg->length - LocalParticle_get_s(part));
#endif
        }
        SRotation_single_particle(part, -cg->sin_z, cg->cos_z);
    }
}


#endif /* XCOLL_CRY_GEOM_H */

#ifndef XOBJ_TYPEDEF_XcollGeometryTestData
#define XOBJ_TYPEDEF_XcollGeometryTestData
typedef   struct XcollGeometryTestData_s * XcollGeometryTestData;
 static inline XcollGeometryTestData XcollGeometryTestData_getp(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  return (XcollGeometryTestData)(( char*) obj+offset);
}
 static inline double XcollGeometryTestData_get__sin_rot_s(const XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryTestData_set__sin_rot_s(XcollGeometryTestData restrict  obj, double value){
  int64_t offset=0;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryTestData_getp__sin_rot_s(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryTestData_get__cos_rot_s(const XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryTestData_set__cos_rot_s(XcollGeometryTestData restrict  obj, double value){
  int64_t offset=0;
  offset+=8;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryTestData_getp__cos_rot_s(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=8;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryTestData_get__shift_x(const XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryTestData_set__shift_x(XcollGeometryTestData restrict  obj, double value){
  int64_t offset=0;
  offset+=16;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryTestData_getp__shift_x(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=16;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryTestData_get__shift_y(const XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryTestData_set__shift_y(XcollGeometryTestData restrict  obj, double value){
  int64_t offset=0;
  offset+=24;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryTestData_getp__shift_y(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=24;
  return ( double*)(( char*) obj+offset);
}
 static inline double XcollGeometryTestData_get__shift_s(const XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return *( double*)(( char*) obj+offset);
}
 static inline void XcollGeometryTestData_set__shift_s(XcollGeometryTestData restrict  obj, double value){
  int64_t offset=0;
  offset+=32;
  *( double*)(( char*) obj+offset)=value;
}
 static inline  double* XcollGeometryTestData_getp__shift_s(XcollGeometryTestData restrict  obj){
  int64_t offset=0;
  offset+=32;
  return ( double*)(( char*) obj+offset);
}
#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_CONSTANTS_H
#define XTRACK_CONSTANTS_H

#if !defined( C_LIGHT )
    #define   C_LIGHT ( 299792458.0 )
#endif /* !defined( C_LIGHT ) */

#if !defined( EPSILON_0 )
    #define   EPSILON_0 (8.854187817620e-12)
#endif /* !defined( EPSILON_0 ) */

#if !defined( PI )
    #define PI (3.1415926535897932384626433832795028841971693993751)
#endif /* !defined( PI ) */

#if !defined( MU_0 )
    #define MU_0 (PI*4.0e-7)
#endif /* !defined( MU_0 ) */

#if !defined( DEG2RAD )
    #define DEG2RAD (0.0174532925199432957692369076848861271344287188854)
#endif /* !defiend( DEG2RAD ) */

#if !defined( RAD2DEG )
    #define RAD2DEG (57.29577951308232087679815481410517033240547246656442)
#endif /* !defiend( RAD2DEG ) */

#if !defined( SQRT_PI )
    #define SQRT_PI (1.7724538509055160272981674833411451827975494561224)
#endif /* !defined( SQRT_PI ) */

#if !defined( QELEM )
    #define QELEM (1.60217662e-19)
#endif /* !defined( QELEM ) */

#if !defined( DBL_MAX )
    #define DBL_MAX (1.7976931348623158e+308)
#endif /* !defined( DBL_MAX ) */

#if !defined( DBL_MIN )
    #define DBL_MIN (2.2250738585072014e-308)
#endif /* !defined( DBL_MIN ) */

#if !defined( DBL_EPSILON )
    #define DBL_EPSILON (2.2204460492503131e-16)
#endif /* !defined( DBL_EPSILON ) */

#endif /* XTRACK_CONSTANTS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_FUNCTIONS_H
#define XTRACK_FUNCTIONS_H

 static inline
void kill_all_particles(LocalParticle* part0, int64_t kill_state) {

    {
//    const int64_t XT_part_block_start_idx = part0->ipart; //only_for_context cpu_openmp
//    const int64_t XT_part_block_end_idx = part0->endpart; //only_for_context cpu_openmp

    const int64_t XT_part_block_start_idx = 0;                                            //only_for_context cpu_serial
    const int64_t XT_part_block_end_idx = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial

    //#pragma omp simd // TODO: currently does not work, needs investigating
    for (int64_t XT_part_block_ii = XT_part_block_start_idx; XT_part_block_ii<XT_part_block_end_idx; XT_part_block_ii++) { //only_for_context cpu_openmp cpu_serial

        LocalParticle lpart = *part0;    //only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;    //only_for_context cpu_serial cpu_openmp
        part->ipart = XT_part_block_ii;  //only_for_context cpu_serial cpu_openmp

//        LocalParticle* part = part0;     //only_for_context opencl cuda

//        if (LocalParticle_get_state(part) > 0) {  //only_for_context cpu_openmp

        LocalParticle_kill_particle(part, kill_state);

//        }  //only_for_context cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
    }

}


 static inline
int8_t assert_tracking(LocalParticle* part, int64_t kill_state){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss.
    if (LocalParticle_get_at_turn(part) < 0){
        LocalParticle_kill_particle(part, kill_state);
        return 0;
    }
    return 1;
}

#endif /* XTRACK_FUNCTIONS_H */

// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //

#ifndef XTRACK_PARTICLE_STATES_H
#define XTRACK_PARTICLE_STATES_H

#define  XT_LOST_ON_APERTURE       0
#define  XT_LOST_ON_LONG_CUT      -2
#define  XT_LOST_ALL_E_IN_SYNRAD -10
#define  RNG_ERR_SEEDS_NOT_SET   -20
#define  RNG_ERR_INVALID_TRACK   -21
#define  RNG_ERR_RUTH_NOT_SET    -22

#endif /* XTRACK_PARTICLE_STATES_H */

#ifndef XTRACK_TRACK_SROTATION_H
#define XTRACK_TRACK_SROTATION_H


 static inline
void SRotation_single_particle(LocalParticle* part, double sin_z, double cos_z){

    double const x  = LocalParticle_get_x(part);
    double const y  = LocalParticle_get_y(part);
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);

    double const x_hat  =  cos_z * x  + sin_z * y;
    double const y_hat  = -sin_z * x  + cos_z * y;
    double const px_hat =  cos_z * px + sin_z * py;
    double const py_hat = -sin_z * px + cos_z * py;

    LocalParticle_set_x(part, x_hat);
    LocalParticle_set_y(part, y_hat);
    LocalParticle_set_px(part, px_hat);
    LocalParticle_set_py(part, py_hat);

}

#endif
// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_DRIFT_H
#define XTRACK_DRIFT_H


 static inline
void Drift_single_particle_expanded(LocalParticle* part, double length){
    double const rpp    = LocalParticle_get_rpp(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const xp     = LocalParticle_get_px(part) * rpp;
    double const yp     = LocalParticle_get_py(part) * rpp;
    double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

    LocalParticle_add_to_x(part, xp * length );
    LocalParticle_add_to_y(part, yp * length );
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length * dzeta );
}


 static inline
void Drift_single_particle_exact(LocalParticle* part, double length){
    double const px = LocalParticle_get_px(part);
    double const py = LocalParticle_get_py(part);
    double const rv0v    = 1./LocalParticle_get_rvv(part);
    double const one_plus_delta = 1. + LocalParticle_get_delta(part);

    double const one_over_pz = 1./sqrt(one_plus_delta*one_plus_delta
                                       - px * px - py * py);
    double const dzeta = 1 - rv0v * one_plus_delta * one_over_pz;

    LocalParticle_add_to_x(part, px * one_over_pz * length);
    LocalParticle_add_to_y(part, py * one_over_pz * length);
    LocalParticle_add_to_zeta(part, dzeta * length);
    LocalParticle_add_to_s(part, length);
}


 static inline
void Drift_single_particle(LocalParticle* part, double length){
    #ifndef XTRACK_USE_EXACT_DRIFTS
        Drift_single_particle_expanded(part, length);
    #else
        Drift_single_particle_exact(part, length);
    #endif
}


#endif /* XTRACK_DRIFT_H */

 static inline
double test_jaw_first(double part_s, double part_x, double part_tan_x, double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side){
    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
    double s = crossing_drift_first(segments, 3, part_s, part_x, part_tan_x);
    destroy_jaw(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_jaw_after_s(double part_s, double part_x, double part_tan_x, double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side, double after_s){
    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
    double s = crossing_drift_after_s(segments, 3, part_s, part_x, part_tan_x, after_s);
    destroy_jaw(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_jaw_vlimit_first(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side, double y_min, double y_max){
    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
    double s = crossing_drift_vlimit_first(segments, 3, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);
    destroy_jaw(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_jaw_vlimit_after_s(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side, double y_min, double y_max, double after_s){
    Segment* segments = create_jaw(s_U, x_U, s_D, x_D, tilt_tan, side);
    double s = crossing_drift_vlimit_after_s(segments, 3, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);
    destroy_jaw(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_polygon_first(double part_s, double part_x, double part_tan_x, double* s_poly, double* x_poly, int8_t num_polys){
    Segment* segments = create_polygon(s_poly, x_poly, num_polys);
    double s = crossing_drift_first(segments, num_polys, part_s, part_x, part_tan_x);
    destroy_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_polygon_after_s(double part_s, double part_x, double part_tan_x, double* s_poly, double* x_poly, int8_t num_polys, double after_s){
    Segment* segments = create_polygon(s_poly, x_poly, num_polys);
    double s = crossing_drift_after_s(segments, num_polys, part_s, part_x, part_tan_x, after_s);
    destroy_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_polygon_vlimit_first(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double* s_poly, double* x_poly, int8_t num_polys, double y_min, double y_max){
    Segment* segments = create_polygon(s_poly, x_poly, num_polys);
    double s = crossing_drift_vlimit_first(segments, num_polys, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);
    destroy_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_polygon_vlimit_after_s(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double* s_poly, double* x_poly, int8_t num_polys, double y_min, double y_max, double after_s){
    Segment* segments = create_polygon(s_poly, x_poly, num_polys);
    double s = crossing_drift_vlimit_after_s(segments, num_polys, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);
    destroy_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_open_polygon_first(double part_s, double part_x, double part_tan_x, double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side){
    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);
    double s = crossing_drift_first(segments, num_polys+1, part_s, part_x, part_tan_x);
    destroy_open_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_open_polygon_after_s(double part_s, double part_x, double part_tan_x, double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side, double after_s){
    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);
    double s = crossing_drift_after_s(segments, num_polys+1, part_s, part_x, part_tan_x, after_s);
    destroy_open_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_open_polygon_vlimit_first(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side, double y_min, double y_max){
    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);
    double s = crossing_drift_vlimit_first(segments, num_polys+1, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);
    destroy_open_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_open_polygon_vlimit_after_s(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side, double y_min, double y_max, double after_s){
    Segment* segments = create_open_polygon(s_poly, x_poly, num_polys, tilt_tan, side);
    double s = crossing_drift_vlimit_after_s(segments, num_polys+1, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);
    destroy_open_polygon(segments, num_polys);  // Important to free memory!!
    return s;
}

 static inline
double test_crystal_first(double part_s, double part_x, double part_tan_x, double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos){
    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);
    double s = crossing_drift_first(segments, 4, part_s, part_x, part_tan_x);
    destroy_crystal(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_crystal_after_s(double part_s, double part_x, double part_tan_x, double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos, double after_s){
    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);
    double s = crossing_drift_after_s(segments, 4, part_s, part_x, part_tan_x, after_s);
    destroy_crystal(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_crystal_vlimit_first(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos, double y_min, double y_max){
    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);
    double s = crossing_drift_vlimit_first(segments, 4, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max);
    destroy_crystal(segments);  // Important to free memory!!
    return s;
}

 static inline
double test_crystal_vlimit_after_s(double part_s, double part_x, double part_tan_x, double part_y, double part_tan_y, double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos, double y_min, double y_max, double after_s){
    Segment* segments = create_crystal(R, width, length, jaw_U, tilt_sin, tilt_cos);
    double s = crossing_drift_vlimit_after_s(segments, 4, part_s, part_x, part_tan_x, part_y, part_tan_y, y_min, y_max, after_s);
    destroy_crystal(segments);  // Important to free memory!!
    return s;
}