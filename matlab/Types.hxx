
#ifndef _TYPES_HXX
#define _TYPES_HXX

#if defined(FLOAT_PVAR)
typedef float PhaseSpaceType; // Type for position and velocity of particle
#elif defined(DOUBLE_PVAR)
typedef double PhaseSpaceType; // Type for position and velocity of particle
#endif

// Enum type for phase space
typedef enum PhaseSpaceDim {
  xdim=0, ydim, zdim, nrdim,
  vxdim=0, vydim, vzdim, nvdim
} PhaseSpaceDim;

#ifdef HDF

typedef float HdfFieldType;
#define DFNT_FIELD DFNT_FLOAT

typedef unsigned short HdfSampleType;
#define DFNT_SAMPLE DFNT_UINT16

typedef unsigned short HdfMultiSampleType;
#define DFNT_MULTISAMPLE DFNT_UINT16

typedef float HdfAxisType;
#define DFNT_AXIS DFNT_FLOAT

#endif


#endif


