SUBROUTINE CLOUD_DETECT_SETUP

!**** Cloud detection setup
!        A. Collard  ECMWF 01/02/06   

!     PURPOSE
!     -------
!        Initialise cloud detection parameters for advanced infrared sounders

!**   INTERFACE
!     ---------
!        cloud_detect_setup is called from defrun

!     METHOD
!     ------
!        Default values are assigned to the cloud detections setup structure.

!     EXTERNALS
!     ---------

!     MODIFICATIONS
!     -------------
!     01/02/06  Original code.     A. Collard
!     27/04/10  Modified the default channel selection for AIRS: Channel
!               1852 to band 3, and 1937 to band 5.  R. Eresmaa
!     10/10/11  Added CRIS TonyMC                  
!     12/02/13  Added setup for imager-assisted cloud detection.  R.Eresmaa
!     T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     19/11/14  C. Lupu  Updates to sensors and platform lists from rttov_consts
!     21/01/15  J. Letertre-Danczak change IASI aerosol setup
!------------------------------------------------------------------------------

USE PARKIND1,  ONLY : JPIM, JPRB
USE YOMHOOK,   ONLY : LHOOK, DR_HOOK
USE YOMLUN,    ONLY : NULOUT, NULUSR3
USE YOMCT0,    ONLY : LECMWF
USE YOMMP0,    ONLY : MYPROC, NPROC 
USE MPL_MODULE,ONLY : MPL_BROADCAST
USE YOMCLDDET, ONLY : S__CLOUD_DETECT_SETUP, &
&                     JP__MIN_SENSOR_INDEX,  &
&                     JP__MAX_SENSOR_INDEX,  &
&                     JP__DIGITAL_FILTER
USE RTTOV_CONST, ONLY: INST_ID_AIRS,          &
&                      INST_ID_IASI,          &
&                      INST_ID_CRIS


IMPLICIT NONE

! Local variables

CHARACTER(LEN=5)   :: CL__INSTRUMENTNAME
CHARACTER(LEN=15)  :: CL__CLOUD_DETECTION_FILE

INTEGER(KIND=JPIM)   :: J, J__SENSOR      ! Loop variables
INTEGER(KIND=JPIM)   :: IOMASTER, ITAG
INTEGER(KIND=JPIM)   :: INIU1, IOS
INTEGER(KIND=JPIM)   :: I__MAXAEROCHANS

REAL(KIND=JPRB)      :: ZHOOK_HANDLE

!-----------------------
! Namelist variables 
!-----------------------

! N.B. Max_Bands must be greater than 5
INTEGER(KIND=JPIM), PARAMETER :: JP__MAX_BANDS    =         8_JPIM
INTEGER(KIND=JPIM), PARAMETER :: JP__MAX_CHANNELS =      8461_JPIM
INTEGER(KIND=JPIM), PARAMETER :: JP__MAX_AEROSOL_CHANS =  200_JPIM

INTEGER(KIND=JPIM)  :: M__SENSOR 
INTEGER(KIND=JPIM)  :: N__FILTER_METHOD 
INTEGER(KIND=JPIM)  :: N__NUM_BANDS 
INTEGER(KIND=JPIM)  :: N__GRADCHKINTERVAL(JP__MAX_BANDS)
INTEGER(KIND=JPIM)  :: N__BAND_SIZE(JP__MAX_BANDS)
INTEGER(KIND=JPIM)  :: N__BANDS(JP__MAX_CHANNELS,JP__MAX_BANDS) 
INTEGER(KIND=JPIM)  :: N__WINDOW_WIDTH(JP__MAX_BANDS)
INTEGER(KIND=JPIM)  :: N__WINDOW_BOUNDS(JP__MAX_BANDS,2)
REAL(KIND=JPRB)     :: R__BT_THRESHOLD(JP__MAX_BANDS) 
REAL(KIND=JPRB)     :: R__GRAD_THRESHOLD(JP__MAX_BANDS)
REAL(KIND=JPRB)     :: R__WINDOW_GRAD_THRESHOLD(JP__MAX_BANDS)
LOGICAL             :: L__DO_QUICK_EXIT  
LOGICAL             :: L__DO_CROSSBAND
INTEGER(KIND=JPIM)  :: N__BANDTOUSE(JP__MAX_BANDS)

! Imager cloud detection
LOGICAL            :: L__DO_IMAGER_CLOUD_DETECTION
INTEGER(KIND=JPIM) :: N__NUM_IMAGER_CHANS
INTEGER(KIND=JPIM) :: N__NUM_IMAGER_CLUSTERS
INTEGER(KIND=JPIM) :: N__IMAGER_CHANS(JP__MAX_BANDS)
REAL(KIND=JPRB)    :: R__STDDEV_THRESHOLD(JP__MAX_BANDS)
REAL(KIND=JPRB)    :: R__COVERAGE_THRESHOLD
REAL(KIND=JPRB)    :: R__FG_DEPARTURE_THRESHOLD

! Aerosol:
LOGICAL            :: L__DO_AEROSOLDETECTION   
INTEGER(KIND=JPIM) :: N__NUM_AEROSOL_TESTS
INTEGER(KIND=JPIM) :: N__NUM_AEROSOL_CHANS(JP__MAX_BANDS)
INTEGER(KIND=JPIM) :: N__AEROSOL_CHANS(JP__MAX_AEROSOL_CHANS,JP__MAX_BANDS)
REAL(KIND=JPRB)    :: R__AEROSOL_ABSCISSAE(JP__MAX_AEROSOL_CHANS,JP__MAX_BANDS)
REAL(KIND=JPRB)    :: R__AEROSOL_THRESHOLD(JP__MAX_BANDS)
REAL(KIND=JPRB)    :: R__AEROSOLMINNORM(JP__MAX_BANDS)

! Buffers used in broadcasting between processing elements.
INTEGER(KIND=JPIM) :: IBUF(10,JP__MAX_BANDS)
REAL(KIND=JPRB)    :: ZBUF(7,JP__MAX_BANDS)
INTEGER(KIND=JPIM) :: III


#include "abor1.intfb.h"
#include "namclddet.nam.h"

!============================================================================

IF (LHOOK) CALL DR_HOOK('CLOUD_DETECT_SETUP',0,ZHOOK_HANDLE)

IF (JP__MAX_BANDS < 5) CALL ABOR1 ('JP__Max_Bands must be greater than 5')

!============================================================================
!   Loop through sensors setting up cloud detection
!============================================================================

SENSORLOOP : DO J__SENSOR = JP__MIN_SENSOR_INDEX, JP__MAX_SENSOR_INDEX

   
   SELECT CASE (J__SENSOR)
      
   CASE(INST_ID_AIRS)
      !====================
      ! Set up AIRS
      !====================

      CL__INSTRUMENTNAME='AIRS'
      CL__CLOUD_DETECTION_FILE = 'AIRS_CLDDET.NL'

      N__FILTER_METHOD = JP__DIGITAL_FILTER
      
      N__NUM_BANDS = 5

      N__BAND_SIZE(:) = 0

      IF (LECMWF) THEN
        N__BAND_SIZE(1:N__NUM_BANDS) =(/141, 36, 54, 23, 65 /)
      ELSE
        N__BAND_SIZE(1:N__NUM_BANDS) =(/137, 36, 53, 22, 63 /)
      ENDIF

      IF (MAXVAL(N__BAND_SIZE(:)) > JP__MAX_CHANNELS) &
&              CALL ABOR1('Too many channels specified in cloud '//&
&                      'detection - increase JP__Max_Channels')
 
      N__BANDS(:,:)= 0 

      IF (LECMWF) THEN
        N__BANDS(1:N__BAND_SIZE(1),1) = &
&      (/ 1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&        22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&        52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&        72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&        83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&       108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&       138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&       165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&       180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&       207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&       232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&       262, 267, 272, 295, 299, 300, 305, 308, 309, 310, &
&       318, 321, 325, 333, 338, 355, 362, 375, 453, 475, &
&       484, 497, 528, 587, 672, 787, 791, 843, 870, 914, &
&       950 /)
      ELSE
        N__BANDS(1:N__BAND_SIZE(1),1) = &
&      (/ 1,   6,   7,  10,  11,  15,  16,  17,  20,  21, &
&        22,  24,  27,  28,  30,  36,  39,  40,  42,  51, &
&        52,  54,  55,  56,  59,  62,  63,  68,  69,  71, &
&        72,  73,  74,  75,  76,  77,  78,  79,  80,  82, &
&        83,  84,  86,  92,  93,  98,  99, 101, 104, 105, &
&       108, 110, 111, 113, 116, 117, 123, 124, 128, 129, &
&       138, 139, 144, 145, 150, 151, 156, 157, 159, 162, &
&       165, 168, 169, 170, 172, 173, 174, 175, 177, 179, &
&       180, 182, 185, 186, 190, 192, 193, 198, 201, 204, &
&       207, 210, 213, 215, 216, 218, 221, 224, 226, 227, &
&       232, 239, 248, 250, 251, 252, 253, 256, 257, 261, &
&       262, 267, 272, 295, 299,      305, 308, 309, 310, &
&       318, 321,      333, 338, 355, 362, 375,      475, &
&       484, 497, 528, 587, 672, 787, 791,      870, 914, &
&       950 /)
      ENDIF

      N__BANDS(1:N__BAND_SIZE(2),2) = &
&    (/ 1003, 1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082,  &
&       1083, 1088, 1090, 1092, 1095, 1104, 1111, 1115, 1116, 1119,  &
&       1120, 1123, 1130, 1138, 1142, 1178, 1199, 1206, 1221, 1237,  &
&       1252, 1260, 1263, 1266, 1278, 1285 /)

      IF (LECMWF) THEN
        N__BANDS(1:N__BAND_SIZE(3),3) = &
&      (/ 1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &
&         1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &
&         1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, &
&         1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
&         1763, 1766, 1771, 1777, 1780, 1783, 1794, 1800, 1803, 1806, &
&         1812, 1826, 1843, 1852 /)

      N__BANDS(1:N__BAND_SIZE(4),4) = &
&      (/ 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 1877, &
&         1881, 1882, 1883, 1884, 1897, 1901, 1911, 1917, 1918, 1921, &
&         1923, 1924, 1928 /)

      N__BANDS(1:N__BAND_SIZE(5),5) = &
&      (/ 1937, 1938, 1939, 1941, 1946, 1947, 1948, 1958, 1971, 1973, &
&         1988, 1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, &
&         2104, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, &
&         2115, 2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, &
&         2134, 2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, &
&         2234, 2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, &
&         2355, 2363, 2370, 2371, 2377  /)  
      ELSE
        N__BANDS(1:N__BAND_SIZE(3),3) = &
&      (/ 1290, 1301, 1304, 1329, 1371, 1382, 1415, 1424, 1449, 1455, &
&         1466, 1471, 1477, 1479, 1488, 1500, 1519, 1520, 1538, 1545, &
&         1565, 1574, 1583, 1593, 1614, 1627, 1636, 1644, 1652, 1669, &
&         1674, 1681, 1694, 1708, 1717, 1723, 1740, 1748, 1751, 1756, &
&         1763, 1766, 1771, 1777, 1780,       1794, 1800, 1803, 1806, &
&         1812, 1826, 1843, 1852 /)

        N__BANDS(1:N__BAND_SIZE(4),4) = &
&      (/ 1865, 1866, 1867, 1868, 1869, 1872, 1873, 1875, 1876, 1877, &
&         1881, 1882, 1883,       1897, 1901, 1911, 1917, 1918, 1921, &
&         1923, 1924, 1928 /)

        N__BANDS(1:N__BAND_SIZE(5),5) = &
&      (/ 1937, 1938, 1939, 1941,       1947, 1948, 1958, 1971, 1973, &
&         1988, 1995, 2084, 2085, 2097, 2098, 2099, 2100, 2101, 2103, &
&               2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, &
&         2115, 2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2128, &
&         2134, 2141, 2145, 2149, 2153, 2164, 2189, 2197, 2209, 2226, &
&         2234, 2280, 2318, 2321, 2325, 2328, 2333, 2339, 2348, 2353, &
&         2355, 2363, 2370, 2371, 2377  /)  
      ENDIF


      N__GRADCHKINTERVAL(:) = 0
      N__GRADCHKINTERVAL(1:N__NUM_BANDS) = (/ 5,5,5,5,5 /)

      N__WINDOW_WIDTH(:) = 0
      N__WINDOW_WIDTH(1:N__NUM_BANDS) = (/ 10,6,8,5,8 /)

      N__WINDOW_BOUNDS(:,:) = 0
      N__WINDOW_BOUNDS(1,1) = 475
      N__WINDOW_BOUNDS(1,2) = 950

      R__BT_THRESHOLD(:) = 0.
      R__BT_THRESHOLD(1:N__NUM_BANDS) = (/ 0.5, 0.5, 0.5, 0.5, 0.5/)

      R__GRAD_THRESHOLD(:) = 0.
      R__GRAD_THRESHOLD(1:N__NUM_BANDS) = (/ 0.02, 0.02, 0.02, 0.02, 0.02 /)

      R__WINDOW_GRAD_THRESHOLD(:) = 0.
      R__WINDOW_GRAD_THRESHOLD(1) = 0.4

      L__DO_QUICK_EXIT = .TRUE.

      ! This is cross-band:
      
      L__DO_CROSSBAND = .TRUE.

      N__BANDTOUSE(:) = 0
      IF (LECMWF) THEN
        N__BANDTOUSE(1:N__NUM_BANDS) = (/ 1,1,1,4,5 /)
      ELSE
        N__BANDTOUSE(1:N__NUM_BANDS) = (/ 1,2,1,4,5 /)
      ENDIF

      ! This is the setup for imager cloud detection

      L__DO_IMAGER_CLOUD_DETECTION = .FALSE.
      N__NUM_IMAGER_CHANS = 0
      N__NUM_IMAGER_CLUSTERS = 0
      N__IMAGER_CHANS(:) = 0
      R__STDDEV_THRESHOLD(:) = 0.0
      R__COVERAGE_THRESHOLD = 0.0
      R__FG_DEPARTURE_THRESHOLD = 0.0

      ! This is aerosol:

      IF (LECMWF) THEN
        L__DO_AEROSOLDETECTION = .TRUE.
      ELSE
        L__DO_AEROSOLDETECTION = .FALSE.
      ENDIF
  
      N__NUM_AEROSOL_TESTS = 1
      N__NUM_AEROSOL_CHANS(:) = 0
      N__NUM_AEROSOL_CHANS(1:N__NUM_AEROSOL_TESTS) = (/ 8 /)
      N__AEROSOL_CHANS(:,:) = 0
      N__AEROSOL_CHANS(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/ 1178, 1199, 1221, 1237, 1252, 1263, 1285, 1290 /)
      R__AEROSOL_ABSCISSAE(:,:) = 0.0
      R__AEROSOL_ABSCISSAE(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/ 1.70747, 1.58687, 1.47237, 1.34777, 1.23624, 0.518735, &
&              0.508225, 0.522959 /)
      R__AEROSOL_THRESHOLD(:) = 0.0
      R__AEROSOL_THRESHOLD(1:N__NUM_AEROSOL_TESTS) = (/ 0.4 /)
      R__AEROSOLMINNORM(:) = 0.0
      R__AEROSOLMINNORM(1:N__NUM_AEROSOL_TESTS) = (/ 2.0 /)

   CASE(INST_ID_IASI)
      !====================
      ! Set up IASI
      !====================

      CL__INSTRUMENTNAME='IASI'
      CL__CLOUD_DETECTION_FILE = 'IASI_CLDDET.NL'

      N__FILTER_METHOD = JP__DIGITAL_FILTER
      
      N__NUM_BANDS = 5

      N__BAND_SIZE(:) = 0

      IF (LECMWF) THEN
        N__BAND_SIZE(1:N__NUM_BANDS) =(/ 193, 15, 116, 4, 15 /)
      ELSE
        N__BAND_SIZE(1:N__NUM_BANDS) =(/ 148, 15, 116, 4, 15 /)
      ENDIF

      IF (MAXVAL(N__BAND_SIZE(:)) > JP__MAX_CHANNELS) &
&              CALL ABOR1('Too many channels specified in cloud '//&
&                      'detection - increase JP__Max_Channels')
 
      N__BANDS(:,:)= 0 

      ! Use the "IASI 366" Subset
      IF (LECMWF) THEN
        N__BANDS(1:N__BAND_SIZE(1),1) = &
&        (/    16,   38,   49,   51,   55,   57,   59,   61,   63,   66, &
&              70,   72,   74,   79,   81,   83,   85,   87,   89,   92, &
&              95,   97,   99,  101,  104,  106,  109,  111,  113,  116, &
&             119,  122,  125,  128,  131,  133,  135,  138,  141,  144, &
&             146,  148,  151,  154,  157,  159,  161,  163,  165,  167, &
&             170,  173,  176,  178,  179,  180,  183,  185,  187,  189, &
&             191,  193,  195,  197,  199,  201,  203,  205,  207,  210, &
&             212,  214,  217,  219,  222,  224,  226,  228,  230,  232, &
&             234,  236,  239,  241,  242,  243,  246,  249,  252,  254, &
&             256,  258,  260,  262,  265,  267,  269,  271,  272,  273, &
&             275,  278,  280,  282,  284,  286,  288,  290,  292,  294, &
&             296,  299,  301,  303,  306,  308,  310,  312,  314,  316, &
&             318,  320,  323,  325,  327,  329,  331,  333,  335,  337, &
&             339,  341,  343,  345,  347,  350,  352,  354,  356,  358, &
&             360,  362,  364,  366,  369,  371,  373,  375,  377,  379, &
&             381,  383,  386,  389,  398,  401,  404,  407,  410,  414, &
&             416,  426,  428,  432,  434,  439,  445,  457,  515,  546, &
&             552,  559,  566,  571,  573,  646,  662,  668,  756,  867, &
&             921, 1027, 1090, 1133, 1191, 1194, 1271, 1805, 1884, 1946, &
&            1991, 2094, 2239 /)
      ELSE
        N__BANDS(1:N__BAND_SIZE(1),1) = &
&        (/    16,   38,   49,   51,   55,   57,   59,   61,   63,   66, &
&              70,   72,   74,   79,   81,   83,   85,   87,   89,   92, &
&              95,   97,   99,  101,  104,  106,  109,  111,  113,  116, &
&             119,  122,  125,  128,  131,  133,  135,  138,  141,  144, &
&             146,  148,  151,  154,  157,  159,  161,  163,        167, &
&             170,  173,  176,        179,  180,        185,  187,       &
&                   193,              199,              205,  207,  210, &
&             212,  214,  217,  219,  222,  224,  226,        230,  232, &
&                   236,  239,        242,  243,  246,  249,  252,  254, &
&                         260,  262,  265,  267,  269,                   &
&             275,        280,  282,                                294, &
&             296,  299,        303,  306,                               &
&                         323,        327,  329,              335,       &
&                               345,  347,  350,        354,  356,       &
&             360,              366,        371,  373,  375,  377,  379, &
&             381,  383,  386,  389,  398,  401,  404,  407,  410,  414, &
&             416,  426,  428,  432,  434,  439,  445,  457,  515,  546, &
&             552,  559,  566,  571,  573,  646,  662,  668,  756,  867, &
&             921, 1027,       1133, 1191, 1194, 1271, 1805, 1884,       &
&            1991, 2094, 2239 /)
      ENDIF

      N__BANDS(1:N__BAND_SIZE(2),2) = &
&      (/ 1479, 1509, 1513, 1521, 1536, 1574, 1579, 1585, 1587, 1626, &
&         1639, 1643, 1652, 1658, 1671  /)

      N__BANDS(1:N__BAND_SIZE(3),3) = &
&      (/ 2119, 2213, 2271, 2321, 2398, 2701, 2741, 2819, 2889, 2907, 2910, &
&         2919, 2939, 2944, 2948, 2951, 2958, 2977, 2985, 2988, 2991, &
&         2993, 3002, 3008, 3014, 3027, 3029, 3036, 3047, 3049, 3053, &
&         3058, 3064, 3069, 3087, 3093, 3098, 3105, 3107, 3110, 3127, &
&         3136, 3151, 3160, 3165, 3168, 3175, 3178, 3207, 3228, 3244, &
&         3248, 3252, 3256, 3263, 3281, 3303, 3309, 3312, 3322, 3375, &
&         3378, 3411, 3438, 3440, 3442, 3444, 3446, 3448, 3450, 3452, &
&         3454, 3458, 3467, 3476, 3484, 3491, 3497, 3499, 3504, 3506, &
&         3509, 3518, 3527, 3555, 3575, 3577, 3580, 3582, 3586, 3589, &
&         3599, 3653, 3658, 3661, 4032, 5368, 5371, 5379, 5381, 5383, &
&         5397, 5399, 5401, 5403, 5405, 5455, 5480, 5483, 5485, 5492, &
&         5502, 5507, 5509, 5517, 5558  /)

      N__BANDS(1:N__BAND_SIZE(4),4) = &
&      (/   5988, 5992, 5994, 6003  /)


      N__BANDS(1:N__BAND_SIZE(5),5) = &
&      (/  6982, 6985, 6987, 6989, 6991, 6993, 6995, 6997, 7267, 7269, &
&          7424, 7426, 7428, 7885, 8007 /)  
   
      N__GRADCHKINTERVAL(:) = 0
      N__GRADCHKINTERVAL(1:N__NUM_BANDS) = (/ 5,5,5,5,5 /)

      N__WINDOW_WIDTH(:) = 0
      N__WINDOW_WIDTH(1:N__NUM_BANDS) = (/ 10,6,8,5,8 /)

      N__WINDOW_BOUNDS(:,:) = 0
      N__WINDOW_BOUNDS(1,1) = 573
      N__WINDOW_BOUNDS(1,2) = 2239

      R__BT_THRESHOLD(:) = 0.
      R__BT_THRESHOLD(1:N__NUM_BANDS) = (/ 0.5, 0.5, 0.5, 0.5, 0.5/)

      R__GRAD_THRESHOLD(:) = 0.
      R__GRAD_THRESHOLD(1:N__NUM_BANDS) = (/ 0.02, 0.02, 0.02, 0.02, 0.02 /)

      R__WINDOW_GRAD_THRESHOLD(:) = 0.
      R__WINDOW_GRAD_THRESHOLD(1) = 0.4

      L__DO_QUICK_EXIT = .TRUE.

      ! This is cross-band:
      
      L__DO_CROSSBAND = .TRUE.

      N__BANDTOUSE(:) = 0
      N__BANDTOUSE(1:N__NUM_BANDS) = (/ 1,1,1,1,1 /)

      ! This is the setup for imager cloud detection

      L__DO_IMAGER_CLOUD_DETECTION = .TRUE.

      N__NUM_IMAGER_CHANS = 2
      N__NUM_IMAGER_CLUSTERS = 7

      N__IMAGER_CHANS(:) = 0
      N__IMAGER_CHANS(1:N__NUM_IMAGER_CHANS) = (/ 2, 3 /)

      R__STDDEV_THRESHOLD(:) = 0.0
      R__STDDEV_THRESHOLD(1:N__NUM_IMAGER_CHANS) = (/ 0.75, 0.80 /)

      R__COVERAGE_THRESHOLD = 0.03
      R__FG_DEPARTURE_THRESHOLD = 1.0

      ! This is aerosol:
      IF (LECMWF) THEN
        L__DO_AEROSOLDETECTION = .TRUE.
      ELSE
        L__DO_AEROSOLDETECTION = .FALSE.
      ENDIF

      N__NUM_AEROSOL_TESTS = 1
      N__NUM_AEROSOL_CHANS(:) = 0
      N__NUM_AEROSOL_CHANS(1:N__NUM_AEROSOL_TESTS) = (/ 5 /)
      N__AEROSOL_CHANS(:,:) = 0
      N__AEROSOL_CHANS(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/ 1340, 1782, 2348, 2356, 7425 /)
      R__AEROSOL_ABSCISSAE(:,:) = 0.0
      R__AEROSOL_ABSCISSAE(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/ 1.70747, 1.58687, 1.47237, 1.34777, 1.23624, 0.518735, &
&              0.508225, 0.522959 /)
      R__AEROSOL_THRESHOLD(:) = 0.0
      R__AEROSOL_THRESHOLD(1:N__NUM_AEROSOL_TESTS) = (/ 0.4 /)
      R__AEROSOLMINNORM(:) = 0.0
      R__AEROSOLMINNORM(1:N__NUM_AEROSOL_TESTS) = (/ 2.0 /)

   CASE(INST_ID_CRIS)
      !====================
      ! Set up CRIS
      !====================

      CL__INSTRUMENTNAME='CRIS'
      CL__CLOUD_DETECTION_FILE = 'CRIS_CLDDET.NL'

      N__FILTER_METHOD = JP__DIGITAL_FILTER
      
      N__NUM_BANDS = 4

      N__BAND_SIZE(:) = 0

      IF (LECMWF) THEN
        N__BAND_SIZE(1:N__NUM_BANDS) =(/ 126, 53, 108, 13 /)
      ELSE
        N__BAND_SIZE(1:N__NUM_BANDS) =(/ 157, 53, 108, 13 /)
      ENDIF

      IF (MAXVAL(N__BAND_SIZE(:)) > JP__MAX_CHANNELS) &
&              CALL ABOR1('Too many channels specified in cloud '//&
&                      'detection - increase JP__Max_Channels')

      N__BANDS(:,:)= 0

      IF (LECMWF) THEN
      ! Use the "CRIS 300" Subset
      N__BANDS(1:N__BAND_SIZE(1),1) = &
&    (/    1,    5,    9,   13,   17,   21,   25,   29,   33,   37, &
&         41,   45,   49,   53,   57,   61,   65,   69,   73,   77, &
&         81,   85,   89,   93,   97,  101,  105,  109,  113,  117, &
&        121,  125,  129,  133,  137,  141,  145,  149,  153,  157, &
&        161,  165,  169,  173,  177,  181,  185,  189,  193,  197, &
&        201,  205,  209,  213,  217,  221,  225,  229,  233,  237, &
&        241,  245,  249,  253,  257,  261,  265,  269,  273,  277, &
&        281,  285,  289,  293,  297,  301,  305,  309,  313,  317, &
&        321,  325,  329,  333,  337,  341,  345,  349,  353,  357, &
&        361,  365,  369,  373,  377,  381,  385,  389,  393,  397, &
&        401,  405,  409,  413,  417,  421,  425,  429,  433,  437, &
&        441,  445,  449,  453,  457,  461,  465,  469,  473,  477, &
&        481,  485,  489,  493,  497,  501 /)

      ELSE

      ! Use the MF "CRIS 331" Subset
      N__BANDS(1:N__BAND_SIZE(1),1) = &
&    (/    1,    3,    5,    7,    9,   11,   13,   15,   17,   19, &
&         21,   23,   25,   27,   29,   31,   33,   35,   37,   39, &
&         41,   43,   45,   47,   49,   51,   53,   55,   57,   59, &
&         61,   63,   65,   67,   69,   71,   73,   75,   77,   79, &
&         81,   83,   85,   87,   89,   91,   93,   95,   97,   99, &
&        101,  103,  105,  107,  109,  111,  113,  117, &
&        121,  125,  127,  129,  131,  133,  137,  141,  145,  147, &
&        149,  153,  157, &
&        161,  165,  169,  173,  177,  181,  185,  189,  193,  197, &
&        201,  205,  209,  213,  217,  221,  225,  229,  233,  237, &
&        241,  245,  249,  253,  257,  261,  265,  269,  273,  277, &
&        281,  285,  289,  293,  297,  301,  305,  309,  313,  317, &
&        321,  325,  329,  333,  337,  341,  345,  349,  353,  357, &
&        361,  365,  369,  373,  377,  381,  385,  389,  393,  397, &
&        401,  405,  409,  413,  417,  421,  425,  429,  433,  437, &
&        441,  445,  449,  453,  457,  461,  465,  469,  473,  477, &
&        481,  485,  489,  493,  497,  501 /)

      ENDIF

      N__BANDS(1:N__BAND_SIZE(2),2) = &
&    (/  505,  509,  513,  517,  521,  525,  529,  533,  537,  541, &
&        545,  549,  553,  557,  561,  565,  569,  573,  577,  581, &
&        585,  589,  593,  597,  601,  605,  609,  613,  617,  621, &
&        625,  629,  633,  637,  641,  645,  649,  653,  657,  661, &
&        665,  669,  673,  677,  681,  685,  689,  693,  697,  701, &
&        705,  709,  713 /)

      N__BANDS(1:N__BAND_SIZE(3),3) = &
&    (/  717,  721,  725,  729,  733,  737,  741,  745,  749,  753, &
&        757,  761,  765,  769,  773,  777,  781,  785,  789,  793, &
&        797,  801,  805,  809,  813,  817,  821,  825,  829,  833, &
&        837,  841,  845,  849,  853,  857,  861,  865,  869,  873, &
&        877,  881,  885,  889,  893,  897,  901,  905,  909,  913, &
&        917,  921,  925,  929,  933,  937,  941,  945,  949,  953, &
&        957,  961,  965,  969,  973,  977,  981,  985,  989,  993, &
&        997, 1001, 1005, 1009, 1013, 1017, 1021, 1025, 1029, 1033, &
&       1037, 1041, 1045, 1049, 1053, 1057, 1061, 1065, 1069, 1073, &
&       1077, 1081, 1085, 1089, 1093, 1097, 1101, 1105, 1109, 1113, &
&       1117, 1121, 1125, 1129, 1133, 1137, 1141, 1145 /)

      N__BANDS(1:N__BAND_SIZE(4),4) = &
&    (/ 1149, 1153, 1157, 1161, 1165, 1169, 1173, 1177, 1181, 1185, &
&       1189, 1193, 1251 /)


      N__GRADCHKINTERVAL(:) = 0
      N__GRADCHKINTERVAL(1:N__NUM_BANDS) = (/ 5,5,5,5 /)

      N__WINDOW_WIDTH(:) = 0
      N__WINDOW_WIDTH(1:N__NUM_BANDS) = (/ 10,6,8,8 /)

      N__WINDOW_BOUNDS(:,:) = 0
      N__WINDOW_BOUNDS(1,1) = 229
      N__WINDOW_BOUNDS(1,2) = 549

      R__BT_THRESHOLD(:) = 0.
      R__BT_THRESHOLD(1:N__NUM_BANDS) = (/ 0.5, 0.5, 0.5, 0.5/)

      R__GRAD_THRESHOLD(:) = 0.
      R__GRAD_THRESHOLD(1:N__NUM_BANDS) = (/ 0.02, 0.02, 0.02, 0.02 /)

      R__WINDOW_GRAD_THRESHOLD(:) = 0.
      R__WINDOW_GRAD_THRESHOLD(1) = 0.4

      L__DO_QUICK_EXIT = .TRUE.

      ! This is cross-band:

      L__DO_CROSSBAND = .TRUE.

      N__BANDTOUSE(:) = 0
      N__BANDTOUSE(1:N__NUM_BANDS) = (/ 1,1,1,1 /)

      ! This is the setup for imager cloud detection

      L__DO_IMAGER_CLOUD_DETECTION = .FALSE.
      N__NUM_IMAGER_CHANS = 0
      N__NUM_IMAGER_CLUSTERS = 0
      N__IMAGER_CHANS(:) = 0
      R__STDDEV_THRESHOLD(:) = 0.0
      R__COVERAGE_THRESHOLD = 0.0
      R__FG_DEPARTURE_THRESHOLD = 0.0

      ! This is aerosol:
      L__DO_AEROSOLDETECTION = .FALSE.

      N__NUM_AEROSOL_TESTS = 1
      N__NUM_AEROSOL_CHANS(:) = 0
      N__NUM_AEROSOL_CHANS(1:N__NUM_AEROSOL_TESTS) = (/ 8 /)
      N__AEROSOL_CHANS(:,:) = 0
      N__AEROSOL_CHANS(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/  704,  707,  709,  714,  728,  730,  732,  734 /)
      R__AEROSOL_ABSCISSAE(:,:) = 0.0
      R__AEROSOL_ABSCISSAE(1:N__NUM_AEROSOL_CHANS(1),1) = &
&           (/ 1.70747, 1.58687, 1.47237, 1.34777, 1.23624, 0.518735, &
&              0.508225, 0.522959 /)
      R__AEROSOL_THRESHOLD(:) = 0.0
      R__AEROSOL_THRESHOLD(1:N__NUM_AEROSOL_TESTS) = (/ 0.4 /)
      R__AEROSOLMINNORM(:) = 0.0
      R__AEROSOLMINNORM(1:N__NUM_AEROSOL_TESTS) = (/ 2.0 /)

   CASE DEFAULT
      CYCLE
   END SELECT

   !------------------------------------------------------------------
   ! Open and read file containing cloud detection setup for the 
   ! current instrument
   !------------------------------------------------------------------
   IOMASTER=1
   MASTERPROCESSOR : IF(MYPROC == IOMASTER) THEN
      INIU1=NULUSR3
      WRITE(NULOUT,*)'READING CLOUD DETECTION FILE FOR ',CL__InstrumentName
      OPEN(INIU1,STATUS='OLD',FORM='FORMATTED',FILE=CL__Cloud_Detection_File, &
&                   IOSTAT=IOS)
      IF (IOS == 0) THEN
         READ(INIU1,NML=NAMCLDDET,IOSTAT=IOS)
         IF (IOS == 0) THEN
            WRITE(NULOUT,*) CL__INSTRUMENTNAME,' CLOUD DETECTION FILE READ OK'
         ELSE
            WRITE(NULOUT,*)'PROBLEM READING '//CL__InstrumentName//&
&                'CLOUD DETECTION FILE: Using Default Values'
         ENDIF
         CLOSE(INIU1)
      ELSE
         WRITE(NULOUT,*)'NO '//CL__InstrumentName//&
&              ' CLOUD DETECTION FILE : Using Default Values'
      ENDIF
      
      IF (MAXVAL(N__BAND_SIZE(:)) > JP__MAX_CHANNELS) &
&              CALL ABOR1('Too many channels specified in cloud '//&
&                      'detection - increase JP__Max_Channels')
      
      
      M__SENSOR = J__SENSOR
      
      !------------------------------------------------------------------
      ! Broadcast values if using multiple PEs
      !------------------------------------------------------------------

      IF (NPROC > 1) THEN
         !   Load buffers used in MPL
         IBUF(1,1)   = M__SENSOR
         IBUF(1,2)   = N__FILTER_METHOD 
         IBUF(1,3)   = N__NUM_BANDS    
         IF (L__DO_QUICK_EXIT) THEN
            IBUF(1,4)   = 1
         ELSE
            IBUF(1,4)   = 0
         ENDIF
         IF (L__DO_CROSSBAND) THEN
            IBUF(1,5)   = 1
         ELSE
            IBUF(1,5)   = 0
         ENDIF
         IF (L__DO_AEROSOLDETECTION) THEN
            IBUF(1,6)   = 1
         ELSE
            IBUF(1,6)   = 0
         ENDIF
         IBUF(1,7) = N__NUM_AEROSOL_TESTS  
         IF (L__DO_IMAGER_CLOUD_DETECTION) THEN
            IBUF(1,8)   = 1
         ELSE
            IBUF(1,8)   = 0
         ENDIF
         IBUF(2,:) = N__BAND_SIZE 
         IBUF(3,:) = N__WINDOW_WIDTH 
         IBUF(4,:) = N__GRADCHKINTERVAL
         IBUF(5,:) = N__BANDTOUSE
         IBUF(6,:) = N__NUM_AEROSOL_CHANS
         IBUF(7,:) = N__WINDOW_BOUNDS(:,1)
         IBUF(8,:) = N__WINDOW_BOUNDS(:,2)
         IBUF(9,1) = N__NUM_IMAGER_CHANS
         IBUF(9,2) = N__NUM_IMAGER_CLUSTERS
         IBUF(10,:) = N__IMAGER_CHANS(:)

         ZBUF(1,:)                = R__BT_THRESHOLD  
         ZBUF(2,:)                = R__GRAD_THRESHOLD 
         ZBUF(3,:)                = R__AEROSOL_THRESHOLD 
         ZBUF(4,:)                = R__AEROSOLMINNORM
         ZBUF(5,:)                = R__WINDOW_GRAD_THRESHOLD 
         ZBUF(6,:)                = R__STDDEV_THRESHOLD(:)
         ZBUF(7,1)                = R__COVERAGE_THRESHOLD
         ZBUF(7,2)                = R__FG_DEPARTURE_THRESHOLD
      ENDIF
      
   ENDIF MASTERPROCESSOR
   
   
   IF(NPROC > 1) THEN
      ITAG=5555
      CALL MPL_BROADCAST(  &
&          IBUF,           &
&          KROOT=IOMASTER, &
&          KTAG=ITAG,      &
&          CDSTRING='CLOUD_DETECT_SETUP:')

      ITAG=5556
      CALL MPL_BROADCAST(  &
&          ZBUF,           &
&          KROOT=IOMASTER, &
&          KTAG=ITAG,      &
&          CDSTRING='CLOUD_DETECT_SETUP:')
      
      !   Unload buffers used in MPL
      M__SENSOR           = IBUF(1,1)
      N__FILTER_METHOD    = IBUF(1,2)
      N__NUM_BANDS        = IBUF(1,3)
      IF (IBUF(1,4) == 1) THEN
         L__DO_QUICK_EXIT = .TRUE.
      ELSE
         L__DO_QUICK_EXIT = .FALSE.
      ENDIF
      IF (IBUF(1,5) == 1) THEN
         L__DO_CROSSBAND = .TRUE.
      ELSE
         L__DO_CROSSBAND = .FALSE.
      ENDIF
      IF (IBUF(1,6) == 1) THEN
         L__DO_AEROSOLDETECTION = .TRUE.
      ELSE
         L__DO_AEROSOLDETECTION = .FALSE.
      ENDIF
      N__NUM_AEROSOL_TESTS  = IBUF(1,7)  
      N__BAND_SIZE          = IBUF(2,:)
      N__WINDOW_WIDTH       = IBUF(3,:)
      N__GRADCHKINTERVAL    = IBUF(4,:)
      N__BANDTOUSE          = IBUF(5,:)  
      N__NUM_AEROSOL_CHANS  = IBUF(6,:)
      N__WINDOW_BOUNDS(:,1) = IBUF(7,:)
      N__WINDOW_BOUNDS(:,2) = IBUF(8,:)
      IF (IBUF(1,8) == 1) THEN
         L__DO_IMAGER_CLOUD_DETECTION = .TRUE.
      ELSE
         L__DO_IMAGER_CLOUD_DETECTION = .FALSE.
      ENDIF
      N__NUM_IMAGER_CHANS    = IBUF(9,1)
      N__NUM_IMAGER_CLUSTERS = IBUF(9,2)
      N__IMAGER_CHANS(:)     = IBUF(10,:)

      R__BT_THRESHOLD      = ZBUF(1,:)
      R__GRAD_THRESHOLD    = ZBUF(2,:)
      R__AEROSOL_THRESHOLD = ZBUF(3,:)
      R__AEROSOLMINNORM    = ZBUF(4,:)
      R__WINDOW_GRAD_THRESHOLD = ZBUF(5,:)
      R__STDDEV_THRESHOLD(:)    = ZBUF(6,:)
      R__COVERAGE_THRESHOLD     = ZBUF(7,1)
      R__FG_DEPARTURE_THRESHOLD = ZBUF(7,2)

      ITAG=5557
      DO III=1,N__NUM_BANDS
      CALL MPL_BROADCAST(  &
&          N__BANDS(1:MAXVAL(N__BAND_SIZE(:)),III), &
&          KROOT=IOMASTER, &
&          KTAG=ITAG,      &
&          CDSTRING='CLOUD_DETECT_SETUP:')
      ENDDO


      IF (L__DO_AEROSOLDETECTION) THEN
         I__MAXAEROCHANS = MAXVAL(N__NUM_AEROSOL_CHANS(:))

         ITAG=5558
         DO III=1,N__NUM_AEROSOL_TESTS
         CALL MPL_BROADCAST(                                             &
&            N__AEROSOL_CHANS(1:I__MAXAEROCHANS,III), &
&            KROOT=IOMASTER,                                             &
&            KTAG=ITAG,                                                  &
&            CDSTRING='CLOUD_DETECT_SETUP:')
         ENDDO

         ITAG=5559
         DO III=1,N__NUM_AEROSOL_TESTS
         CALL MPL_BROADCAST(                                                 &
&            R__AEROSOL_ABSCISSAE(1:I__MAXAEROCHANS,III), &
&            KROOT=IOMASTER,                                                 &
&            KTAG=ITAG,                                                      &
&            CDSTRING='CLOUD_DETECT_SETUP:')
         ENDDO
      ENDIF


   ENDIF

   !------------------------------------------------------------------
   ! Set up the S__Cloud_Detect_Setup structure for current sensor
   !------------------------------------------------------------------

   S__CLOUD_DETECT_SETUP(J__SENSOR) % M__SENSOR = M__SENSOR
   
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__FILTER_METHOD = N__FILTER_METHOD
   
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__NUM_BANDS = N__NUM_BANDS
   
   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BAND_SIZE(N__NUM_BANDS) )
   
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BAND_SIZE(:) = &
&             N__BAND_SIZE(1:N__NUM_BANDS)
   
   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BANDS & 
&            (MAXVAL(N__BAND_SIZE(:)), N__NUM_BANDS))
      
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BANDS(:,:) = 0

   DO J = 1, N__NUM_BANDS
      S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BANDS(1:N__BAND_SIZE(J),J) = &
&          N__BANDS(1:N__BAND_SIZE(J),J)
      IF (MYPROC == 1) THEN
         WRITE(NULOUT,*) 'Sensor ',M__Sensor,' Band ',J,' has',&
&          N__Band_Size(J),' Cloud Detection Channels :'
         WRITE(NULOUT,*) N__BANDS(1:N__BAND_SIZE(J),J)
      ENDIF
   ENDDO
      
   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % N__WINDOW_WIDTH(N__NUM_BANDS) )
      
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__WINDOW_WIDTH(:) = &
&        N__WINDOW_WIDTH(1:N__NUM_BANDS) 
      
      
   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % R__BT_THRESHOLD(N__NUM_BANDS) )
   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__BT_THRESHOLD(:) = &
&        R__BT_THRESHOLD(1:N__NUM_BANDS) 
      
   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % R__GRAD_THRESHOLD(N__NUM_BANDS))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__GRAD_THRESHOLD(:) = &
&        R__GRAD_THRESHOLD(1:N__NUM_BANDS) 
      
   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        R__WINDOW_GRAD_THRESHOLD(N__NUM_BANDS))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__WINDOW_GRAD_THRESHOLD(:) = &
&        R__WINDOW_GRAD_THRESHOLD(1:N__NUM_BANDS) 
      
   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % N__GRADCHKINTERVAL(N__NUM_BANDS))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__GRADCHKINTERVAL(:) = &
&        N__GRADCHKINTERVAL(1:N__NUM_BANDS) 
      
   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % N__WINDOW_BOUNDS(N__NUM_BANDS,2))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__WINDOW_BOUNDS(:,:) = &
&        N__WINDOW_BOUNDS(1:N__NUM_BANDS,:) 
      
   S__CLOUD_DETECT_SETUP(J__SENSOR) % L__DO_QUICK_EXIT = L__DO_QUICK_EXIT
   
   !-------------
   ! Cross Band
   !-------------

   S__CLOUD_DETECT_SETUP(J__SENSOR) % L__DO_CROSSBAND = L__DO_CROSSBAND
      
   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BANDTOUSE(N__NUM_BANDS) )
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__BANDTOUSE(:) = &
&        N__BANDTOUSE(1:N__NUM_BANDS) 

   !-------------
   ! Imager cloud detection
   !-------------

   S__CLOUD_DETECT_SETUP(J__SENSOR) % L__DO_IMAGER_CLOUD_DETECTION = &
&        L__DO_IMAGER_CLOUD_DETECTION

   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__NUM_IMAGER_CHANS = &
&        N__NUM_IMAGER_CHANS

   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__NUM_IMAGER_CLUSTERS = &
&        N__NUM_IMAGER_CLUSTERS

   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        N__IMAGER_CHANS(N__NUM_IMAGER_CHANS))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__IMAGER_CHANS(:) = &
&        N__IMAGER_CHANS(1:N__NUM_IMAGER_CHANS)

   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        R__STDDEV_THRESHOLD(N__NUM_IMAGER_CHANS))
   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__STDDEV_THRESHOLD(:) = &
&        R__STDDEV_THRESHOLD(1:N__NUM_IMAGER_CHANS)

   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__COVERAGE_THRESHOLD = &
&        R__COVERAGE_THRESHOLD

   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__FG_DEPARTURE_THRESHOLD = &
&        R__FG_DEPARTURE_THRESHOLD

   !-------------
   ! Aerosol
   !-------------

   S__CLOUD_DETECT_SETUP(J__SENSOR) % L__DO_AEROSOLDETECTION = &
&        L__DO_AEROSOLDETECTION

   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__NUM_AEROSOL_TESTS = &
&        N__NUM_AEROSOL_TESTS

   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        N__NUM_AEROSOL_CHANS(N__NUM_AEROSOL_TESTS) )

   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__NUM_AEROSOL_CHANS(:) = &
&        N__NUM_AEROSOL_CHANS(1:N__NUM_AEROSOL_TESTS) 

   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % N__AEROSOL_CHANS & 
&            (MAXVAL(N__NUM_AEROSOL_CHANS(:)), N__NUM_AEROSOL_TESTS))

   S__CLOUD_DETECT_SETUP(J__SENSOR) % N__AEROSOL_CHANS(:,:) = 0
   DO J = 1, N__NUM_AEROSOL_TESTS
      S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&          N__AEROSOL_CHANS(1:N__NUM_AEROSOL_CHANS(J),J) = &
&          N__AEROSOL_CHANS(1:N__NUM_AEROSOL_CHANS(J),J)
   ENDDO

   ALLOCATE(S__CLOUD_DETECT_SETUP(J__SENSOR) % R__AEROSOL_ABSCISSAE & 
&            (MAXVAL(N__NUM_AEROSOL_CHANS(:)), N__NUM_AEROSOL_TESTS))

   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__AEROSOL_ABSCISSAE(:,:) = 0.0
   DO J = 1, N__NUM_AEROSOL_TESTS
      S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&          R__AEROSOL_ABSCISSAE(1:N__NUM_AEROSOL_CHANS(J),J) = &
&          R__AEROSOL_ABSCISSAE(1:N__NUM_AEROSOL_CHANS(J),J)
   ENDDO

   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        R__AEROSOL_THRESHOLD(N__NUM_AEROSOL_TESTS) )

   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__AEROSOL_THRESHOLD(:) = &
&        R__AEROSOL_THRESHOLD(1:N__NUM_AEROSOL_TESTS) 

   ALLOCATE( S__CLOUD_DETECT_SETUP(J__SENSOR) % &
&        R__AEROSOLMINNORM(N__NUM_AEROSOL_TESTS) )

   S__CLOUD_DETECT_SETUP(J__SENSOR) % R__AEROSOLMINNORM(:) = &
&        R__AEROSOLMINNORM(1:N__NUM_AEROSOL_TESTS) 

ENDDO SENSORLOOP

IF (LHOOK) CALL DR_HOOK('CLOUD_DETECT_SETUP',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUD_DETECT_SETUP
