typedef struct _LaserList  {
   int mode;
   double k0;
   double amplitude;   //unit is a0.
   double rU;
   double rD;
   double flat;

   double rayleighLength;
   double beamWaist;
   double focus;
   int direction;

   struct _LaserList *next;
} LaserList;

typedef struct _Laser  {
   double real;
   double img;
} Laser;
