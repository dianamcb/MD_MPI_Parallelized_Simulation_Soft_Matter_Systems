# include <cstdlib>
# include <cstdio>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
void compute (double w, int np, int nd, double pos[], double vel[], double mass, double f[], double &pot, double &kin, float tau, float temp_0, double &temp, double &press, int step );
double cpu_time ( );
double dist (double w, int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
void timestamp ( );
void update ( int np, int nd, double pos[], double vel[], double f[], double acc[], double mass, double &w, double dt, float tau, float tau_0, double beta, float temp_0, float press_0, double temp, double press, int step );
void print_pos_xyz ( int np, int nd, double pos[], int step, int actu );
void periodical_conditions ( int np, int nd, double w, double pos[] );
int sign ( double a );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MD.
//
//  Discussion:
//
//    MD implements a simple molecular dynamics simulation.
//
//    The velocity Verlet time integration scheme is used. 
//
//    The particles interact with a central pair potential.
//
//    This program is based on a FORTRAN90 program by Bill Magro.
//
//  Usage:
//
//    md nd np step_num dt
//
//    where
//
//    * nd is the spatial dimension (2 or 3);
//    * np is the number of particles (500, for instance);
//    * step_num is the number of time steps (500, for instance).
//    * dt is the time step (0.1 for instance)
//
//
{
  double *acc;
  double ctime;
  double dt;
  double e0;
  double *force;
  double kinetic;
  double mass = 1;
  float sigma = 1;//Unit of length
  float epsilon = 1;//Unit of energy
  int nd = 3;
  int np;
  double *pos;
  double potential;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;
  int actu;
  bool os;
  float diam;
  double w;
  double temperature; //Termostato.
  double pressure; //Barostato.
  float temperature_0; //Termostato.
  float pressure_0; //Barostato.
  float tau; //Termostato.
  float tau_p; //Barostato.
  float beta; //Barostato.
  
  system( "mkdir data" );
  
  timestamp ( );
  cout << "\n";
  cout << "MD\n";
  cout << "  C++ version\n";
  cout << "  A molecular dynamics program.\n";
//
//  Get the number of particles.
//
  if ( 1 < argc )
  {
    np = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter NP, the number of particles (500, for instance).\n";
    cin >> np;
  }
//
//  Get the number of time steps.
//
  if ( 2 < argc )
  {
    step_num = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter STEP_NUM, the number of time steps (500 or 1000, for instance).\n";
    cin >> step_num;
  }
//
//  Get the time step.
//
  if ( 3 < argc )
  {
    dt = atof ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter DT, the time step size (0.1, for instance).\n";
    cin >> dt;
  }
//
//  Get the cell dimentions.
//
  if ( 4 < argc )
  {
    dt = atof ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter W, the cell width.\n";
	cin >> w;
  }
//
//  Get the number of steps to print positions.
//  Consider that steps most be a multiple of actu.
//
  if ( 5 < argc )
  {
    dt = atof ( argv[5] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter ACTU, the number of steps to pass through before printing positions.\n";
    cin >> actu;
  }
//
//  Get the coupling factor.
//
  if ( 6 < argc )
  {
    dt = atof ( argv[6] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter TAU, the coupling factor.\n"; //Termostato.
    cin >> tau;
  }
//
//  Get the temperature of the heat bath.
//
  if ( 7 < argc )
  {
    dt = atof ( argv[7] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter TEMPERATURE_0, the temperature of the heat bath.\n"; //Termostato.
    cin >> temperature_0;
  }
//
//  Get the time constant for the coupling.
//
  if ( 8 < argc )
  {
    dt = atof ( argv[8] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter TAU_P, the time constant for the coupling to pressure bath.\n"; //Barostato.
    cin >> tau_p;
  }
//
//  Get the pressure of the pressure bath.
//
  if ( 9 < argc )
  {
    dt = atof ( argv[9] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter PRESSURE_0, the pressure of the bath.\n"; //Barostato.
    cin >> pressure_0;
  }
//
//  Get the isotermal compressibility.
//
  if ( 10 < argc )
  {
    dt = atof ( argv[10] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter BETA, isotermal compressibility.\n"; //Barostato.
    cin >> beta;
  }
//
//  Report.
//
  cout << "\n";
  cout << "  ND, the spatial dimension, is " << nd << "\n";
  cout << "  NP, the number of particles in the simulation is " << np << "\n";
  cout << "  STEP_NUM, the number of time steps, is " << step_num << "\n";
  cout << "  DT, the size of each time step, is " << dt << "\n";
  cout << "  W, the cell width, is " << w << "\n";
  cout << "  ACTU, the number of steps to pass through, is " << actu << "\n";
  cout << "  TAU, the coupling factor, is " << tau << "\n"; //Termostato.
  cout << "  TAU_P, the time constant for the coupling, is " << tau_p << "\n"; //Barostato.
  cout << "  BETA, isotermal compressibility, is " << beta << "\n"; //Barostato.
  cout << "  TEMPERATURE_0, the temperature of the heat bath, is " << temperature_0 << "\n"; //Termostato.
  cout << "  PRESSURE_0, the pressure of the bath, is " << pressure_0 << "\n"; //Barostato.
  cout << "  The density is " << np * 1.0 / pow( w, 3 ) << "\n";
//
//  Allocate memory.
//
  acc = new double[nd*np];
  force = new double[nd*np];
  pos = new double[nd*np];
  vel = new double[nd*np];
//
//    This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//
  cout << "\n";
  cout << "  At each step, we report the potential and kinetic energies.\n";
  cout << "  The sum of these energies should be a constant.\n";
  cout << "  As an accuracy check, we also print the relative error\n";
  cout << "  in the total energy.\n";
  cout << "\n";
  cout << "      Step      Potential       Kinetic        Temperature        Pressure        (P+K-E0)/E0\n"; //Termostato. Barostato.
  cout << "                Energy P        Energy K       T                  p               Relative Energy Error\n";
  cout << "\n";

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;
  
  ctime = cpu_time ( );

  for ( step = 0; step <= step_num; step++ )
  {
    if ( step == 0 ) initialize ( np, nd, pos, vel, acc );
    else update ( np, nd, pos, vel, force, acc, mass, w, dt, tau, tau_p, beta, temperature_0, pressure_0, temperature, pressure, step );
    
	periodical_conditions ( np, nd, w, pos );
	
    compute (w, np, nd, pos, vel, mass, force, potential, kinetic, tau, temperature_0, temperature, pressure, step );
	
	if ( step%actu == 0 ) print_pos_xyz(np, nd, pos, step, actu);
	
    if ( step == 0 ) e0 = potential + kinetic;

    if ( step == step_print )
    {
      cout << "  " << setw(8) << step
           << "  " << setw(14) << potential
           << "  " << setw(14) << kinetic
           << "  " << setw(14) << temperature //Termostato.
           << "  " << setw(14) << pressure //Barostato.
           << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 
           << "  " << setw(14) << ( potential + kinetic) << "\n";
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }
  }
//
//  Report timing.
//
  ctime = cpu_time ( ) - ctime;
  cout << "\n";
  cout << "  Elapsed cpu time " << ctime << " seconds.\n";
//
//  Free memory.
//
  delete [] acc;
  delete [] force;
  delete [] pos;
  delete [] vel;
//
//  Terminate.
//
  cout << "\n";
  cout << "MD\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
  
  return 0;
}
//****************************************************************************80

void compute ( double w, int np, int nd, double pos[], double vel[], double mass,
  double f[], double &pot, double &kin, float tau, float temp_0, 
  double &temp,  double &press, int step )

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces, energies and temperature.
//
//  Discussion:
//
//    The computation of forces and energies is fully parallel.
//
//    The potential function V(X) is a harmonic well which smoothly
//    saturates to a maximum value at PI/2:
//
//      v(x) = ( sin ( min ( x, PI2 ) ) )**2
//
//    The derivative of the potential is:
//
//      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
//            = sin ( 2.0 * min ( x, PI2 ) )
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double POS[ND*NP], the position of each particle.
//
//    Input, double VEL[ND*NP], the velocity of each particle.
//
//    Input, double MASS, the mass of each particle.
//
//    Input, double W, the cell dimensions.
//
//    Output, double F[ND*NP], the forces.
//
//    Output, double &POT, the total potential energy.
//
//    Output, double &KIN, the total kinetic energy.
//
//    Output, double &PRESS, the system pressure.
//
//    Input, float TAU, the coupling factor. //Termostato.
//
//    Input, float TEMP_0, the temperature of the heat bath. //Termostato.
//
//    Output, double &TEMP, the temperature of the system. //Termostato.
//
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];
  double vir;

  pot = 0.0;
  kin = 0.0;
  vir = 0.0;
  press = 0.0;

//
//  Compute the kinetic energy and the temperature. //Termostato.
//
  for ( k = 0; k < np; k++ )
  {
  	for ( i = 0; i < nd; i++ )
  	{
  	  kin += vel[i+k*nd] * vel[i+k*nd];
	}
  }
  kin *= 0.5 * mass;
  temp = kin * 2 / ( 3.0 * np );

//
//  Compute the potential energy and forces.
//
  for ( k = 0; k < np; k++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist (w, nd, pos+k*nd, pos+j*nd, rij );
//
//  Attribute half of the potential energy to particle J.
//
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pot += 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < nd; i++ )
        {
          if ( step < 1 )
		  {
		    f[i+k*nd] -= rij[i] * sin ( 2.0 * d2 ) / d;
	      }
	      else
	      {
	        f[i+k*nd] -= rij[i] * sin ( 2.0 * d2 ) / d - ( temp_0 / temp - 1 ) * mass * vel[i+k*nd] / ( 2 * tau ); //Termostato.
	      }
	      if ( k < j ) vir += rij[i] * f[i+k*nd]; //Barostato.
        }
      }
    }
  }

//
//  Compute pressure.
//
  vir *= -0.5; //Barostato.
  if (step > 0) press = 2 * ( kin - vir ) / ( 3.0 * w * w * w ); //Barostato.
  
  return;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double dist (double w, int nd, double r1[], double r2[], double dr[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIST computes the displacement (and its norm) between two particles.
//
//
//  Parameters:
//    
//    Input, double W, the cell lenght.
//    
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double R1[ND], R2[ND], the positions.
//
//    Output, double DR[ND], the displacement vector.
//
//    Output, double D, the Euclidean norm of the displacement.
//
{
  double d;
  int i;
  double r2i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {	  
    r2i = r2[i];
	if ( fabs ( r1[i] - r2[i] ) > w*1.0/2 ) r2i += w*sign ( r1[i] - r2[i] );//Mínima imagen
    dr[i] = r1[i] - r2i;
    d += dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
//****************************************************************************80

void initialize ( int np, int nd, double pos[], double vel[], double acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes the positions, velocities, and accelerations.
//
//    
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Output, double POS[ND*NP], the positions.
//
//    Output, double VEL[ND*NP], the velocities.
//
//    Output, double ACC[ND*NP], the accelerations.
//
{
  int i;
  int j;
  int seed;
//
//  Set the positions.
//
  seed = 123456789;
  
  r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed, pos );
//
//  Set the velocities.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
//
//  Set the accelerations.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }
  
  return;
}
//****************************************************************************80

void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  References:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

  void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double &w, double dt, float tau, 
  float tau_p, double beta, float temp_0, float press_0, double temp, 
  double press, int step )

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE updates positions, velocities and accelerations.
//
//  Discussion:
//
//    The time integration is fully parallel.
//
//    A velocity Verlet algorithm is used for the updating.
//
//    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
//    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
//    a(t+dt) = f(t) / m
//
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input/output, double POS[ND*NP], the positions.
//
//    Input/output, double VEL[ND*NP], the velocities.
//
//    Input, double F[ND*NP], the forces.
//
//    Input/output, double ACC[ND*NP], the accelerations.
//
//    Input, double MASS, the mass of each particle.
//
//    Output, double &W, the cell dimensions.
//
//    Input, double DT, the time step.
//
//    Input, float TAU, the coupling factor. //Termostato.
//
//    Input, float TAU_P, the time constant coupling. //Barostato.
//
//    Input, double BETA, the isotermal compressibility. //Barostato.
//
//    Input, float TEMP_0, the temperature of the heat bath. //Termostato.
//
//    Input, float PRESS_0, the pressure of the bath. //Barostato.
//
//    Input, double TEMP, the temperature of the system. //Termostato.
//
//    Input, double PRESS, the system pressure. //Barostato.
//
{
  int i;
  int j;
  double rmass;
  double lambda;
  double alpha; //Barostato.
  double mu; //Barostato.

  rmass = 1.0 / mass;
  
  if ( step > 1 )
  {
    lambda = sqrt( 1 + dt / tau * ( temp_0 / temp - 1 ) ); //Termostato.
    alpha = -beta * ( press_0 - press )/( 3 * tau_p ); //Barostato.
    mu = 1 + dt * alpha; //Barostato.
    w *= mu; //Barostato.
  }

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      if ( step <= 1 )
	  {
	  	pos[i+j*nd] += vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt; 
	  	vel[i+j*nd] += 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
	  } 
      else
	  {
		pos[i+j*nd] += ( vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt ) * mu; //Barostato.
		vel[i+j*nd] += 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] ) * lambda; //Termostato.
	  } 
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }
  
  return;
}
//****************************************************************************80

  void print_pos_xyz( int np, int nd, double pos[], int step, int actu ) 

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_POS_XYZ print the position of each particle in xyz format.
//
//	Parameters:
//
//		Input/output, int NP, the number of particles.
//
//		Input, int ND, the number of spatial dimentions.
//
//		Input/putput, double POS[NP*ND], the positions.
//
//		Input, int STEP, the time step.
//
//		Input, int ACTU, the number of steps to pass through.
{
  FILE *dat;
  char nombre[20];

  sprintf( nombre, "data/position.xyz" );	
  if( step == 0 ) dat = fopen ( nombre, "w" );
  else dat = fopen ( nombre, "a" );
  
  fprintf ( dat, "%i\n", np );
  fprintf ( dat, "xyz\n" );  
  for ( int j = 1; j <= np; j++ )
  {
    for ( int i = 1; i <= nd; i += nd )
    {
	  if ( nd == 2 ) fprintf ( dat, "C\t%e\t%e\t0\n", pos[i+j*nd], pos[(i+1)+j*nd] );
	  else fprintf ( dat, "C\t%e\t%e\t%e\n", pos[i+j*nd], pos[(i+1)+j*nd], pos[(i+2)+j*nd] );
    }
  }
  
  fclose ( dat );
  
  return;
}
//****************************************************************************80

  int sign ( double a ) 

//****************************************************************************80
//
//  Purpose:
//
//    SIGN is the common sign function.
//
//	Parameters:
//
//	  Input, double a.
{
  if ( a>0 ) return ( 1 );
  else if ( a<0 ) return ( -1 );
  else return ( 0 );
}
//****************************************************************************80

  void periodical_conditions ( int np, int nd, double w, double pos[] ) 

//****************************************************************************80
//
//  Purpose:
//
//    PERIODICAL_CONDITIONS establishes the particles displacements in a 
//    simulation cell.
//
//	Parameters:
//
//	  Input, int NP, the number of particles.
//
//	  Input, int ND, the number of spatial dimentions.
//
//	  Input, double W, the cell lenght.
//
//	  Input/output, double POS[NP*ND], the positions.
{
  for ( int j = 1; j <= np; j++ )
  {
    for ( int i = 1; i <= nd; i += nd )
    {
      if ( ( pos[i+j*nd] < 0 ) || ( pos[i+j*nd] >= w ) )
      {
        pos[i+j*nd] += w*sign( -pos[i+j*nd] );
      }
      if ( ( pos[(i+1)+j*nd] < 0 ) || ( pos[(i+1)+j*nd] >= w ) )
      {
        pos[(i+1)+j*nd] += w*sign ( -pos[(i+1)+j*nd] );
      }
      if ( ( pos[(i+2)+j*nd] < 0 ) || ( pos[(i+2)+j*nd] >= w ) )
      {
        pos[(i+2)+j*nd] += w*sign ( -pos[(i+2)+j*nd] );
      }
    }
  }
  
  return;
}
