#include <fstream>
#include <cstring>
#include <iostream>
#include <limits>
#include <cmath>

// fltk includes
#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <FL/Fl_Hor_Value_Slider.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Choice.H>

using namespace std;

const int N = 100;

float* DRAWME;
float U[N+1];
float U_prev[N+1];
float U_new[N+1];

static float sim_dx;
static float sim_dt;
static float sim_speed;

void init_wave()
{
  DRAWME = &U[0];
  // initial wave condition
  //U[0] = 0; U[N] = 0;
  for ( int i = 0; i < N+1; ++i ) {
      if ( ( i * sim_dx > 0.2 )
           && ( i * sim_dx < 0.35 ) )
      U[i] = 1.0f;
    else
      U[i] = 0.2; //numeric_limits<float>::epsilon(); //0.0f;
    U_prev[i] = U[i];
  }

}


/*! Solve the wave equation. Numerically using explicit temporal
  differencing, and centered spatial differencing (ETCS)

  wave equation is:
  u_tt = c^2 * u_xx
*/
void wave1d_etcs( float wavespeed, float dx, float dt )
{
  /* 
     u_tt = ( c^2 / dx^2 ) * ( u_x+1 - 2 u_x + u_x-1 )
     C = c^2 / dx^2

     u_t+1 = dt^2 * C * ( u_x+1 - 2 u_x + u_x-1 ) + 2 u_x - u_x,t-1
  */
  // solve
  float C = dt * dt * wavespeed * wavespeed / ( dx * dx );
  for ( int i = 1; i < N; i++ )
    U_new[i] = C * ( U[i+1] - 2 * U[i] + U[i-1] ) + 2 * U[i] - U_prev[i];
  
  // update solution
  for ( int i = 1; i < N; i++ ) {
    U_prev[i] = U[i];
    U[i] = U_new[i];
  }

}


/*! Solve the wave equation. Numerically using explicit temporal
  differencing, and centered spatial differencing (ETCS)

  wave equation is:
  u_tt = c^2 * u_xx
*/
void wave1d_lax( float wavespeed, float dx, float dt )
{
  /* 
     u_tt = ( c^2 / dx^2 ) * ( u_x+1 - 2 u_x + u_x-1 )
     C = c^2 / dx^2

     u_t+1 = dt^2 * C * ( u_x+1 - 2 u_x + u_x-1 ) + 2 u_x - u_x,t-1
  */

  for( int i = 1; i < N; i++ )
    U_new[i] = U[i+i] + U[i-1] - U_prev[i];
}

/*! Solve the advection equation. Numerically using explicit temporal
  differencing, and centered spatial differencing (ETCS)

  advection equation is:

  u_t = -C * u_x
*/
void advect1d_etcs( float wavespeed, float dx, float dt )
{
  float C = wavespeed * dt / ( dx * 2.0f );
  for ( int i = 1; i < N; i++ )
    U_new[i] = U[i] - C * ( U[i+1] - U[i-1] );

  // update
  for ( int i = 1; i < N; i++ )
    U[i] = U_new[i];
}


/*! Solve the advection equation. Numerically using Lax scheme

advection equation is:

u_t = -C * u_x
*/
void advect1d_lax( float wavespeed, float dx, float dt )
{
  float C = wavespeed * dt / ( dx * 2.0f );
  for ( int i = 1; i < N; i++ )
    U_new[i] = 0.5 * ( U[i+1] + U[i-1] ) - C * ( U[i+1] - U[i-1] );

  // update
  for ( int i = 1; i < N; i++ )
    U[i] = U_new[i];
}


/*! Solve the advection equation. Numerically using upwind differencing

advection equation is:

u_t = -C * u_x
*/
void advect1d_upwind( float wavespeed, float dx, float dt )
{
  float C = wavespeed * dt / dx;
  for ( int i = 1; i < N; i++ )
    U_new[i] = U[i] - C * ( U[i] - U[i-1] );

  // update
  for ( int i = 1; i < N; i++ )
    U[i] = U_new[i];
}


float U_step1[N+1];
float U_step2[N+2];
void advect1d_bfecc( float wavespeed, float dx, float dt )
{
  /* 1. forward step
     2. backwards step
     3. error estimate
     4. corrected final step
  */

  const float C = wavespeed * dt / dx;

  // #1
  for( int i = 1; i < N; ++i )
    U_step1[i] = ( 1 - C ) * U[i] + C * U_prev[i-1];

  // #2
  for( int i = 1; i < N; ++i )
    U_step2[i] = ( 1 - C ) * U_step1[i] + C * U_step1[i+1];

  // #3
  for( int i = 1; i < N; ++i )
    U_step2[i] = U[i] - 0.5 * ( U_step2[i] - U[i] );

  for( int i = 1; i < N; ++i )
    U_prev[i] = U[i];
  
  // #4
  for( int i = 1; i < N; ++i )
    U[i] = ( 1 - C ) * U_step2[i] + C * U_step2[i-1];
}


void advect1d_bfecc_mod( float wavespeed, float dx, float dt )
{
  /* 1. forward step
     2. backwards step
     3. error estimate
     4. corrected final step
  */

  const float C = wavespeed * dt / dx;

  // #1
  for( int i = 1; i < N; ++i )
    U_step1[i] = ( 1 - C ) * U[i] + C * U_prev[i-1];

  // #2
  for( int i = 1; i < N; ++i )
    U_step2[i] = ( 1 - C ) * U_step1[i] + C * U_step1[i+1];

  // #3
  for( int i = 1; i < N; ++i )
    U_step2[i] = U[i] - 0.5 * ( U_step2[i] - U[i] );

  for( int i = 1; i < N; ++i )
    U_prev[i] = U[i];
  
  // #4
  float trial, min, max;
  for( int i = 1; i < N; ++i ) {
      if( U[i] < U[i-1] ) {
          min = U[i];
          max = U[i-1];
      }
      else {
          min = U[i-1];
          max = U[i];
      }
      trial = ( 1 - C ) * U_step2[i] + C * U_step2[i-1];
      trial = std::max( min, std::min( trial, max ) );
      U_new[i] = trial;
  }

  for( int i = 1; i < N; ++i )
    U[i] = U_new[i];
}


void advect1d_maccormack( float wavespeed, float dx, float dt )
{
  /* 1. forward step
     2. backwards step
     3. error estimate
     4. corrected final step
  */

  const float C = wavespeed * dt / dx;

  // #1
  for( int i = 1; i < N; ++i )
    U_step1[i] = ( 1 - C ) * U[i] + C * U_prev[i-1];

  // #2
  for( int i = 1; i < N; ++i )
    U_step2[i] = ( 1 - C ) * U_step1[i] + C * U_step1[i+1];

  for( int i = 1; i < N; ++i )
    U_prev[i] = U[i];
  
  // #4
  for( int i = 1; i < N; ++i )
    U[i] = U_step1[i] - 0.5 * ( U_step2[i] - U[i] );
}


void advect1d_maccormack_clamped( float wavespeed, float dx, float dt )
{
  /* 1. forward step
     2. backwards step
     3. error estimate
     4. corrected final step
  */

  const float C = wavespeed * dt / dx;

  // #1
  for( int i = 1; i < N; ++i )
    U_step1[i] = ( 1 - C ) * U[i] + C * U_prev[i-1];

  // #2
  for( int i = 1; i < N; ++i )
    U_step2[i] = ( 1 - C ) * U_step1[i] + C * U_step1[i+1];

  for( int i = 1; i < N; ++i )
    U_prev[i] = U[i];
  
  // #4
  float trial, min, max;
  for( int i = 1; i < N; ++i ) {
      trial = U_step1[i] - 0.5 * ( U_step2[i] - U[i] );
      
      if( U[i] < U[i-1] ) {
          min = U[i];
          max = U[i-1];
      }
      else {
          min = U[i-1];
          max = U[i];
      }

      trial = std::max( min, std::min( trial, max ) );

      U_new[i] = trial;
  }

  memcpy( &U[0], &U_new[0], sizeof(float) * N );
  
}


float PHI[N+1];

/*!
  x1 = x0 
 */
void advect_steinhoff( float wavespeed, float dx, float dt )
{
  float stepsf = fabs( ceil( wavespeed * dt / dx ) );
  unsigned int steps = stepsf;
  dt /= stepsf;
  
  float C = wavespeed * dt / ( dx * 2.0f );
  float dxsqr_inv = 1.0f / (dx * dx);
  
  float mu = 0.15; //1e-10;
  float beta = 1e-10;
  float eps = 0.325; //mu + beta * dx * dx;

  cout << "Taking: " << stepsf << " steps\r" << flush;

  for( unsigned int step = 0; step < steps; ++step ) {
      // compute the sort of non-linear mean
      for( int i = 1; i < N; ++i )
        PHI[i] = 3.f/ ( 1/U[i+1] + 1/U[i] + 1/U[i-1] );

#define RHO(i) (mu * U[i] - eps * PHI[i])
  
      for ( int i = 1; i < N; i++ )
        U_new[i] = U[i] - C * ( U[i+1] - U[i-1] )
                   + ( RHO(i+1) - 2.f * RHO(i) + RHO(i-1) ) * dxsqr_inv;

#undef RHO
  
      // update
      memcpy( &U[0], &U_new[0], sizeof(float) * N );
  }

  //DRAWME = &PHI[0];
}

namespace {
  enum{ WaveETCS = 0,
        WaveLax,
        AdvectETCS,
        AdvectLax,
        AdvectUpwind,
        AdvectBfecc,
        AdvectBfeccMod,
        AdvectMacCormack,
        AdvectMacCormack_Mod,
        AdvectSteinhoff,
  };

  static unsigned int SolverType = AdvectETCS;
}

class SignalWin : public Fl_Gl_Window
{
  public:
    SignalWin( const int samples, float* data,
               int x, int y, int w, int h,
               const char* label = "SignalWin" )
        : Fl_Gl_Window( x, y, w, h, label ),
          my_numsamples( samples ),
          my_data( data )
      { }
    
    void draw()
      {
        if ( ! valid() ) {
            initGL();
            valid( 1 );
        }
        glClear( GL_COLOR_BUFFER_BIT );

        // draw a connected segment for each sample
        float g = 0;
        switch( SolverType ) {
          case WaveETCS:
              wave1d_etcs( sim_speed, sim_dx, sim_dt );
              break;
          case WaveLax:
              wave1d_lax( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectETCS:
              advect1d_etcs( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectLax:
              advect1d_lax( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectUpwind:
              advect1d_upwind( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectBfecc:
              advect1d_bfecc( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectBfeccMod:
              advect1d_bfecc_mod( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectMacCormack:
              advect1d_maccormack( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectMacCormack_Mod:
              advect1d_maccormack_clamped( sim_speed, sim_dx, sim_dt );
              break;
          case AdvectSteinhoff:
              advect_steinhoff( sim_speed, sim_dx, sim_dt );
              break;
        };

    
        glColor3f( 1.0, g, 0.0 );
        glBegin( GL_LINE_STRIP );
        for ( int i = 0; i < my_numsamples; i++ )
          glVertex2f( float(i), DRAWME[i] );
          //glVertex2f( float(i), my_data[i] );
        glEnd();

      }
    
    void initGL()
      {
        glColor3f( 1.0, 0.0, 0.0 );
        glClearColor( 0., 0., 0., 0. );
        
        glViewport ( 0, 0, w(), h() );
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();
        
        glOrtho( -1, my_numsamples + 1, // left, right
                 -1.1, 1.1,             // bottom, top
                 -1., 1. );             // near, far
  
        glMatrixMode (GL_MODELVIEW);
      }

    int handle( int e )
      {
        switch( e ) {
          case FL_SHORTCUT:
              if ( Fl::event_key() == 'c' ) {
                  init_wave();
                  return 1;
              }
        };

        return Fl_Gl_Window::handle( e );
      }
    
  private:
    int my_numsamples;
    const float* my_data;
};

void reinit_sim()
{
  // initial wave condition
  U[0] = 0; U[N] = 0;
  for ( int i = 1; i < N; ++i ) {
    if ( i * sim_dx < 0.2 )
      U[i] = 1.0f;
    else
      U[i] = 0.0f;
    U_prev[i] = U[i];
  }
}

void simtimer_cb( void* drawable )
{
  //advect1d_upwind( sim_speed, sim_dx, sim_dt );
  static_cast<SignalWin*>( drawable )->redraw();
  
  //Fl::repeat_timeout( 1.0/24.0, simtimer_cb, drawable );
  Fl::repeat_timeout( 1.0/60.0, simtimer_cb, drawable );
  //Fl::repeat_timeout( 0, simtimer_cb, drawable );
}

void wavespeed_cb( Fl_Widget* sld, void* dataDraw )
{
  sim_speed = float( static_cast<Fl_Hor_Value_Slider*>( sld )->value() );
}

void timestep_cb( Fl_Widget* sld, void* dataDraw )
{
  sim_dt = float( static_cast<Fl_Hor_Value_Slider*>( sld )->value() );
}


void solver_choice_cb( Fl_Widget* w, void* v )
{
  init_wave();
  SolverType = (unsigned int) v;
}

Fl_Menu_Item solver_pulldown[] = {
    { "Wave - ETCS", 0, solver_choice_cb, (void*)WaveETCS },
    { "Wave - Lax",  0, solver_choice_cb, (void*)WaveLax },
    { "Advection - ETCS", 0, solver_choice_cb, (void*)AdvectETCS },
    { "Advection - Lax", 0, solver_choice_cb, (void*)AdvectLax },
    { "Advection - Upwind", 0, solver_choice_cb, (void*)AdvectUpwind },
    { "Advection - BFECC", 0, solver_choice_cb, (void*)AdvectBfecc },
    { "Advection - BFECC Clamped", 0, solver_choice_cb, (void*)AdvectBfeccMod },
    { "Advection - MacCormack", 0, solver_choice_cb, (void*)AdvectMacCormack },
    { "Advection - MacCormack Clamped", 0, solver_choice_cb, (void*)AdvectMacCormack_Mod },
    { "Advection - Steinhoff", 0, solver_choice_cb, (void*)AdvectSteinhoff },
    { 0 }
};

int main( int argc, char** argv )
{
  sim_dx = 1.0f / N;
  sim_dt = 0.001;
  sim_speed  = 0.2; // wave speed m/s

  init_wave();
  // setup data display
  Fl_Window mainWin( 410, 475 );
  SignalWin dataDraw( N+1, U, 5, 5, 400, 400, "1D Wave" );
  int ypos = 410;

  Fl::add_timeout( 1.0/24.0, simtimer_cb, &dataDraw );

  Fl_Choice choiceMenu( 90, ypos, 150, 20, "Solvers:" );
  choiceMenu.menu( solver_pulldown );
  choiceMenu.value( SolverType );
  ypos += 25;
  
  Fl_Hor_Value_Slider speedSld( 90, ypos, 310, 15, "Wavespeed:" );
  speedSld.bounds( 0, 0.5 );
  speedSld.step( 1e-2 );
  speedSld.value( sim_speed );
  speedSld.align( FL_ALIGN_LEFT );
  speedSld.callback( wavespeed_cb, &dataDraw );
  ypos += 20;

  Fl_Hor_Value_Slider timestepSld( 90, ypos, 310, 15, "Timestep:" );
  timestepSld.bounds( 0, 0.01 );
  timestepSld.step( 1e-5 );
  timestepSld.value( sim_dt );
  timestepSld.align( FL_ALIGN_LEFT );
  timestepSld.callback( timestep_cb, &dataDraw );

  mainWin.end();
  mainWin.show( argc, argv );

  Fl::run();
  return 0;
}


