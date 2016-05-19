/*************************************************************************************
			       DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
					   		     IMPERIAL COLLEGE LONDON 

 				      EE 3.19: Real Time Digital Signal Processing
					       Dr Paul Mitcheson and Daniel Harvey

					       Eusebius Ngemera and Prahnav Sharma

				        		  LAB 5: IIR filter

 				            ********* I N T I O. C **********

	Performs digital filtering using an included IIR filter from a text file

 *************************************************************************************
 				Updated for use on 6713 DSK by Danny Harvey: May-Aug 2006
				Updated for CCS V4 Sept 10
				Modified to do digital filtering using an IIR filter
 ************************************************************************************/
/*
 *	You should modify the code so that interrupts are used to service the 
 *  audio port.
 */
/**************************** Pre-processor statements ******************************/

#include <stdlib.h>
//  Included so program can make use of DSP/BIOS configuration tool.  
#include "dsp_bios_cfg.h"

/* The file dsk6713.h must be included in every program that uses the BSL.  This 
   example also includes dsk6713_aic23.h because it uses the 
   AIC23 codec module (audio interface). */
#include "dsk6713.h"
#include "dsk6713_aic23.h"

// math library (trig functions)
#include <math.h>

// Some functions to help with writing/reading the audio ports when using interrupts.
#include <helper_functions_ISR.h>

#include "iir_coef.txt"
// include MATLAB-generated IIR Filter coefficients file
// #define M (length of a[]) if different from N
// #define N (order+1)
// double a[] = {...}
// double b[] = {...}

/******************************** Global variables **********************************/
double x[N];
double y[M];

// dynamically-allocated array
double* z;

/******************************* Global declarations ********************************/

/* Audio port configuration settings: these values set registers in the AIC23 audio 
   interface to configure it. See TI doc SLWS106D 3-3 to 3-10 for more info. */
DSK6713_AIC23_Config Config = { \
			 /**********************************************************************/
			 /*   REGISTER	            FUNCTION			      SETTINGS         */ 
			 /**********************************************************************/\
    0x0017,  /* 0 LEFTINVOL  Left line input channel volume  0dB                   */\
    0x0017,  /* 1 RIGHTINVOL Right line input channel volume 0dB                   */\
    0x01f9,  /* 2 LEFTHPVOL  Left channel headphone volume   0dB                   */\
    0x01f9,  /* 3 RIGHTHPVOL Right channel headphone volume  0dB                   */\
    0x0011,  /* 4 ANAPATH    Analog audio path control       DAC on, Mic boost 20dB*/\
    0x0000,  /* 5 DIGPATH    Digital audio path control      All Filters off       */\
    0x0000,  /* 6 DPOWERDOWN Power down control              All Hardware on       */\
    0x0043,  /* 7 DIGIF      Digital audio interface format  16 bit                */\
    0x008d,  /* 8 SAMPLERATE Sample rate control             8 KHZ                 */\
    0x0001   /* 9 DIGACT     Digital interface activation    On                    */\
			 /**********************************************************************/
};


// Codec handle:- a variable used to identify audio interface  
DSK6713_AIC23_CodecHandle H_Codec;

 /******************************* Function prototypes ********************************/
void init_hardware(void);     
void init_HWI(void);
void ISR_AIC(void);
// IIR implementations
double simple_IIR(short new_sample);
double DF2_IIR(short new_input);
double DF2T_IIR(short new_input);
/********************************** Main routine ************************************/
void main(){      

	// initialize board and the audio port
 	init_hardware();

 	// dynamic memory allocation, with zero initiliasing
	z = (double*)calloc(N, sizeof(double));
	memset(z, 0, sizeof(z));
	
 	// initialize hardware interrupts
	init_HWI();

	// loop indefinitely, waiting for interrupts  					
	while(1) 
	{};
}
        
/********************************** init_hardware() **********************************/  
void init_hardware()
{
    // Initialize the board support library, must be called first 
    DSK6713_init();
    
    // Start the AIC23 codec using the settings defined above in config 
    H_Codec = DSK6713_AIC23_openCodec(0, &Config);

	/* Function below sets the number of bits in word used by MSBSP (serial port) for 
	receives from AIC23 (audio port). We are using a 32 bit packet containing two 
	16 bit numbers hence 32BIT is set for  receive */
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);
	
	/* Configures interrupt to activate on each consecutive available 32 bits 
	from Audio port hence an interrupt is generated for each L & R sample pair */	
	MCBSP_FSETS(SPCR1, RINTM, FRM);

	/* These commands do the same thing as above but applied to data transfers to  
	the audio port */
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);	
	MCBSP_FSETS(SPCR1, XINTM, FRM);	
	
}

/********************************** init_HWI() **************************************/  
void init_HWI(void)
{
	IRQ_globalDisable();			// Globally disables interrupts
	IRQ_nmiEnable();				// Enables the NMI interrupt (used by the debugger)
	IRQ_map(IRQ_EVT_RINT1,4);		// Maps an event to a physical interrupt
	IRQ_enable(IRQ_EVT_RINT1);		// Enables the event
	IRQ_globalEnable();				// Globally enables interrupts

} 

/*************************** INTERRUPT SERVICE ROUTINE ******************************/  
void ISR_AIC(void)
{
	// reads mono-channel input as type Int16 and outputs the filtered output
	// as type Int16 in the mono-channel output
	 
	mono_write_16Bit( (short) DF2T_IIR( mono_read_16Bit() ));
	// in-line coding is used to reduce time spent in the Interrupt
}

/************************** IIR Filter Implementations *****************************/  
// IIR assuming a and b are of different lengths, a[M] and b[N]
double simple_IIR(short new_sample)
{
	int i;					// for-loop iteration index
	double sum = 0.0;		// data type matches IIR coefficients
	
	// shifting x, the past inputs
	for (i = N-1; i > 0; i--)
	{
		// move data one index up, discarding last (oldest) input
		x[i] = x[i-1];
	}
	// prepend new sample, to be the first element
	x[0] = new_sample;
	
	// Convolution between b coefficients and past inputs
	for (i = 0; i < N; i++)
	{ 
		sum += b[i] * x[i];
	}
	
	// Convolution between a coefficients and past outputs
	// a[0] = 1, convolution starts from a[1]
	for (i = 1; i < M; i++)
	{ 
		sum -= a[i] * y[i];
	}
	
	// shifting y, the past outputs
	for (i = M-1; i > 0; i--)
	{
		// move data one index up, discarding last (oldest) output
		y[i] = y[i-1];
	}
	
	// prepend final output value to y buffer
	// y index matches a coefficient array
	y[1] = sum;
	
	return sum;
}

// Direct-Form 2 IIR filtering
double DF2_IIR(short new_input)
{
	int i;
	double output = 0.0;
	
	// z[0] is the top, intermediate node
	z[0] = new_input;
	
	// (accumulated) convolution of the left branch
	for (i = 1; i < N; i++)
	{
		z[0] -= a[i] * z[i];
	}
	
	output = b[0] * z[0];
	
	// (accumulated) convolution of the right branch
	for (i = 1; i < N; i++)
	{ 
		output += b[i] * z[i];
	}
	
	// shifting z down
	for (i = N-1; i > 0; i--)
	{
		z[i] = z[i-1];
	}
	
	return output;
}

// Direct-Form 2 transposed IIR filtering
double DF2T_IIR(short new_input)
{
	int i;
	
	double output = (new_input * b[0]) + z[0];
	
	// iterate from right to left
	// last summation does not have z[N-1], and is an edge case
	for (i = 1; i < N-1; i++)
	{
		z[i-1] = z[i] + (new_input * b[i]) - (output * a[i]);
	}
	z[N-2] = (new_input * b[N-1]) - (output * a[N-1]);
	
	return output;
}
\end{minted}



