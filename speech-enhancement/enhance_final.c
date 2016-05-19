/*************************************************************************************
			       DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
					   		     IMPERIAL COLLEGE LONDON 

 				      EE 3.19: Real Time Digital Signal Processing
					       Dr Paul Mitcheson and Daniel Harvey
					       Prahnav Sharma and Eusebius Ngemera

				        		 PROJECT: Frame Processing

 				            ********* ENHANCE. C **********
							 Shell for speech enhancement 

  		Utilises overlap-add frame processing (interrupt driven) on the DSK. 
  		Enhances speech in noisy signal through several enhancements

 *************************************************************************************
 				             By Danny Harvey: 21 July 2006
							 Updated for use on CCS v4 Sept 2010
							 Added core funtionality to skeleton Mar 2016
 ************************************************************************************/
//  library required when using calloc
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

/* Some functions to help with Complex algebra and FFT. */
#include "cmplx.h"      
#include "fft_functions.h"  

// Some functions to help with writing/reading the audio ports when using interrupts.
#include <helper_functions_ISR.h>

#define WINCONST 0.85185			/* 0.46/0.54 for Hamming window */
#define FSAMP 8000.0		/* sample frequency, ensure this matches Config for AIC */
#define FFTLEN 256					/* fft length = frame length 256/8000 = 32 ms*/
#define NFREQ (1+FFTLEN/2)			/* number of frequency bins from a real FFT */
#define OVERSAMP 4					/* oversampling ratio (2 or 4) */  
#define FRAMEINC (FFTLEN/OVERSAMP)	/* Frame increment */
#define CIRCBUF (FFTLEN+FRAMEINC)	/* length of I/O buffers */

#define OUTGAIN 50000.0				/* Output gain for DAC */
#define INGAIN  (1.0/16000.0)		/* Input gain for ADC  */
// PI defined here for use in your code 
#define PI 3.141592653589793
#define TFRAME FRAMEINC/FSAMP       /* time between calculation of each frame */

#define min(a,b) (((a) < (b)) ? (a):(b))
#define max(a,b) (((a) > (b)) ? (a):(b))
#define cmplx_min(a,b) (((cabs(a)) < cabs(b)) ? (a):(b))
#define cmplx_max(a,b) (((cabs(a)) > cabs(b)) ? (a):(b))
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
    0x008d,  /* 8 SAMPLERATE Sample rate control        8 KHZ-ensure matches FSAMP */\
    0x0001   /* 9 DIGACT     Digital interface activation    On                    */\
			 /**********************************************************************/
};

// Codec handle:- a variable used to identify audio interface  
DSK6713_AIC23_CodecHandle H_Codec;

float *inbuffer, *outbuffer;   		/* Input/output circular buffers */
float *inframe, *outframe;          /* Input and output frames */
float *inwin, *outwin;              /* Input and output windows */

// float arrays for G(w) and |X(w)|
float *mag, *original_mag;

// output Y(w)=X(w)*|G(w)|

// complex arrays to hold FFT's and IFFT's
// enhancement 8 involves history of |Y(w)|
complex *complex_array, *complex_array_now, *complex_array_prev;

// Bins of minimums in frequency-domain magnitude
float *M1, *M2, *M3, *M4;

float *noise; 		// noise estimate magnitude, |N(w)|

// Low-pass filtered versions
float *lpf_input, *lpf_noise;

float *noise_mag;	// simply a pointer to the noise buffer to use

float ingain, outgain;				/* ADC and DAC gains */ 
float cpufrac; 						/* Fraction of CPU time used */
volatile int io_ptr=0;              /* Input/ouput pointer for circular buffers */
volatile int frame_ptr=0;           /* Frame pointer */

int detection_period = 10;	// period for detection of minimum noise
int frame_count = 0;		// keep track of how many frames processed for current Minimum bin

// Parameters
float lambda = 0.01;
float alpha = 4;
float tau = 0.030;			// milliseconds
float k_pole;				// Low-pass filter pole, exp(-TFRAME/tau)

// enhancement 6
float alpha_high = 3;		// alpha scaling for low SNR
float iSNR_threshold = 6;	// threshold for alpha scaling; simply |N(w)|/|X(w)|

// enhancement 8
float *SNR_now;				// array of the previous SNR values
float noise_threshold = 1;	// |X(w)|/|N(w)| threshold 
int i_threshold = 3;

// Extra enhancement; e_speech_gain
// FFT integer indices
int low_freq = 15;
int high_freq = 35;
int cut_off = 41;
float speech_gain = 4.0;

// enhancement switches
int allpass = 0;
int e1 = 1;
int e2 = 0;
int e3 = 1;
int e4 = 1; // 0 to 4
int e5 = 0; // 0 to 4
int e6 = 1;
int e8 = 0;
int e_speech_gain = 0;

 /******************************* Function prototypes *******************************/
void init_hardware(void);    	/* Initialize codec */ 
void init_HWI(void);            /* Initialize hardware interrupts */
void ISR_AIC(void);             /* Interrupt service routine for codec */
void process_frame(void);       /* Frame processing routine */

void update_minimums(float* input);	// update M bins and estimate noise
float complex_mag(complex input);	// custom implementation of cabs(.)
/********************************** Main routine ************************************/
void main()
{      
  	int k; // used in various for loops
  	
  	k_pole = exp(-TFRAME/tau);	// Low-Pass Filter pole
  
	/*  Initialize and zero fill arrays */  

	inbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Input array */
    outbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Output array */
	inframe		= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    outframe	= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    inwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Input window */
    outwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Output window */
    
    mag			= (float *) calloc(FFTLEN, sizeof(float));
    original_mag= (float *) calloc(FFTLEN, sizeof(float));

	complex_array = (complex *) calloc(FFTLEN, sizeof(complex));
	complex_array_prev = (complex *) calloc(FFTLEN, sizeof(complex));
	complex_array_now = (complex *) calloc(FFTLEN, sizeof(complex));
	
	// Minimum bins
	M1			= (float *) calloc(FFTLEN, sizeof(float));
	M2			= (float *) calloc(FFTLEN, sizeof(float));
	M3			= (float *) calloc(FFTLEN, sizeof(float));
	M4			= (float *) calloc(FFTLEN, sizeof(float));
	noise		= (float *) calloc(FFTLEN, sizeof(float));
	
	lpf_input 	= (float *) calloc(FFTLEN, sizeof(float));
	lpf_noise	= (float *) calloc(FFTLEN, sizeof(float));
	SNR_now 	= (float *) calloc(FFTLEN, sizeof(float));
	
	noise_mag = noise;
	
	/* initialize board and the audio port */
  	init_hardware();
  
  	/* initialize hardware interrupts */
  	init_HWI();    
  
	/* initialize algorithm constants */  
                       
  	for (k=0; k<FFTLEN; k++)
	{                           
		inwin[k] = sqrt((1.0-WINCONST*cos(PI*(2*k+1)/FFTLEN))/OVERSAMP);
		outwin[k] = inwin[k]; 
	}

  	ingain=INGAIN;
  	outgain=OUTGAIN;        
 							
  	/* main loop, wait for interrupt */  
  	while(1) 	process_frame();
}
    
/********************************** init_hardware() *********************************/  
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

	/* These commands do the same thing as above but applied to data transfers to the 
	audio port */
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
        
/******************************** process_frame() ***********************************/  
void process_frame(void)
{
	int k, m, kk; 
	int io_ptr0;
	float factor, lambda_factor;
	complex *temp_ca;

	/* work out fraction of available CPU time used by algorithm */    
	cpufrac = ((float) (io_ptr & (FRAMEINC - 1)))/FRAMEINC;  
		
	/* wait until io_ptr is at the start of the current frame */ 	
	while((io_ptr/FRAMEINC) != frame_ptr); 
	
	/* then increment the framecount (wrapping if required) */ 
	if (++frame_ptr >= (CIRCBUF/FRAMEINC)) frame_ptr=0;
 	
 	/* save a pointer to the position in the I/O buffers (inbuffer/outbuffer) where the 
 	data should be read (inbuffer) and saved (outbuffer) for the purpose of processing */
 	io_ptr0=frame_ptr * FRAMEINC;
	
	/* copy input data from inbuffer into inframe (starting from the pointer position) */ 
	 
	m=io_ptr0;
    for (k=0;k<FFTLEN;k++)
	{                           
		inframe[k] = inbuffer[m] * inwin[k];
		complex_array[k] = cmplx(inframe[k],0.0);	// copy inframe into a complex array
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	} 
	
	/************************* PROCESSING OF FRAME  HERE **************************/
		
	// put FFT in complex_array
	fft(FFTLEN, complex_array);
	
	for (k = 0; k<FFTLEN; k++)
	{
		if (e8 == 1){
			// store previous SNR values
			SNR_now[k] = noise_mag[k]/original_mag[k];
		}
		original_mag[k] = cabs(complex_array[k]); 
	}
	
	/*********** enhancement 1 *****************/
	// Low-pass filter the input frame magnitude
	if (e1 == 1){
		for (k = 0; k < FFTLEN; k++){
			lpf_input[k] = (1-k_pole)*original_mag[k] + (k_pole*lpf_input[k]);
		}
		update_minimums(lpf_input);
	}
	/*********** enhancement 2 *****************/
	// Low-pass filter the input frame power
	else if (e2 == 1){
		for (k = 0; k < FFTLEN; k++){
			lpf_input[k] = sqrt((1-k_pole)*original_mag[k]*original_mag[k] + (k_pole*lpf_input[k]*lpf_input[k]));
		}
		update_minimums(lpf_input);
	}
	else{
		// default case
		update_minimums(original_mag);
	}
	
	/****************** enhancement 3 *****************/
	// Low-pass filter the noise estimate magnitude
	if (e3 == 1){
		for (k = 0; k < FFTLEN; k++){
			lpf_noise[k] = (1-k_pole)*noise[k] + (k_pole*lpf_noise[k]);
		}
		noise_mag = lpf_noise;
	}
	// Low-pass filter the noise estimate power
	else if (e3 == 2){
		for (k = 0; k < FFTLEN; k++){
			lpf_noise[k] = sqrt((1-k_pole)*noise[k]*noise[k] + (k_pole*lpf_noise[k]*lpf_noise[k]));
		}
		noise_mag = lpf_noise;
	}
	else{
		noise_mag = noise;
	}
	
	// noise subtraction for-loop
	for (k = 0; k < FFTLEN; k++)
	{
		// calculate G(w)
		// G(w) = max(lamba_factor, lambda)
		/************** enhancement 4 ****************/
		if (e5 == 0){
			switch(e4){		
				case 2:
					factor = 1 - noise_mag[k]/original_mag[k];
					lambda_factor = lambda * noise_mag[k]/original_mag[k];
					break;
				
				case 3:
					lambda_factor = lambda * lpf_input[k]/original_mag[k];
					factor = 1 - noise_mag[k]/original_mag[k];
					break;
				
				case 4:
					lambda_factor = lambda * noise_mag[k]/lpf_input[k];
					factor = 1 - noise_mag[k]/lpf_input[k];
					break;
				
				case 5:
					lambda_factor = lambda;
					factor = 1 - noise_mag[k]/lpf_input[k];
					break;
				
				default: //1
					lambda_factor = lambda;
					factor = 1 - noise_mag[k]/original_mag[k];
			}
		}
		/**************** enhancement 5 ****************/
		else{
			switch(e5){	
				case 2:
					lambda_factor = lambda * sqrt(noise_mag[k]*noise_mag[k]/ (original_mag[k]*original_mag[k]));
					factor = sqrt(1 - noise_mag[k]*noise_mag[k]/(original_mag[k]*original_mag[k]));
					break;
				
				case 3:
					lambda_factor = lambda * sqrt(lpf_input[k]*lpf_input[k]/ (original_mag[k]*original_mag[k]));
					factor = sqrt(1 - noise_mag[k]*noise_mag[k]/ (original_mag[k]*original_mag[k]));
					break;
				
				case 4:
					lambda_factor = lambda * sqrt(noise_mag[k]*noise_mag[k]/ (lpf_input[k]*lpf_input[k]));
					factor = sqrt(1 - noise_mag[k]*noise_mag[k]/ (lpf_input[k]* lpf_input[k]));
					break;
				
				case 5:
					lambda_factor = lambda;
					factor = sqrt(1 - noise_mag[k]*noise_mag[k]/ (lpf_input[k]*lpf_input[k]));
					break;
				
				default: //1
					lambda_factor = lambda;
					factor = sqrt(1 - noise_mag[k]*noise_mag[k]/ (original_mag[k]*original_mag[k]));
			}
		}
		
		mag[k] = max(lambda_factor, factor);
		// output Y(w)=X(w)*|G(w)|
		complex_array[k].r = mag[k] * complex_array[k].r;
		complex_array[k].i = mag[k] * complex_array[k].i;
		
		/************ Enhancement8 ***************/
		// complex_array is the processed FFT that will be delayed to output
		// complex_array_now is the FFT to be output now
		// complex_array_prev is the FFT of the previous output
		
		if (e8 == 1){
			if ((SNR_now[k]) > noise_threshold){
				// complex_array_now = complex_with_smallest_mag(complex_array, complex_array_now, complex_array_prev)
				// three adjacent frames
				
				if (complex_mag(complex_array_now[k]) > complex_mag(complex_array[k])){
					complex_array_now[k] = complex_array[k];
				}
				if (complex_mag(complex_array_now[k]) > complex_mag(complex_array_prev[k])){
					complex_array_now[k] = complex_array_prev[k];
				}
			} 
		}
		
		/****** Extra enhancement: Frequency-domain filtering *******/
		if (e_speech_gain == 1){
			
			// mirror the right side to the left
			if (k > (FFTLEN/2)){
				kk = FFTLEN - k;
			}
			else{
				kk = k;
			}
			
			// Pass-band filter: amplify the speech frequency range
			if (kk >= low_freq && kk <= high_freq){
				complex_array[k].r = complex_array[k].r * speech_gain;
				complex_array[k].i = complex_array[k].i * speech_gain;
			}

			// Sharp, low-pass filtering
			else if (kk >= cut_off){
				complex_array[k] =  cmplx(0.0, 0.0);
			}
		}
	}
	
	// Enhancement 8: rotate array pointers of Y(w) history
	if (e8 == 1){
		temp_ca = complex_array;
		complex_array = complex_array_prev;
		complex_array_prev = complex_array_now;
		complex_array_now = temp_ca;
		
		ifft(FFTLEN, complex_array_prev);
	}
	else{
		ifft(FFTLEN, complex_array);
	}
	

    for (k=0; k<FFTLEN; k++)
	{                           
		if (allpass == 1){
			// copy input straight into output
			outframe[k] = inframe[k];
		}
		else if (e8 == 1){
			// enhancement 8: one frame delay
			outframe[k] = complex_array_prev[k].r; // _now before c
		}
		else{
			outframe[k] = complex_array[k].r;
		}
	}
	
	/********************************************************************************/
	
    /* multiply outframe by output window and overlap-add into output buffer */  
                           
	m=io_ptr0;
    
    for (k=0;k<(FFTLEN-FRAMEINC);k++) 
	{    										/* this loop adds into outbuffer */                       
	  	outbuffer[m] = outbuffer[m]+outframe[k]*outwin[k];   
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	}         
    for (;k<FFTLEN;k++) 
	{                           
		outbuffer[m] = outframe[k]*outwin[k];   /* this loop over-writes outbuffer */        
	    m++;
	}	                                   
}        
/*************************** INTERRUPT SERVICE ROUTINE  *****************************/

// Map this to the appropriate interrupt in the CDB file
   
void ISR_AIC(void)
{       
	short sample;
	/* Read and write the ADC and DAC using inbuffer and outbuffer */
	
	sample = mono_read_16Bit();
	inbuffer[io_ptr] = ((float)sample)*ingain;
		/* write new output data */
	mono_write_16Bit((int)(outbuffer[io_ptr]*outgain)); 
	
	/* update io_ptr and check for buffer wraparound */    
	
	if (++io_ptr >= CIRCBUF) io_ptr=0;
}

/************************************************************************************/

void update_minimums(float* input)
{
	int k;
	float *temp;
	float SNR, min_noise;
	
	if (frame_count == 0){
		// first frame input for this M bin (every 2.5s)
		for (k = 0; k < FFTLEN; k++){
			M1[k] = input[k];
		}
	}
	else{
		for (k = 0; k < FFTLEN; k++){
			// M1 = min(M1, input_frame_mag) for each f
			if (input[k] < M1[k]){
				M1[k] = input[k];
			}
		}
	}
	
	frame_count++;
	
	// if 10s of minumums collected
	if (frame_count >= ((int) (detection_period/(4*TFRAME))) ){
		frame_count = 0;
		
		// rotate pointers to Minimum bins
		temp = M4;
		M4 = M3;
		M3 = M2;
		M2 = M1;
		M1 = temp;
	
		// M1 is made to always be the current minimum bin
		
		k_pole = exp(-TFRAME/tau);	// Low-Pass Filter pole
		
		for (k = 0; k < FFTLEN; k++){			
			min_noise = min(M1[k], min(M2[k],min(M3[k],M4[k])));
			
			if (e6 == 1){
				SNR = input[k]/ (min_noise); // (computationally) simplified SNR value
				if (SNR < iSNR_threshold){
					// increase estimate of the noise by alpha_high scaling
					noise[k] = alpha * alpha_high * min_noise;
				}
				else{
					noise[k] = alpha * min_noise;
				}
			}
			else{
				noise[k] = alpha * min_noise;
			}

		}
	}
}

float complex_mag(complex c){
	return sqrt(c.r * c.r + c.i * c.i);
}
