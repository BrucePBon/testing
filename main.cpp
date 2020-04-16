#include <iostream>
#include <chrono>
#include <fftw3.h>

int main( int argc, char * argv [] ) 
{
	int IA0  = 10; 
	int max  = 60;
	int reps = 100; 

	fftw_complex* in1; 
	fftw_complex* in2; 
	fftw_plan plan1;  // inplace  fft on in1
	fftw_plan plan2;  // inplace  fft on in2
	fftw_plan iplan;  // inplace ifft on in1

	int IA; 
	for ( IA=IA0; IA <= max; IA++ ) 
	{
		double t = 0.0; // time will be accumulted here
		int n = IA*IA; 
		double corr[n]; 

		in1 = (fftw_complex*) fftw_malloc ( sizeof (fftw_complex) * n );
		in2 = (fftw_complex*) fftw_malloc ( sizeof (fftw_complex) * n );
		for ( int i = 0; i < n; i++ ) { // Initializing to 0.0's. Maybe not necessary
			in1[i][0] = 0.0; 
			in1[i][1] = 0.0; 
			in2[i][0] = 0.0; 
			in2[i][1] = 0.0; 
		}

		plan1 = fftw_plan_dft_2d( IA, IA, in1, in1, FFTW_FORWARD , FFTW_MEASURE );
		plan2 = fftw_plan_dft_2d( IA, IA, in2, in2, FFTW_FORWARD , FFTW_MEASURE );
		iplan = fftw_plan_dft_2d( IA, IA, in1, in1, FFTW_BACKWARD, FFTW_MEASURE );

		for (int r=1; r <= reps; r++ ) 
		{
			// Here starts the execution of cross-correlation
    		auto start = std::chrono::high_resolution_clock::now(); 

			fftw_execute( plan1 ); 
			fftw_execute( plan2 ); 

			int i; 
			double a, b, c, d; 
			// replacing in1 with: conj.( in1 ) .* in2
			for ( i=0; i < n; i++ ) {
				a = in1[i][0]; 
				b = in1[i][1]; 
				c = in2[i][0]; 
				d = in2[i][1]; 
				in1[i][0] = a * c + b * d; 
				in1[i][1] = a * d - b * c;
			}

			fftw_execute(iplan); 

			// storing the real part of in1 into corr
			for ( i=0; i < n; i++ ) {
				corr[i] = in1[i][0]; 
			}

			// that's all for cross-correlation
			auto stop     = std::chrono::high_resolution_clock::now(); 
    		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);    
    		t = t + duration.count(); 
		}	

	t = t/reps; 

	std::cout << IA << ": " << t << "\n"; 
	}

	fftw_free( in1 ); 
	fftw_free( in2 ); 
	fftw_destroy_plan( plan1 ); 
	fftw_destroy_plan( plan2 ); 
	fftw_destroy_plan( iplan );
	
	return false; 
}