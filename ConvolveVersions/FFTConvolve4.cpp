#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdint.h>
#include <unistd.h>
#include <string>

#define SWAP(a,b)	tempr=(a);(a)=(b);(b)=tempr

using namespace std;

// Define classes
class AudioHeader
{
    
    public:
        uint8_t         riff[4], wave[4], fmt[4], subChunk2ID[4];
        uint32_t        chunkSize, subChunk1Size, samplesPS, bytesPS, subChunk2Size;
        uint16_t        audioFormat, numChannel, blockAlign, bitsPerSample;
        
        uint8_t*		extraData;
        
        double*         dataA;
        AudioHeader(string);
        
};

//Declare Functions
void four1 (double*, int, int);
void convolve(double*, int, double*, int, double*, int);
void writeOutput(string, double*, int, int);
void printHeader(AudioHeader);


AudioHeader::AudioHeader(string fileName){
    ifstream input;
    input.open(fileName, ios::binary);
    
    input.read((char*)&riff,4);
    input.read((char*)&chunkSize,4);
    input.read((char*)&wave, 4);
    input.read((char*)&fmt, 4);
    input.read((char*)&subChunk1Size, 4);
    input.read((char*)&audioFormat,2);
    input.read((char*)&numChannel, 2);
    input.read((char*)&samplesPS, 4);
    input.read((char*)&bytesPS, 4);
    input.read((char*)&blockAlign, 2);
    input.read((char*)&bitsPerSample, 2);

    extraData = new uint8_t[subChunk1Size-16];
    input.read((char*)&extraData, subChunk1Size-16);
    input.read((char*)&subChunk2ID, 4);
    input.read((char*)&subChunk2Size, 4);
    
    dataA = new double[subChunk2Size / (bitsPerSample/8)];
    
    int16_t sample;
    
    // Go through 2 bytes at a time
    for(int i =0; i < subChunk2Size / (bitsPerSample/8); i++){
        input.read((char*)&sample, 2);
        
        double limit = (double) sample / (double)INT16_MAX;
        
        // Check if beyond max negative range, if so set to max negative limit
        if(limit < -1.0){
            limit = (-1.0);
        }
        
        dataA[i] = limit;
        
    }
    
    // close inputfile stream
    input.close();
    
}


// Code from given handouts
void four1(double data[], int nn, int isign){
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;
	
	n = nn << 1;
	j = 1;
	
	for (i = 1; i < n; i += 2){
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m){
			j -= m;
			m >>= 1;
		}
		
		j += m;
	}
	
	mmax = 2;
	while (n > mmax){
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2){
			for (i = m; i <= n; i += istep){
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}	
}






void convolve(double signalInput[], int numInput, double impulseResponse[], int numImpulse, double output[], int numResult){
    
    int i, maxLength;
    
    double maxRange = 0.0;

    
    cout << "reached convolve function" << endl;

	cout << "NumInput = " << numInput << " numImpulse = " << numImpulse << endl;
	
	if(numInput < numImpulse){
		maxLength = numImpulse;
	}else{
		maxLength = numInput;
	}	
	
//	maxLength = max(numInput,numImpulse);
	
	int numSquared = 1;
		
	
	// find n in nn
	while(numSquared < maxLength){
		numSquared = numSquared << 1;
	} 

//	int nHolder = numSquared;
    
 //   numSquared = numSquared << 1;

	// Create double arrays of nn size and output array for complex multi 
	double* inputPadded = new double[numSquared];
	double* impulsePadded = new double[numSquared];
		
	// Zero padding new arrays inputPadded and impulsePadded

	for(i=0; i < numSquared/2;i++){
		inputPadded[i] = 0;
        inputPadded[i+1] = 0;
   		impulsePadded[i] = 0;
        impulsePadded[i+1] = 0;
	}
	
	// Put audio data into new arrays inputPadded
	for (i = 0; i < numInput; i++){
		inputPadded[i<<1] = signalInput[i];
	}
	// Put audio data into new array impulsePadded
	for (i = 0; i < numImpulse; i++){
		impulsePadded[i<<1] = impulseResponse[i];
	}
	
	cout << "start four1 algorithm" << endl;
		
	four1(inputPadded - 1, numSquared, 1);
	four1(impulsePadded - 1, numSquared, 1);
	
	cout << "end four1 algorithm" << endl;
	

	double* oData = new double[numSquared];


	// Complex multiplcation
	for(i=0; i < numSquared * 2; i+=2){
		// Real = (R * R) - (I * I)
		oData[i] = (inputPadded[i] * impulsePadded[i]) - (inputPadded[i+1] * impulsePadded[i+1]); 
		
		// Imaginary = (R * I) + (R * I) 
		oData[i+1] = (inputPadded[i] * impulsePadded[i+1]) + (inputPadded[i+1] * impulsePadded[i]);
	}
	
	// repass into four1 with inverse
	four1(oData-1, numSquared, -1);
	
	// Pass real values back into array
	for(i=0; i < numResult; i++){
		output[i] = oData[i<<1];
		
		// deal with out of range
		if((abs(oData[i<<1])) > maxRange){
			maxRange = abs(oData[i<<1]);
		}
	}
	
	for(i = 0; i < numResult; i++){
        output[i] = output[i] / maxRange;
    }

}


void writeOutput(string oFile, double oData[], int oLength, int sampleRate){
    ofstream oStream;
    oStream.open(oFile, ios::binary | ios::out);
    
    short audioFormat = 1, numChannel = 1, bitsPS = 16;
    int byteRate = numChannel * sampleRate * (bitsPS / 8);
    int blockAlign = numChannel * (bitsPS/8);
    int dLength = oLength * 2;
    int chunkSize = dLength + 36;
    
    //writing RIFF chunk descriptor
    char* chunkID = new char[4]{'R', 'I', 'F', 'F'};
    oStream.write(chunkID, 4);
    oStream.write((char*)&chunkSize, 4);
    
    char* format = new char[4]{'W','A','V','E'};
    oStream.write(format, 4);
    
    int subChunk1Size = 16;
    //Writing fmt sub-chunk
    char* subChunk1ID = new char[4]{'f','m','t',' '};
    oStream.write(subChunk1ID, 4); 
    oStream.write((char*)&subChunk1Size, 4);
    oStream.write((char*)&audioFormat, 2);
    oStream.write((char*)&numChannel, 2);
    oStream.write((char*)&sampleRate, 4);
    oStream.write((char*)&byteRate, 4);
    oStream.write((char*)&blockAlign, 2);
    oStream.write((char*)&bitsPS, 2);
    
    // Writing data sub-chunk
    char* subChunk2ID = new char[4]{'d','a','t','a'};
    oStream.write(subChunk2ID, 4);
    oStream.write((char*)&dLength, 4);
    
    // Writing Data
    int16_t outputData;
    for(int i = 0; i < oLength; i++){
        outputData = (int16_t) (oData[i] * INT16_MAX);
        oStream.write((char*)&outputData, 2);
    }
    
    oStream.close();
    
}





//Prints the wave header
void printHeader(AudioHeader header){
    cout << "--- RIFF Chunk Descriptor---" << endl;
    cout << "RIFF header                :" << header.riff[0] << header.riff[1] << header.riff[2] << header.riff[3] << endl;
    cout << "Chunk Size                 :" << header.chunkSize << endl;
    cout << "Format                     :" << header.fmt[0] << header.fmt[1] << header.fmt[2] << header.fmt[3] << endl;

    cout << "------ fmt sub-chunk  ------" << endl;
    cout << "subChunk1ID                :" << header.fmt[0] << header.fmt[1] << header.fmt[2] << header.fmt [3] << endl;
    cout << "subChunk1Size              :" << header.subChunk1Size << endl;
    cout << "audioFormat                :" << header.audioFormat << endl;
    cout << "numChannel                 :" << header.numChannel << endl;
    cout << "sampleRate                 :" << header.samplesPS << endl;
    cout << "byteRate                   :" << header.bytesPS << endl;
    cout << "blockAlign                 :" << header.blockAlign << endl;
    cout << "bitsPerSample              :" << header.bitsPerSample << endl;
    
    cout << "------ data sub-chunk ------" << endl;
    cout << "subChunk2ID                :" << header.subChunk2ID[0] << header.subChunk2ID[1] << header.subChunk2ID[2] << header.subChunk2ID[3] << endl;
    cout << "subChunk2Size              :" << header.subChunk2Size << endl;
}


int main(int argc, char * const argv[]){
    
    if (argc != 4){
        cout << "Wrong number of arguments provided. Closing program" << endl;
        exit(-1);
    }
    
    cout << "reading audio file: " << argv[1] << endl;
    AudioHeader iAudio(argv[1]);
    
    cout << "reading IR file: " << argv[2] << endl;
    AudioHeader irAudio(argv[2]);
    
    cout << "iAudio header information:" << endl;
    printHeader(iAudio);
	cout << endl;
    
    cout << "irAudio header information:" << endl;
    printHeader(irAudio);
    
    cout << "finished reading wave files" << endl;
    
    int iLength = iAudio.subChunk2Size / 2;
    int irLength = irAudio.subChunk2Size / 2;
    int rLength = iLength + irLength - 1;
    
    double* result = new double[rLength];
    
	cout << "convolving" << endl;
    convolve(iAudio.dataA, iLength, irAudio.dataA, irLength, result, rLength);
    cout << "finished convolving" << endl;
    
    cout << "writing wave file" << endl;
    writeOutput(argv[3], result, rLength, iAudio.samplesPS);
    cout << "finished writing wave file" << endl;
    
    return 0;    

}
