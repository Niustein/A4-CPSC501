#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <stdint.h>
#include <unistd.h>
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
int getFileSize(FILE *wavFile);
void printHeader(AudioHeader header);
void convolve(double signalInput[], int numInput, double impulseResponse[], int numImpulse, double output[], int numResult);
void writeOutput(string oFile, double oData[], int oLength, int sampleRate);


int finalLength;



// find the file size
int getFileSize(FILE* wavFile)
{
    int fileSize = 0;
    fseek(wavFile, 0, SEEK_END);

    fileSize = ftell(wavFile);

    fseek(wavFile, 0, SEEK_SET);
    return fileSize;
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
    
    cout << "test" << endl;
}

//Code from given handouts
void convolve(double signalInput[], int numInput, double impulseResponse[], int numImpulse, double output[], int numResult){
    
    int i, j;
    double maxRange = 0;
    
    cout << "reached convolve function" << endl;
    
    if (numResult != (numInput + numImpulse -1)){
        cout << "numResult is the wrong size" << endl;
        return;
    }

    // Clear output buffer
    for(i =0; i < numResult; i++){
        output[i] = 0.0;
    }
    
    cout << "cleared buffer" << endl;
    
    for(i = 0; i < numInput; i++){
        for(j = 0; j < numImpulse; j++){
            output[i+j] += signalInput[i] * impulseResponse[j];
                        
            // deal with out of range
            if((abs(output[i+j])) > maxRange){
                maxRange = abs(output[i+j]);
			}
		}
    }
    
    cout << "end convolve step" << endl;
    
    for(i = 0; i < numResult; i++){
        output[i] = output[i] / maxRange;
    }
    
    cout << "end convolve function" << endl;
    
}

void writeOutput(string oFile, double oData[], int oLength, int sampleRate){
    ofstream oStream;
    oStream.open(oFile, ios::binary | ios::trunc);
    
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
    for(int i =0; i < (subChunk2Size / (bitsPerSample/8)); i++){
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

int main(int argc, char * const argv[]){
    
    if (argc != 4){
        printf("Wrong number of arguments provided. Closing program");
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
    
    AudioHeader testAudio(argv[3]);
    
    printHeader(testAudio);
    
    return 0;    

}
