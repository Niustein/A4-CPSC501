#include <iostream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>
#include stdint.h>
#include <unistd.h>
using namespace std;

// Define classes
class AudioHeader
{
    
    public:
        uint8_t         riff[4], wave[4], fmt[4], subChunk2ID[4];
        uint32_t        chunkSize, subChunk1Size, samplesPS, bytesPS, subChunk2Size;
        uint16_t        audioFormat, numChannel, blockAlign, bitsPerSample;
        
        double*         dataArray;
        Header(string);
        
}

//Declare Functions
int getFileSize(FILE *inFile);
void printHeader(AudioHeader header);

int finalLength;

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

    extraData = new char[subChunk1Size-16];
    input.read(extraData, subChunk1Size-16);
    input.read((char*)&subChunk2ID, 4);
    input.read((char*)&subChunk2Size, 4);
    
    dataA = new double[subChunk2Size / (bitsPerSample/8)];
    
    int16_t sample;
    
    for(int i =0; i < (dataChunkSize / (bitsPerSample/8); i++){
        input.read((char*)&sample, 2);
        
        double convert = (double) sample / (double)INT16_MAX;
        
        // Check if beyond max range, if so set to max limit
        if(convert < -1.0){
            convert = (-1.0);
        }
        
        dataA[i] = convert;
        
    }
    
    // close inputfile stream
    input.close();
    
}

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
    cout << "--- RIFF Chunk Descriptor---" << endl
    cout << "RIFF header                :" << header.riff[0] << header.riff[1] << header.riff[2] << header.riff[3] << endl;
    cout << "Chunk Size                 :" << header.chunkSize << endL
    cout << "Format                     :" << header.fmt[0] << header.fmt[1] << header.fmt[2] << header.fmt[3] << endl;

    cout << "------ fmt sub-chunk  ------" << endl
    cout << "subChunk1ID                :" << "'" << header.fmt[0] << header.fmt[1] << header.fmt[2] << header.fmt [3] << "'" << endl;
    cout << "subChunk1Size              :" << header.subChunk1Size << endl;
    cout << "audioFormat                :" << header.audioFormat << endl;
    cout << "numChannels                :" << header.numChannel << endl;
    cout << "sampleRate                 :" << header.samplesPS << endl;
    cout << "byteRate                   :" << header.bytesPS << endl;
    cout << "blockAlign                 :" << header.blockAlign << endl;
    cout << "bitsPerSample              :" << header.bitesPerSample << endl;
    
    cout << "------ data sub-chunk ------" << endl
    cout << "subChunk2ID                :" << header.subChunkID2[0] << header.subChunkID2[1] << header.subChunkID2[2] << header.subChunkID2[3] << endl;
    cout << "subChunk2Size              :" << header.subChunk2Size << endl;
}

//Code from given handouts
void convolve(float signalInput[], int numInput, float impulseResponse[], int numImpulse, float output[], int numResult ){
    
    int i, j;
    double maxRange;
    
    printf("reached convolve function");
    
    if (numResult != (numInput + numImpulse -1){
        printf("numResult is the wrong size");
        return;
    }
    
    // Clear output buffer
    for(i =0; i < numResult; i++){
        output[i] = 0.0;
    }
    
    for(i = 0; i < numInput; i++){
        for(j = 0; j < numInpulse; j++){
            ouput[i+j] += signalInput[i] * impulseResponse[j];
            
            // deal with out of range
            if((abs(output[i+j])) > maxRange):
                maxRange = abs(y[n+m]);
        }
    }
    
    for(i = 0; i < finalLength; i++){
        output[i] = output[i] / maxRange]
    }
    
}

void writeOutput(String oFile, double oData[], int oLength, int sampleRate){
    ofstream oStream;
    oStream.open(oFile, ios::binary | ios::trunc);
    
    
    int subChunk1Size = 16
    short audioFormat = 1, numChannels = 1, bitsPS = 16;
    int byteRate = numChannels * sampleRate * (bitsPS / 8);
    int blockAlign = numChannels * (bitsPS/8)
    
    //writing RIFF chunk descriptor
    char* chunkID = new char[4]{"R", "I", "F", "F"};
    oStream.write(chunkID, 4);
    
    int dLength = oLength * 2;
    int chunkSize = dLength + 36;
    oStream.write((char*)&chunkSize, 4);
    
    char* format = new char[4]{"W","A","V","E"}
    oStream.write(format, 4);
    
    //Writing fmt sub-chunk
    char* subChunk1ID = new char[4]{"f","m","t"," "}
    oStream.write(subChunk1ID, 4); 
    oStream.write((char*)&subChunk1Size, 4);
    oStream.write((char*)&audioFormat, 2);
    oStream.write((char*)&numChannels, 2);
    oStream.write((char*)&sampleRate, 4);
    oStream.write((char*)&byteRate, 4);
    oStream.write((char*)&blockAlign, 2);
    oStream.write((char*)&bitsPS, 2);
    
    // Writing data sub-chunk
    char* subChunk2ID = new char[4]{"d","a","t","a"};
    oStream.write(subChunk2ID, 4);
    oStream.write((char*)&dLength, 4);
    
    // Writing Data
    int16_t oData;
    for(int i = 0; i < oLength; i++){
        oData = (int16_t) (oData[i] * INT16_MAX);
        oStream.write((char*)&oData, 2);
    }
    
    oStream.close();
    
}

int main(int argc, char * const argv[]){
    
    if (argc != 4){
        printf("Wrong number of arguments provided. Closing program");
        exit(-1);
    }
    
    printf("reading audio file: ", argv[1]);
    AudioHeader iAudio(argv[1]);
    
    printf("reading IR file: ", argv[2]);
    AudioHeader irAudio(argv[2]);
    
    PrintHeader(iAudio);
    printf("");
    
    printHeader(irAudio);
    printf("");
    
    printf("finished reading wave files");
    
    int iLength = iAudio.subChunk2Size / 2;
    int irLength = irAudio.subChunk2Size / 2
    int rLength = iLength + irLength - 1
    
    result = new double[rLength];
    
    printf("convulving");
    convulve(iAudio.dataArray, iLength, irAudio.dataArray, irLength, result, rLength);
    printf("finished convulving");
    
    printf("writing wave file");
    writeOutput(argv[3], result, rLength, iAudio.samplesPS);
    printf("finished writing wave file");
    
    return 0;    

}