//==========================================================================

typedef struct {
    float a0, a1, a2, b1, b2;
} LRCoefficients;


struct Filter{
    
    LRCoefficients hpfCoeffs;
    LRCoefficients lpfCoeffs;
    
    
    void hpfLRCoeffs(float f_crossover, float fs)
    {
        float theta = 2 * M_PI * f_crossover / fs;
        float Wc = M_PI * f_crossover;
        float k = Wc / tan(theta);
        float d = pow(k, 2.0) + pow(Wc, 2.0) + 2.0 * k * Wc;
        
        hpfCoeffs.a0 = pow(Wc, 2.0) / d;
        hpfCoeffs.a1 = -2.0 * pow(Wc, 2.0) / d;
        hpfCoeffs.a2 = hpfCoeffs.a0;
        hpfCoeffs.b1 = (-2.0 * pow(k, 2.0) + 2.0 * pow(Wc, 2.0)) / d;
        hpfCoeffs.b2 = (-2.0 * k * Wc + pow(k, 2.0) + pow(Wc, 2.0)) / d;
    }

    void lpfLRCoeffs(float f_crossover, float fs)
    {
        float theta = 2 * M_PI * f_crossover / fs;
        float Wc = M_PI * f_crossover;
        float k = Wc / tan(theta);
        float d = pow(k, 2.0) + pow(Wc, 2.0) + 2.0 * k * Wc;
        
        lpfCoeffs.a0 = pow(Wc, 2.0) / d;
        lpfCoeffs.a1 = 2.0 * pow(Wc, 2.0) / d;
        lpfCoeffs.a2 = lpfCoeffs.a0;
        lpfCoeffs.b1 = (-2.0 * pow(k, 2.0) + 2.0 * pow(Wc, 2.0)) / d;
        lpfCoeffs.b2 = (-2.0 * k * Wc + pow(k, 2.0) + pow(Wc, 2.0)) / d;
    }

    float lowpass_filter(float input, float *state1, float *state2, float a0, float a1, float a2, float b1, float b2) {
        float output = a0 * input + a1 * (*state1) + a2 * (*state2);
        *state2 = *state1;
        *state1 = input - b1 * (*state1) - b2 * (*state2);
        return output;
    }

    float highpass_filter(float input, float *state1, float *state2, float a0, float a1, float a2, float b1, float b2) {
        
        float output = a0 * input + a1 * (*state1) + a2 * (*state2);
        *state2 = *state1;
        *state1 = input - b1 * (*state1) - b2 * (*state2);
        return output * (-1);
    }
};


struct TwoBandLinkwitzRileyFilter
{
    
    float* high_states_1;
    float* high_states_2;
    float* low_states_1;
    float* low_states_2;
    float* outputSamples;
    float* low_outputs;
    float* high_outputs;
    float* dist_lows;
    
    float f_crossover = 0.f; // Crossover frequency
    float initValue = 0.f;
    Filter* filters;

    float** inBuffer;
    float** outBuffer;
    float** highBand;
    float** lowBand;

    int _nMaxChannels;
    int _nMaxBlockSize;
    float _fGain_01;
    float fs;

    
    void Init() {
        inBuffer =  NULL;
        outBuffer =  NULL;
        lowBand = NULL;
        highBand = NULL;
        _nMaxChannels = 1;
        _nMaxBlockSize = 1;
        _fGain_01 = 1;
        
        high_states_1 =  NULL;
        high_states_2 =  NULL;
        low_states_1 =  NULL;
        low_states_2 =  NULL;
        outputSamples =  NULL;
        low_outputs =  NULL;
        high_outputs =  NULL;
        dist_lows =  NULL;
        
    } //RRS: All initializations needed for your DSP, memory allocations are allowed inside

    //RRS: Memory allocations are allowed inside
    void SetMaxBlockSize(int a_nMaxBlockSize)
    {
        if (_nMaxBlockSize != a_nMaxBlockSize)
        {
            _nMaxBlockSize = a_nMaxBlockSize;

            _ReAllocInternalBuffers(_nMaxChannels);
        }
    }
    
    void SetCrossoverFrequency(float a_nCrossoverFreq)
    {
        f_crossover = a_nCrossoverFreq;
        
        for(int channel = 0; channel < _nMaxChannels; channel++)
        {
            filters[channel].hpfLRCoeffs(f_crossover, fs);
            filters[channel].lpfLRCoeffs(f_crossover, fs);
        }
    }
    
    void SetMaxChannels(int a_nMaxChannels)
    {
        if (_nMaxChannels != a_nMaxChannels)
            _ReAllocInternalBuffers(a_nMaxChannels);
        
    }

    //RRS: Sample rate is not constant, so you have to reinitialize your sample rate dependent params (such as filters coeffs) on this call from our framework
    void SetSampleRate(float a_fSampleRate_Hz) {
        fs = a_fSampleRate_Hz;
        
        for(int channel = 0; channel < _nMaxChannels; channel++)
        {
            filters[channel].hpfLRCoeffs(f_crossover, fs);
            filters[channel].lpfLRCoeffs(f_crossover, fs);
        }
    } //RRS: Memory allocations are allowed inside

    void SetGain(float a_fGain_01) { _fGain_01 = a_fGain_01; } //RRS: Assertion: No memory allocations are allowed inside!
    void SetSomeParam1(float a_fSomeParam1Value) {} //RRS: Assertion: No memory allocations are allowed inside!
    void SetSomeParam2(float a_fSomeParam2Value) {} //RRS: Assertion: No memory allocations are allowed inside!
        
    void Release() { _ReleaseInternalBuffers();
        //RRS: All previously allocated memory can be deallocated here
    }
    
    float ProcessSampleLow(float* readData, Filter channelFilter, int channel, int index){
        auto sample = channelFilter.lowpass_filter(readData[index], &low_states_1[channel], &low_states_2[channel], channelFilter.lpfCoeffs.a0, channelFilter.lpfCoeffs.a1, channelFilter.lpfCoeffs.a2, channelFilter.lpfCoeffs.b1, channelFilter.lpfCoeffs.b2);
        
        return sample;
    }

    float ProcessSampleHigh(float* readData, Filter channelFilter, int channel, int index){
        auto sample = channelFilter.highpass_filter(readData[index], &high_states_1[channel], &high_states_2[channel], channelFilter.hpfCoeffs.a0, channelFilter.hpfCoeffs.a1, channelFilter.hpfCoeffs.a2, channelFilter.hpfCoeffs.b1, channelFilter.hpfCoeffs.b2);
        
        return sample;
    }
    
    float* ProcessHighBand(float* inStream, float** band, Filter channelFilter, int a_nChannels, int a_nSampleCount, int channel)
    {
        for (int i = 0; i < a_nSampleCount; ++i)
        {
            high_outputs[channel] = channelFilter.highpass_filter(inStream[i], &high_states_1[channel], &high_states_2[channel], channelFilter.hpfCoeffs.a0, channelFilter.hpfCoeffs.a1, channelFilter.hpfCoeffs.a2, channelFilter.hpfCoeffs.b1, channelFilter.hpfCoeffs.b2);
            
            band[channel][i] = high_outputs[channel];
            
        }
        return band[channel];
    }
    
    float* ProcessLowBand(float* inStream, float** band, Filter channelFilter, int a_nChannels, int a_nSampleCount, int channel)
    {
        for (int i = 0; i < a_nSampleCount; ++i)
        {
            low_outputs[channel] = channelFilter.lowpass_filter(inStream[i], &low_states_1[channel], &low_states_2[channel], channelFilter.lpfCoeffs.a0, channelFilter.lpfCoeffs.a1, channelFilter.lpfCoeffs.a2, channelFilter.lpfCoeffs.b1, channelFilter.lpfCoeffs.b2);
            
            band[channel][i] = low_outputs[channel];
        }
        
        return band[channel];
    }
        
        
    //RRS: Assertion: No memory allocations are allowed inside!
    void Process(float** a_vAudioBlocksInPlace, int a_nChannels, int a_nSampleCount)
    {
        //RRS: Assertion: a_nChannels less or equal _nMaxChannels set in SetMaxChannels()
        //RRS: Assertion: a_nSampleCount less or equal _nMaxBlockSize set in SetMaxBlockSize()
        
        //RRS: Here we copy the input audio data into internal buffers, apply Gain and then copy the processed data back to the input buffers (in-place).
        //RRS: There is no need for this copy operation (and these internal buffers): it's just a demonstration of how to work with the internal buffers if you need them.

        
        for (int channel = 0; channel < a_nChannels; ++channel)
        {
            memcpy(inBuffer[channel], a_vAudioBlocksInPlace[channel], a_nSampleCount * sizeof(float));
            
            float* readData = inBuffer[channel];
            float* writeData = outBuffer[channel];
            auto monoFilter = filters[channel];
            
            // Process audio samples
            highBand[channel] = ProcessHighBand(readData, highBand, monoFilter, a_nChannels, a_nSampleCount, channel);
            lowBand[channel] = ProcessLowBand(readData, lowBand, monoFilter, a_nChannels, a_nSampleCount, channel);

            for (int i = 0; i < a_nSampleCount; ++i)
            {
                writeData[i] = (highBand[channel][i] + lowBand[channel][i]) * .707f;  // Sum Signals
            }
                        
            memcpy(a_vAudioBlocksInPlace[channel], outBuffer[channel], a_nSampleCount * sizeof(float));
        }
    }

    void _ReleaseInternalBuffers()
    {
        if (inBuffer)
        {
            for (int n = 0; n < _nMaxChannels; ++n)
            {
                delete[] inBuffer[n];
            }
            
            delete[] inBuffer;
        }
        
        if (outBuffer)
        {
            for (int n = 0; n < _nMaxChannels; ++n)
            {
                delete[] outBuffer[n];
            }

            delete[] outBuffer;
        }
        
        if (highBand)
        {
            for (int n = 0; n < _nMaxChannels; ++n)
            {
                delete[] highBand[n];
            }
            
            delete[] highBand;
        }
        
        if (lowBand)
        {
            for (int n = 0; n < _nMaxChannels; ++n)
            {
                delete[] lowBand[n];
            }

            delete[] lowBand;
        }
        
        if(filters)
        {
            delete[] filters;
        }
        
        
        inBuffer =  NULL;
        outBuffer =  NULL;
        highBand =  NULL;
        lowBand =  NULL;
        filters = NULL;
    }

    void _ReAllocInternalBuffers(int a_nNewMaxChannels)
    {
        _ReleaseInternalBuffers();
        
        inBuffer = new float*[_nMaxChannels = a_nNewMaxChannels];
        outBuffer = new float*[_nMaxChannels = a_nNewMaxChannels];
        highBand = new float*[_nMaxChannels = a_nNewMaxChannels];
        lowBand = new float*[_nMaxChannels = a_nNewMaxChannels];
        filters = (Filter*) malloc(_nMaxChannels * sizeof(Filter));

        high_states_1 = (float *) malloc(_nMaxChannels);
        high_states_2 = (float *) malloc(_nMaxChannels);
        low_states_1 = (float *) malloc(_nMaxChannels);
        low_states_2 = (float *) malloc(_nMaxChannels);
        outputSamples = (float *) malloc(_nMaxChannels);
        low_outputs = (float *) malloc(_nMaxChannels);
        high_outputs = (float *) malloc(_nMaxChannels);
        dist_lows = (float *) malloc(_nMaxChannels);
        
        for (int n = 0; n < _nMaxChannels; ++n)
        {
            inBuffer[n] = new float[_nMaxBlockSize];
            outBuffer[n] = new float[_nMaxBlockSize];
            highBand[n] = new float[_nMaxBlockSize];
            lowBand[n] = new float[_nMaxBlockSize];
        }
    }
};

