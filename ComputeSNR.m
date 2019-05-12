function [SNR] = ComputeSNR(signal, fs, annotatedOnsets, annotatedOffsets)
		% ComputeSNR2 - Estimates the signal-to-noise ratio (SNR) given a set of manual annotations (in seconds).
    	% 
	% This algorithm estimates the signal-to-noise ratio in recordings that contains manually annotations of the desired and noise segments. The assumption that the noise is completely uncorrelated from the desired signal is made (e.g. additive noise). Recordings without information regarding the sound intensity or the voltage-sound-intensity relationship in the recording devices may benefit.
	% 
	% The basis formula (eq. 1) is defined like so:
	% `SNR = 10 * log10((P(desired)−P(noise))/P(noise))`, (eq. 1)
	% 
	% where:
	% * P(s) is the power of a finite period signal and 
	% * P(s) = 1/N * (|s(1)|^2+|s(2)|^2+...+|s(N)|^2)
	% 
	% The desired signal is constructed by appending the entries of the mixed signal specified by the annotated segments (onset[i]-offset[i] pair), and the noise signal is constructed by appending the entries outside the annotations. Good annotations should make the noise signal contain only silence or not-desired segments. Turn flag DEBUG=true to save the segmented signals in two files called "desired.wav" and "noise.wav". This also prints a graph with the annotations and the oscillogram of the saved *.wav files.
	% 
	% References:
	% 1. Papadopoulos, P., & Tsiartas, A. (2014). A supervised signal-to-noise ratio estimation of desired signals. IEEE International Conference on Acoustics desired and Signal Processing, 8287–8291
	% 2. Beritelli, F., Casale, S., Grasso, R., & Spadaccini, A. (2010). Performance Evaluation of SNR Estimation Methods in Forensic Speaker Recognition. In 2010 Fourth International Conference on Emerging Security Information, Systems and Technologies (pp. 88–92)
	% 3. SNR estimation in segmented speech signals. URL: http://dsp.stackexchange.com/questions/35041/snr-estimation-in-segmented-desired-signals (last consulted on May 12th 2019)
   	% 
	% Syntax: 
	% * ComputeSNR2(signal, fs, annotatedOnsets, annotatedOffsets)
	% 
	% Inputs:
	% * signal - an array containing the mixed signals: noise and desired
	% * fs - sampling rate
	% * annotatedOnsets - array with manually annotated onsets of the desired signal
	% * annotatedOffsets - array with manually annotated offsets of the desired signal
	% 
	% Outputs:
	% * SNR - Signal-to-noise ratio estimate in dB
	% 
	% Notes:
	% * Other m-files required: none
	% * Subfunctions: none
	% * MAT-files required: none
	% 
	% ---
	% 
	% Author: Juan Manuel Fonseca-Solís
	% San José, Costa Rica
	% Email: juanma2268@gmail.com
	% Oct 2016; Last revision: 26-Oct-2016


	DEBUG = true;

	signalLength = numel(signal);	
	desired = [];	% desired signal
	noise = [];	% noise signal
	dt = 1/fs;
	signalDuration = dt*(signalLength-1); % total signal duration in seconds

	% build arrays for the onsets and length of the desired signal
	annotationsLength = numel(annotatedOnsets);
	segmentOnsets = zeros(annotationsLength, 1);
	segmentLengths = zeros(annotationsLength, 1);	
	
	for iAnnotations=1:annotationsLength
		
		% calculate onset and offset indexes for the current segment
		iOnset = floor(annotatedOnsets(iAnnotations)/signalDuration * signalLength);	
		iOnset = ifelse(0==iOnset,1,iOnset); % case when annotations is t=0		
		iOffset = floor(annotatedOffsets(iAnnotations)/signalDuration * signalLength); 		
		
		if DEBUG		
			fprintf("Onset: %i - Offset: %i\n",iOnset,iOffset)		
		endif	

		% sum energy for partially valid segments (i.e. only offset is outside) and valid segments
		isOnsetValid = 1<=iOnset<=signalLength;
		isOffsetValid = iOnset<=iOffset<=signalLength;		
		isSegmentValid = isOffsetValid && isOffsetValid;				
		
		if(isSegmentValid)
			segmentOnsets(iAnnotations) = iOnset;	

			% avoid segments to overlap
			if(iAnnotations<annotationsLength)
				iNextOnset = floor(annotatedOnsets(iAnnotations+1)/signalDuration * signalLength);
				if(iNextOnset<iOffset)
					iOffset = iNextOnset;
				endif
			endif		
		
			segmentLengths(iAnnotations) = iOffset-iOnset;

		elseif(isOnsetValid)		
			% at least the onset is valid, count the remaining until the end of the signal
			segmentOnsets(iAnnotations) = iOnset;
			segmentLengths(iAnnotations) = signalLength-iOnset;
		else
			% the segment is invalid, store size = 0
			fprintf("ERROR: found an invalid segment: %i-%i\n",iOnset,iOffset)
			segmentLengths(iAnnotations) = 0;
		endif
		
	end
	
	if DEBUG
		figure
		subplot(311)
		stem(segmentOnsets,ones(annotationsLength,1),'b--','linewidth',2)
		hold on	
		stem(segmentOnsets+segmentLengths,ones(annotationsLength,1),'r--','linewidth',2)
		title('Annotations (onsets in blue, offsets in red)')
	endif	

	% calcultate the desired and noise lengths
	desiredLength = sum(segmentLengths);
	noiseLength = signalLength - desiredLength;	
	desired = zeros(desiredLength, 1);
	noise = zeros(noiseLength, 1);
	
	% build up the desired and noise signals
	iAnnotations = idesired = iNoise = 1;
	iOffset = 0; % trick to pick the first segment at the start
	for iSignal=1:signalLength
		
		isSegmentOld = iOffset < iSignal;
		if(isSegmentOld && iAnnotations<=annotationsLength) % jump to next segment
			iOnset = segmentOnsets(iAnnotations);
			iOffset = segmentOnsets(iAnnotations) + segmentLengths(iAnnotations);
			iAnnotations += 1;
		endif
		
		isEntrydesired = iOnset <= iSignal && iSignal < iOffset;
		if(isEntrydesired)	
			desired(idesired) = signal(iSignal);
			idesired += 1;
		else
			noise(iNoise) = signal(iSignal);
			iNoise += 1;	
		endif

	endfor

	if DEBUG
		subplot(312)
		plot(desired)
		title('desired.wav')
		subplot(313)
		plot(noise)
		title('noise.wav')
	endif

	% make the signals zero-mean
	desired = desired - mean(desired);
	noise = noise - mean(noise);

	% compute the powers of desired and noise	
	desiredPower = noisePower = 0;	
	for(i=1:desiredLength)
		desiredPower += abs(desired(i))^2; % abs in case of complex signal
	endfor
	for(i=1:noiseLength)
		noisePower += abs(noise(i))^2; % abs in case of complex signal
	endfor
	desiredPower /= desiredLength;
	noisePower /= noiseLength;
	
	% SNR calculation
	SNR = (desiredPower - noisePower)/noisePower;	
	SNR = 10 * log10(SNR);

	if DEBUG
		wavwrite(desired,fs,'desired.wav');
		wavwrite(noise,fs,'noise.wav');
	endif

end

%!test
%! pkg load signal
%! annotatedOnsets = [0.798292];
%! annotatedOffsets = [1.213205]; 
%! [signal, fs]= wavread("wav/253387__hauntedswkids__hello.wav");
%! SNR = ComputeSNR(signal, fs, annotatedOnsets, annotatedOffsets);	% 21.831104
%! fprintf("SNR = %f\n", SNR)
%!
%!test
%! pkg load signal
%! annotatedOnsets = [0.798292];
%! annotatedOffsets = [1.213205]; 
%! [signal, fs]= wavread("wav/253387__hauntedswkids__hello_noisy.wav");	% -3.559825
%! SNR = ComputeSNR(signal, fs, annotatedOnsets, annotatedOffsets);
%! fprintf("SNR = %f\n", SNR)
%!


