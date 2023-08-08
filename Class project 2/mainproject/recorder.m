fs = 16000;
recObj = audiorecorder(fs,16,1);%16000Hz,16bit, single channel
disp('Start speaking.')
recordblocking(recObj, 2);
myRecording = getaudiodata(recObj);
disp('End of Recording.');
filename = 'myvoice.wav'; 
audiowrite(filename,myRecording,fs);