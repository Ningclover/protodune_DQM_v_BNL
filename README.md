# protodune_DQM_BNL

This is a demo for protoDUNE DQM.

Initial framework is from Gabriela Vitti Stenico 

## Usage
```
lar -n1 -c standard_reco_stage1_protodunehd_keepup.fcl <input.hdf5>
root -l -b -q 'convert_artroot_DQM_2.C("reco_stage1_artroot.root","<Path to put output plots>")'
```

## convert_artroot_DQM_2.C 

The input file is artROOT file from larsoft with tag *reco_stage1*.

Now it will output plots:

- Waveforms: all waveform, before and after noise filter; waveform in seperate APA and plane.
- Baseline;
- RMS;
- Noise spectrum vs frequency;
- Correlation channel to channel;

Also output datas in a runNo*.dat file:

- Averaged RMS after NF for each event
- Averaged peak of noise spectra after NF for each event


## Draw_summary.cc 

Can be used to draw the summary plots, with input from the *.dat file generated from *convert_artroot_DQM_2.C*.


## convert_artroot_DQM_3.C 

The input file is artROOT file from larsoft with tag *reco_stage2*.

Now it will output plots:

- decon Waveforms: all waveform, before and after noise filter; waveform in seperate APA and plane.
- Hit information: Number of hits found in each channel; Total ADC of all hits in each channel



## load_metadata_my.py

Generate a DQM page online


I am keep running a DQM page on dunegpvm11. 

If you are interested, you can see it by:

  ```
  ssh -N -f -L localhost:8050:localhost:8050 <username>@dunegpvm11.fnal.gov
  
  http://127.0.0.1:8050/
  ```

## processdata_loop.sh

A script to process a list of artROOT file on *reco_stage1*

*list.dat*: should contain a list of reco_stage1 artROOT files.

*outpath*: where you want to store the output plots.
