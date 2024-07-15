# protodune_DQM_BNL

This is a demo for protoDUNE DQM.

Initial framework is from Gabriela Vitti Stenico 

## Usage
```
lar -n1 -c standard_reco_stage1_protodunehd_keepup.fcl <input.hdf5>
root -l -b -q 'convert_artroot_DQM_2.C("reco_stage1_artroot.root")'
```

## convert_artroot_DQM_2.C 

The input file is artROOT file from larsoft with tag *reco_stage1*.

Now it will output plots:

- Waveforms: all waveform, before and after noise filter; waveform in seperate APA and plane.
- Baseline;
- RMS;
- Niose spectrum vs frequency;
- Correlation channel to channel;

## load_metadata_my.py

Generate a DQM page online


I am keep running a DQM page on dunegpvm11. 

If you are interested, you can see it by:

  
  ssh -N -f -L localhost:8050:localhost:8050 \<username\>@dunegpvm11.fnal.gov
  
  
  http://127.0.0.1:8050/
  
