# WashU Epigenome Browser

[http://epigenomegateway.wustl.edu/browser/](http://epigenomegateway.wustl.edu/browser/)

- Only partial support for UCSC hubs (http://wiki.wubrowse.org/V29_(2_of_4):_displaying_track_hubs_from_UCSC_Genome_Browser), which is insufficient for our data
  - Doesn't support bigBED
  - For complete support, would need to make a datahub (JSON file)
  - [http://wiki.wubrowse.org/Datahub](http://wiki.wubrowse.org/Datahub)
- Renders bigWig as barplot :(
- They support a novel 'MethylC' track type (http://wiki.wubrowse.org/MethylC_track)
  - mCG, mCHH, and mCHG, optionally stranded
	- Would plot rawMeth
	- Need to create additional bigWig from BSseq objects
- The Matplot rendering is pretty good (http://wiki.wubrowse.org/Matplot)
	- But seems to interpolate through zero, which looks ugly for mC (mCG, especially)
	- No transparency makes hard to visualise when one sample is low and the other is high
	  - [Example of tracks being hard to see (NAcc is low while BA9 is high); Fig 3c of our paper](http://epigenomegateway.wustl.edu/browser/t/174832858.pdf)
- They do offer a 'service' to publish data as a hub (http://wiki.wubrowse.org/V22_(2_of_3):_public_datahubs)
