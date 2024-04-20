(root_env) ubuntu@ubuntu2204:prototype$ root 'project.cxx(1000,1)'
   ------------------------------------------------------------------
  | Welcome to ROOT 6.30/04                        https://root.cern |
  | (c) 1995-2024, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxarm64 on Feb 07 2024, 00:33:26                    |
  | From heads/master@tags/v6-30-04                                  |
  | With                                                             |
  | Try '.help'/'.?', '.demo', '.license', '.credits', '.quit'/'.q'  |
   ------------------------------------------------------------------

root [0] 
Processing project.cxx(1000,1)...
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
Number of signal events: 2000
Number of background events: 1000
(int) 0
root [1] .x project.cxx(1000,10)
Number of signal events: 11000
Number of background events: 10000
(int) 0
root [2] .x project.cxx(100000,10)
Number of signal events: 1.1e+06
Number of background events: 1e+06
(int) 0
root [3] .x project.cxx(1000,10)
Number of signal events: 11000
Number of background events: 10000
(int) 0
root [4] .x project.cxx(10000,10)
Number of signal events: 110000
Number of background events: 100000
(int) 0
root [5] .q
