OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38116954) q[0];
sx q[0];
rz(-1.086735) q[0];
sx q[0];
rz(0.83597393) q[0];
rz(2.6354191) q[2];
sx q[2];
rz(-1.0867599) q[2];
sx q[2];
rz(1.6137705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.3904966) q[1];
sx q[1];
rz(3.1027604) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19102328) q[3];
sx q[3];
rz(-1.028562) q[3];
sx q[3];
rz(-0.99381002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(-0.65223637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95603847) q[0];
sx q[0];
rz(-1.599405) q[0];
sx q[0];
rz(1.5831328) q[0];
x q[1];
rz(0.89276887) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(-2.6004651) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80174175) q[1];
sx q[1];
rz(-2.7738214) q[1];
sx q[1];
rz(0.66373177) q[1];
rz(0.90918031) q[3];
sx q[3];
rz(-1.0717857) q[3];
sx q[3];
rz(3.1391075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7499381) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(-1.2664938) q[0];
rz(-1.0565874) q[2];
sx q[2];
rz(-0.075767013) q[2];
sx q[2];
rz(2.5839992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.459356) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(-0.45046803) q[1];
rz(-pi) q[2];
rz(0.0098185929) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(-2.1992418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.13658) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.2269292) q[0];
rz(-pi) q[1];
rz(-0.58156275) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(0.046422596) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1069146) q[1];
sx q[1];
rz(-0.87768302) q[1];
sx q[1];
rz(-2.9702859) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9070204) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(0.58805874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(0.23434815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317961) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(1.1774506) q[0];
rz(1.3946103) q[2];
sx q[2];
rz(-2.641045) q[2];
sx q[2];
rz(-1.8793775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8824132) q[1];
sx q[1];
rz(-2.3098574) q[1];
sx q[1];
rz(-1.773136) q[1];
x q[2];
rz(-1.8946394) q[3];
sx q[3];
rz(-1.6456592) q[3];
sx q[3];
rz(-0.74136855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.9979427) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8457984) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(-3.0315115) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.3442163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.027187849) q[1];
sx q[1];
rz(-1.6611551) q[1];
sx q[1];
rz(1.1272217) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2915217) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(-2.6078893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.9082665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9702643) q[0];
sx q[0];
rz(-1.5818705) q[0];
sx q[0];
rz(1.1867255) q[0];
rz(-pi) q[1];
rz(-2.2912824) q[2];
sx q[2];
rz(-2.0680973) q[2];
sx q[2];
rz(2.8720299) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8804633) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-0.1111828) q[1];
rz(-pi) q[2];
rz(2.535378) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(-1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(-2.1255778) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163517) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.258237) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60511741) q[0];
sx q[0];
rz(-1.1400756) q[0];
sx q[0];
rz(0.14819781) q[0];
rz(-pi) q[1];
rz(2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.9922436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59626034) q[1];
sx q[1];
rz(-0.55570554) q[1];
sx q[1];
rz(-2.943379) q[1];
rz(-pi) q[2];
rz(2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-3.1138611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160307) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-0.82657878) q[0];
rz(-pi) q[1];
rz(2.2698195) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(-2.2733462) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7293538) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(-2.9669697) q[1];
x q[2];
rz(-1.2003044) q[3];
sx q[3];
rz(-0.071404608) q[3];
sx q[3];
rz(0.84884531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(-2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.1429943) q[0];
rz(-pi) q[1];
rz(-0.72980482) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(0.44621106) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32138667) q[1];
sx q[1];
rz(-2.4452129) q[1];
sx q[1];
rz(0.022298261) q[1];
rz(2.0503644) q[3];
sx q[3];
rz(-1.3551095) q[3];
sx q[3];
rz(1.8857764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(1.1889585) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(-2.0541035) q[3];
sx q[3];
rz(-1.8046422) q[3];
sx q[3];
rz(1.3840152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];