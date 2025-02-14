OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4580392) q[0];
sx q[0];
rz(-0.91705051) q[0];
sx q[0];
rz(0.24669692) q[0];
rz(3.9858272) q[1];
sx q[1];
rz(1.2935473) q[1];
sx q[1];
rz(10.727439) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5555252) q[0];
sx q[0];
rz(-1.2159776) q[0];
sx q[0];
rz(0.3094425) q[0];
rz(-pi) q[1];
rz(2.6370722) q[2];
sx q[2];
rz(-2.2103643) q[2];
sx q[2];
rz(0.96790403) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0820391) q[1];
sx q[1];
rz(-2.5055725) q[1];
sx q[1];
rz(2.9530557) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4250303) q[3];
sx q[3];
rz(-1.4633388) q[3];
sx q[3];
rz(-1.859457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.3641597) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(-1.9470661) q[2];
rz(-1.901769) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164923) q[0];
sx q[0];
rz(-1.1273552) q[0];
sx q[0];
rz(-0.66705739) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(-2.0858696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919843) q[0];
sx q[0];
rz(-2.4469564) q[0];
sx q[0];
rz(1.1437835) q[0];
rz(-pi) q[1];
rz(1.6237359) q[2];
sx q[2];
rz(-0.38489562) q[2];
sx q[2];
rz(0.062907779) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70433863) q[1];
sx q[1];
rz(-1.2341712) q[1];
sx q[1];
rz(-3.047154) q[1];
x q[2];
rz(-0.74852826) q[3];
sx q[3];
rz(-1.0077884) q[3];
sx q[3];
rz(0.23458086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51862741) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(2.5248027) q[2];
rz(-2.8957497) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(0.1804633) q[0];
rz(0.47897419) q[1];
sx q[1];
rz(-0.45061794) q[1];
sx q[1];
rz(1.3165547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0947803) q[0];
sx q[0];
rz(-2.1553596) q[0];
sx q[0];
rz(0.78740904) q[0];
rz(-pi) q[1];
rz(-0.20256217) q[2];
sx q[2];
rz(-1.1656467) q[2];
sx q[2];
rz(-0.7140401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41539792) q[1];
sx q[1];
rz(-1.2737927) q[1];
sx q[1];
rz(2.4261977) q[1];
rz(-2.1785424) q[3];
sx q[3];
rz(-0.96745771) q[3];
sx q[3];
rz(-0.36336366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58671826) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(2.6678616) q[2];
rz(-2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(-1.0152869) q[0];
rz(0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-0.36605787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81824077) q[0];
sx q[0];
rz(-2.2466772) q[0];
sx q[0];
rz(0.9267207) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7160077) q[2];
sx q[2];
rz(-1.337707) q[2];
sx q[2];
rz(2.9200302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28747044) q[1];
sx q[1];
rz(-2.4123998) q[1];
sx q[1];
rz(-2.6602938) q[1];
x q[2];
rz(-0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(-0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4801243) q[2];
sx q[2];
rz(-1.3084359) q[2];
sx q[2];
rz(-2.4465731) q[2];
rz(0.85092893) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224391) q[0];
sx q[0];
rz(-0.16534403) q[0];
sx q[0];
rz(0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(0.42627898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3836455) q[0];
sx q[0];
rz(-1.8611127) q[0];
sx q[0];
rz(-2.2955672) q[0];
x q[1];
rz(-2.2441909) q[2];
sx q[2];
rz(-2.7798712) q[2];
sx q[2];
rz(2.7182686) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0002928) q[1];
sx q[1];
rz(-1.1361377) q[1];
sx q[1];
rz(3.0463957) q[1];
rz(-pi) q[2];
rz(-0.3235281) q[3];
sx q[3];
rz(-1.1051205) q[3];
sx q[3];
rz(-2.9126008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2686501) q[2];
sx q[2];
rz(-2.908417) q[2];
sx q[2];
rz(-3.1128913) q[2];
rz(0.21102333) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(-0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(2.9826214) q[0];
rz(0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(-1.5370625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66012525) q[0];
sx q[0];
rz(-0.9430389) q[0];
sx q[0];
rz(-2.5808236) q[0];
rz(0.36559029) q[2];
sx q[2];
rz(-1.2477008) q[2];
sx q[2];
rz(1.1981724) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3517396) q[1];
sx q[1];
rz(-2.5765214) q[1];
sx q[1];
rz(1.6050102) q[1];
rz(-1.8810684) q[3];
sx q[3];
rz(-1.2168443) q[3];
sx q[3];
rz(1.9913395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0752461) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(0.016544841) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(0.56104463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(-2.2313927) q[0];
rz(0.74288145) q[1];
sx q[1];
rz(-1.0123092) q[1];
sx q[1];
rz(1.2339309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260516) q[0];
sx q[0];
rz(-1.431246) q[0];
sx q[0];
rz(-1.6427137) q[0];
rz(-pi) q[1];
rz(1.0077613) q[2];
sx q[2];
rz(-2.6673311) q[2];
sx q[2];
rz(2.8104818) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9561718) q[1];
sx q[1];
rz(-0.90683504) q[1];
sx q[1];
rz(-2.0645622) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5911932) q[3];
sx q[3];
rz(-0.66662753) q[3];
sx q[3];
rz(2.3638099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10741216) q[2];
sx q[2];
rz(-1.567652) q[2];
sx q[2];
rz(-0.91599715) q[2];
rz(2.8902174) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(-0.34534064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6420355) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(-0.04059759) q[0];
rz(1.1740855) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(1.0601128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1112633) q[0];
sx q[0];
rz(-1.666774) q[0];
sx q[0];
rz(2.7594271) q[0];
rz(1.7383582) q[2];
sx q[2];
rz(-3.02387) q[2];
sx q[2];
rz(-0.37230834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9684978) q[1];
sx q[1];
rz(-1.2588333) q[1];
sx q[1];
rz(-0.63271823) q[1];
rz(-1.0193018) q[3];
sx q[3];
rz(-2.1658033) q[3];
sx q[3];
rz(-0.76535705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64330953) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(-2.3504284) q[2];
rz(1.6061973) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(0.78996381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20752792) q[0];
sx q[0];
rz(-2.9495033) q[0];
sx q[0];
rz(-0.15765634) q[0];
rz(-0.12570307) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(2.8368565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9466772) q[0];
sx q[0];
rz(-1.4875571) q[0];
sx q[0];
rz(-1.2437245) q[0];
rz(-1.9900121) q[2];
sx q[2];
rz(-0.97330026) q[2];
sx q[2];
rz(-2.1066163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9090635) q[1];
sx q[1];
rz(-1.705085) q[1];
sx q[1];
rz(1.568455) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35924201) q[3];
sx q[3];
rz(-0.83165681) q[3];
sx q[3];
rz(-0.11906448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2711082) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(-1.452272) q[2];
rz(-0.74635402) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-3.0284184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8993503) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(-2.6089456) q[0];
rz(-0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(-2.9740082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78412752) q[0];
sx q[0];
rz(-1.7622158) q[0];
sx q[0];
rz(0.1249868) q[0];
x q[1];
rz(0.26248502) q[2];
sx q[2];
rz(-1.0392611) q[2];
sx q[2];
rz(-2.558311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1780562) q[1];
sx q[1];
rz(-1.2884226) q[1];
sx q[1];
rz(-2.8368188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7277023) q[3];
sx q[3];
rz(-0.87855708) q[3];
sx q[3];
rz(0.62769753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8297537) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(2.6746542) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5655831) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(0.78454984) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(0.85781893) q[2];
sx q[2];
rz(-2.1265278) q[2];
sx q[2];
rz(2.3139755) q[2];
rz(-0.9355024) q[3];
sx q[3];
rz(-2.7059434) q[3];
sx q[3];
rz(-0.18999204) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
