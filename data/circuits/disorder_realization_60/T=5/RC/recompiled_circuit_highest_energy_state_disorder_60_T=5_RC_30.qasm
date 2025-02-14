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
rz(5.3661348) q[0];
sx q[0];
rz(9.6714749) q[0];
rz(0.84423455) q[1];
sx q[1];
rz(-1.2935473) q[1];
sx q[1];
rz(-1.8389314) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58606746) q[0];
sx q[0];
rz(-1.2159776) q[0];
sx q[0];
rz(0.3094425) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6370722) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(-2.1736886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0595536) q[1];
sx q[1];
rz(-2.5055725) q[1];
sx q[1];
rz(0.18853699) q[1];
rz(-pi) q[2];
rz(0.10860034) q[3];
sx q[3];
rz(-1.425877) q[3];
sx q[3];
rz(-2.8371881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.777433) q[2];
sx q[2];
rz(-1.0182764) q[2];
sx q[2];
rz(1.9470661) q[2];
rz(1.901769) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(-1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0251004) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-2.4745353) q[0];
rz(2.8705813) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8264814) q[0];
sx q[0];
rz(-1.3024863) q[0];
sx q[0];
rz(-2.2196191) q[0];
x q[1];
rz(-1.1863883) q[2];
sx q[2];
rz(-1.5509275) q[2];
sx q[2];
rz(-1.5569614) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42506524) q[1];
sx q[1];
rz(-2.7924573) q[1];
sx q[1];
rz(-1.8340111) q[1];
rz(0.74852826) q[3];
sx q[3];
rz(-1.0077884) q[3];
sx q[3];
rz(2.9070118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(2.8957497) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(-1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(-0.47897419) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(-1.825038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265357) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(-2.3903484) q[0];
x q[1];
rz(2.9390305) q[2];
sx q[2];
rz(-1.1656467) q[2];
sx q[2];
rz(-0.7140401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90637246) q[1];
sx q[1];
rz(-0.89284578) q[1];
sx q[1];
rz(-1.9560019) q[1];
x q[2];
rz(-2.4495921) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5548744) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(0.47373104) q[2];
rz(-2.6321865) q[3];
sx q[3];
rz(-1.3207057) q[3];
sx q[3];
rz(1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944775) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(1.0152869) q[0];
rz(-0.21233755) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(-2.7755348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.191712) q[0];
sx q[0];
rz(-1.0832583) q[0];
sx q[0];
rz(-0.78678188) q[0];
rz(-pi) q[1];
rz(-2.9061142) q[2];
sx q[2];
rz(-1.7120541) q[2];
sx q[2];
rz(1.3154674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8541222) q[1];
sx q[1];
rz(-2.4123998) q[1];
sx q[1];
rz(-2.6602938) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96694209) q[3];
sx q[3];
rz(-1.8181385) q[3];
sx q[3];
rz(0.027146904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(0.69501957) q[2];
rz(2.2906637) q[3];
sx q[3];
rz(-1.7780108) q[3];
sx q[3];
rz(-1.883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(0.31546053) q[0];
rz(0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-0.42627898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.016092) q[0];
sx q[0];
rz(-2.3707485) q[0];
sx q[0];
rz(1.9941814) q[0];
x q[1];
rz(0.23172863) q[2];
sx q[2];
rz(-1.2905057) q[2];
sx q[2];
rz(-2.0120399) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1412999) q[1];
sx q[1];
rz(-2.0054549) q[1];
sx q[1];
rz(-3.0463957) q[1];
x q[2];
rz(1.0067389) q[3];
sx q[3];
rz(-0.56017471) q[3];
sx q[3];
rz(-0.87040802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2686501) q[2];
sx q[2];
rz(-0.23317569) q[2];
sx q[2];
rz(0.028701393) q[2];
rz(2.9305693) q[3];
sx q[3];
rz(-0.70091453) q[3];
sx q[3];
rz(2.4215462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215304) q[0];
sx q[0];
rz(-1.2380607) q[0];
sx q[0];
rz(-2.9826214) q[0];
rz(-2.4863415) q[1];
sx q[1];
rz(-2.7924004) q[1];
sx q[1];
rz(1.5370625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55729729) q[0];
sx q[0];
rz(-2.0157776) q[0];
sx q[0];
rz(2.2792982) q[0];
rz(-pi) q[1];
rz(-1.9150429) q[2];
sx q[2];
rz(-1.9166528) q[2];
sx q[2];
rz(0.49357061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3517396) q[1];
sx q[1];
rz(-0.56507128) q[1];
sx q[1];
rz(1.5365824) q[1];
x q[2];
rz(0.37015583) q[3];
sx q[3];
rz(-1.8612544) q[3];
sx q[3];
rz(-0.53122093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0663466) q[2];
sx q[2];
rz(-1.9074351) q[2];
sx q[2];
rz(-2.9998903) q[2];
rz(3.1250478) q[3];
sx q[3];
rz(-2.8576272) q[3];
sx q[3];
rz(-0.56104463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(0.91019994) q[0];
rz(-2.3987112) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(1.9076617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54523477) q[0];
sx q[0];
rz(-1.4995793) q[0];
sx q[0];
rz(-0.13990732) q[0];
x q[1];
rz(-1.1612438) q[2];
sx q[2];
rz(-1.3245593) q[2];
sx q[2];
rz(0.72803942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2382395) q[1];
sx q[1];
rz(-0.80447865) q[1];
sx q[1];
rz(-2.5969905) q[1];
rz(-pi) q[2];
rz(-0.55039946) q[3];
sx q[3];
rz(-0.66662753) q[3];
sx q[3];
rz(-0.77778274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10741216) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(2.2255955) q[2];
rz(-0.25137526) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995572) q[0];
sx q[0];
rz(-2.636203) q[0];
sx q[0];
rz(-3.1009951) q[0];
rz(1.1740855) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(-2.0814799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745251) q[0];
sx q[0];
rz(-0.39345783) q[0];
sx q[0];
rz(2.8889546) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1218704) q[2];
sx q[2];
rz(-1.4547299) q[2];
sx q[2];
rz(-0.20360064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9684978) q[1];
sx q[1];
rz(-1.2588333) q[1];
sx q[1];
rz(0.63271823) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1222909) q[3];
sx q[3];
rz(-0.97578933) q[3];
sx q[3];
rz(0.76535705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4982831) q[2];
sx q[2];
rz(-2.0556367) q[2];
sx q[2];
rz(2.3504284) q[2];
rz(1.5353954) q[3];
sx q[3];
rz(-2.199506) q[3];
sx q[3];
rz(-0.78996381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(-2.9839363) q[0];
rz(-0.12570307) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(-2.8368565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7939111) q[0];
sx q[0];
rz(-1.2448989) q[0];
sx q[0];
rz(-0.087875701) q[0];
rz(-pi) q[1];
rz(1.1515806) q[2];
sx q[2];
rz(-2.1682924) q[2];
sx q[2];
rz(2.1066163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.803639) q[1];
sx q[1];
rz(-1.5684761) q[1];
sx q[1];
rz(3.0073037) q[1];
rz(-pi) q[2];
rz(-0.79873884) q[3];
sx q[3];
rz(-1.307957) q[3];
sx q[3];
rz(-1.2039487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8704845) q[2];
sx q[2];
rz(-1.5799589) q[2];
sx q[2];
rz(1.452272) q[2];
rz(-2.3952386) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(-0.11317429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8993503) q[0];
sx q[0];
rz(-0.32651383) q[0];
sx q[0];
rz(0.53264701) q[0];
rz(2.7592754) q[1];
sx q[1];
rz(-0.89235726) q[1];
sx q[1];
rz(2.9740082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3310247) q[0];
sx q[0];
rz(-1.6934886) q[0];
sx q[0];
rz(-1.3779089) q[0];
rz(-pi) q[1];
rz(-2.1176867) q[2];
sx q[2];
rz(-1.7963855) q[2];
sx q[2];
rz(1.1228665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6217612) q[1];
sx q[1];
rz(-1.863136) q[1];
sx q[1];
rz(-1.8660493) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3066809) q[3];
sx q[3];
rz(-1.8855699) q[3];
sx q[3];
rz(-1.2164468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8297537) q[2];
sx q[2];
rz(-1.6646174) q[2];
sx q[2];
rz(2.6746542) q[2];
rz(-1.1030819) q[3];
sx q[3];
rz(-1.9481877) q[3];
sx q[3];
rz(-1.9095437) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(-2.3570428) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(-0.81132728) q[2];
sx q[2];
rz(-0.87292508) q[2];
sx q[2];
rz(1.2909918) q[2];
rz(0.26950027) q[3];
sx q[3];
rz(-1.2242347) q[3];
sx q[3];
rz(0.4927529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
