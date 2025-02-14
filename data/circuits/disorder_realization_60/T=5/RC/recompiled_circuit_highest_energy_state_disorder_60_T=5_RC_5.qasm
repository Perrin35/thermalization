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
rz(3.9858272) q[1];
sx q[1];
rz(1.2935473) q[1];
sx q[1];
rz(10.727439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87410754) q[0];
sx q[0];
rz(-1.8603854) q[0];
sx q[0];
rz(1.9417614) q[0];
rz(2.2751426) q[2];
sx q[2];
rz(-1.1725008) q[2];
sx q[2];
rz(0.9212538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0595536) q[1];
sx q[1];
rz(-2.5055725) q[1];
sx q[1];
rz(-2.9530557) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10860034) q[3];
sx q[3];
rz(-1.7157156) q[3];
sx q[3];
rz(-0.30440457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.777433) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(1.9470661) q[2];
rz(1.2398237) q[3];
sx q[3];
rz(-1.4907962) q[3];
sx q[3];
rz(-1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(0.27101135) q[1];
sx q[1];
rz(-1.4457116) q[1];
sx q[1];
rz(1.0557231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5496084) q[0];
sx q[0];
rz(-0.69463629) q[0];
sx q[0];
rz(-1.9978092) q[0];
x q[1];
rz(-1.6237359) q[2];
sx q[2];
rz(-0.38489562) q[2];
sx q[2];
rz(-0.062907779) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89773399) q[1];
sx q[1];
rz(-1.4816726) q[1];
sx q[1];
rz(-1.9088163) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74774489) q[3];
sx q[3];
rz(-2.2390331) q[3];
sx q[3];
rz(1.2838319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51862741) q[2];
sx q[2];
rz(-2.3845606) q[2];
sx q[2];
rz(-0.61678994) q[2];
rz(-2.8957497) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(1.8776548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(-0.1804633) q[0];
rz(-2.6626185) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(1.825038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0468123) q[0];
sx q[0];
rz(-2.1553596) q[0];
sx q[0];
rz(0.78740904) q[0];
rz(-pi) q[1];
rz(1.1321849) q[2];
sx q[2];
rz(-0.45044611) q[2];
sx q[2];
rz(0.23368719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7261947) q[1];
sx q[1];
rz(-1.8677999) q[1];
sx q[1];
rz(-2.4261977) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4495921) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(0.52317515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5548744) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(2.6678616) q[2];
rz(0.50940618) q[3];
sx q[3];
rz(-1.3207057) q[3];
sx q[3];
rz(1.9853076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2471152) q[0];
sx q[0];
rz(-0.04016567) q[0];
sx q[0];
rz(2.1263057) q[0];
rz(-0.21233755) q[1];
sx q[1];
rz(-1.5539853) q[1];
sx q[1];
rz(-0.36605787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057945874) q[0];
sx q[0];
rz(-2.2444631) q[0];
sx q[0];
rz(-2.4988453) q[0];
rz(-pi) q[1];
rz(0.23547844) q[2];
sx q[2];
rz(-1.4295385) q[2];
sx q[2];
rz(-1.3154674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6547114) q[1];
sx q[1];
rz(-1.8843448) q[1];
sx q[1];
rz(0.66982072) q[1];
rz(1.9892392) q[3];
sx q[3];
rz(-2.4949346) q[3];
sx q[3];
rz(1.2025976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66146835) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(-0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(0.42627898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5770059) q[0];
sx q[0];
rz(-2.2590911) q[0];
sx q[0];
rz(2.7618963) q[0];
rz(-0.89740173) q[2];
sx q[2];
rz(-0.36172141) q[2];
sx q[2];
rz(2.7182686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0002928) q[1];
sx q[1];
rz(-2.0054549) q[1];
sx q[1];
rz(-0.095196916) q[1];
rz(1.0067389) q[3];
sx q[3];
rz(-2.5814179) q[3];
sx q[3];
rz(-2.2711846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8729426) q[2];
sx q[2];
rz(-0.23317569) q[2];
sx q[2];
rz(-0.028701393) q[2];
rz(-2.9305693) q[3];
sx q[3];
rz(-0.70091453) q[3];
sx q[3];
rz(-2.4215462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92006224) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(0.15897121) q[0];
rz(0.65525118) q[1];
sx q[1];
rz(-0.34919229) q[1];
sx q[1];
rz(1.6045301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66012525) q[0];
sx q[0];
rz(-2.1985538) q[0];
sx q[0];
rz(-2.5808236) q[0];
x q[1];
rz(2.7760024) q[2];
sx q[2];
rz(-1.8938918) q[2];
sx q[2];
rz(-1.9434203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3517396) q[1];
sx q[1];
rz(-2.5765214) q[1];
sx q[1];
rz(1.5365824) q[1];
rz(-pi) q[2];
rz(0.37015583) q[3];
sx q[3];
rz(-1.8612544) q[3];
sx q[3];
rz(-0.53122093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0663466) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(-0.1417024) q[2];
rz(3.1250478) q[3];
sx q[3];
rz(-0.2839655) q[3];
sx q[3];
rz(0.56104463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905106) q[0];
sx q[0];
rz(-0.53618479) q[0];
sx q[0];
rz(-0.91019994) q[0];
rz(-0.74288145) q[1];
sx q[1];
rz(-2.1292834) q[1];
sx q[1];
rz(-1.9076617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.493422) q[0];
sx q[0];
rz(-0.1568846) q[0];
sx q[0];
rz(2.6687401) q[0];
x q[1];
rz(1.1612438) q[2];
sx q[2];
rz(-1.8170333) q[2];
sx q[2];
rz(0.72803942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2382395) q[1];
sx q[1];
rz(-0.80447865) q[1];
sx q[1];
rz(2.5969905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59071006) q[3];
sx q[3];
rz(-1.241467) q[3];
sx q[3];
rz(1.2423837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0341805) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-2.2255955) q[2];
rz(-0.25137526) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6420355) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(-3.1009951) q[0];
rz(1.9675072) q[1];
sx q[1];
rz(-0.65316713) q[1];
sx q[1];
rz(-2.0814799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030329372) q[0];
sx q[0];
rz(-1.666774) q[0];
sx q[0];
rz(-2.7594271) q[0];
rz(-1.7383582) q[2];
sx q[2];
rz(-3.02387) q[2];
sx q[2];
rz(-2.7692843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0012944) q[1];
sx q[1];
rz(-2.4457275) q[1];
sx q[1];
rz(-0.49927478) q[1];
rz(-pi) q[2];
rz(-1.0193018) q[3];
sx q[3];
rz(-0.97578933) q[3];
sx q[3];
rz(-2.3762356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(-2.3504284) q[2];
rz(-1.6061973) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(-2.3516288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20752792) q[0];
sx q[0];
rz(-2.9495033) q[0];
sx q[0];
rz(-0.15765634) q[0];
rz(-3.0158896) q[1];
sx q[1];
rz(-1.1748284) q[1];
sx q[1];
rz(-0.30473614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19491542) q[0];
sx q[0];
rz(-1.6540356) q[0];
sx q[0];
rz(-1.8978682) q[0];
x q[1];
rz(1.9900121) q[2];
sx q[2];
rz(-0.97330026) q[2];
sx q[2];
rz(-1.0349764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9090635) q[1];
sx q[1];
rz(-1.705085) q[1];
sx q[1];
rz(-1.5731377) q[1];
x q[2];
rz(2.3428538) q[3];
sx q[3];
rz(-1.8336356) q[3];
sx q[3];
rz(-1.9376439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8704845) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(-1.452272) q[2];
rz(-0.74635402) q[3];
sx q[3];
rz(-1.4515667) q[3];
sx q[3];
rz(0.11317429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.2492354) q[1];
sx q[1];
rz(-2.9740082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3310247) q[0];
sx q[0];
rz(-1.448104) q[0];
sx q[0];
rz(1.7636838) q[0];
rz(-1.0239059) q[2];
sx q[2];
rz(-1.3452072) q[2];
sx q[2];
rz(1.1228665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3321489) q[1];
sx q[1];
rz(-2.7291307) q[1];
sx q[1];
rz(2.3730538) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7277023) q[3];
sx q[3];
rz(-0.87855708) q[3];
sx q[3];
rz(-2.5138951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8297537) q[2];
sx q[2];
rz(-1.6646174) q[2];
sx q[2];
rz(-0.46693841) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.1934049) q[3];
sx q[3];
rz(-1.232049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5760096) q[0];
sx q[0];
rz(-2.3143815) q[0];
sx q[0];
rz(-2.6630493) q[0];
rz(0.78454984) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(-0.85781893) q[2];
sx q[2];
rz(-1.0150649) q[2];
sx q[2];
rz(-0.82761717) q[2];
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
