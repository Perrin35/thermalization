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
rz(3.5340478) q[0];
sx q[0];
rz(3.4253149) q[0];
sx q[0];
rz(14.520159) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(5.0168283) q[1];
sx q[1];
rz(14.556806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2806704) q[0];
sx q[0];
rz(-0.09083561) q[0];
sx q[0];
rz(-1.8152186) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63090905) q[2];
sx q[2];
rz(-2.563347) q[2];
sx q[2];
rz(-1.9672245) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6952327) q[1];
sx q[1];
rz(-0.44539136) q[1];
sx q[1];
rz(2.8421418) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2972732) q[3];
sx q[3];
rz(-1.0680388) q[3];
sx q[3];
rz(1.1532602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2426131) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3886609) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(-3.062124) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-0.4963378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74018986) q[0];
sx q[0];
rz(-1.8215792) q[0];
sx q[0];
rz(-2.0987505) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9251516) q[2];
sx q[2];
rz(-0.34578824) q[2];
sx q[2];
rz(-1.3160637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1524758) q[1];
sx q[1];
rz(-0.82738872) q[1];
sx q[1];
rz(-2.6914196) q[1];
rz(-pi) q[2];
rz(-3.1383508) q[3];
sx q[3];
rz(-1.7326983) q[3];
sx q[3];
rz(1.1461794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(0.17862865) q[2];
rz(-1.6102839) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71132201) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(-1.3871644) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(1.656146) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324155) q[0];
sx q[0];
rz(-2.2542291) q[0];
sx q[0];
rz(1.528549) q[0];
x q[1];
rz(1.9958271) q[2];
sx q[2];
rz(-2.3299937) q[2];
sx q[2];
rz(1.5852864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8456257) q[1];
sx q[1];
rz(-2.6085771) q[1];
sx q[1];
rz(-2.19341) q[1];
rz(-pi) q[2];
rz(-0.36538045) q[3];
sx q[3];
rz(-0.93910223) q[3];
sx q[3];
rz(-1.9393348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2014655) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(2.6489769) q[2];
rz(2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(-2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1241322) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(0.19805743) q[0];
rz(2.4824704) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448387) q[0];
sx q[0];
rz(-1.0403883) q[0];
sx q[0];
rz(-1.5317937) q[0];
rz(-pi) q[1];
rz(-2.9022759) q[2];
sx q[2];
rz(-1.2916358) q[2];
sx q[2];
rz(2.5365732) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9777364) q[1];
sx q[1];
rz(-1.2173554) q[1];
sx q[1];
rz(2.5720235) q[1];
rz(3.1048685) q[3];
sx q[3];
rz(-2.1154599) q[3];
sx q[3];
rz(2.7721318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78843242) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(1.2737459) q[2];
rz(1.5431917) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(0.25632349) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0565599) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(3.0388167) q[0];
rz(3.00434) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(1.7030254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6213087) q[0];
sx q[0];
rz(-2.1946593) q[0];
sx q[0];
rz(-0.72786967) q[0];
rz(-pi) q[1];
rz(1.6167783) q[2];
sx q[2];
rz(-2.0334938) q[2];
sx q[2];
rz(1.3563987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66617353) q[1];
sx q[1];
rz(-1.3120781) q[1];
sx q[1];
rz(-2.9889876) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6827312) q[3];
sx q[3];
rz(-1.4482968) q[3];
sx q[3];
rz(-1.732633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(-2.4414818) q[2];
rz(-0.89811283) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20823088) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5640278) q[0];
rz(2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(0.97445625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12688836) q[0];
sx q[0];
rz(-3.1278353) q[0];
sx q[0];
rz(0.62928726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48759385) q[2];
sx q[2];
rz(-1.7238443) q[2];
sx q[2];
rz(2.7014554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8169957) q[1];
sx q[1];
rz(-1.5600509) q[1];
sx q[1];
rz(0.50159494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10147126) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(1.056788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4727891) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(-2.2514553) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(0.11372456) q[0];
rz(1.744489) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(-3.1166335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.002288) q[0];
sx q[0];
rz(-1.6253259) q[0];
sx q[0];
rz(1.789458) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8005902) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(0.18358809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.125122) q[1];
sx q[1];
rz(-0.85713398) q[1];
sx q[1];
rz(0.086507052) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9437441) q[3];
sx q[3];
rz(-1.6258844) q[3];
sx q[3];
rz(2.4428989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52594841) q[2];
sx q[2];
rz(-1.4333466) q[2];
sx q[2];
rz(2.136266) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.7361448) q[3];
sx q[3];
rz(2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-1.5263851) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(-1.1108387) q[0];
rz(-0.088518294) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(-1.7875338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8691886) q[0];
sx q[0];
rz(-1.4581632) q[0];
sx q[0];
rz(-2.3478481) q[0];
x q[1];
rz(0.1536151) q[2];
sx q[2];
rz(-1.3014463) q[2];
sx q[2];
rz(0.17698032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14446324) q[1];
sx q[1];
rz(-1.6465997) q[1];
sx q[1];
rz(-0.24667011) q[1];
rz(-pi) q[2];
x q[2];
rz(2.380221) q[3];
sx q[3];
rz(-2.0132445) q[3];
sx q[3];
rz(0.82514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-2.3044846) q[2];
sx q[2];
rz(0.38541547) q[2];
rz(2.391583) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(-3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239922) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(-2.3093834) q[0];
rz(1.0819134) q[1];
sx q[1];
rz(-0.81580201) q[1];
sx q[1];
rz(1.5218511) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8502626) q[0];
sx q[0];
rz(-1.9614152) q[0];
sx q[0];
rz(-2.0807092) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3004901) q[2];
sx q[2];
rz(-0.37174598) q[2];
sx q[2];
rz(-1.8018307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.202692) q[1];
sx q[1];
rz(-1.0098411) q[1];
sx q[1];
rz(2.2748018) q[1];
x q[2];
rz(-0.91401871) q[3];
sx q[3];
rz(-1.3784215) q[3];
sx q[3];
rz(-0.094948204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(-1.291412) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(-3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.130403) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(2.3977872) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(-1.1415175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7615712) q[0];
sx q[0];
rz(-2.0788045) q[0];
sx q[0];
rz(-2.6045799) q[0];
rz(-pi) q[1];
rz(-1.1820993) q[2];
sx q[2];
rz(-0.98983374) q[2];
sx q[2];
rz(-0.81755359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1591349) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(-2.3680658) q[1];
rz(-pi) q[2];
rz(-1.7520245) q[3];
sx q[3];
rz(-0.93398213) q[3];
sx q[3];
rz(-0.94094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49012524) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(0.30533314) q[2];
rz(-1.7134282) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74844985) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(2.582386) q[1];
sx q[1];
rz(-1.3312085) q[1];
sx q[1];
rz(-2.8937173) q[1];
rz(2.7087173) q[2];
sx q[2];
rz(-2.4268711) q[2];
sx q[2];
rz(-2.101244) q[2];
rz(-1.4356267) q[3];
sx q[3];
rz(-1.6886938) q[3];
sx q[3];
rz(2.0791048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
