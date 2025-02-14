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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5555252) q[0];
sx q[0];
rz(-1.925615) q[0];
sx q[0];
rz(2.8321502) q[0];
rz(2.6370722) q[2];
sx q[2];
rz(-0.93122831) q[2];
sx q[2];
rz(2.1736886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9907551) q[1];
sx q[1];
rz(-2.1937943) q[1];
sx q[1];
rz(1.4332818) q[1];
rz(-0.10860034) q[3];
sx q[3];
rz(-1.7157156) q[3];
sx q[3];
rz(0.30440457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3641597) q[2];
sx q[2];
rz(-2.1233163) q[2];
sx q[2];
rz(1.9470661) q[2];
rz(-1.2398237) q[3];
sx q[3];
rz(-1.6507964) q[3];
sx q[3];
rz(-1.7279846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251004) q[0];
sx q[0];
rz(-2.0142374) q[0];
sx q[0];
rz(-0.66705739) q[0];
rz(-0.27101135) q[1];
sx q[1];
rz(-1.695881) q[1];
sx q[1];
rz(1.0557231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8264814) q[0];
sx q[0];
rz(-1.8391063) q[0];
sx q[0];
rz(0.92197355) q[0];
rz(-pi) q[1];
rz(-0.021432545) q[2];
sx q[2];
rz(-1.9551245) q[2];
sx q[2];
rz(-0.0057980428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89773399) q[1];
sx q[1];
rz(-1.6599201) q[1];
sx q[1];
rz(1.2327764) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2818774) q[3];
sx q[3];
rz(-2.1840349) q[3];
sx q[3];
rz(2.2656253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6229652) q[2];
sx q[2];
rz(-0.75703207) q[2];
sx q[2];
rz(0.61678994) q[2];
rz(0.24584298) q[3];
sx q[3];
rz(-1.5893693) q[3];
sx q[3];
rz(-1.2639379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.371405) q[0];
sx q[0];
rz(-2.4076732) q[0];
sx q[0];
rz(2.9611294) q[0];
rz(-2.6626185) q[1];
sx q[1];
rz(-2.6909747) q[1];
sx q[1];
rz(1.825038) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265357) q[0];
sx q[0];
rz(-0.94158544) q[0];
sx q[0];
rz(2.3903484) q[0];
rz(-pi) q[1];
rz(0.20256217) q[2];
sx q[2];
rz(-1.1656467) q[2];
sx q[2];
rz(-2.4275526) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2352202) q[1];
sx q[1];
rz(-2.2487469) q[1];
sx q[1];
rz(-1.9560019) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69200055) q[3];
sx q[3];
rz(-0.82847906) q[3];
sx q[3];
rz(2.6184175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5548744) q[2];
sx q[2];
rz(-0.46963936) q[2];
sx q[2];
rz(0.47373104) q[2];
rz(2.6321865) q[3];
sx q[3];
rz(-1.820887) q[3];
sx q[3];
rz(-1.1562851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944775) q[0];
sx q[0];
rz(-3.101427) q[0];
sx q[0];
rz(2.1263057) q[0];
rz(-2.9292551) q[1];
sx q[1];
rz(-1.5876074) q[1];
sx q[1];
rz(2.7755348) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057945874) q[0];
sx q[0];
rz(-0.89712954) q[0];
sx q[0];
rz(-0.64274733) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54738657) q[2];
sx q[2];
rz(-2.8676708) q[2];
sx q[2];
rz(-0.78597921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8541222) q[1];
sx q[1];
rz(-2.4123998) q[1];
sx q[1];
rz(0.4812989) q[1];
rz(-pi) q[2];
rz(-0.29764953) q[3];
sx q[3];
rz(-2.1538056) q[3];
sx q[3];
rz(1.7109554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66146835) q[2];
sx q[2];
rz(-1.8331567) q[2];
sx q[2];
rz(0.69501957) q[2];
rz(0.85092893) q[3];
sx q[3];
rz(-1.3635819) q[3];
sx q[3];
rz(1.2581576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224391) q[0];
sx q[0];
rz(-2.9762486) q[0];
sx q[0];
rz(0.31546053) q[0];
rz(-0.28469616) q[1];
sx q[1];
rz(-1.5226786) q[1];
sx q[1];
rz(-2.7153137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3836455) q[0];
sx q[0];
rz(-1.28048) q[0];
sx q[0];
rz(2.2955672) q[0];
rz(-pi) q[1];
rz(-1.8583723) q[2];
sx q[2];
rz(-1.3482665) q[2];
sx q[2];
rz(2.6351647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1412999) q[1];
sx q[1];
rz(-1.1361377) q[1];
sx q[1];
rz(-0.095196916) q[1];
x q[2];
rz(0.3235281) q[3];
sx q[3];
rz(-1.1051205) q[3];
sx q[3];
rz(2.9126008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8729426) q[2];
sx q[2];
rz(-0.23317569) q[2];
sx q[2];
rz(-3.1128913) q[2];
rz(2.9305693) q[3];
sx q[3];
rz(-2.4406781) q[3];
sx q[3];
rz(0.72004643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92006224) q[0];
sx q[0];
rz(-1.9035319) q[0];
sx q[0];
rz(2.9826214) q[0];
rz(-2.4863415) q[1];
sx q[1];
rz(-2.7924004) q[1];
sx q[1];
rz(-1.6045301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6626016) q[0];
sx q[0];
rz(-2.3259386) q[0];
sx q[0];
rz(-2.2032477) q[0];
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
rz(-0.75204471) q[1];
sx q[1];
rz(-1.5524781) q[1];
sx q[1];
rz(-2.1356029) q[1];
x q[2];
rz(-1.2605242) q[3];
sx q[3];
rz(-1.2168443) q[3];
sx q[3];
rz(-1.9913395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0752461) q[2];
sx q[2];
rz(-1.2341576) q[2];
sx q[2];
rz(-0.1417024) q[2];
rz(-0.016544841) q[3];
sx q[3];
rz(-0.2839655) q[3];
sx q[3];
rz(0.56104463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905106) q[0];
sx q[0];
rz(-2.6054079) q[0];
sx q[0];
rz(2.2313927) q[0];
rz(0.74288145) q[1];
sx q[1];
rz(-1.0123092) q[1];
sx q[1];
rz(1.2339309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155411) q[0];
sx q[0];
rz(-1.7103467) q[0];
sx q[0];
rz(-1.6427137) q[0];
rz(-pi) q[1];
rz(-2.1338314) q[2];
sx q[2];
rz(-2.6673311) q[2];
sx q[2];
rz(-0.33111085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0651111) q[1];
sx q[1];
rz(-1.9533159) q[1];
sx q[1];
rz(-0.72648813) q[1];
rz(-0.55039946) q[3];
sx q[3];
rz(-2.4749651) q[3];
sx q[3];
rz(-2.3638099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10741216) q[2];
sx q[2];
rz(-1.5739406) q[2];
sx q[2];
rz(-2.2255955) q[2];
rz(-2.8902174) q[3];
sx q[3];
rz(-0.97094691) q[3];
sx q[3];
rz(-2.796252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4995572) q[0];
sx q[0];
rz(-0.50538969) q[0];
sx q[0];
rz(3.1009951) q[0];
rz(1.1740855) q[1];
sx q[1];
rz(-2.4884255) q[1];
sx q[1];
rz(1.0601128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1112633) q[0];
sx q[0];
rz(-1.4748186) q[0];
sx q[0];
rz(-2.7594271) q[0];
x q[1];
rz(-0.019722299) q[2];
sx q[2];
rz(-1.6868627) q[2];
sx q[2];
rz(-0.20360064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5225142) q[1];
sx q[1];
rz(-0.97303094) q[1];
sx q[1];
rz(-1.9512216) q[1];
rz(-0.67146639) q[3];
sx q[3];
rz(-1.1219624) q[3];
sx q[3];
rz(-1.1374813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4982831) q[2];
sx q[2];
rz(-1.085956) q[2];
sx q[2];
rz(-0.79116428) q[2];
rz(-1.6061973) q[3];
sx q[3];
rz(-0.94208661) q[3];
sx q[3];
rz(0.78996381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9340647) q[0];
sx q[0];
rz(-0.19208935) q[0];
sx q[0];
rz(2.9839363) q[0];
rz(0.12570307) q[1];
sx q[1];
rz(-1.9667642) q[1];
sx q[1];
rz(-2.8368565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7939111) q[0];
sx q[0];
rz(-1.2448989) q[0];
sx q[0];
rz(-3.053717) q[0];
x q[1];
rz(2.6025099) q[2];
sx q[2];
rz(-2.4266908) q[2];
sx q[2];
rz(1.7049005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9090635) q[1];
sx q[1];
rz(-1.4365077) q[1];
sx q[1];
rz(-1.568455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35924201) q[3];
sx q[3];
rz(-2.3099358) q[3];
sx q[3];
rz(3.0225282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8704845) q[2];
sx q[2];
rz(-1.5616337) q[2];
sx q[2];
rz(1.452272) q[2];
rz(2.3952386) q[3];
sx q[3];
rz(-1.690026) q[3];
sx q[3];
rz(-0.11317429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8993503) q[0];
sx q[0];
rz(-2.8150788) q[0];
sx q[0];
rz(0.53264701) q[0];
rz(-0.38231725) q[1];
sx q[1];
rz(-2.2492354) q[1];
sx q[1];
rz(-2.9740082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3574651) q[0];
sx q[0];
rz(-1.7622158) q[0];
sx q[0];
rz(-3.0166059) q[0];
x q[1];
rz(-1.1551935) q[2];
sx q[2];
rz(-0.58718449) q[2];
sx q[2];
rz(-0.09584643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96353647) q[1];
sx q[1];
rz(-1.85317) q[1];
sx q[1];
rz(-2.8368188) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0224277) q[3];
sx q[3];
rz(-2.3529624) q[3];
sx q[3];
rz(-3.1166704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31183895) q[2];
sx q[2];
rz(-1.4769752) q[2];
sx q[2];
rz(-2.6746542) q[2];
rz(2.0385108) q[3];
sx q[3];
rz(-1.1934049) q[3];
sx q[3];
rz(1.9095437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5655831) q[0];
sx q[0];
rz(-0.82721114) q[0];
sx q[0];
rz(0.4785434) q[0];
rz(-2.3570428) q[1];
sx q[1];
rz(-2.7922834) q[1];
sx q[1];
rz(-2.2193411) q[1];
rz(2.2837737) q[2];
sx q[2];
rz(-1.0150649) q[2];
sx q[2];
rz(-0.82761717) q[2];
rz(-2.2060903) q[3];
sx q[3];
rz(-0.43564921) q[3];
sx q[3];
rz(2.9516006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
