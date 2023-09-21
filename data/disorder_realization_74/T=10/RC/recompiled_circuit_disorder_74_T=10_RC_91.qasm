OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6362808) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(1.0907537) q[0];
rz(-pi) q[1];
rz(-3.1172328) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(-2.9565405) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9413486) q[1];
sx q[1];
rz(-1.779088) q[1];
sx q[1];
rz(-1.8336481) q[1];
x q[2];
rz(1.8208002) q[3];
sx q[3];
rz(-2.8651926) q[3];
sx q[3];
rz(-1.1652511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267589) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(-2.360789) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.3134726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-2.8333089) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(0.75454933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62994176) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(-2.6444142) q[1];
x q[2];
rz(1.4649269) q[3];
sx q[3];
rz(-1.4956116) q[3];
sx q[3];
rz(-2.4014846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69497067) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(0.95412811) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057782877) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(0.35481915) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(0.16608873) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27667339) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(-2.9898781) q[0];
rz(-2.5327352) q[2];
sx q[2];
rz(-1.9945842) q[2];
sx q[2];
rz(2.599803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6925466) q[1];
sx q[1];
rz(-1.0967626) q[1];
sx q[1];
rz(-2.0708592) q[1];
rz(-pi) q[2];
rz(2.8852799) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(-3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.3004998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(-0.77484432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64530173) q[0];
sx q[0];
rz(-3.0662363) q[0];
sx q[0];
rz(-0.2199748) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11544322) q[2];
sx q[2];
rz(-2.7713159) q[2];
sx q[2];
rz(1.8625129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0171417) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(3.0567398) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1713486) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(-1.1495429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.501804) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-0.10770527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62892249) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-0.098350071) q[0];
rz(-pi) q[1];
rz(0.55107112) q[2];
sx q[2];
rz(-1.0017918) q[2];
sx q[2];
rz(-1.0702733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9039893) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(0.49828766) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1343469) q[3];
sx q[3];
rz(-1.8661024) q[3];
sx q[3];
rz(-1.9293279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(3.1048408) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(2.1246134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74422979) q[0];
sx q[0];
rz(-1.1083535) q[0];
sx q[0];
rz(-0.37325333) q[0];
x q[1];
rz(-1.594627) q[2];
sx q[2];
rz(-1.3696635) q[2];
sx q[2];
rz(2.3383274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.348154) q[1];
sx q[1];
rz(-0.70421709) q[1];
sx q[1];
rz(-0.80146316) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.059504) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54414526) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(-2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6780753) q[0];
sx q[0];
rz(-3.0419555) q[0];
sx q[0];
rz(-2.7479322) q[0];
x q[1];
rz(0.90791038) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(2.8223035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68901686) q[1];
sx q[1];
rz(-2.2170482) q[1];
sx q[1];
rz(2.8313682) q[1];
x q[2];
rz(-0.089902417) q[3];
sx q[3];
rz(-1.9780759) q[3];
sx q[3];
rz(0.75059964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(-2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(-1.1720852) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09181) q[0];
sx q[0];
rz(-0.70735332) q[0];
sx q[0];
rz(-2.011767) q[0];
rz(-pi) q[1];
rz(1.2572844) q[2];
sx q[2];
rz(-2.0743309) q[2];
sx q[2];
rz(1.3877649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2891149) q[1];
sx q[1];
rz(-2.5534938) q[1];
sx q[1];
rz(-2.8647793) q[1];
rz(-1.8141361) q[3];
sx q[3];
rz(-1.0662603) q[3];
sx q[3];
rz(2.6017021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1428712) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(0.55465737) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(-2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2471837) q[0];
sx q[0];
rz(-1.7315454) q[0];
sx q[0];
rz(-2.2649691) q[0];
rz(1.4417068) q[2];
sx q[2];
rz(-0.61868762) q[2];
sx q[2];
rz(0.29030756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.29533169) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(2.1178513) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6375332) q[3];
sx q[3];
rz(-0.70422322) q[3];
sx q[3];
rz(-3.0610386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(2.7049086) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.7601097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94986445) q[0];
sx q[0];
rz(-2.0725595) q[0];
sx q[0];
rz(-1.5230595) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6673615) q[2];
sx q[2];
rz(-1.4719163) q[2];
sx q[2];
rz(0.41254166) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30365147) q[1];
sx q[1];
rz(-0.90180574) q[1];
sx q[1];
rz(-0.53955697) q[1];
x q[2];
rz(2.8548106) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(-1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.3960719) q[2];
sx q[2];
rz(-2.8833485) q[2];
sx q[2];
rz(-1.9042518) q[2];
rz(-2.6315401) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
