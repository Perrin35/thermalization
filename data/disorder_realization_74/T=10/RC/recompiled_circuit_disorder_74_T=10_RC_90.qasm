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
rz(-0.98288012) q[0];
sx q[0];
rz(-2.13184) q[0];
rz(-0.17619625) q[1];
sx q[1];
rz(-2.2390525) q[1];
sx q[1];
rz(-1.2892105) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50531189) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(-2.050839) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.870954) q[2];
sx q[2];
rz(-0.025279609) q[2];
sx q[2];
rz(1.6563005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28469052) q[1];
sx q[1];
rz(-0.33387091) q[1];
sx q[1];
rz(0.88792172) q[1];
rz(-pi) q[2];
rz(0.07006499) q[3];
sx q[3];
rz(-1.8383887) q[3];
sx q[3];
rz(-2.23578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267589) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(1.82812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03598729) q[0];
sx q[0];
rz(-2.8337038) q[0];
sx q[0];
rz(0.76460989) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0560527) q[2];
sx q[2];
rz(-1.2962356) q[2];
sx q[2];
rz(0.95825125) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0158851) q[1];
sx q[1];
rz(-1.7310377) q[1];
sx q[1];
rz(1.6586152) q[1];
rz(-pi) q[2];
rz(2.1900858) q[3];
sx q[3];
rz(-3.0118239) q[3];
sx q[3];
rz(-1.4459923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-2.1874645) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(-2.4296956) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-0.16608873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27667339) q[0];
sx q[0];
rz(-1.2827946) q[0];
sx q[0];
rz(2.9898781) q[0];
rz(-pi) q[1];
rz(0.60885749) q[2];
sx q[2];
rz(-1.1470084) q[2];
sx q[2];
rz(-2.599803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7158311) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(2.3900044) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8852799) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.8410929) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(3.0680455) q[0];
rz(2.7735633) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(-2.7421943) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91037726) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(-1.7675722) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99153783) q[3];
sx q[3];
rz(-0.30617985) q[3];
sx q[3];
rz(1.3907719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0094770771) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(2.5715128) q[2];
rz(-1.8360957) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-0.10770527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62892249) q[0];
sx q[0];
rz(-1.218601) q[0];
sx q[0];
rz(-3.0432426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2568251) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(-1.9212854) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2376033) q[1];
sx q[1];
rz(-2.206982) q[1];
sx q[1];
rz(0.49828766) q[1];
rz(-pi) q[2];
rz(0.34541901) q[3];
sx q[3];
rz(-2.1072227) q[3];
sx q[3];
rz(-2.9649343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(-1.1389114) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74422979) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(0.37325333) q[0];
rz(1.594627) q[2];
sx q[2];
rz(-1.3696635) q[2];
sx q[2];
rz(-2.3383274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0000227) q[1];
sx q[1];
rz(-1.1035898) q[1];
sx q[1];
rz(-2.1187374) q[1];
x q[2];
rz(-1.059504) q[3];
sx q[3];
rz(-2.0492616) q[3];
sx q[3];
rz(-2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54414526) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(-2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25062659) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(-1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.2197781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0680925) q[0];
sx q[0];
rz(-1.4788027) q[0];
sx q[0];
rz(1.6091225) q[0];
rz(-pi) q[1];
rz(-1.9856521) q[2];
sx q[2];
rz(-1.875669) q[2];
sx q[2];
rz(0.65326234) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9417604) q[1];
sx q[1];
rz(-2.4344749) q[1];
sx q[1];
rz(1.9553528) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1620429) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(0.85588928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(-3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(2.3855456) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-2.8102002) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(0.40922871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17657875) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(0.9126419) q[0];
rz(-1.8843083) q[2];
sx q[2];
rz(-2.0743309) q[2];
sx q[2];
rz(-1.7538278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2891149) q[1];
sx q[1];
rz(-2.5534938) q[1];
sx q[1];
rz(-0.27681338) q[1];
rz(-1.3274566) q[3];
sx q[3];
rz(-2.0753324) q[3];
sx q[3];
rz(2.6017021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(-0.64152843) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.6329637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2752339) q[0];
sx q[0];
rz(-0.70952053) q[0];
sx q[0];
rz(1.8190246) q[0];
rz(0.091392322) q[2];
sx q[2];
rz(-2.1835727) q[2];
sx q[2];
rz(2.6932655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.846261) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(-2.1178513) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86767254) q[3];
sx q[3];
rz(-1.613986) q[3];
sx q[3];
rz(-1.5411351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-0.43668401) q[2];
rz(-1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(-0.25282192) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(1.3814829) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94986445) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.5230595) q[0];
rz(-1.4597458) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(-1.2088838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5188462) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-2.3153789) q[1];
x q[2];
rz(2.8548106) q[3];
sx q[3];
rz(-0.87773318) q[3];
sx q[3];
rz(1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.272841) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(-3.0957072) q[2];
sx q[2];
rz(-1.3165717) q[2];
sx q[2];
rz(-1.7236621) q[2];
rz(-1.4894555) q[3];
sx q[3];
rz(-2.079439) q[3];
sx q[3];
rz(0.54872201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
