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
rz(1.0097526) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50531189) q[0];
sx q[0];
rz(-1.2828151) q[0];
sx q[0];
rz(2.050839) q[0];
rz(-pi) q[1];
rz(1.5775561) q[2];
sx q[2];
rz(-1.5464371) q[2];
sx q[2];
rz(-1.3855795) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28469052) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(0.88792172) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3025769) q[3];
sx q[3];
rz(-1.6383639) q[3];
sx q[3];
rz(-2.4951631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-3.0060449) q[2];
rz(-0.1401976) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267589) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(-0.78080368) q[0];
rz(-0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.3134726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.388893) q[0];
sx q[0];
rz(-1.3503195) q[0];
sx q[0];
rz(1.3541143) q[0];
x q[1];
rz(-2.0560527) q[2];
sx q[2];
rz(-1.8453571) q[2];
sx q[2];
rz(-2.1833414) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5116509) q[1];
sx q[1];
rz(-0.18254666) q[1];
sx q[1];
rz(-2.6444142) q[1];
x q[2];
rz(2.1900858) q[3];
sx q[3];
rz(-0.12976876) q[3];
sx q[3];
rz(1.4459923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(2.4296956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(0.69586786) q[0];
rz(0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-0.16608873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9244708) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(1.0990418) q[0];
x q[1];
rz(0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(-0.49612507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7158311) q[1];
sx q[1];
rz(-0.67486963) q[1];
sx q[1];
rz(-0.75158822) q[1];
rz(-0.2563128) q[3];
sx q[3];
rz(-2.267572) q[3];
sx q[3];
rz(-3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-0.29378763) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64530173) q[0];
sx q[0];
rz(-3.0662363) q[0];
sx q[0];
rz(2.9216179) q[0];
x q[1];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.5291011) q[2];
sx q[2];
rz(0.3993984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0171417) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(-3.0567398) q[1];
x q[2];
rz(0.99153783) q[3];
sx q[3];
rz(-2.8354128) q[3];
sx q[3];
rz(-1.3907719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(-2.5715128) q[2];
rz(1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(1.6397887) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775979) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-0.10770527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340794) q[0];
sx q[0];
rz(-0.36511746) q[0];
sx q[0];
rz(1.3097197) q[0];
rz(-pi) q[1];
rz(-0.55107112) q[2];
sx q[2];
rz(-2.1398009) q[2];
sx q[2];
rz(2.0713194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4957461) q[1];
sx q[1];
rz(-1.9653814) q[1];
sx q[1];
rz(-0.87162019) q[1];
rz(-pi) q[2];
rz(0.34541901) q[3];
sx q[3];
rz(-2.1072227) q[3];
sx q[3];
rz(0.1766583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(-1.5104729) q[0];
rz(2.0026813) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(1.0169792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65360083) q[0];
sx q[0];
rz(-1.9032318) q[0];
sx q[0];
rz(2.0623273) q[0];
x q[1];
rz(-1.5469656) q[2];
sx q[2];
rz(-1.7719291) q[2];
sx q[2];
rz(-0.80326524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79343866) q[1];
sx q[1];
rz(-2.4373756) q[1];
sx q[1];
rz(-0.80146316) q[1];
x q[2];
rz(1.059504) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(0.90261501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(-2.5069359) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(0.44657782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.73262501) q[1];
sx q[1];
rz(1.9218146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0680925) q[0];
sx q[0];
rz(-1.4788027) q[0];
sx q[0];
rz(1.5324701) q[0];
x q[1];
rz(-1.1559406) q[2];
sx q[2];
rz(-1.875669) q[2];
sx q[2];
rz(2.4883303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4525758) q[1];
sx q[1];
rz(-0.92454443) q[1];
sx q[1];
rz(-2.8313682) q[1];
rz(-3.0516902) q[3];
sx q[3];
rz(-1.1635167) q[3];
sx q[3];
rz(-2.390993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-0.66292271) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-2.1989609) q[0];
sx q[0];
rz(-2.7917042) q[0];
rz(-pi) q[1];
rz(0.52493237) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(-0.33821019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2891149) q[1];
sx q[1];
rz(-2.5534938) q[1];
sx q[1];
rz(-2.8647793) q[1];
x q[2];
rz(-0.41142923) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-3.0761449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(-0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2239969) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.6329637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89440896) q[0];
sx q[0];
rz(-1.4100473) q[0];
sx q[0];
rz(-0.87662351) q[0];
rz(0.95605085) q[2];
sx q[2];
rz(-1.4960669) q[2];
sx q[2];
rz(1.9664619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.846261) q[1];
sx q[1];
rz(-1.8896855) q[1];
sx q[1];
rz(-2.1178513) q[1];
x q[2];
rz(-0.86767254) q[3];
sx q[3];
rz(-1.5276067) q[3];
sx q[3];
rz(1.6004576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9987954) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(-1.3302594) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(2.57634) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.5230595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21391602) q[2];
sx q[2];
rz(-0.48366085) q[2];
sx q[2];
rz(-2.1733401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6773721) q[1];
sx q[1];
rz(-0.83253011) q[1];
sx q[1];
rz(2.1470451) q[1];
rz(-pi) q[2];
rz(1.2425209) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(-1.5433943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(1.8252774) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(-1.6521372) q[3];
sx q[3];
rz(-1.0621536) q[3];
sx q[3];
rz(-2.5928706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
