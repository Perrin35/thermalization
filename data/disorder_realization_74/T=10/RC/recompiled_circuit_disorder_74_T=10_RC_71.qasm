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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754817) q[0];
sx q[0];
rz(-0.55395836) q[0];
sx q[0];
rz(-2.1411116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024359811) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(0.18505219) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42613712) q[1];
sx q[1];
rz(-1.313756) q[1];
sx q[1];
rz(-0.21547683) q[1];
rz(-pi) q[2];
rz(-1.3207924) q[3];
sx q[3];
rz(-0.27640009) q[3];
sx q[3];
rz(-1.9763415) q[3];
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
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(2.8813598) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.27539) q[0];
sx q[0];
rz(-1.35944) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-pi) q[1];
rz(2.0560527) q[2];
sx q[2];
rz(-1.2962356) q[2];
sx q[2];
rz(0.95825125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62994176) q[1];
sx q[1];
rz(-2.959046) q[1];
sx q[1];
rz(-0.49717848) q[1];
rz(3.0659862) q[3];
sx q[3];
rz(-1.4652271) q[3];
sx q[3];
rz(0.82270634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(2.4296956) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9244708) q[0];
sx q[0];
rz(-0.32453254) q[0];
sx q[0];
rz(-2.0425509) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4737349) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(-2.6454676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42576158) q[1];
sx q[1];
rz(-2.466723) q[1];
sx q[1];
rz(-2.3900044) q[1];
rz(-pi) q[2];
rz(2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0825901) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(-0.77484432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7168717) q[0];
sx q[0];
rz(-1.6443335) q[0];
sx q[0];
rz(1.5872692) q[0];
rz(-0.36802937) q[2];
sx q[2];
rz(-1.6124915) q[2];
sx q[2];
rz(2.7421943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2312154) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(1.7675722) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3121698) q[3];
sx q[3];
rz(-1.4050409) q[3];
sx q[3];
rz(0.37763077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0094770771) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-2.9602125) q[0];
rz(0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(3.0338874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5126702) q[0];
sx q[0];
rz(-1.9229917) q[0];
sx q[0];
rz(-0.098350071) q[0];
rz(-pi) q[1];
rz(-0.88476752) q[2];
sx q[2];
rz(-0.77026412) q[2];
sx q[2];
rz(1.9212854) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6458466) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-0.87162019) q[1];
rz(-pi) q[2];
rz(2.7961736) q[3];
sx q[3];
rz(-2.1072227) q[3];
sx q[3];
rz(2.9649343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(0.28953141) q[2];
rz(3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(0.010820476) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058379563) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(1.0169792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4879918) q[0];
sx q[0];
rz(-1.9032318) q[0];
sx q[0];
rz(-2.0623273) q[0];
rz(-0.20118841) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(-2.3692998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.444126) q[1];
sx q[1];
rz(-2.0545469) q[1];
sx q[1];
rz(0.53375785) q[1];
rz(-pi) q[2];
x q[2];
rz(1.059504) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(-2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(-2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-0.84386688) q[0];
rz(1.6183052) q[1];
sx q[1];
rz(-0.73262501) q[1];
sx q[1];
rz(1.9218146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0735002) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(-1.5324701) q[0];
x q[1];
rz(-1.1559406) q[2];
sx q[2];
rz(-1.875669) q[2];
sx q[2];
rz(2.4883303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.8169889) q[1];
sx q[1];
rz(-0.90087955) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7759833) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(0.52683559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(-0.037242446) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1036296) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(-2.7323639) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(2.7917042) q[0];
rz(-pi) q[1];
rz(0.52493237) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(2.8033825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2891149) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(0.27681338) q[1];
rz(-0.51729047) q[3];
sx q[3];
rz(-1.7833157) q[3];
sx q[3];
rz(-1.9912491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1428712) q[2];
sx q[2];
rz(-1.2850392) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(1.508629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976534) q[0];
sx q[0];
rz(-2.2542852) q[0];
sx q[0];
rz(-0.20792122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.091392322) q[2];
sx q[2];
rz(-2.1835727) q[2];
sx q[2];
rz(0.44832715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.75084) q[1];
sx q[1];
rz(-0.62493338) q[1];
sx q[1];
rz(2.1363439) q[1];
rz(-pi) q[2];
rz(-0.056592654) q[3];
sx q[3];
rz(-0.86846272) q[3];
sx q[3];
rz(3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14279723) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(1.8113332) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(-1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.7601097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488659) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(3.0548274) q[0];
rz(-0.21391602) q[2];
sx q[2];
rz(-0.48366085) q[2];
sx q[2];
rz(-2.1733401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46422058) q[1];
sx q[1];
rz(-2.3090625) q[1];
sx q[1];
rz(2.1470451) q[1];
rz(0.28678203) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(-2.5449469) q[2];
rz(2.6560442) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-3.1141282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116466) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(-1.3163153) q[2];
sx q[2];
rz(-1.6152059) q[2];
sx q[2];
rz(-0.16441328) q[2];
rz(-2.9968895) q[3];
sx q[3];
rz(-0.51454138) q[3];
sx q[3];
rz(0.71458057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];