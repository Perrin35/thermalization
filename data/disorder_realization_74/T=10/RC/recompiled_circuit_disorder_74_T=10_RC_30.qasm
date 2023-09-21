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
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122983) q[0];
sx q[0];
rz(-2.0295143) q[0];
sx q[0];
rz(2.8192768) q[0];
rz(-pi) q[1];
rz(-2.870954) q[2];
sx q[2];
rz(-3.116313) q[2];
sx q[2];
rz(1.4852922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2002441) q[1];
sx q[1];
rz(-1.779088) q[1];
sx q[1];
rz(1.3079446) q[1];
rz(-3.0715277) q[3];
sx q[3];
rz(-1.8383887) q[3];
sx q[3];
rz(0.90581264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5499128) q[2];
sx q[2];
rz(-0.16273558) q[2];
sx q[2];
rz(-0.13554779) q[2];
rz(-3.0013951) q[3];
sx q[3];
rz(-2.086816) q[3];
sx q[3];
rz(2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(0.78080368) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-2.5550911) q[1];
sx q[1];
rz(-1.3134726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.388893) q[0];
sx q[0];
rz(-1.7912731) q[0];
sx q[0];
rz(1.3541143) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0275231) q[2];
sx q[2];
rz(-0.55210219) q[2];
sx q[2];
rz(-0.13763025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(-2.980742) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95150681) q[3];
sx q[3];
rz(-0.12976876) q[3];
sx q[3];
rz(1.6956003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.6080674) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3375219) q[0];
sx q[0];
rz(-1.7162168) q[0];
sx q[0];
rz(-1.8619596) q[0];
rz(1.0679929) q[2];
sx q[2];
rz(-1.0223801) q[2];
sx q[2];
rz(1.8333679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6925466) q[1];
sx q[1];
rz(-1.0967626) q[1];
sx q[1];
rz(-1.0707335) q[1];
rz(-pi) q[2];
rz(1.2766248) q[3];
sx q[3];
rz(-0.73495451) q[3];
sx q[3];
rz(2.8923349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(0.98177838) q[2];
rz(-2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(2.847805) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-0.77484432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(0.073547151) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7735633) q[2];
sx q[2];
rz(-1.6124915) q[2];
sx q[2];
rz(-0.3993984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6619819) q[1];
sx q[1];
rz(-1.4927215) q[1];
sx q[1];
rz(1.1675203) q[1];
rz(-0.1713486) q[3];
sx q[3];
rz(-1.315794) q[3];
sx q[3];
rz(-1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-1.501804) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(3.0338874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2337423) q[0];
sx q[0];
rz(-1.478501) q[0];
sx q[0];
rz(1.924563) q[0];
x q[1];
rz(-2.2147419) q[2];
sx q[2];
rz(-1.1139718) q[2];
sx q[2];
rz(0.18075519) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2376033) q[1];
sx q[1];
rz(-0.9346107) q[1];
sx q[1];
rz(-0.49828766) q[1];
rz(2.0884573) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(-2.351458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(-3.1048408) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058379563) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(-1.1389114) q[1];
sx q[1];
rz(-1.6612256) q[1];
sx q[1];
rz(-2.1246134) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3973629) q[0];
sx q[0];
rz(-2.0332391) q[0];
sx q[0];
rz(0.37325333) q[0];
rz(1.594627) q[2];
sx q[2];
rz(-1.7719291) q[2];
sx q[2];
rz(-0.80326524) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.348154) q[1];
sx q[1];
rz(-2.4373756) q[1];
sx q[1];
rz(-2.3401295) q[1];
x q[2];
rz(2.6050657) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(-2.7262053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54414526) q[2];
sx q[2];
rz(-1.0015254) q[2];
sx q[2];
rz(2.5069359) q[2];
rz(-1.0774111) q[3];
sx q[3];
rz(-1.0297188) q[3];
sx q[3];
rz(2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0680925) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(-1.6091225) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9856521) q[2];
sx q[2];
rz(-1.2659237) q[2];
sx q[2];
rz(2.4883303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0724832) q[1];
sx q[1];
rz(-1.8169889) q[1];
sx q[1];
rz(-2.2407131) q[1];
rz(-1.3656093) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(-2.6147571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62961489) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(2.3855456) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-0.66292271) q[1];
sx q[1];
rz(2.7323639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5361355) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(2.7917042) q[0];
rz(0.52493237) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(2.8033825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85247773) q[1];
sx q[1];
rz(-0.58809885) q[1];
sx q[1];
rz(-0.27681338) q[1];
rz(-2.6243022) q[3];
sx q[3];
rz(-1.358277) q[3];
sx q[3];
rz(-1.9912491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(0.085111246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(-3.1316485) q[0];
rz(2.0857281) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.6329637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89440896) q[0];
sx q[0];
rz(-1.4100473) q[0];
sx q[0];
rz(-0.87662351) q[0];
x q[1];
rz(-3.0502003) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(0.44832715) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0548045) q[1];
sx q[1];
rz(-2.0874223) q[1];
sx q[1];
rz(-0.36887849) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.085) q[3];
sx q[3];
rz(-0.86846272) q[3];
sx q[3];
rz(-3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9987954) q[2];
sx q[2];
rz(-1.5974416) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-2.1222474) q[1];
sx q[1];
rz(-1.3814829) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0488659) q[0];
sx q[0];
rz(-2.6377569) q[0];
sx q[0];
rz(-0.086765246) q[0];
rz(1.4597458) q[2];
sx q[2];
rz(-1.0990708) q[2];
sx q[2];
rz(-1.2088838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8379412) q[1];
sx q[1];
rz(-2.2397869) q[1];
sx q[1];
rz(-0.53955697) q[1];
x q[2];
rz(0.28678203) q[3];
sx q[3];
rz(-2.2638595) q[3];
sx q[3];
rz(1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.4505724) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(-0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52994603) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(1.8687517) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(1.3163153) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(0.5100526) q[3];
sx q[3];
rz(-1.4997713) q[3];
sx q[3];
rz(-0.98239519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];