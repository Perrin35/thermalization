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
rz(1.8523822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122983) q[0];
sx q[0];
rz(-1.1120783) q[0];
sx q[0];
rz(0.32231583) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1172328) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(-2.9565405) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8569021) q[1];
sx q[1];
rz(-2.8077217) q[1];
sx q[1];
rz(0.88792172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3207924) q[3];
sx q[3];
rz(-0.27640009) q[3];
sx q[3];
rz(1.9763415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(-3.0060449) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(2.8285817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.27539) q[0];
sx q[0];
rz(-1.7821527) q[0];
sx q[0];
rz(-2.9160121) q[0];
rz(-pi) q[1];
rz(2.0560527) q[2];
sx q[2];
rz(-1.8453571) q[2];
sx q[2];
rz(-0.95825125) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7105512) q[1];
sx q[1];
rz(-1.4841054) q[1];
sx q[1];
rz(0.1608506) q[1];
rz(-pi) q[2];
rz(1.6766657) q[3];
sx q[3];
rz(-1.6459811) q[3];
sx q[3];
rz(-2.4014846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69497067) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-0.95412811) q[2];
rz(1.4387087) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838098) q[0];
sx q[0];
rz(-2.5397781) q[0];
sx q[0];
rz(-0.69586786) q[0];
rz(2.7867735) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-0.16608873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8649193) q[0];
sx q[0];
rz(-1.858798) q[0];
sx q[0];
rz(-2.9898781) q[0];
x q[1];
rz(0.66785779) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(2.6454676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3661763) q[1];
sx q[1];
rz(-1.1300547) q[1];
sx q[1];
rz(2.612545) q[1];
rz(-pi) q[2];
rz(2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(-0.13845201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-0.91527462) q[2];
sx q[2];
rz(2.1598143) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(2.847805) q[0];
rz(0.44149533) q[1];
sx q[1];
rz(-1.0499294) q[1];
sx q[1];
rz(-2.3667483) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64530173) q[0];
sx q[0];
rz(-0.075356396) q[0];
sx q[0];
rz(0.2199748) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11544322) q[2];
sx q[2];
rz(-0.37027678) q[2];
sx q[2];
rz(1.2790797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0171417) q[1];
sx q[1];
rz(-1.1688197) q[1];
sx q[1];
rz(-0.084852858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1713486) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(-1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(2.5715128) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.7385087) q[0];
sx q[0];
rz(0.18138012) q[0];
rz(-0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5126702) q[0];
sx q[0];
rz(-1.218601) q[0];
sx q[0];
rz(-0.098350071) q[0];
rz(2.2147419) q[2];
sx q[2];
rz(-2.0276208) q[2];
sx q[2];
rz(-2.9608375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6377423) q[1];
sx q[1];
rz(-2.3554194) q[1];
sx q[1];
rz(-2.1450858) q[1];
rz(2.0884573) q[3];
sx q[3];
rz(-2.5128551) q[3];
sx q[3];
rz(-0.79013463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-0.28953141) q[2];
rz(3.1048408) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(-0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.0832131) q[0];
sx q[0];
rz(-2.932974) q[0];
sx q[0];
rz(1.6311197) q[0];
rz(1.1389114) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6769584) q[0];
sx q[0];
rz(-2.5559253) q[0];
sx q[0];
rz(2.2023489) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11634155) q[2];
sx q[2];
rz(-2.939072) q[2];
sx q[2];
rz(-0.6845189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6974666) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(0.53375785) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3853112) q[3];
sx q[3];
rz(-0.68538266) q[3];
sx q[3];
rz(1.7862198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54414526) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(1.0774111) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(2.6950148) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(1.6183052) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.2197781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6780753) q[0];
sx q[0];
rz(-3.0419555) q[0];
sx q[0];
rz(0.3936605) q[0];
x q[1];
rz(-2.8104066) q[2];
sx q[2];
rz(-1.965431) q[2];
sx q[2];
rz(-1.0489724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.3246037) q[1];
sx q[1];
rz(-2.2407131) q[1];
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
x q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(-0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-0.75604701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(2.8102002) q[0];
rz(-1.1720852) q[1];
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
rz(-1.0497826) q[0];
sx q[0];
rz(-2.4342393) q[0];
sx q[0];
rz(-2.011767) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6166603) q[2];
sx q[2];
rz(-1.8443174) q[2];
sx q[2];
rz(-2.8033825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48620263) q[1];
sx q[1];
rz(-1.7230002) q[1];
sx q[1];
rz(-0.57032077) q[1];
x q[2];
rz(-2.6243022) q[3];
sx q[3];
rz(-1.7833157) q[3];
sx q[3];
rz(-1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.5107091) q[0];
sx q[0];
rz(3.1316485) q[0];
rz(1.0558646) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.508629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89440896) q[0];
sx q[0];
rz(-1.4100473) q[0];
sx q[0];
rz(0.87662351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95605085) q[2];
sx q[2];
rz(-1.6455257) q[2];
sx q[2];
rz(1.1751307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.75084) q[1];
sx q[1];
rz(-0.62493338) q[1];
sx q[1];
rz(-1.0052488) q[1];
x q[2];
rz(-1.6375332) q[3];
sx q[3];
rz(-0.70422322) q[3];
sx q[3];
rz(3.0610386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(-1.8113332) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-1.0288303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10903877) q[0];
sx q[0];
rz(-2.3333896) q[0];
sx q[0];
rz(-2.57634) q[0];
rz(2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.7601097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0927267) q[0];
sx q[0];
rz(-0.5038358) q[0];
sx q[0];
rz(-3.0548274) q[0];
x q[1];
rz(2.9276766) q[2];
sx q[2];
rz(-2.6579318) q[2];
sx q[2];
rz(2.1733401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5188462) q[1];
sx q[1];
rz(-1.1559876) q[1];
sx q[1];
rz(2.3153789) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8990717) q[3];
sx q[3];
rz(-2.4007113) q[3];
sx q[3];
rz(1.5981984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4505724) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(0.59664574) q[2];
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
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116466) q[0];
sx q[0];
rz(-1.9367138) q[0];
sx q[0];
rz(1.1116897) q[0];
rz(-1.8687517) q[1];
sx q[1];
rz(-1.0796937) q[1];
sx q[1];
rz(-0.83273522) q[1];
rz(1.3163153) q[2];
sx q[2];
rz(-1.5263867) q[2];
sx q[2];
rz(2.9771794) q[2];
rz(-0.14470312) q[3];
sx q[3];
rz(-2.6270513) q[3];
sx q[3];
rz(-2.4270121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
