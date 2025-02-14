OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5069198) q[0];
sx q[0];
rz(-1.4022175) q[0];
sx q[0];
rz(0.75048598) q[0];
rz(0.063272417) q[1];
sx q[1];
rz(7.1019389) q[1];
sx q[1];
rz(11.547312) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3844168) q[0];
sx q[0];
rz(-1.0121167) q[0];
sx q[0];
rz(-1.7930828) q[0];
x q[1];
rz(-1.8797714) q[2];
sx q[2];
rz(-0.7337386) q[2];
sx q[2];
rz(0.66622615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77324919) q[1];
sx q[1];
rz(-1.4344402) q[1];
sx q[1];
rz(-2.5565867) q[1];
rz(2.597911) q[3];
sx q[3];
rz(-1.8988109) q[3];
sx q[3];
rz(1.8193434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2941306) q[2];
sx q[2];
rz(-0.99404603) q[2];
sx q[2];
rz(-1.6729986) q[2];
rz(2.9351506) q[3];
sx q[3];
rz(-1.3279746) q[3];
sx q[3];
rz(2.4895721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9701397) q[0];
sx q[0];
rz(-1.7296706) q[0];
sx q[0];
rz(-0.10301244) q[0];
rz(0.39762321) q[1];
sx q[1];
rz(-2.7151974) q[1];
sx q[1];
rz(0.85223371) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9270902) q[0];
sx q[0];
rz(-0.53686857) q[0];
sx q[0];
rz(-1.4528627) q[0];
rz(-pi) q[1];
rz(-0.28362922) q[2];
sx q[2];
rz(-1.8967702) q[2];
sx q[2];
rz(-0.21421283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.045984775) q[1];
sx q[1];
rz(-1.5846445) q[1];
sx q[1];
rz(2.0944089) q[1];
rz(-pi) q[2];
rz(1.6144595) q[3];
sx q[3];
rz(-1.763673) q[3];
sx q[3];
rz(-1.5571644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6323866) q[2];
sx q[2];
rz(-2.1337324) q[2];
sx q[2];
rz(0.81727916) q[2];
rz(2.9362074) q[3];
sx q[3];
rz(-2.9257264) q[3];
sx q[3];
rz(-0.5425905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83194724) q[0];
sx q[0];
rz(-2.0537856) q[0];
sx q[0];
rz(2.621424) q[0];
rz(-2.4179516) q[1];
sx q[1];
rz(-1.8579204) q[1];
sx q[1];
rz(2.9626194) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13143927) q[0];
sx q[0];
rz(-0.69973937) q[0];
sx q[0];
rz(-2.418451) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51836022) q[2];
sx q[2];
rz(-0.39769618) q[2];
sx q[2];
rz(-0.30425554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5120257) q[1];
sx q[1];
rz(-0.24228748) q[1];
sx q[1];
rz(-2.971388) q[1];
x q[2];
rz(-0.34647219) q[3];
sx q[3];
rz(-2.4541509) q[3];
sx q[3];
rz(0.32780743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2554864) q[2];
sx q[2];
rz(-2.8612558) q[2];
sx q[2];
rz(0.96957668) q[2];
rz(0.3420091) q[3];
sx q[3];
rz(-1.0229144) q[3];
sx q[3];
rz(-1.3268933) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2372357) q[0];
sx q[0];
rz(-2.3396753) q[0];
sx q[0];
rz(-3.0661769) q[0];
rz(0.30741179) q[1];
sx q[1];
rz(-1.9785898) q[1];
sx q[1];
rz(-2.7154162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3097174) q[0];
sx q[0];
rz(-0.85098905) q[0];
sx q[0];
rz(2.4144774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6226042) q[2];
sx q[2];
rz(-0.52757571) q[2];
sx q[2];
rz(2.0466387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28171912) q[1];
sx q[1];
rz(-1.3744756) q[1];
sx q[1];
rz(2.1540889) q[1];
rz(-1.5463065) q[3];
sx q[3];
rz(-1.1653641) q[3];
sx q[3];
rz(1.0129077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9093466) q[2];
sx q[2];
rz(-1.1541977) q[2];
sx q[2];
rz(2.2171891) q[2];
rz(1.1030446) q[3];
sx q[3];
rz(-1.1140991) q[3];
sx q[3];
rz(1.8614785) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0019317146) q[0];
sx q[0];
rz(-2.6502471) q[0];
sx q[0];
rz(0.77144462) q[0];
rz(-1.8292684) q[1];
sx q[1];
rz(-1.7701365) q[1];
sx q[1];
rz(-1.2244474) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1405914) q[0];
sx q[0];
rz(-3.0016905) q[0];
sx q[0];
rz(1.7732267) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61842711) q[2];
sx q[2];
rz(-1.569241) q[2];
sx q[2];
rz(0.33829442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33358303) q[1];
sx q[1];
rz(-1.4794655) q[1];
sx q[1];
rz(-1.9568981) q[1];
rz(-pi) q[2];
rz(-1.1125074) q[3];
sx q[3];
rz(-0.84064129) q[3];
sx q[3];
rz(-2.2793913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15765849) q[2];
sx q[2];
rz(-2.5658786) q[2];
sx q[2];
rz(2.0470587) q[2];
rz(-0.076722773) q[3];
sx q[3];
rz(-2.6321497) q[3];
sx q[3];
rz(0.69242394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7102966) q[0];
sx q[0];
rz(-2.4831979) q[0];
sx q[0];
rz(2.9265213) q[0];
rz(2.1521425) q[1];
sx q[1];
rz(-0.96885252) q[1];
sx q[1];
rz(-1.6542124) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7152255) q[0];
sx q[0];
rz(-2.0837113) q[0];
sx q[0];
rz(-1.3839534) q[0];
x q[1];
rz(-1.8645668) q[2];
sx q[2];
rz(-1.1907628) q[2];
sx q[2];
rz(0.88949163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.9729197) q[1];
sx q[1];
rz(-0.56747371) q[1];
sx q[1];
rz(0.85166906) q[1];
rz(2.0607674) q[3];
sx q[3];
rz(-2.1741011) q[3];
sx q[3];
rz(-1.5321466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6626849) q[2];
sx q[2];
rz(-1.7272793) q[2];
sx q[2];
rz(-1.5642081) q[2];
rz(-2.1936737) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(-1.4028367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3941512) q[0];
sx q[0];
rz(-1.489137) q[0];
sx q[0];
rz(-0.58962756) q[0];
rz(1.6472752) q[1];
sx q[1];
rz(-1.7604156) q[1];
sx q[1];
rz(-0.078701198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32959461) q[0];
sx q[0];
rz(-2.643739) q[0];
sx q[0];
rz(-0.6377687) q[0];
rz(-pi) q[1];
rz(2.6769889) q[2];
sx q[2];
rz(-1.7853933) q[2];
sx q[2];
rz(0.047296798) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4884686) q[1];
sx q[1];
rz(-0.32652509) q[1];
sx q[1];
rz(-2.4420183) q[1];
x q[2];
rz(2.614903) q[3];
sx q[3];
rz(-1.1136276) q[3];
sx q[3];
rz(1.4032193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2176723) q[2];
sx q[2];
rz(-0.87409002) q[2];
sx q[2];
rz(-0.11202845) q[2];
rz(0.40515408) q[3];
sx q[3];
rz(-1.743789) q[3];
sx q[3];
rz(-2.9356094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21796313) q[0];
sx q[0];
rz(-2.1681652) q[0];
sx q[0];
rz(-1.4189036) q[0];
rz(-2.0619242) q[1];
sx q[1];
rz(-2.7579312) q[1];
sx q[1];
rz(-1.7344281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0100219) q[0];
sx q[0];
rz(-2.6162806) q[0];
sx q[0];
rz(-0.66862392) q[0];
rz(-pi) q[1];
rz(0.43145572) q[2];
sx q[2];
rz(-1.4979216) q[2];
sx q[2];
rz(-0.88957722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8720334) q[1];
sx q[1];
rz(-2.5692794) q[1];
sx q[1];
rz(0.78929269) q[1];
rz(-1.8457165) q[3];
sx q[3];
rz(-1.8105728) q[3];
sx q[3];
rz(-1.6473242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5088523) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(0.73842326) q[2];
rz(-2.9449055) q[3];
sx q[3];
rz(-2.5236712) q[3];
sx q[3];
rz(0.64016199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7554373) q[0];
sx q[0];
rz(-0.25179395) q[0];
sx q[0];
rz(-0.039903076) q[0];
rz(0.2208605) q[1];
sx q[1];
rz(-1.5233636) q[1];
sx q[1];
rz(1.9853282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98469767) q[0];
sx q[0];
rz(-1.3067208) q[0];
sx q[0];
rz(2.9913802) q[0];
x q[1];
rz(-2.4100163) q[2];
sx q[2];
rz(-2.8551364) q[2];
sx q[2];
rz(-2.4713304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6917518) q[1];
sx q[1];
rz(-0.66206703) q[1];
sx q[1];
rz(-1.0074703) q[1];
rz(-pi) q[2];
rz(-3.0262104) q[3];
sx q[3];
rz(-1.4658864) q[3];
sx q[3];
rz(-2.2656588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3952737) q[2];
sx q[2];
rz(-0.36397448) q[2];
sx q[2];
rz(-3.086997) q[2];
rz(-0.78628457) q[3];
sx q[3];
rz(-1.9449284) q[3];
sx q[3];
rz(-1.981885) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38128024) q[0];
sx q[0];
rz(-2.3535643) q[0];
sx q[0];
rz(0.81083167) q[0];
rz(2.6122818) q[1];
sx q[1];
rz(-1.4365173) q[1];
sx q[1];
rz(-1.1108105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4585574) q[0];
sx q[0];
rz(-2.1330569) q[0];
sx q[0];
rz(-2.0319776) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1322242) q[2];
sx q[2];
rz(-1.9931442) q[2];
sx q[2];
rz(-0.37882159) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6880092) q[1];
sx q[1];
rz(-0.78026672) q[1];
sx q[1];
rz(0.23201402) q[1];
rz(-pi) q[2];
rz(1.9930028) q[3];
sx q[3];
rz(-0.65311382) q[3];
sx q[3];
rz(-0.70142581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79798737) q[2];
sx q[2];
rz(-0.56861773) q[2];
sx q[2];
rz(-0.24151754) q[2];
rz(-2.8999515) q[3];
sx q[3];
rz(-1.4915978) q[3];
sx q[3];
rz(-0.1514761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3342313) q[0];
sx q[0];
rz(-1.7932899) q[0];
sx q[0];
rz(0.67247969) q[0];
rz(2.6895651) q[1];
sx q[1];
rz(-0.94257911) q[1];
sx q[1];
rz(-0.51530757) q[1];
rz(1.6967026) q[2];
sx q[2];
rz(-0.44582146) q[2];
sx q[2];
rz(1.3498342) q[2];
rz(2.0911471) q[3];
sx q[3];
rz(-1.1523714) q[3];
sx q[3];
rz(1.8377792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
