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
rz(1.9361629) q[0];
sx q[0];
rz(-0.054447629) q[0];
sx q[0];
rz(0.86583889) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(-0.058252637) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6893495) q[0];
sx q[0];
rz(-1.5238683) q[0];
sx q[0];
rz(-1.4062792) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0664332) q[2];
sx q[2];
rz(-1.6358536) q[2];
sx q[2];
rz(3.0954454) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4076947) q[1];
sx q[1];
rz(-1.0419894) q[1];
sx q[1];
rz(2.4771792) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9179814) q[3];
sx q[3];
rz(-1.2652186) q[3];
sx q[3];
rz(-1.7479727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3129348) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(-0.9642967) q[2];
rz(-3.0621081) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013414772) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(-1.5403904) q[0];
rz(0.84114289) q[1];
sx q[1];
rz(-1.7336188) q[1];
sx q[1];
rz(-0.85743633) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6591561) q[0];
sx q[0];
rz(-2.0093587) q[0];
sx q[0];
rz(-0.0060122251) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93480627) q[2];
sx q[2];
rz(-1.9611214) q[2];
sx q[2];
rz(1.1771415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7341566) q[1];
sx q[1];
rz(-0.70824558) q[1];
sx q[1];
rz(-2.4328961) q[1];
rz(-2.8174761) q[3];
sx q[3];
rz(-0.29359152) q[3];
sx q[3];
rz(-2.8568134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2555344) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(1.9319755) q[2];
rz(-1.3813193) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7968314) q[0];
sx q[0];
rz(-1.826257) q[0];
sx q[0];
rz(1.8094081) q[0];
rz(1.4765129) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(0.61980334) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42821845) q[0];
sx q[0];
rz(-3.0626903) q[0];
sx q[0];
rz(2.9419241) q[0];
rz(-1.8263426) q[2];
sx q[2];
rz(-1.4076621) q[2];
sx q[2];
rz(2.852885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9732144) q[1];
sx q[1];
rz(-1.4992979) q[1];
sx q[1];
rz(2.8236103) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26472802) q[3];
sx q[3];
rz(-0.94784289) q[3];
sx q[3];
rz(-2.5271551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0009813112) q[2];
sx q[2];
rz(-2.8162214) q[2];
sx q[2];
rz(-0.78835431) q[2];
rz(-1.8654478) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(0.86281323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514423) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(1.066712) q[0];
rz(-1.7881296) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(1.1563168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3301834) q[0];
sx q[0];
rz(-1.5895939) q[0];
sx q[0];
rz(1.7123332) q[0];
x q[1];
rz(1.5300445) q[2];
sx q[2];
rz(-1.5446071) q[2];
sx q[2];
rz(1.9350236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9741648) q[1];
sx q[1];
rz(-1.2819703) q[1];
sx q[1];
rz(1.9603189) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.636999) q[3];
sx q[3];
rz(-0.90747031) q[3];
sx q[3];
rz(-2.5732793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(3.0827674) q[2];
rz(-2.5022653) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-0.85734573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095379742) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(1.5785716) q[1];
sx q[1];
rz(-2.450921) q[1];
sx q[1];
rz(2.6522327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90411579) q[0];
sx q[0];
rz(-1.4129708) q[0];
sx q[0];
rz(-2.514181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7362664) q[2];
sx q[2];
rz(-2.2203682) q[2];
sx q[2];
rz(-2.6774466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64397821) q[1];
sx q[1];
rz(-2.3757739) q[1];
sx q[1];
rz(0.52180565) q[1];
rz(-0.33133026) q[3];
sx q[3];
rz(-3.1090571) q[3];
sx q[3];
rz(-2.0842541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2995149) q[2];
sx q[2];
rz(-1.0966417) q[2];
sx q[2];
rz(-3.1411324) q[2];
rz(0.67356235) q[3];
sx q[3];
rz(-0.898415) q[3];
sx q[3];
rz(-2.4323997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8507268) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(-1.3036183) q[0];
rz(-0.29779008) q[1];
sx q[1];
rz(-0.89556634) q[1];
sx q[1];
rz(0.29388014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760982) q[0];
sx q[0];
rz(-0.84653234) q[0];
sx q[0];
rz(0.80418555) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7104425) q[2];
sx q[2];
rz(-2.2047055) q[2];
sx q[2];
rz(2.0946787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3658947) q[1];
sx q[1];
rz(-2.0112717) q[1];
sx q[1];
rz(-1.0431402) q[1];
x q[2];
rz(0.63346432) q[3];
sx q[3];
rz(-2.7567299) q[3];
sx q[3];
rz(-0.72609392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6728354) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(2.920816) q[2];
rz(1.8686434) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(-1.1481736) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918168) q[0];
sx q[0];
rz(-2.5204372) q[0];
sx q[0];
rz(0.57449269) q[0];
rz(-0.54221398) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(2.7508459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6764561) q[0];
sx q[0];
rz(-2.1628597) q[0];
sx q[0];
rz(1.3964064) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4847267) q[2];
sx q[2];
rz(-0.54605267) q[2];
sx q[2];
rz(-2.4459185) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73394164) q[1];
sx q[1];
rz(-1.9875257) q[1];
sx q[1];
rz(1.5768361) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0454759) q[3];
sx q[3];
rz(-0.36727723) q[3];
sx q[3];
rz(-2.8678081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4603525) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(2.528842) q[2];
rz(-2.1593306) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(-2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14923444) q[0];
sx q[0];
rz(-1.5748011) q[0];
sx q[0];
rz(-3.0862578) q[0];
rz(-1.5845567) q[1];
sx q[1];
rz(-2.0535856) q[1];
sx q[1];
rz(-2.7016644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.330141) q[0];
sx q[0];
rz(-0.52882551) q[0];
sx q[0];
rz(-0.50409601) q[0];
x q[1];
rz(-1.7753168) q[2];
sx q[2];
rz(-2.0530982) q[2];
sx q[2];
rz(2.4647146) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7288558) q[1];
sx q[1];
rz(-1.2555033) q[1];
sx q[1];
rz(2.5607177) q[1];
rz(0.28195076) q[3];
sx q[3];
rz(-2.1190756) q[3];
sx q[3];
rz(-1.5632526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3736734) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(0.31201735) q[2];
rz(0.9196552) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(0.21617226) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31174082) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(3.1134636) q[0];
rz(-1.6955388) q[1];
sx q[1];
rz(-1.9606934) q[1];
sx q[1];
rz(-2.6820954) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.979769) q[0];
sx q[0];
rz(-1.5422675) q[0];
sx q[0];
rz(1.5351487) q[0];
rz(0.33513432) q[2];
sx q[2];
rz(-0.85079934) q[2];
sx q[2];
rz(2.013226) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0952009) q[1];
sx q[1];
rz(-1.4262137) q[1];
sx q[1];
rz(-0.88674366) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9528804) q[3];
sx q[3];
rz(-1.700125) q[3];
sx q[3];
rz(0.36256177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10099899) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(2.8864268) q[2];
rz(-0.98662871) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3808909) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(1.799452) q[0];
rz(-3.1396719) q[1];
sx q[1];
rz(-1.0425967) q[1];
sx q[1];
rz(1.049918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0107897) q[0];
sx q[0];
rz(-1.3585886) q[0];
sx q[0];
rz(2.9899068) q[0];
x q[1];
rz(0.0012108525) q[2];
sx q[2];
rz(-1.7467611) q[2];
sx q[2];
rz(2.3893285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0319977) q[1];
sx q[1];
rz(-1.2032615) q[1];
sx q[1];
rz(-1.3369186) q[1];
x q[2];
rz(2.9258419) q[3];
sx q[3];
rz(-0.85236824) q[3];
sx q[3];
rz(-1.488648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2125825) q[2];
sx q[2];
rz(-1.9276103) q[2];
sx q[2];
rz(0.49599656) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(1.1869441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37954189) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(0.63915359) q[1];
sx q[1];
rz(-1.879138) q[1];
sx q[1];
rz(-2.383147) q[1];
rz(-1.8201309) q[2];
sx q[2];
rz(-2.032866) q[2];
sx q[2];
rz(0.55351071) q[2];
rz(2.2517754) q[3];
sx q[3];
rz(-2.2636236) q[3];
sx q[3];
rz(-1.3286535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
