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
rz(1.2122756) q[0];
sx q[0];
rz(-2.0097998) q[0];
sx q[0];
rz(1.3728859) q[0];
rz(-1.1276487) q[1];
sx q[1];
rz(-1.375066) q[1];
sx q[1];
rz(-2.3553203) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776913) q[0];
sx q[0];
rz(-0.4253952) q[0];
sx q[0];
rz(-1.0546272) q[0];
rz(-0.89531548) q[2];
sx q[2];
rz(-2.3626872) q[2];
sx q[2];
rz(-1.9112183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3481625) q[1];
sx q[1];
rz(-1.5409267) q[1];
sx q[1];
rz(1.2814327) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7336596) q[3];
sx q[3];
rz(-1.0108447) q[3];
sx q[3];
rz(-1.5247281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0218411) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(-1.8915668) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982872) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(3.0702718) q[0];
rz(1.9885063) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(-0.52871314) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77069717) q[0];
sx q[0];
rz(-0.57790725) q[0];
sx q[0];
rz(2.175287) q[0];
rz(2.9527093) q[2];
sx q[2];
rz(-1.802465) q[2];
sx q[2];
rz(2.846039) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15859998) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(0.92293585) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22457476) q[3];
sx q[3];
rz(-1.16515) q[3];
sx q[3];
rz(-0.39164513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(3.1129692) q[2];
rz(-0.35401595) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(-2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4383168) q[0];
sx q[0];
rz(-1.0030614) q[0];
sx q[0];
rz(0.89135528) q[0];
rz(0.50093961) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(3.0827789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0100493) q[0];
sx q[0];
rz(-1.5230522) q[0];
sx q[0];
rz(3.0555807) q[0];
rz(-1.0436613) q[2];
sx q[2];
rz(-1.322515) q[2];
sx q[2];
rz(2.6206526) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2477639) q[1];
sx q[1];
rz(-0.86834748) q[1];
sx q[1];
rz(2.48163) q[1];
x q[2];
rz(-1.5988878) q[3];
sx q[3];
rz(-1.3904316) q[3];
sx q[3];
rz(1.9244058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7848876) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(-0.2365665) q[2];
rz(-0.92075721) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22223602) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(0.87523571) q[0];
rz(1.5454166) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(0.045305591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047642853) q[0];
sx q[0];
rz(-1.4685306) q[0];
sx q[0];
rz(-0.065667466) q[0];
rz(-pi) q[1];
rz(-2.1854109) q[2];
sx q[2];
rz(-0.76757694) q[2];
sx q[2];
rz(2.0337348) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9074752) q[1];
sx q[1];
rz(-1.4839412) q[1];
sx q[1];
rz(0.90076561) q[1];
rz(-pi) q[2];
rz(-1.8604467) q[3];
sx q[3];
rz(-0.60263205) q[3];
sx q[3];
rz(-1.8528191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28216013) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(-1.947594) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0328261) q[0];
sx q[0];
rz(-2.325401) q[0];
sx q[0];
rz(-2.4591675) q[0];
rz(1.8531063) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(-1.8103745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2829943) q[0];
sx q[0];
rz(-0.72688519) q[0];
sx q[0];
rz(-1.5954533) q[0];
rz(-2.8814828) q[2];
sx q[2];
rz(-1.9528198) q[2];
sx q[2];
rz(1.1191238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9115548) q[1];
sx q[1];
rz(-1.7111519) q[1];
sx q[1];
rz(2.7799003) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0430452) q[3];
sx q[3];
rz(-1.1748136) q[3];
sx q[3];
rz(0.92786232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(-2.4746573) q[2];
rz(-0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126467) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(-2.4378648) q[0];
rz(-2.9723889) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-2.9827548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.279351) q[0];
sx q[0];
rz(-0.18687525) q[0];
sx q[0];
rz(0.86743768) q[0];
rz(-2.4689697) q[2];
sx q[2];
rz(-0.19572283) q[2];
sx q[2];
rz(-1.4162404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.006268) q[1];
sx q[1];
rz(-2.6416831) q[1];
sx q[1];
rz(0.27854021) q[1];
rz(2.5037836) q[3];
sx q[3];
rz(-1.7553698) q[3];
sx q[3];
rz(-2.0710433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77330294) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(0.85757315) q[2];
rz(-1.9722021) q[3];
sx q[3];
rz(-1.2811067) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9661949) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(-2.4229557) q[0];
rz(-2.8596558) q[1];
sx q[1];
rz(-1.3749296) q[1];
sx q[1];
rz(-0.66863543) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91808337) q[0];
sx q[0];
rz(-1.232172) q[0];
sx q[0];
rz(2.2302367) q[0];
rz(-pi) q[1];
rz(-2.1988772) q[2];
sx q[2];
rz(-1.2758024) q[2];
sx q[2];
rz(1.4843954) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9397647) q[1];
sx q[1];
rz(-1.2039021) q[1];
sx q[1];
rz(-2.7903627) q[1];
rz(1.5626146) q[3];
sx q[3];
rz(-2.2020209) q[3];
sx q[3];
rz(0.62731987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4947027) q[2];
sx q[2];
rz(-1.1855482) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(0.63452619) q[3];
sx q[3];
rz(-0.98792881) q[3];
sx q[3];
rz(-0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46102872) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(-2.4155937) q[0];
rz(-1.1622102) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(1.1964218) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57121074) q[0];
sx q[0];
rz(-1.397126) q[0];
sx q[0];
rz(0.028868361) q[0];
rz(-pi) q[1];
x q[1];
rz(3.093946) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(-2.7757211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5550674) q[1];
sx q[1];
rz(-1.0690926) q[1];
sx q[1];
rz(1.0341767) q[1];
rz(-1.9184018) q[3];
sx q[3];
rz(-0.93033997) q[3];
sx q[3];
rz(0.15061041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90765816) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-2.2080803) q[2];
rz(-2.524611) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(-0.84701076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-0.053255178) q[0];
rz(-0.28757295) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(-0.25921777) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5968551) q[0];
sx q[0];
rz(-1.6238295) q[0];
sx q[0];
rz(-3.1174942) q[0];
x q[1];
rz(2.0854112) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(-0.76902232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1146015) q[1];
sx q[1];
rz(-0.82102126) q[1];
sx q[1];
rz(-3.0575086) q[1];
rz(-pi) q[2];
rz(0.15786981) q[3];
sx q[3];
rz(-2.3693858) q[3];
sx q[3];
rz(-1.7534353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(0.1813691) q[2];
rz(0.3950611) q[3];
sx q[3];
rz(-0.22160465) q[3];
sx q[3];
rz(0.980353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(2.7594866) q[0];
rz(0.29640472) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.872725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110376) q[0];
sx q[0];
rz(-1.0998304) q[0];
sx q[0];
rz(2.6593326) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13258719) q[2];
sx q[2];
rz(-0.90182038) q[2];
sx q[2];
rz(-0.35945177) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8812534) q[1];
sx q[1];
rz(-1.9104382) q[1];
sx q[1];
rz(3.1189671) q[1];
rz(-1.6927229) q[3];
sx q[3];
rz(-1.1200179) q[3];
sx q[3];
rz(-1.6850231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84539139) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(0.38983795) q[2];
rz(-1.2975533) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525477) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(2.6806954) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(-0.78430243) q[3];
sx q[3];
rz(-1.3275066) q[3];
sx q[3];
rz(2.2307997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
