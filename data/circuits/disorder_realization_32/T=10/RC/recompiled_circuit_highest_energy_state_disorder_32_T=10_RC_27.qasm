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
rz(2.0139439) q[1];
sx q[1];
rz(4.5166587) q[1];
sx q[1];
rz(5.4969129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776913) q[0];
sx q[0];
rz(-2.7161975) q[0];
sx q[0];
rz(2.0869654) q[0];
rz(-pi) q[1];
rz(-0.9144056) q[2];
sx q[2];
rz(-2.0255651) q[2];
sx q[2];
rz(0.17780534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8189803) q[1];
sx q[1];
rz(-0.29085813) q[1];
sx q[1];
rz(-1.4664654) q[1];
rz(-2.7336596) q[3];
sx q[3];
rz(-2.130748) q[3];
sx q[3];
rz(-1.6168646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11975153) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(1.2500259) q[2];
rz(0.73578468) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.5031987) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(2.6128795) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3246177) q[0];
sx q[0];
rz(-1.2551089) q[0];
sx q[0];
rz(-1.0782918) q[0];
rz(-pi) q[1];
rz(0.89854449) q[2];
sx q[2];
rz(-2.8437584) q[2];
sx q[2];
rz(0.98984776) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0866249) q[1];
sx q[1];
rz(-1.0180915) q[1];
sx q[1];
rz(-2.5235559) q[1];
rz(1.0924322) q[3];
sx q[3];
rz(-2.6809664) q[3];
sx q[3];
rz(0.13308976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66136393) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(-3.1129692) q[2];
rz(2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(-0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(0.89135528) q[0];
rz(2.640653) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1315434) q[0];
sx q[0];
rz(-1.5230522) q[0];
sx q[0];
rz(3.0555807) q[0];
rz(-pi) q[1];
rz(1.1040014) q[2];
sx q[2];
rz(-0.57764556) q[2];
sx q[2];
rz(1.6922598) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2477639) q[1];
sx q[1];
rz(-2.2732452) q[1];
sx q[1];
rz(-0.6599627) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15283382) q[3];
sx q[3];
rz(-0.18251576) q[3];
sx q[3];
rz(1.0618126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3567051) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9193566) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-0.87523571) q[0];
rz(-1.596176) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-0.045305591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5298669) q[0];
sx q[0];
rz(-1.5054724) q[0];
sx q[0];
rz(-1.4683112) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50778054) q[2];
sx q[2];
rz(-0.96755257) q[2];
sx q[2];
rz(-0.3321307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6958624) q[1];
sx q[1];
rz(-0.67477422) q[1];
sx q[1];
rz(-1.710102) q[1];
rz(-pi) q[2];
rz(-1.8604467) q[3];
sx q[3];
rz(-0.60263205) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8594325) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(0.68242514) q[0];
rz(-1.2884864) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.8103745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73062452) q[0];
sx q[0];
rz(-1.5544116) q[0];
sx q[0];
rz(2.2975305) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8814828) q[2];
sx q[2];
rz(-1.9528198) q[2];
sx q[2];
rz(-1.1191238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.694931) q[1];
sx q[1];
rz(-0.38685054) q[1];
sx q[1];
rz(2.7617161) q[1];
rz(-pi) q[2];
rz(2.6909037) q[3];
sx q[3];
rz(-1.0875877) q[3];
sx q[3];
rz(2.2775028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1149301) q[2];
sx q[2];
rz(-0.59938359) q[2];
sx q[2];
rz(-0.66693532) q[2];
rz(-0.55495787) q[3];
sx q[3];
rz(-2.4542377) q[3];
sx q[3];
rz(0.2989029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126467) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(2.4378648) q[0];
rz(2.9723889) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-0.15883787) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1554398) q[0];
sx q[0];
rz(-1.4503398) q[0];
sx q[0];
rz(-1.7140165) q[0];
rz(2.4689697) q[2];
sx q[2];
rz(-2.9458698) q[2];
sx q[2];
rz(-1.4162404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68137121) q[1];
sx q[1];
rz(-1.4386144) q[1];
sx q[1];
rz(-2.6580407) q[1];
rz(2.8377526) q[3];
sx q[3];
rz(-0.66038495) q[3];
sx q[3];
rz(0.25750289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77330294) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(-2.2840195) q[2];
rz(1.9722021) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17539772) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(0.71863693) q[0];
rz(-2.8596558) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(0.66863543) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40067264) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(-0.4192062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1988772) q[2];
sx q[2];
rz(-1.8657902) q[2];
sx q[2];
rz(1.4843954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.201828) q[1];
sx q[1];
rz(-1.2039021) q[1];
sx q[1];
rz(0.35122996) q[1];
rz(0.63124048) q[3];
sx q[3];
rz(-1.5774014) q[3];
sx q[3];
rz(2.1932878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(1.0364214) q[2];
rz(0.63452619) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.6805639) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(0.72599894) q[0];
rz(1.1622102) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(1.1964218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0045753) q[0];
sx q[0];
rz(-1.5423623) q[0];
sx q[0];
rz(1.397055) q[0];
rz(-3.093946) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(2.7757211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4359522) q[1];
sx q[1];
rz(-2.0355823) q[1];
sx q[1];
rz(2.57354) q[1];
rz(2.7128588) q[3];
sx q[3];
rz(-0.7168684) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2339345) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-2.2080803) q[2];
rz(-0.61698169) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(-2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1652949) q[0];
sx q[0];
rz(-0.10460654) q[0];
sx q[0];
rz(-0.053255178) q[0];
rz(0.28757295) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(-2.8823749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142562) q[0];
sx q[0];
rz(-1.5948609) q[0];
sx q[0];
rz(-1.5177478) q[0];
x q[1];
rz(1.0561814) q[2];
sx q[2];
rz(-2.8764203) q[2];
sx q[2];
rz(2.3725703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6280243) q[1];
sx q[1];
rz(-1.6322989) q[1];
sx q[1];
rz(0.81925591) q[1];
x q[2];
rz(1.7227371) q[3];
sx q[3];
rz(-2.3309618) q[3];
sx q[3];
rz(-1.6068271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(2.7465316) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(0.980353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(0.3821061) q[0];
rz(-2.8451879) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.872725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5871966) q[0];
sx q[0];
rz(-2.4807839) q[0];
sx q[0];
rz(-0.83202485) q[0];
rz(-pi) q[1];
rz(1.7364794) q[2];
sx q[2];
rz(-2.4615917) q[2];
sx q[2];
rz(0.14762893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8812534) q[1];
sx q[1];
rz(-1.9104382) q[1];
sx q[1];
rz(-0.022625523) q[1];
rz(-pi) q[2];
rz(-1.6927229) q[3];
sx q[3];
rz(-1.1200179) q[3];
sx q[3];
rz(1.4565695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2962013) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(-0.38983795) q[2];
rz(-1.2975533) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98904499) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(1.7851495) q[1];
sx q[1];
rz(-1.0379797) q[1];
sx q[1];
rz(-1.5300068) q[1];
rz(2.6806954) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(1.233558) q[3];
sx q[3];
rz(-2.3261286) q[3];
sx q[3];
rz(0.89589768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
