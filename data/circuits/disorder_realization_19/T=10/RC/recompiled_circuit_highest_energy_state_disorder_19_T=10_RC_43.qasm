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
rz(0.28111464) q[0];
sx q[0];
rz(-1.3507564) q[0];
sx q[0];
rz(-0.95066345) q[0];
rz(-2.9188393) q[1];
sx q[1];
rz(-0.25382257) q[1];
sx q[1];
rz(-2.4737127) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0216103) q[0];
sx q[0];
rz(-1.2111145) q[0];
sx q[0];
rz(1.5403454) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.736343) q[2];
sx q[2];
rz(-1.5551463) q[2];
sx q[2];
rz(-2.1652255) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3006134) q[1];
sx q[1];
rz(-1.6924227) q[1];
sx q[1];
rz(3.0517314) q[1];
rz(-pi) q[2];
rz(1.5667772) q[3];
sx q[3];
rz(-1.360047) q[3];
sx q[3];
rz(0.8523418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5369947) q[2];
sx q[2];
rz(-0.94702417) q[2];
sx q[2];
rz(-0.14744082) q[2];
rz(-1.7668746) q[3];
sx q[3];
rz(-1.3876029) q[3];
sx q[3];
rz(1.714777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073567063) q[0];
sx q[0];
rz(-2.1275213) q[0];
sx q[0];
rz(2.9378939) q[0];
rz(-0.112946) q[1];
sx q[1];
rz(-2.4599383) q[1];
sx q[1];
rz(-0.019498652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9984922) q[0];
sx q[0];
rz(-1.2930627) q[0];
sx q[0];
rz(2.2953675) q[0];
rz(-1.2157529) q[2];
sx q[2];
rz(-0.40242919) q[2];
sx q[2];
rz(-2.7192164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7884036) q[1];
sx q[1];
rz(-2.5063597) q[1];
sx q[1];
rz(0.53643463) q[1];
rz(0.56657378) q[3];
sx q[3];
rz(-0.972803) q[3];
sx q[3];
rz(1.2815042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6273177) q[2];
sx q[2];
rz(-1.5760199) q[2];
sx q[2];
rz(-2.9493098) q[2];
rz(0.2187885) q[3];
sx q[3];
rz(-1.8407121) q[3];
sx q[3];
rz(-1.3505664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4718276) q[0];
sx q[0];
rz(-0.28194031) q[0];
sx q[0];
rz(0.48926735) q[0];
rz(1.6911814) q[1];
sx q[1];
rz(-1.9905118) q[1];
sx q[1];
rz(-0.53389186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035504) q[0];
sx q[0];
rz(-2.8189604) q[0];
sx q[0];
rz(-3.0586092) q[0];
x q[1];
rz(-3.0025173) q[2];
sx q[2];
rz(-1.5267045) q[2];
sx q[2];
rz(1.2033552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.008153) q[1];
sx q[1];
rz(-2.714909) q[1];
sx q[1];
rz(1.2231136) q[1];
rz(-pi) q[2];
rz(1.7524285) q[3];
sx q[3];
rz(-2.2881621) q[3];
sx q[3];
rz(2.7434012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.78407946) q[2];
sx q[2];
rz(-2.1046905) q[2];
sx q[2];
rz(0.64884031) q[2];
rz(0.31323788) q[3];
sx q[3];
rz(-2.1514838) q[3];
sx q[3];
rz(-1.8296957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.677815) q[0];
sx q[0];
rz(-2.2721993) q[0];
sx q[0];
rz(-2.38548) q[0];
rz(0.4568049) q[1];
sx q[1];
rz(-1.1242547) q[1];
sx q[1];
rz(0.66064984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27460262) q[0];
sx q[0];
rz(-1.1701692) q[0];
sx q[0];
rz(-1.7652049) q[0];
rz(-0.22539745) q[2];
sx q[2];
rz(-1.9323914) q[2];
sx q[2];
rz(2.0394675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3676388) q[1];
sx q[1];
rz(-2.3970897) q[1];
sx q[1];
rz(1.8402694) q[1];
rz(-pi) q[2];
rz(-0.9432299) q[3];
sx q[3];
rz(-0.5411202) q[3];
sx q[3];
rz(-2.5515722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75055355) q[2];
sx q[2];
rz(-1.6497352) q[2];
sx q[2];
rz(0.98461119) q[2];
rz(1.6330968) q[3];
sx q[3];
rz(-2.0506141) q[3];
sx q[3];
rz(2.7896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9462747) q[0];
sx q[0];
rz(-2.1066505) q[0];
sx q[0];
rz(-2.4786733) q[0];
rz(1.4074869) q[1];
sx q[1];
rz(-0.88277849) q[1];
sx q[1];
rz(-0.9443121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51296455) q[0];
sx q[0];
rz(-1.4127534) q[0];
sx q[0];
rz(0.66642828) q[0];
rz(2.216249) q[2];
sx q[2];
rz(-2.1133377) q[2];
sx q[2];
rz(-1.134657) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8710526) q[1];
sx q[1];
rz(-0.72977282) q[1];
sx q[1];
rz(0.60128984) q[1];
x q[2];
rz(-1.8098894) q[3];
sx q[3];
rz(-1.0464365) q[3];
sx q[3];
rz(0.6839377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.044109) q[2];
sx q[2];
rz(-0.58610836) q[2];
sx q[2];
rz(-0.07494542) q[2];
rz(-2.1465541) q[3];
sx q[3];
rz(-1.1040265) q[3];
sx q[3];
rz(1.1546571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5465882) q[0];
sx q[0];
rz(-2.0363448) q[0];
sx q[0];
rz(0.30584359) q[0];
rz(2.2587237) q[1];
sx q[1];
rz(-0.63778937) q[1];
sx q[1];
rz(2.2960704) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99793154) q[0];
sx q[0];
rz(-0.73135644) q[0];
sx q[0];
rz(1.9173918) q[0];
rz(1.8550625) q[2];
sx q[2];
rz(-1.5719116) q[2];
sx q[2];
rz(-1.0363611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9041992) q[1];
sx q[1];
rz(-0.9353726) q[1];
sx q[1];
rz(2.2255484) q[1];
x q[2];
rz(-1.7843397) q[3];
sx q[3];
rz(-2.41666) q[3];
sx q[3];
rz(1.8740579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66158295) q[2];
sx q[2];
rz(-0.84605399) q[2];
sx q[2];
rz(2.0704849) q[2];
rz(-0.24460159) q[3];
sx q[3];
rz(-0.94541234) q[3];
sx q[3];
rz(-1.6036114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5474608) q[0];
sx q[0];
rz(-0.56532156) q[0];
sx q[0];
rz(-0.93455899) q[0];
rz(0.67340988) q[1];
sx q[1];
rz(-0.62973657) q[1];
sx q[1];
rz(-3.0455132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6357416) q[0];
sx q[0];
rz(-1.4560501) q[0];
sx q[0];
rz(3.1155354) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3548278) q[2];
sx q[2];
rz(-0.25100916) q[2];
sx q[2];
rz(2.2089843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2404732) q[1];
sx q[1];
rz(-0.34209004) q[1];
sx q[1];
rz(2.6035477) q[1];
rz(-0.9137093) q[3];
sx q[3];
rz(-1.6471225) q[3];
sx q[3];
rz(-2.9655991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38976321) q[2];
sx q[2];
rz(-1.7182173) q[2];
sx q[2];
rz(0.54330379) q[2];
rz(-0.027103847) q[3];
sx q[3];
rz(-2.6684941) q[3];
sx q[3];
rz(1.0075587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84697023) q[0];
sx q[0];
rz(-0.68055081) q[0];
sx q[0];
rz(-2.3934225) q[0];
rz(-1.6698042) q[1];
sx q[1];
rz(-1.7282149) q[1];
sx q[1];
rz(0.69563037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75891209) q[0];
sx q[0];
rz(-1.6864538) q[0];
sx q[0];
rz(1.8291436) q[0];
x q[1];
rz(2.8606121) q[2];
sx q[2];
rz(-2.4644445) q[2];
sx q[2];
rz(2.3317762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66113093) q[1];
sx q[1];
rz(-0.5782402) q[1];
sx q[1];
rz(0.64774143) q[1];
x q[2];
rz(-1.4616667) q[3];
sx q[3];
rz(-1.4617451) q[3];
sx q[3];
rz(-0.29666479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8814322) q[2];
sx q[2];
rz(-2.3023534) q[2];
sx q[2];
rz(3.0879171) q[2];
rz(-1.3850877) q[3];
sx q[3];
rz(-1.799492) q[3];
sx q[3];
rz(-0.16689859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3593035) q[0];
sx q[0];
rz(-1.8817236) q[0];
sx q[0];
rz(2.2480929) q[0];
rz(2.9529849) q[1];
sx q[1];
rz(-2.392417) q[1];
sx q[1];
rz(2.9395054) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30464393) q[0];
sx q[0];
rz(-0.82425398) q[0];
sx q[0];
rz(-3.1008284) q[0];
rz(-3.0918504) q[2];
sx q[2];
rz(-2.0936077) q[2];
sx q[2];
rz(-2.3573014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1043865) q[1];
sx q[1];
rz(-1.361711) q[1];
sx q[1];
rz(1.0250183) q[1];
rz(-1.8340183) q[3];
sx q[3];
rz(-1.9303476) q[3];
sx q[3];
rz(1.6040126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14822745) q[2];
sx q[2];
rz(-1.8348285) q[2];
sx q[2];
rz(0.6692878) q[2];
rz(2.3547442) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(2.440786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65449077) q[0];
sx q[0];
rz(-0.50906068) q[0];
sx q[0];
rz(1.2855726) q[0];
rz(0.14699832) q[1];
sx q[1];
rz(-1.1451984) q[1];
sx q[1];
rz(-2.1910892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5122172) q[0];
sx q[0];
rz(-0.75429082) q[0];
sx q[0];
rz(0.74441433) q[0];
x q[1];
rz(-0.062531907) q[2];
sx q[2];
rz(-2.1075296) q[2];
sx q[2];
rz(-1.5944531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78792469) q[1];
sx q[1];
rz(-1.0312925) q[1];
sx q[1];
rz(0.9307534) q[1];
rz(-pi) q[2];
rz(-0.70159961) q[3];
sx q[3];
rz(-2.9605013) q[3];
sx q[3];
rz(2.4871684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46108437) q[2];
sx q[2];
rz(-0.66404873) q[2];
sx q[2];
rz(2.7456679) q[2];
rz(2.8094273) q[3];
sx q[3];
rz(-1.1484523) q[3];
sx q[3];
rz(2.240326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.0300776) q[0];
sx q[0];
rz(-2.5319396) q[0];
sx q[0];
rz(1.4421705) q[0];
rz(-3.1051927) q[1];
sx q[1];
rz(-1.8489238) q[1];
sx q[1];
rz(1.363149) q[1];
rz(-1.2716952) q[2];
sx q[2];
rz(-0.80812412) q[2];
sx q[2];
rz(-3.0199188) q[2];
rz(3.0721138) q[3];
sx q[3];
rz(-2.4705349) q[3];
sx q[3];
rz(0.51746838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
