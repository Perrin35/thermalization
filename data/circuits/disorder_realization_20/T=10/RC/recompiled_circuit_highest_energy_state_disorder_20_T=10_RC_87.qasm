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
rz(-2.1751997) q[0];
sx q[0];
rz(4.8973358) q[0];
sx q[0];
rz(10.020221) q[0];
rz(1.2915986) q[1];
sx q[1];
rz(-2.5505677) q[1];
sx q[1];
rz(-0.12330595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57014556) q[0];
sx q[0];
rz(-0.71200221) q[0];
sx q[0];
rz(0.90306905) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78211981) q[2];
sx q[2];
rz(-2.189866) q[2];
sx q[2];
rz(-1.182488) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1327269) q[1];
sx q[1];
rz(-1.4767329) q[1];
sx q[1];
rz(1.7857185) q[1];
x q[2];
rz(-2.7601065) q[3];
sx q[3];
rz(-2.7499928) q[3];
sx q[3];
rz(-2.7158383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3143602) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(-3.0296791) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.7258485) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(-0.14990212) q[0];
rz(-2.8556178) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(1.1289977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492831) q[0];
sx q[0];
rz(-0.2899) q[0];
sx q[0];
rz(1.4647746) q[0];
rz(-pi) q[1];
rz(3.1270119) q[2];
sx q[2];
rz(-1.5111856) q[2];
sx q[2];
rz(0.73734847) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2098238) q[1];
sx q[1];
rz(-1.9289513) q[1];
sx q[1];
rz(-1.7478554) q[1];
rz(2.3150857) q[3];
sx q[3];
rz(-0.77946589) q[3];
sx q[3];
rz(-2.731088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4631606) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(1.1506259) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(-3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(2.7599957) q[0];
rz(-1.2941788) q[1];
sx q[1];
rz(-2.0353863) q[1];
sx q[1];
rz(0.19827422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851882) q[0];
sx q[0];
rz(-1.745178) q[0];
sx q[0];
rz(1.3402433) q[0];
rz(-pi) q[1];
rz(-1.6724104) q[2];
sx q[2];
rz(-2.4696015) q[2];
sx q[2];
rz(0.72173126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1041996) q[1];
sx q[1];
rz(-1.3379119) q[1];
sx q[1];
rz(1.7892455) q[1];
x q[2];
rz(0.020129344) q[3];
sx q[3];
rz(-1.745589) q[3];
sx q[3];
rz(0.59390261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76564378) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(-0.51885968) q[2];
rz(1.8910003) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(-2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64312235) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(-0.44922391) q[0];
rz(-1.6953702) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(2.5208688) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0064471) q[0];
sx q[0];
rz(-2.2250611) q[0];
sx q[0];
rz(1.3885137) q[0];
x q[1];
rz(2.5081257) q[2];
sx q[2];
rz(-1.4711498) q[2];
sx q[2];
rz(1.5657305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6606969) q[1];
sx q[1];
rz(-2.6891853) q[1];
sx q[1];
rz(-0.47641944) q[1];
rz(-pi) q[2];
rz(-1.659538) q[3];
sx q[3];
rz(-0.54026287) q[3];
sx q[3];
rz(-1.0080573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(2.980496) q[2];
rz(2.7437239) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.38254) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(1.8880728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81957626) q[0];
sx q[0];
rz(-1.6068234) q[0];
sx q[0];
rz(-3.1234804) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3316111) q[2];
sx q[2];
rz(-1.7415541) q[2];
sx q[2];
rz(2.8851938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1977928) q[1];
sx q[1];
rz(-1.6159004) q[1];
sx q[1];
rz(-0.79774858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7944466) q[3];
sx q[3];
rz(-1.5724564) q[3];
sx q[3];
rz(-1.1356789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72838655) q[2];
sx q[2];
rz(-1.9325958) q[2];
sx q[2];
rz(-1.8750635) q[2];
rz(-2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(-2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(0.025064502) q[0];
rz(-1.2114245) q[1];
sx q[1];
rz(-1.2592659) q[1];
sx q[1];
rz(2.8866344) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1392649) q[0];
sx q[0];
rz(-0.049590913) q[0];
sx q[0];
rz(-3.0399486) q[0];
x q[1];
rz(1.0827052) q[2];
sx q[2];
rz(-1.0399264) q[2];
sx q[2];
rz(-1.7386029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7595664) q[1];
sx q[1];
rz(-2.3236378) q[1];
sx q[1];
rz(-0.76063971) q[1];
x q[2];
rz(-0.21277748) q[3];
sx q[3];
rz(-1.0085953) q[3];
sx q[3];
rz(2.8977545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-1.0422372) q[2];
sx q[2];
rz(2.6825405) q[2];
rz(-1.6445271) q[3];
sx q[3];
rz(-3.0247757) q[3];
sx q[3];
rz(0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(-0.086061867) q[0];
rz(0.91880265) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.930621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051832599) q[0];
sx q[0];
rz(-0.38516949) q[0];
sx q[0];
rz(-0.56398192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7007799) q[2];
sx q[2];
rz(-2.0939504) q[2];
sx q[2];
rz(-0.8280821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1001491) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(0.99975296) q[1];
rz(-pi) q[2];
rz(1.9350697) q[3];
sx q[3];
rz(-1.2281872) q[3];
sx q[3];
rz(0.87876696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3226037) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(2.4244579) q[2];
rz(-2.1142193) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-2.7023442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913919) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(0.014658654) q[0];
rz(0.75153366) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(-1.470648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58107046) q[0];
sx q[0];
rz(-1.7092686) q[0];
sx q[0];
rz(2.951589) q[0];
x q[1];
rz(-2.884955) q[2];
sx q[2];
rz(-1.4932639) q[2];
sx q[2];
rz(2.3687349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5353706) q[1];
sx q[1];
rz(-2.722123) q[1];
sx q[1];
rz(-1.2910045) q[1];
x q[2];
rz(3.1378391) q[3];
sx q[3];
rz(-1.225718) q[3];
sx q[3];
rz(2.0451289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8810001) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(-0.27349681) q[2];
rz(1.1547487) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1028033) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(-2.0661085) q[0];
rz(3.0134046) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(-0.34559616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3227279) q[0];
sx q[0];
rz(-0.68536192) q[0];
sx q[0];
rz(-1.8249874) q[0];
x q[1];
rz(1.0021474) q[2];
sx q[2];
rz(-1.2969742) q[2];
sx q[2];
rz(-1.3783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3010878) q[1];
sx q[1];
rz(-1.5078168) q[1];
sx q[1];
rz(2.5999864) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4050499) q[3];
sx q[3];
rz(-2.1624544) q[3];
sx q[3];
rz(1.4480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(-1.5209939) q[2];
rz(1.3711551) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(-2.7267406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82776752) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(0.23571043) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.3689573) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536606) q[0];
sx q[0];
rz(-2.2412123) q[0];
sx q[0];
rz(-2.5817342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4659285) q[2];
sx q[2];
rz(-2.0620637) q[2];
sx q[2];
rz(1.5671687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5342218) q[1];
sx q[1];
rz(-1.0080308) q[1];
sx q[1];
rz(-2.9396179) q[1];
x q[2];
rz(-0.24457358) q[3];
sx q[3];
rz(-0.37062708) q[3];
sx q[3];
rz(-1.2115492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77091757) q[2];
sx q[2];
rz(-1.8149899) q[2];
sx q[2];
rz(0.053000432) q[2];
rz(-0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(-0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5175405) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(-1.7908295) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(2.9861957) q[2];
sx q[2];
rz(-0.70320319) q[2];
sx q[2];
rz(2.7413766) q[2];
rz(0.97196058) q[3];
sx q[3];
rz(-2.2836015) q[3];
sx q[3];
rz(-2.9290269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
