OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96707764) q[0];
sx q[0];
rz(-1.4580026) q[0];
sx q[0];
rz(-2.8414677) q[0];
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(-1.6405029) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97273443) q[0];
sx q[0];
rz(-1.5746813) q[0];
sx q[0];
rz(-0.015104276) q[0];
rz(-pi) q[1];
rz(-1.1741224) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(-1.8044745) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38084322) q[1];
sx q[1];
rz(-2.2135206) q[1];
sx q[1];
rz(1.8442783) q[1];
rz(-pi) q[2];
rz(-3.1367125) q[3];
sx q[3];
rz(-2.3327851) q[3];
sx q[3];
rz(0.67833662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7370721) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(-2.4386151) q[2];
rz(-0.56973714) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1341781) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(-1.2331569) q[0];
rz(-0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(-2.0436683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5220338) q[0];
sx q[0];
rz(-0.21169835) q[0];
sx q[0];
rz(-1.8926527) q[0];
rz(-pi) q[1];
rz(0.64867257) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(-1.5209188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.040702601) q[1];
sx q[1];
rz(-1.4933673) q[1];
sx q[1];
rz(3.1162683) q[1];
x q[2];
rz(3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7132831) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(1.606288) q[2];
rz(0.71896583) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035901) q[0];
sx q[0];
rz(-2.328023) q[0];
sx q[0];
rz(-2.549262) q[0];
rz(1.2447641) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62632769) q[0];
sx q[0];
rz(-2.0538616) q[0];
sx q[0];
rz(-1.417571) q[0];
x q[1];
rz(-1.4136366) q[2];
sx q[2];
rz(-2.4470235) q[2];
sx q[2];
rz(2.5099011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1864724) q[1];
sx q[1];
rz(-1.8545517) q[1];
sx q[1];
rz(-2.1452745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0823194) q[3];
sx q[3];
rz(-1.3311738) q[3];
sx q[3];
rz(0.75877305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.738872) q[2];
sx q[2];
rz(-1.0307743) q[2];
sx q[2];
rz(1.2325475) q[2];
rz(2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(3.1109911) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.4786913) q[1];
sx q[1];
rz(-2.0487002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83704889) q[0];
sx q[0];
rz(-2.2185159) q[0];
sx q[0];
rz(2.442028) q[0];
rz(-pi) q[1];
rz(-2.3762114) q[2];
sx q[2];
rz(-1.5281406) q[2];
sx q[2];
rz(2.9431107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2177201) q[1];
sx q[1];
rz(-1.7752547) q[1];
sx q[1];
rz(-2.2108498) q[1];
x q[2];
rz(-2.2956156) q[3];
sx q[3];
rz(-2.1095157) q[3];
sx q[3];
rz(-1.908345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0577724) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(2.8687381) q[2];
rz(1.3911635) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(-2.1896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6735753) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(1.2982298) q[0];
rz(-2.4507554) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(-1.5265436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218791) q[0];
sx q[0];
rz(-1.0391616) q[0];
sx q[0];
rz(0.74703947) q[0];
rz(-pi) q[1];
rz(-0.66730209) q[2];
sx q[2];
rz(-2.7030253) q[2];
sx q[2];
rz(1.9090609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5461108) q[1];
sx q[1];
rz(-2.271436) q[1];
sx q[1];
rz(3.107772) q[1];
rz(0.36187343) q[3];
sx q[3];
rz(-1.2905777) q[3];
sx q[3];
rz(-1.1079196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7319298) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(-2.9742677) q[2];
rz(-1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(-2.4840241) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.14199) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(-1.863119) q[0];
rz(-1.7059884) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631133) q[0];
sx q[0];
rz(-2.7857384) q[0];
sx q[0];
rz(0.87219091) q[0];
rz(2.2907545) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(-1.6825324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3842114) q[1];
sx q[1];
rz(-2.7486292) q[1];
sx q[1];
rz(2.1860366) q[1];
x q[2];
rz(0.2578854) q[3];
sx q[3];
rz(-2.3688859) q[3];
sx q[3];
rz(1.7125285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(2.2806878) q[2];
rz(-1.143645) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025573108) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(2.386911) q[0];
rz(-2.2017551) q[1];
sx q[1];
rz(-2.266423) q[1];
sx q[1];
rz(-1.8584724) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59903501) q[0];
sx q[0];
rz(-1.577566) q[0];
sx q[0];
rz(-1.8775131) q[0];
x q[1];
rz(0.36410113) q[2];
sx q[2];
rz(-2.0023097) q[2];
sx q[2];
rz(1.3442775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93604532) q[1];
sx q[1];
rz(-0.48921555) q[1];
sx q[1];
rz(0.65566109) q[1];
rz(-pi) q[2];
rz(-0.96638443) q[3];
sx q[3];
rz(-2.3570286) q[3];
sx q[3];
rz(2.0469637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(-2.0580573) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.98899406) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(-0.75491828) q[0];
rz(-2.6424291) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(-2.4430433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4071199) q[0];
sx q[0];
rz(-1.880058) q[0];
sx q[0];
rz(0.63068181) q[0];
rz(0.28251799) q[2];
sx q[2];
rz(-2.660224) q[2];
sx q[2];
rz(-0.11061874) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6900151) q[1];
sx q[1];
rz(-1.5507586) q[1];
sx q[1];
rz(2.103014) q[1];
x q[2];
rz(-1.1517161) q[3];
sx q[3];
rz(-1.569783) q[3];
sx q[3];
rz(-2.5496497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2373206) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(0.80238706) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9629843) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-1.1075903) q[0];
rz(1.3811318) q[1];
sx q[1];
rz(-1.7778722) q[1];
sx q[1];
rz(-1.6345056) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3699781) q[0];
sx q[0];
rz(-0.42015892) q[0];
sx q[0];
rz(-2.7922947) q[0];
x q[1];
rz(-2.2198943) q[2];
sx q[2];
rz(-0.49921152) q[2];
sx q[2];
rz(1.2865024) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.448882) q[1];
sx q[1];
rz(-1.4254942) q[1];
sx q[1];
rz(2.7753745) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0595257) q[3];
sx q[3];
rz(-0.82857271) q[3];
sx q[3];
rz(0.16211331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(1.2723119) q[2];
rz(1.6614527) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(-0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965974) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(-3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(0.12558118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.784362) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(-3.1249376) q[0];
x q[1];
rz(2.0415885) q[2];
sx q[2];
rz(-0.94816899) q[2];
sx q[2];
rz(-2.560844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6148147) q[1];
sx q[1];
rz(-1.1452951) q[1];
sx q[1];
rz(2.6761901) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9405648) q[3];
sx q[3];
rz(-1.8327692) q[3];
sx q[3];
rz(-2.6233545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(0.8030836) q[2];
rz(-2.965029) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17978996) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(1.5079458) q[2];
sx q[2];
rz(-1.812915) q[2];
sx q[2];
rz(-1.6356638) q[2];
rz(0.17038067) q[3];
sx q[3];
rz(-0.34210404) q[3];
sx q[3];
rz(0.31755527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
