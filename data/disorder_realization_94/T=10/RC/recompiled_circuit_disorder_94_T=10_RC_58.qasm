OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(4.6946445) q[0];
sx q[0];
rz(10.932218) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98696729) q[0];
sx q[0];
rz(-1.1095424) q[0];
sx q[0];
rz(-2.2485562) q[0];
rz(-2.7841714) q[2];
sx q[2];
rz(-1.3961785) q[2];
sx q[2];
rz(1.071196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33949172) q[1];
sx q[1];
rz(-1.9994945) q[1];
sx q[1];
rz(2.2080253) q[1];
x q[2];
rz(0.30652133) q[3];
sx q[3];
rz(-1.3984826) q[3];
sx q[3];
rz(2.4337208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-0.33828503) q[2];
rz(-1.4398549) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-3.1112444) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47027662) q[0];
sx q[0];
rz(-2.94194) q[0];
sx q[0];
rz(-1.5623564) q[0];
rz(-pi) q[1];
rz(3.121071) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(3.0128535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46088947) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(2.1417888) q[1];
x q[2];
rz(-1.1094692) q[3];
sx q[3];
rz(-2.9328049) q[3];
sx q[3];
rz(1.3267645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(1.3267481) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-2.8258064) q[0];
rz(2.2029927) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2825851) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(-0.99143272) q[0];
rz(2.2601068) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(1.880868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4430868) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(0.94633533) q[1];
x q[2];
rz(-0.39450816) q[3];
sx q[3];
rz(-1.9022226) q[3];
sx q[3];
rz(3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0079460572) q[0];
sx q[0];
rz(-0.61230731) q[0];
sx q[0];
rz(2.4164819) q[0];
x q[1];
rz(-2.6949276) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(2.5926673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69869631) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(-0.15028468) q[1];
rz(0.96890038) q[3];
sx q[3];
rz(-0.36877353) q[3];
sx q[3];
rz(-0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3209155) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(0.25462338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4655612) q[0];
sx q[0];
rz(-1.6499632) q[0];
sx q[0];
rz(-0.013750793) q[0];
rz(-pi) q[1];
rz(2.1483634) q[2];
sx q[2];
rz(-1.5585871) q[2];
sx q[2];
rz(0.73405594) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9533206) q[1];
sx q[1];
rz(-2.4362262) q[1];
sx q[1];
rz(-0.377368) q[1];
rz(-2.0650495) q[3];
sx q[3];
rz(-2.3688865) q[3];
sx q[3];
rz(-0.62437526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(0.20656955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26039133) q[0];
sx q[0];
rz(-2.3058389) q[0];
sx q[0];
rz(-0.3221237) q[0];
rz(1.0429522) q[2];
sx q[2];
rz(-1.2252508) q[2];
sx q[2];
rz(0.31574677) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6321322) q[1];
sx q[1];
rz(-1.7473979) q[1];
sx q[1];
rz(-2.280974) q[1];
rz(2.2264678) q[3];
sx q[3];
rz(-1.8135999) q[3];
sx q[3];
rz(-1.1806928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(1.3279703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-2.6866331) q[0];
sx q[0];
rz(1.8332464) q[0];
x q[1];
rz(2.1778657) q[2];
sx q[2];
rz(-2.4237195) q[2];
sx q[2];
rz(-1.3340064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2139637) q[1];
sx q[1];
rz(-0.8968401) q[1];
sx q[1];
rz(-0.2143292) q[1];
rz(-1.9973203) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(1.1748479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(1.4303738) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5247105) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(1.3046706) q[2];
sx q[2];
rz(-1.5106491) q[2];
sx q[2];
rz(-1.6595449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9329405) q[1];
sx q[1];
rz(-1.3622074) q[1];
sx q[1];
rz(-1.5044466) q[1];
rz(-0.37787921) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(0.32285238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-0.075909464) q[2];
rz(-0.54801303) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.7653718) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8575681) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(0.88733034) q[0];
rz(2.9218036) q[2];
sx q[2];
rz(-0.79384365) q[2];
sx q[2];
rz(1.1500051) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0162504) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(2.7793105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.057007313) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.9407335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.8359258) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(2.1059039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10782345) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(-2.4643722) q[0];
rz(-pi) q[1];
rz(-2.5752441) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(-1.4373506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4020821) q[1];
sx q[1];
rz(-1.2716736) q[1];
sx q[1];
rz(1.8226536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8519248) q[3];
sx q[3];
rz(-0.77054671) q[3];
sx q[3];
rz(-0.25911261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(2.776896) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(0.097461854) q[2];
sx q[2];
rz(-2.7290191) q[2];
sx q[2];
rz(1.4636427) q[2];
rz(0.012398331) q[3];
sx q[3];
rz(-2.5231902) q[3];
sx q[3];
rz(-1.4004422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
