OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(-0.32615647) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.125995) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(0.21685812) q[0];
x q[1];
rz(0.60249451) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(0.22533016) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5988679) q[1];
sx q[1];
rz(-1.8017144) q[1];
sx q[1];
rz(-2.8938107) q[1];
x q[2];
rz(-0.78674973) q[3];
sx q[3];
rz(-1.5610715) q[3];
sx q[3];
rz(1.7973289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-2.2564783) q[2];
rz(0.72201133) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(2.0425178) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(-0.87759334) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7769593) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(2.4448256) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36466937) q[2];
sx q[2];
rz(-1.4450577) q[2];
sx q[2];
rz(-1.0276577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3493537) q[1];
sx q[1];
rz(-1.6041479) q[1];
sx q[1];
rz(-1.4637714) q[1];
x q[2];
rz(2.8072076) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(-0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(-2.8515942) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(0.07382948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66729546) q[0];
sx q[0];
rz(-1.4814261) q[0];
sx q[0];
rz(-1.6325634) q[0];
x q[1];
rz(-1.9129842) q[2];
sx q[2];
rz(-1.295919) q[2];
sx q[2];
rz(0.16155044) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9170345) q[1];
sx q[1];
rz(-0.18637603) q[1];
sx q[1];
rz(-1.8774377) q[1];
rz(1.6894433) q[3];
sx q[3];
rz(-1.2216179) q[3];
sx q[3];
rz(-2.196764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(1.1289319) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9639503) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(-1.9529163) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(1.7046938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68543363) q[0];
sx q[0];
rz(-1.6628656) q[0];
sx q[0];
rz(2.409056) q[0];
rz(1.2595348) q[2];
sx q[2];
rz(-1.2584104) q[2];
sx q[2];
rz(2.600008) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9096133) q[1];
sx q[1];
rz(-0.59514272) q[1];
sx q[1];
rz(2.2243648) q[1];
x q[2];
rz(-1.7014245) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(3.0076722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(-1.1222703) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(0.22398082) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1236498) q[0];
sx q[0];
rz(-1.1153478) q[0];
sx q[0];
rz(-2.7148867) q[0];
x q[1];
rz(2.3773642) q[2];
sx q[2];
rz(-0.35596213) q[2];
sx q[2];
rz(-0.9133577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1662235) q[1];
sx q[1];
rz(-1.6227286) q[1];
sx q[1];
rz(3.1293392) q[1];
x q[2];
rz(-0.76977323) q[3];
sx q[3];
rz(-0.7036182) q[3];
sx q[3];
rz(-2.9911656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0405154) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(2.6082883) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.8213173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1307756) q[0];
sx q[0];
rz(-0.59033075) q[0];
sx q[0];
rz(-0.70461313) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27310246) q[2];
sx q[2];
rz(-2.5540076) q[2];
sx q[2];
rz(-2.7762129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49530416) q[1];
sx q[1];
rz(-2.1332281) q[1];
sx q[1];
rz(-0.56519392) q[1];
rz(2.8517013) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0423353) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(0.60097224) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(-1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(3.1069966) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.3909891) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334675) q[0];
sx q[0];
rz(-1.2398749) q[0];
sx q[0];
rz(-0.56751859) q[0];
rz(-pi) q[1];
rz(2.7752635) q[2];
sx q[2];
rz(-1.804367) q[2];
sx q[2];
rz(-0.69586588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(0.89819737) q[1];
rz(-2.3659412) q[3];
sx q[3];
rz(-1.5175879) q[3];
sx q[3];
rz(-2.720811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81327072) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(2.7461046) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(1.4377726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4576365) q[0];
sx q[0];
rz(-2.5914765) q[0];
sx q[0];
rz(1.9253299) q[0];
rz(-pi) q[1];
rz(0.17022325) q[2];
sx q[2];
rz(-1.7049689) q[2];
sx q[2];
rz(1.2302878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7319665) q[1];
sx q[1];
rz(-2.1831174) q[1];
sx q[1];
rz(0.97169505) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1215454) q[3];
sx q[3];
rz(-2.0181977) q[3];
sx q[3];
rz(-2.945154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64547223) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(2.2802165) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(-0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20539595) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(2.4049092) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(2.0827983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4207142) q[0];
sx q[0];
rz(-1.8433365) q[0];
sx q[0];
rz(-2.1331302) q[0];
rz(2.7339897) q[2];
sx q[2];
rz(-2.0098445) q[2];
sx q[2];
rz(-1.1328732) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79566369) q[1];
sx q[1];
rz(-1.9469896) q[1];
sx q[1];
rz(-0.13161195) q[1];
x q[2];
rz(-1.1684253) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(-2.7311834) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.8695628) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(-1.4321009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390022) q[0];
sx q[0];
rz(-1.0800835) q[0];
sx q[0];
rz(1.9309994) q[0];
rz(-2.6851995) q[2];
sx q[2];
rz(-0.49955873) q[2];
sx q[2];
rz(-1.6299786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4869625) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(-0.55200465) q[1];
rz(-pi) q[2];
rz(1.7581975) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(0.24500971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6802406) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(2.0001901) q[2];
rz(1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(2.6721568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(1.776406) q[2];
sx q[2];
rz(-1.3484577) q[2];
sx q[2];
rz(1.2586013) q[2];
rz(-2.7183919) q[3];
sx q[3];
rz(-1.7508218) q[3];
sx q[3];
rz(-1.4765061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
