OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(4.6255339) q[0];
sx q[0];
rz(12.892527) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.125995) q[0];
sx q[0];
rz(-1.3552357) q[0];
sx q[0];
rz(0.21685812) q[0];
rz(-pi) q[1];
rz(-0.60249451) q[2];
sx q[2];
rz(-1.7817111) q[2];
sx q[2];
rz(0.22533016) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3781158) q[1];
sx q[1];
rz(-0.33707481) q[1];
sx q[1];
rz(0.76428767) q[1];
x q[2];
rz(1.5845675) q[3];
sx q[3];
rz(-2.3574986) q[3];
sx q[3];
rz(-0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2797543) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(-3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279813) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-2.2639993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646334) q[0];
sx q[0];
rz(-1.7517125) q[0];
sx q[0];
rz(-0.69676708) q[0];
rz(0.36466937) q[2];
sx q[2];
rz(-1.4450577) q[2];
sx q[2];
rz(-2.113935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5223915) q[1];
sx q[1];
rz(-3.0295105) q[1];
sx q[1];
rz(-1.8735318) q[1];
rz(-0.33438501) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(-0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.4484423) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(0.07382948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66729546) q[0];
sx q[0];
rz(-1.4814261) q[0];
sx q[0];
rz(1.6325634) q[0];
rz(1.9129842) q[2];
sx q[2];
rz(-1.295919) q[2];
sx q[2];
rz(-0.16155044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5362629) q[1];
sx q[1];
rz(-1.3932091) q[1];
sx q[1];
rz(-3.0847343) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4521493) q[3];
sx q[3];
rz(-1.2216179) q[3];
sx q[3];
rz(0.94482869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(0.64669615) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98722092) q[0];
sx q[0];
rz(-0.73723307) q[0];
sx q[0];
rz(-0.1371951) q[0];
rz(-pi) q[1];
rz(0.32707733) q[2];
sx q[2];
rz(-1.8665258) q[2];
sx q[2];
rz(2.0138274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77364591) q[1];
sx q[1];
rz(-1.2229496) q[1];
sx q[1];
rz(1.0775954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8021637) q[3];
sx q[3];
rz(-0.18076104) q[3];
sx q[3];
rz(-0.67644955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58549762) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(2.7686152) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81985695) q[0];
sx q[0];
rz(-0.61367354) q[0];
sx q[0];
rz(-0.86921285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3773642) q[2];
sx q[2];
rz(-2.7856305) q[2];
sx q[2];
rz(2.228235) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1662235) q[1];
sx q[1];
rz(-1.6227286) q[1];
sx q[1];
rz(3.1293392) q[1];
rz(1.0373711) q[3];
sx q[3];
rz(-2.0539527) q[3];
sx q[3];
rz(-2.0869568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0405154) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(2.6082883) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(-2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(2.2816351) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(-1.3202753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1307756) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-0.70461313) q[0];
x q[1];
rz(1.3930416) q[2];
sx q[2];
rz(-2.1338935) q[2];
sx q[2];
rz(-2.4515738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37636138) q[1];
sx q[1];
rz(-2.3665641) q[1];
sx q[1];
rz(2.2750957) q[1];
rz(1.201186) q[3];
sx q[3];
rz(-0.69023057) q[3];
sx q[3];
rz(1.6604916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0423353) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(-0.48505923) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(-3.1069966) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(1.3909891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334675) q[0];
sx q[0];
rz(-1.9017178) q[0];
sx q[0];
rz(-2.5740741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.820307) q[2];
sx q[2];
rz(-1.9267285) q[2];
sx q[2];
rz(-2.3552259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78087378) q[1];
sx q[1];
rz(-2.2370173) q[1];
sx q[1];
rz(2.9800376) q[1];
x q[2];
rz(-2.3659412) q[3];
sx q[3];
rz(-1.6240048) q[3];
sx q[3];
rz(-0.42078161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3283219) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(-1.8528806) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.46273461) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(1.4920374) q[0];
rz(0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.4377726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8673082) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67280025) q[2];
sx q[2];
rz(-0.2163419) q[2];
sx q[2];
rz(-2.1397482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2807323) q[1];
sx q[1];
rz(-0.8287462) q[1];
sx q[1];
rz(2.4651205) q[1];
rz(-1.5290456) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(2.8988422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64547223) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(2.9628741) q[2];
rz(-2.2802165) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20539595) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(-1.0937011) q[0];
rz(-2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(2.0827983) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537542) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(2.0536325) q[0];
x q[1];
rz(0.40760298) q[2];
sx q[2];
rz(-1.1317481) q[2];
sx q[2];
rz(2.0087194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4150548) q[1];
sx q[1];
rz(-1.448436) q[1];
sx q[1];
rz(-1.1916257) q[1];
rz(-pi) q[2];
rz(1.1684253) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8686691) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(-2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.8695628) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(-0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(-1.2105932) q[0];
rz(-1.3347689) q[2];
sx q[2];
rz(-2.0152976) q[2];
sx q[2];
rz(-2.1399463) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0309279) q[1];
sx q[1];
rz(-1.0193766) q[1];
sx q[1];
rz(-1.621978) q[1];
rz(-pi) q[2];
rz(-2.0026783) q[3];
sx q[3];
rz(-1.6500041) q[3];
sx q[3];
rz(-1.495804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(-1.6067778) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5948982) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(-2.7728511) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(2.914631) q[2];
sx q[2];
rz(-1.3703176) q[2];
sx q[2];
rz(-0.3581518) q[2];
rz(2.7244444) q[3];
sx q[3];
rz(-0.45776164) q[3];
sx q[3];
rz(-0.2840857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
