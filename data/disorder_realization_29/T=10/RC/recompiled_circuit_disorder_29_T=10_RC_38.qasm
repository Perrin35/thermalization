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
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(-1.7915373) q[1];
sx q[1];
rz(-1.5426558) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60228981) q[0];
sx q[0];
rz(-1.7825583) q[0];
sx q[0];
rz(-1.7913626) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5390981) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(-2.9162625) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5988679) q[1];
sx q[1];
rz(-1.8017144) q[1];
sx q[1];
rz(2.8938107) q[1];
rz(-pi) q[2];
rz(-1.5845675) q[3];
sx q[3];
rz(-0.78409401) q[3];
sx q[3];
rz(-0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8618384) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-2.2564783) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136114) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(1.0990748) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-2.2639993) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1359522) q[0];
sx q[0];
rz(-0.71605039) q[0];
sx q[0];
rz(-0.27766772) q[0];
rz(-pi) q[1];
rz(1.4363102) q[2];
sx q[2];
rz(-1.9324537) q[2];
sx q[2];
rz(-2.550617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9237325) q[1];
sx q[1];
rz(-1.6777615) q[1];
sx q[1];
rz(3.1080493) q[1];
rz(-pi) q[2];
rz(-2.8072076) q[3];
sx q[3];
rz(-1.9352203) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(-1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(-0.6668123) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(-0.07382948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66729546) q[0];
sx q[0];
rz(-1.6601666) q[0];
sx q[0];
rz(-1.5090293) q[0];
rz(-0.87186558) q[2];
sx q[2];
rz(-2.7060894) q[2];
sx q[2];
rz(1.0812024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.044588305) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(-1.3929277) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35145268) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(-0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6510216) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9639503) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98722092) q[0];
sx q[0];
rz(-2.4043596) q[0];
sx q[0];
rz(0.1371951) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32707733) q[2];
sx q[2];
rz(-1.8665258) q[2];
sx q[2];
rz(-2.0138274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23197933) q[1];
sx q[1];
rz(-0.59514272) q[1];
sx q[1];
rz(0.91722782) q[1];
rz(1.7014245) q[3];
sx q[3];
rz(-1.6960973) q[3];
sx q[3];
rz(3.0076722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.556095) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.8251098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3554768) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(-2.064408) q[0];
rz(-pi) q[1];
rz(-1.8226132) q[2];
sx q[2];
rz(-1.8251112) q[2];
sx q[2];
rz(-0.11676678) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74354467) q[1];
sx q[1];
rz(-0.053357031) q[1];
sx q[1];
rz(1.3392901) q[1];
x q[2];
rz(-0.76977323) q[3];
sx q[3];
rz(-2.4379745) q[3];
sx q[3];
rz(-0.15042703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(2.1300952) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(-0.85995752) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.3202753) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1307756) q[0];
sx q[0];
rz(-0.59033075) q[0];
sx q[0];
rz(2.4369795) q[0];
x q[1];
rz(1.3930416) q[2];
sx q[2];
rz(-1.0076992) q[2];
sx q[2];
rz(2.4515738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49530416) q[1];
sx q[1];
rz(-1.0083645) q[1];
sx q[1];
rz(-2.5763987) q[1];
x q[2];
rz(-0.28989132) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(-0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.3909891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334675) q[0];
sx q[0];
rz(-1.2398749) q[0];
sx q[0];
rz(0.56751859) q[0];
rz(0.58629845) q[2];
sx q[2];
rz(-2.7099897) q[2];
sx q[2];
rz(-1.4177711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.890306) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(-0.89819737) q[1];
rz(-pi) q[2];
rz(0.77565149) q[3];
sx q[3];
rz(-1.5175879) q[3];
sx q[3];
rz(0.42078161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(-2.7461046) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46273461) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(1.4920374) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.7038201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2742845) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(-pi) q[1];
rz(-1.4346801) q[2];
sx q[2];
rz(-1.4021177) q[2];
sx q[2];
rz(-0.31751925) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2807323) q[1];
sx q[1];
rz(-2.3128465) q[1];
sx q[1];
rz(-2.4651205) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.612547) q[3];
sx q[3];
rz(-0.44781993) q[3];
sx q[3];
rz(-0.24275045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64547223) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20539595) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(-1.0937011) q[0];
rz(-0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-2.0827983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7208784) q[0];
sx q[0];
rz(-1.2982561) q[0];
sx q[0];
rz(2.1331302) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86973127) q[2];
sx q[2];
rz(-2.5517002) q[2];
sx q[2];
rz(-1.9260977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79566369) q[1];
sx q[1];
rz(-1.1946031) q[1];
sx q[1];
rz(0.13161195) q[1];
rz(-1.9731673) q[3];
sx q[3];
rz(-1.1306922) q[3];
sx q[3];
rz(0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-0.69512853) q[2];
sx q[2];
rz(-1.4808562) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.4321009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(-1.2105932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6851995) q[2];
sx q[2];
rz(-2.6420339) q[2];
sx q[2];
rz(1.511614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0309279) q[1];
sx q[1];
rz(-1.0193766) q[1];
sx q[1];
rz(-1.621978) q[1];
rz(-pi) q[2];
rz(1.3833952) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(-0.24500971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6802406) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(-1.1414026) q[2];
rz(1.6067778) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(2.6721568) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948982) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(2.4070807) q[2];
sx q[2];
rz(-0.30167689) q[2];
sx q[2];
rz(-2.6405356) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-1.9867298) q[3];
sx q[3];
rz(0.17476535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
