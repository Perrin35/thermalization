OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7862406) q[0];
sx q[0];
rz(-2.8347637) q[0];
sx q[0];
rz(2.8428349) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(-0.63036418) q[1];
sx q[1];
rz(-1.6685553) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0399892) q[0];
sx q[0];
rz(-0.77058661) q[0];
sx q[0];
rz(-1.6165074) q[0];
rz(-1.2700419) q[2];
sx q[2];
rz(-1.5879121) q[2];
sx q[2];
rz(-0.37671664) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70570488) q[1];
sx q[1];
rz(-2.0667564) q[1];
sx q[1];
rz(-0.75261799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24362917) q[3];
sx q[3];
rz(-0.28377658) q[3];
sx q[3];
rz(0.42967969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95639688) q[2];
sx q[2];
rz(-0.38072017) q[2];
sx q[2];
rz(-1.3107276) q[2];
rz(-0.091536097) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(0.53102791) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(1.1473038) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(-2.4400585) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6525878) q[0];
sx q[0];
rz(-0.93237823) q[0];
sx q[0];
rz(2.5030563) q[0];
rz(2.0364136) q[2];
sx q[2];
rz(-1.1934416) q[2];
sx q[2];
rz(1.5355664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5614723) q[1];
sx q[1];
rz(-2.0144723) q[1];
sx q[1];
rz(2.2444112) q[1];
rz(-pi) q[2];
rz(-0.40579943) q[3];
sx q[3];
rz(-2.4630952) q[3];
sx q[3];
rz(-3.0140227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39701617) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(-1.5783295) q[2];
rz(-0.27966106) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(-2.7320812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3438943) q[0];
sx q[0];
rz(-2.7293623) q[0];
sx q[0];
rz(-1.718234) q[0];
rz(-0.1238981) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.5843676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9950969) q[0];
sx q[0];
rz(-2.2687758) q[0];
sx q[0];
rz(-0.3905889) q[0];
x q[1];
rz(1.9818457) q[2];
sx q[2];
rz(-2.374994) q[2];
sx q[2];
rz(-1.3992978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1184956) q[1];
sx q[1];
rz(-1.4916991) q[1];
sx q[1];
rz(1.3239856) q[1];
rz(-2.4361472) q[3];
sx q[3];
rz(-1.0407606) q[3];
sx q[3];
rz(-0.3769266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.065217) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(-2.2067113) q[2];
rz(2.8377418) q[3];
sx q[3];
rz(-2.5287718) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3728751) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(2.8745162) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-2.5649773) q[1];
sx q[1];
rz(-2.3042302) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1601279) q[0];
sx q[0];
rz(-2.5257887) q[0];
sx q[0];
rz(1.0976904) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5576511) q[2];
sx q[2];
rz(-0.9204694) q[2];
sx q[2];
rz(-0.038829858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4761047) q[1];
sx q[1];
rz(-1.3766411) q[1];
sx q[1];
rz(0.23745115) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.046810026) q[3];
sx q[3];
rz(-1.6322281) q[3];
sx q[3];
rz(2.2991179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5460633) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(-1.1692125) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(-2.9235212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021407481) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(2.2688493) q[0];
rz(-1.1574289) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(0.016955888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30045045) q[0];
sx q[0];
rz(-0.73472856) q[0];
sx q[0];
rz(0.018201306) q[0];
rz(-pi) q[1];
rz(1.5808231) q[2];
sx q[2];
rz(-2.1705049) q[2];
sx q[2];
rz(-1.6184652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53345799) q[1];
sx q[1];
rz(-2.3124586) q[1];
sx q[1];
rz(-1.9749179) q[1];
rz(-pi) q[2];
rz(-1.8618552) q[3];
sx q[3];
rz(-1.8292871) q[3];
sx q[3];
rz(-0.84421009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0852647) q[2];
sx q[2];
rz(-1.4339829) q[2];
sx q[2];
rz(2.8223574) q[2];
rz(2.4625835) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(2.3398248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.7769258) q[0];
sx q[0];
rz(-0.32061446) q[0];
sx q[0];
rz(-0.6426245) q[0];
rz(1.7721666) q[1];
sx q[1];
rz(-2.5182928) q[1];
sx q[1];
rz(-0.97698897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5394655) q[0];
sx q[0];
rz(-1.4475679) q[0];
sx q[0];
rz(-1.4719523) q[0];
rz(2.0989292) q[2];
sx q[2];
rz(-1.1234049) q[2];
sx q[2];
rz(1.3716979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20583892) q[1];
sx q[1];
rz(-2.4495857) q[1];
sx q[1];
rz(-1.5490973) q[1];
rz(1.4897519) q[3];
sx q[3];
rz(-0.5872927) q[3];
sx q[3];
rz(-2.0964266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52648181) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(-1.9132445) q[2];
rz(-0.86722106) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(0.76798463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277471) q[0];
sx q[0];
rz(-2.9242046) q[0];
sx q[0];
rz(-2.9687498) q[0];
rz(2.3364283) q[1];
sx q[1];
rz(-2.5925437) q[1];
sx q[1];
rz(3.0056312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522536) q[0];
sx q[0];
rz(-0.87119188) q[0];
sx q[0];
rz(-0.1564349) q[0];
x q[1];
rz(-2.3195761) q[2];
sx q[2];
rz(-1.4944585) q[2];
sx q[2];
rz(-0.99886307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.48366212) q[1];
sx q[1];
rz(-1.6854981) q[1];
sx q[1];
rz(1.5516267) q[1];
rz(1.7053267) q[3];
sx q[3];
rz(-1.2438275) q[3];
sx q[3];
rz(3.0643058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2148296) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(-2.9570441) q[2];
rz(2.9522225) q[3];
sx q[3];
rz(-2.6034077) q[3];
sx q[3];
rz(0.93809938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1592584) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(0.92639297) q[0];
rz(3.1023846) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(-2.1047986) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85425692) q[0];
sx q[0];
rz(-1.2764036) q[0];
sx q[0];
rz(-2.9578985) q[0];
rz(0.37929566) q[2];
sx q[2];
rz(-0.46078983) q[2];
sx q[2];
rz(0.57578218) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9766624) q[1];
sx q[1];
rz(-1.9769985) q[1];
sx q[1];
rz(-1.3066533) q[1];
rz(1.3831656) q[3];
sx q[3];
rz(-1.1867282) q[3];
sx q[3];
rz(1.7373067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2490273) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(-2.3426775) q[2];
rz(0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(-2.5911205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35140458) q[0];
sx q[0];
rz(-0.0049954448) q[0];
sx q[0];
rz(1.5193526) q[0];
rz(-1.5937357) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(-0.69040745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6192157) q[0];
sx q[0];
rz(-1.3478312) q[0];
sx q[0];
rz(2.7821543) q[0];
x q[1];
rz(-3.0418441) q[2];
sx q[2];
rz(-0.80173758) q[2];
sx q[2];
rz(0.24101098) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4286297) q[1];
sx q[1];
rz(-0.74364122) q[1];
sx q[1];
rz(2.3548467) q[1];
rz(1.810174) q[3];
sx q[3];
rz(-0.55174151) q[3];
sx q[3];
rz(2.5249979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-0.60882336) q[2];
sx q[2];
rz(2.1717066) q[2];
rz(-2.8841833) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.6064706) q[0];
sx q[0];
rz(-2.4115925) q[0];
sx q[0];
rz(3.0560793) q[0];
rz(0.27587786) q[1];
sx q[1];
rz(-2.52067) q[1];
sx q[1];
rz(1.1891018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4258041) q[0];
sx q[0];
rz(-0.5606519) q[0];
sx q[0];
rz(2.406757) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8315808) q[2];
sx q[2];
rz(-1.5793206) q[2];
sx q[2];
rz(3.0954158) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.117529) q[1];
sx q[1];
rz(-1.6895362) q[1];
sx q[1];
rz(1.1681144) q[1];
rz(-pi) q[2];
rz(-2.889685) q[3];
sx q[3];
rz(-2.7828476) q[3];
sx q[3];
rz(-0.22049604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4756061) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(2.6574668) q[2];
rz(-3.0020946) q[3];
sx q[3];
rz(-0.77553427) q[3];
sx q[3];
rz(-0.16105306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(-1.4576661) q[1];
sx q[1];
rz(-1.4977581) q[1];
sx q[1];
rz(-1.3118634) q[1];
rz(2.8756782) q[2];
sx q[2];
rz(-0.66638808) q[2];
sx q[2];
rz(-0.93925977) q[2];
rz(3.1112989) q[3];
sx q[3];
rz(-2.3843063) q[3];
sx q[3];
rz(1.6142308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
