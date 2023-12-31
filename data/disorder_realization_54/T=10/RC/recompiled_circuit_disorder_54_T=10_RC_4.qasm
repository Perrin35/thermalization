OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(2.4324774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37492) q[0];
sx q[0];
rz(-1.8639038) q[0];
sx q[0];
rz(2.2105182) q[0];
rz(-1.66967) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(-1.0613943) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4830556) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(-0.87941054) q[1];
x q[2];
rz(0.96592824) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(-1.0498429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(-0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(2.2252749) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5641862) q[0];
sx q[0];
rz(-2.2803377) q[0];
sx q[0];
rz(2.2400411) q[0];
x q[1];
rz(3.0099478) q[2];
sx q[2];
rz(-0.83877124) q[2];
sx q[2];
rz(1.6694348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.037109) q[1];
sx q[1];
rz(-2.6625405) q[1];
sx q[1];
rz(0.80161174) q[1];
rz(-pi) q[2];
rz(-1.6244435) q[3];
sx q[3];
rz(-2.7430153) q[3];
sx q[3];
rz(-0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(-3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(0.78330529) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(2.6170513) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4315955) q[0];
sx q[0];
rz(-1.4072197) q[0];
sx q[0];
rz(-3.1113935) q[0];
x q[1];
rz(0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7150535) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(-3.0798562) q[1];
rz(-pi) q[2];
rz(-2.9933661) q[3];
sx q[3];
rz(-1.9865415) q[3];
sx q[3];
rz(-1.1566597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(2.7041919) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(2.451992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3084429) q[0];
sx q[0];
rz(-1.8652548) q[0];
sx q[0];
rz(-1.4008646) q[0];
x q[1];
rz(-1.7961411) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5213485) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(2.0116624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.425399) q[3];
sx q[3];
rz(-2.0510011) q[3];
sx q[3];
rz(2.8499545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(0.62292567) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.4978283) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935055) q[0];
sx q[0];
rz(-1.3110647) q[0];
sx q[0];
rz(2.3555059) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3439212) q[2];
sx q[2];
rz(-1.6399709) q[2];
sx q[2];
rz(0.56874146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.076482) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-2.4379424) q[1];
rz(-pi) q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65537611) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(2.2359713) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435796) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.4354338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84076607) q[0];
sx q[0];
rz(-0.78221417) q[0];
sx q[0];
rz(1.3602815) q[0];
rz(-pi) q[1];
rz(0.058850364) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(-2.5318052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1286436) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(0.13601555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4581497) q[3];
sx q[3];
rz(-0.23922353) q[3];
sx q[3];
rz(-2.3000172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(2.215109) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(0.57762161) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(2.1320027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59238626) q[0];
sx q[0];
rz(-1.9027332) q[0];
sx q[0];
rz(0.50360002) q[0];
rz(-pi) q[1];
rz(-0.52092123) q[2];
sx q[2];
rz(-0.43653566) q[2];
sx q[2];
rz(-1.6233363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13751444) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(2.532258) q[1];
rz(-0.76241242) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0499095) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(-1.9419149) q[2];
rz(-0.47232929) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(2.902466) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(-0.4447287) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037576588) q[0];
sx q[0];
rz(-0.74066478) q[0];
sx q[0];
rz(-0.91233493) q[0];
x q[1];
rz(2.4773981) q[2];
sx q[2];
rz(-1.0401298) q[2];
sx q[2];
rz(-1.2150089) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0900314) q[1];
sx q[1];
rz(-1.4434442) q[1];
sx q[1];
rz(-0.8831555) q[1];
rz(0.41110699) q[3];
sx q[3];
rz(-1.6931705) q[3];
sx q[3];
rz(-1.3947595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-0.81531173) q[2];
rz(1.0347962) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(2.8544193) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(2.8093991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9438292) q[0];
sx q[0];
rz(-1.257886) q[0];
sx q[0];
rz(0.80425941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5154778) q[2];
sx q[2];
rz(-0.69960591) q[2];
sx q[2];
rz(0.31457065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5933983) q[1];
sx q[1];
rz(-1.3151004) q[1];
sx q[1];
rz(0.42773186) q[1];
x q[2];
rz(-0.62854564) q[3];
sx q[3];
rz(-0.62918951) q[3];
sx q[3];
rz(2.7043846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(-2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-2.8870781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6063639) q[0];
sx q[0];
rz(-1.6441206) q[0];
sx q[0];
rz(-1.3396157) q[0];
x q[1];
rz(0.43912402) q[2];
sx q[2];
rz(-2.6115978) q[2];
sx q[2];
rz(0.32967552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75603112) q[1];
sx q[1];
rz(-2.8061211) q[1];
sx q[1];
rz(-0.24753333) q[1];
rz(-1.47154) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(-1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(1.1709471) q[2];
sx q[2];
rz(-2.5794537) q[2];
sx q[2];
rz(-0.07515547) q[2];
rz(2.9110254) q[3];
sx q[3];
rz(-2.1721526) q[3];
sx q[3];
rz(-0.9823907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
