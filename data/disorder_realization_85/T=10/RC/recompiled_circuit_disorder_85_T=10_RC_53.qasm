OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(2.5420904) q[1];
sx q[1];
rz(11.190344) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(-2.2493275) q[0];
rz(-pi) q[1];
rz(-2.0618093) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(0.71066463) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7264759) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(-0.16023689) q[1];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(-1.1799312) q[0];
rz(-2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(0.72431272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(-2.2749167) q[0];
x q[1];
rz(-2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(-0.85180887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38072941) q[1];
sx q[1];
rz(-2.7298096) q[1];
sx q[1];
rz(-2.1058583) q[1];
rz(-pi) q[2];
rz(2.4507387) q[3];
sx q[3];
rz(-2.8375531) q[3];
sx q[3];
rz(3.0512878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(0.36188564) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1299382) q[0];
sx q[0];
rz(-1.9217976) q[0];
sx q[0];
rz(2.8610693) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8457882) q[2];
sx q[2];
rz(-2.25053) q[2];
sx q[2];
rz(2.5909397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0137274) q[1];
sx q[1];
rz(-1.8121108) q[1];
sx q[1];
rz(1.024854) q[1];
rz(-pi) q[2];
rz(-3.1261256) q[3];
sx q[3];
rz(-1.8387715) q[3];
sx q[3];
rz(-2.1249352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(-0.82675654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137705) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(-0.69112372) q[0];
x q[1];
rz(0.76765676) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.3354288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(-2.4946458) q[1];
rz(-pi) q[2];
rz(-1.2590253) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(-1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(3.0692696) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(-2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(2.856423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.973424) q[0];
sx q[0];
rz(-2.0452721) q[0];
sx q[0];
rz(-0.04767496) q[0];
rz(-1.9589013) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-0.91458048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4080216) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(-1.8838521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(-0.13084403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.8834615) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.3670115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-2.0220387) q[0];
rz(-2.5622257) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(-0.52057779) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.651598) q[1];
sx q[1];
rz(1.8606436) q[1];
rz(-pi) q[2];
rz(-1.4671765) q[3];
sx q[3];
rz(-1.2195671) q[3];
sx q[3];
rz(0.87408376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0872333) q[0];
sx q[0];
rz(-1.5477675) q[0];
sx q[0];
rz(1.5646294) q[0];
rz(-pi) q[1];
x q[1];
rz(1.621775) q[2];
sx q[2];
rz(-1.4547252) q[2];
sx q[2];
rz(-0.34250868) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2611744) q[1];
sx q[1];
rz(-2.2194127) q[1];
sx q[1];
rz(2.3748114) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8740158) q[3];
sx q[3];
rz(-1.056864) q[3];
sx q[3];
rz(-0.60790387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(0.47510535) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(2.8444667) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32983366) q[0];
sx q[0];
rz(-2.9114897) q[0];
sx q[0];
rz(1.1128725) q[0];
rz(-pi) q[1];
rz(-1.7392776) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(-0.70665765) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4375953) q[1];
sx q[1];
rz(-0.62347177) q[1];
sx q[1];
rz(0.61203476) q[1];
rz(2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(-2.3074647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(-2.6240255) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(0.10841766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71589564) q[0];
sx q[0];
rz(-1.513039) q[0];
sx q[0];
rz(-0.9066559) q[0];
rz(-1.7494781) q[2];
sx q[2];
rz(-1.6166501) q[2];
sx q[2];
rz(-0.97464857) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0260967) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(2.2309169) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7897723) q[3];
sx q[3];
rz(-0.7162381) q[3];
sx q[3];
rz(-2.4718474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1291528) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(0.33995315) q[2];
rz(-2.7231976) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(-2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(0.055158786) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.508651) q[0];
sx q[0];
rz(-1.635449) q[0];
sx q[0];
rz(2.7642194) q[0];
rz(-pi) q[1];
rz(-2.5818985) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(2.6924804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2436279) q[1];
sx q[1];
rz(-0.98476714) q[1];
sx q[1];
rz(0.90378739) q[1];
rz(-pi) q[2];
rz(-1.1147898) q[3];
sx q[3];
rz(-1.6037233) q[3];
sx q[3];
rz(-1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(-2.6514163) q[2];
rz(0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(2.8293777) q[2];
sx q[2];
rz(-1.6630465) q[2];
sx q[2];
rz(1.9065471) q[2];
rz(0.24003868) q[3];
sx q[3];
rz(-1.5272899) q[3];
sx q[3];
rz(1.8361113) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
