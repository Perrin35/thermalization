OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1546254) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(-0.89303645) q[0];
x q[1];
rz(2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-2.0703966) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33949172) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(0.93356737) q[1];
x q[2];
rz(-1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(0.030348226) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.617584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.671316) q[0];
sx q[0];
rz(-2.94194) q[0];
sx q[0];
rz(1.5623564) q[0];
rz(-pi) q[1];
rz(0.020521684) q[2];
sx q[2];
rz(-2.0251209) q[2];
sx q[2];
rz(3.0128535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6807032) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(0.99980385) q[1];
x q[2];
rz(-1.7582943) q[3];
sx q[3];
rz(-1.4783995) q[3];
sx q[3];
rz(-2.932991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7559738) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(-0.31578627) q[0];
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
rz(-2.2825851) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(-2.1501599) q[0];
rz(-pi) q[1];
rz(-0.88148586) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(1.880868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6985059) q[1];
sx q[1];
rz(-1.3830739) q[1];
sx q[1];
rz(-0.94633533) q[1];
rz(-pi) q[2];
rz(-2.7470845) q[3];
sx q[3];
rz(-1.2393701) q[3];
sx q[3];
rz(3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(2.4345051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0079460572) q[0];
sx q[0];
rz(-0.61230731) q[0];
sx q[0];
rz(2.4164819) q[0];
rz(1.445304) q[2];
sx q[2];
rz(-1.1271994) q[2];
sx q[2];
rz(-1.0759629) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9634339) q[1];
sx q[1];
rz(-1.6903094) q[1];
sx q[1];
rz(-0.91576373) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1726923) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(-0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3209155) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(-1.0955411) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(-0.25462338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4655612) q[0];
sx q[0];
rz(-1.4916294) q[0];
sx q[0];
rz(-0.013750793) q[0];
x q[1];
rz(-0.014572797) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(-2.2968959) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4659781) q[1];
sx q[1];
rz(-1.8120159) q[1];
sx q[1];
rz(-2.4720008) q[1];
rz(-pi) q[2];
rz(2.280064) q[3];
sx q[3];
rz(-1.9083175) q[3];
sx q[3];
rz(2.5634114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(-1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(-2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(2.9350231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20139192) q[0];
sx q[0];
rz(-0.79027806) q[0];
sx q[0];
rz(1.2339562) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1913387) q[2];
sx q[2];
rz(-2.5197919) q[2];
sx q[2];
rz(1.7813462) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8787074) q[1];
sx q[1];
rz(-2.4135114) q[1];
sx q[1];
rz(1.8379777) q[1];
rz(-1.1849095) q[3];
sx q[3];
rz(-2.4486802) q[3];
sx q[3];
rz(-3.0544359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(2.858813) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(1.8136224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049411557) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(2.0121502) q[0];
rz(-pi) q[1];
rz(-0.46220025) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(2.5525023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59203397) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(1.3105427) q[1];
rz(-0.51289576) q[3];
sx q[3];
rz(-1.9478056) q[3];
sx q[3];
rz(0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(0.95091933) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2598341) q[0];
sx q[0];
rz(-1.330733) q[0];
sx q[0];
rz(1.865922) q[0];
x q[1];
rz(-3.0792564) q[2];
sx q[2];
rz(-1.3051635) q[2];
sx q[2];
rz(3.0692284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9329405) q[1];
sx q[1];
rz(-1.7793852) q[1];
sx q[1];
rz(1.5044466) q[1];
rz(-pi) q[2];
rz(2.7637134) q[3];
sx q[3];
rz(-1.9482908) q[3];
sx q[3];
rz(2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(-0.075909464) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-2.4818647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2840246) q[0];
sx q[0];
rz(-1.0463456) q[0];
sx q[0];
rz(0.88733034) q[0];
x q[1];
rz(2.9218036) q[2];
sx q[2];
rz(-2.347749) q[2];
sx q[2];
rz(-1.1500051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1253423) q[1];
sx q[1];
rz(-2.2012735) q[1];
sx q[1];
rz(-0.36228212) q[1];
rz(1.6660059) q[3];
sx q[3];
rz(-0.54064893) q[3];
sx q[3];
rz(1.3117865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(0.025645105) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(2.1059039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64338387) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(-1.2163175) q[0];
x q[1];
rz(0.56634855) q[2];
sx q[2];
rz(-2.1841335) q[2];
sx q[2];
rz(1.4373506) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90696883) q[1];
sx q[1];
rz(-1.8112438) q[1];
sx q[1];
rz(0.30827) q[1];
rz(-pi) q[2];
rz(-0.74929897) q[3];
sx q[3];
rz(-1.7710847) q[3];
sx q[3];
rz(1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62853652) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.1516085) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7619027) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-1.6133616) q[2];
sx q[2];
rz(-1.9812937) q[2];
sx q[2];
rz(-1.5716256) q[2];
rz(3.1291943) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
