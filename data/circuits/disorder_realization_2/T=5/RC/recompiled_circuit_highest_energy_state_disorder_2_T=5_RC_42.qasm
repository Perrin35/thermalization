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
rz(-1.4827363) q[0];
sx q[0];
rz(-2.1602477) q[0];
sx q[0];
rz(2.0438097) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(-0.0069590574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46411447) q[0];
sx q[0];
rz(-0.65271806) q[0];
sx q[0];
rz(1.663289) q[0];
x q[1];
rz(0.9240146) q[2];
sx q[2];
rz(-1.3818463) q[2];
sx q[2];
rz(1.3786157) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0557441) q[1];
sx q[1];
rz(-1.6043449) q[1];
sx q[1];
rz(-1.335024) q[1];
x q[2];
rz(0.29421774) q[3];
sx q[3];
rz(-0.3539043) q[3];
sx q[3];
rz(0.82415165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2385345) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(-1.6713589) q[2];
rz(-0.12601958) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(2.457705) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11494342) q[0];
sx q[0];
rz(-0.4489972) q[0];
sx q[0];
rz(-0.63585109) q[0];
rz(0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(-1.4453452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46585402) q[0];
sx q[0];
rz(-2.2096118) q[0];
sx q[0];
rz(2.5771228) q[0];
x q[1];
rz(-2.0436297) q[2];
sx q[2];
rz(-0.84542984) q[2];
sx q[2];
rz(-0.54976094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0750097) q[1];
sx q[1];
rz(-2.420418) q[1];
sx q[1];
rz(-2.9020792) q[1];
rz(-pi) q[2];
rz(3.0731455) q[3];
sx q[3];
rz(-1.7020149) q[3];
sx q[3];
rz(2.69773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14629743) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(-0.13776097) q[2];
rz(-2.6855101) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(-0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.544203) q[0];
sx q[0];
rz(-1.5578288) q[0];
sx q[0];
rz(-2.5624045) q[0];
rz(0.10313615) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(-0.78027049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4057377) q[0];
sx q[0];
rz(-3.0445703) q[0];
sx q[0];
rz(2.9461622) q[0];
rz(0.68124007) q[2];
sx q[2];
rz(-1.4461755) q[2];
sx q[2];
rz(-1.7523426) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6578778) q[1];
sx q[1];
rz(-0.45770744) q[1];
sx q[1];
rz(-1.0545516) q[1];
rz(-1.3079726) q[3];
sx q[3];
rz(-2.4287831) q[3];
sx q[3];
rz(1.3808911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4812193) q[2];
sx q[2];
rz(-1.6796835) q[2];
sx q[2];
rz(-0.090864651) q[2];
rz(-1.6615435) q[3];
sx q[3];
rz(-1.3250947) q[3];
sx q[3];
rz(0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12260967) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(0.51938272) q[0];
rz(-1.9620365) q[1];
sx q[1];
rz(-2.4603381) q[1];
sx q[1];
rz(-0.048351668) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8285168) q[0];
sx q[0];
rz(-2.3017602) q[0];
sx q[0];
rz(2.2832995) q[0];
x q[1];
rz(0.3860571) q[2];
sx q[2];
rz(-2.0252844) q[2];
sx q[2];
rz(-0.41815652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1241403) q[1];
sx q[1];
rz(-2.0582948) q[1];
sx q[1];
rz(1.5392787) q[1];
rz(-pi) q[2];
rz(2.8818733) q[3];
sx q[3];
rz(-0.81213299) q[3];
sx q[3];
rz(-2.557723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68253303) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(1.691386) q[2];
rz(-1.1059443) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76097101) q[0];
sx q[0];
rz(-0.80046099) q[0];
sx q[0];
rz(-2.3642484) q[0];
rz(-2.6457973) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(2.9487603) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7161761) q[0];
sx q[0];
rz(-1.2649853) q[0];
sx q[0];
rz(2.5772463) q[0];
rz(-pi) q[1];
rz(-0.76754153) q[2];
sx q[2];
rz(-1.7450252) q[2];
sx q[2];
rz(-1.7310639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4905688) q[1];
sx q[1];
rz(-1.2747163) q[1];
sx q[1];
rz(0.091551642) q[1];
rz(-pi) q[2];
rz(-1.939658) q[3];
sx q[3];
rz(-2.7141124) q[3];
sx q[3];
rz(-2.1457246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(-2.3053816) q[2];
rz(0.95544514) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(2.6098765) q[3];
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
rz(-0.29109508) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(0.033705458) q[0];
rz(2.2233502) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(0.69923002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13051899) q[0];
sx q[0];
rz(-2.1222181) q[0];
sx q[0];
rz(0.60451492) q[0];
rz(1.5169473) q[2];
sx q[2];
rz(-1.757818) q[2];
sx q[2];
rz(2.5973926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3311309) q[1];
sx q[1];
rz(-1.7469566) q[1];
sx q[1];
rz(-2.6947652) q[1];
x q[2];
rz(2.8630303) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(2.183941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0596727) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(-0.89422798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28602257) q[0];
sx q[0];
rz(-2.6595071) q[0];
sx q[0];
rz(1.9739738) q[0];
rz(2.7883912) q[1];
sx q[1];
rz(-2.628852) q[1];
sx q[1];
rz(0.28164992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5198811) q[0];
sx q[0];
rz(-2.6586652) q[0];
sx q[0];
rz(1.3819225) q[0];
rz(-pi) q[1];
rz(3.0007576) q[2];
sx q[2];
rz(-2.4895034) q[2];
sx q[2];
rz(1.2229133) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0593917) q[1];
sx q[1];
rz(-1.7052461) q[1];
sx q[1];
rz(1.4217292) q[1];
x q[2];
rz(2.6572833) q[3];
sx q[3];
rz(-2.0748027) q[3];
sx q[3];
rz(-0.33226099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1436651) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(1.3227051) q[2];
rz(-1.6392684) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(-3.0433906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1439576) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(2.8676046) q[0];
rz(0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(0.53057539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218567) q[0];
sx q[0];
rz(-1.6101735) q[0];
sx q[0];
rz(-1.5669797) q[0];
rz(-pi) q[1];
rz(-3.0184348) q[2];
sx q[2];
rz(-0.55610699) q[2];
sx q[2];
rz(-1.3431988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0116327) q[1];
sx q[1];
rz(-1.7372565) q[1];
sx q[1];
rz(-2.7618206) q[1];
rz(1.382393) q[3];
sx q[3];
rz(-2.419988) q[3];
sx q[3];
rz(2.5327275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8990367) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(0.25699082) q[2];
rz(1.1993923) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.5523819) q[0];
sx q[0];
rz(-2.1871545) q[0];
sx q[0];
rz(-0.5300262) q[0];
rz(1.5026622) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(-2.2030305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8535271) q[0];
sx q[0];
rz(-0.99366436) q[0];
sx q[0];
rz(2.1398628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6907127) q[2];
sx q[2];
rz(-1.7133822) q[2];
sx q[2];
rz(0.8220807) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5512498) q[1];
sx q[1];
rz(-1.9971202) q[1];
sx q[1];
rz(2.9077282) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1562002) q[3];
sx q[3];
rz(-1.6847576) q[3];
sx q[3];
rz(-2.0843907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3878801) q[2];
sx q[2];
rz(-0.91999274) q[2];
sx q[2];
rz(2.06125) q[2];
rz(-0.8148109) q[3];
sx q[3];
rz(-1.9459414) q[3];
sx q[3];
rz(-1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.998488) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(-2.0527573) q[0];
rz(3.0740652) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(2.3823104) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6367594) q[0];
sx q[0];
rz(-2.8593316) q[0];
sx q[0];
rz(-2.0105848) q[0];
rz(0.49726059) q[2];
sx q[2];
rz(-0.25575519) q[2];
sx q[2];
rz(-1.349468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85222178) q[1];
sx q[1];
rz(-1.3557649) q[1];
sx q[1];
rz(2.380886) q[1];
x q[2];
rz(-0.57709007) q[3];
sx q[3];
rz(-2.7365757) q[3];
sx q[3];
rz(1.7187723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33409432) q[2];
sx q[2];
rz(-2.0457902) q[2];
sx q[2];
rz(2.5661772) q[2];
rz(-2.5790162) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0811049) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(2.0545215) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(-1.7123863) q[2];
sx q[2];
rz(-1.7238486) q[2];
sx q[2];
rz(0.42949745) q[2];
rz(-2.6230276) q[3];
sx q[3];
rz(-2.4760078) q[3];
sx q[3];
rz(2.7005394) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
