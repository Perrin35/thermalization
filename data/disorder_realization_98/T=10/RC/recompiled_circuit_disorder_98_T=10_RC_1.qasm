OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7604231) q[0];
sx q[0];
rz(-1.086735) q[0];
sx q[0];
rz(0.83597393) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3157273) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(0.65537383) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.522361) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(-0.038832263) q[1];
rz(-1.876086) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-0.65223637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61440496) q[0];
sx q[0];
rz(-1.5584649) q[0];
sx q[0];
rz(-3.1129818) q[0];
rz(-pi) q[1];
rz(-2.9035283) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(0.1711947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0374239) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(1.3377405) q[1];
x q[2];
rz(0.60450508) q[3];
sx q[3];
rz(-2.1406056) q[3];
sx q[3];
rz(-1.9250211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91857598) q[0];
sx q[0];
rz(-2.4239459) q[0];
sx q[0];
rz(2.7735604) q[0];
rz(-pi) q[1];
rz(1.6367958) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(-2.6413692) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50074358) q[1];
sx q[1];
rz(-1.3355681) q[1];
sx q[1];
rz(-2.6231223) q[1];
x q[2];
rz(-3.1317741) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-0.90399495) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7061493) q[0];
sx q[0];
rz(-0.94067803) q[0];
sx q[0];
rz(-2.8773727) q[0];
rz(1.8965917) q[2];
sx q[2];
rz(-2.1278283) q[2];
sx q[2];
rz(-1.7011124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1069146) q[1];
sx q[1];
rz(-2.2639096) q[1];
sx q[1];
rz(-0.17130674) q[1];
x q[2];
rz(2.9070204) q[3];
sx q[3];
rz(-1.4751225) q[3];
sx q[3];
rz(2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2146384) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7974632) q[0];
sx q[0];
rz(-1.9088233) q[0];
sx q[0];
rz(0.56030886) q[0];
rz(-pi) q[1];
rz(2.064803) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(-2.6780724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.036451642) q[1];
sx q[1];
rz(-0.76117939) q[1];
sx q[1];
rz(2.9245604) q[1];
rz(-pi) q[2];
rz(1.8022728) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(-1.0486697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.9979427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8457984) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(0.11008115) q[0];
rz(-pi) q[1];
rz(-2.4684858) q[2];
sx q[2];
rz(-1.3654725) q[2];
sx q[2];
rz(-2.7538607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1144048) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(-1.1272217) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2915217) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(-0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.2333262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(1.5412488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2912824) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(0.26956272) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26112939) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(0.1111828) q[1];
rz(-pi) q[2];
rz(-0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(-2.1255778) q[2];
rz(-0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381666) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(-1.1358791) q[0];
x q[1];
rz(-0.7779185) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.9922436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.313169) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(-1.4491175) q[1];
x q[2];
rz(2.2860252) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(0.95203979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68226472) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(3.1138611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-2.3150139) q[0];
x q[1];
rz(-2.8743923) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(0.53182488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4122388) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(-0.17462294) q[1];
rz(-3.1157007) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(-2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-0.85723248) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(2.8840816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.9985984) q[0];
rz(1.0742513) q[2];
sx q[2];
rz(-2.0601344) q[2];
sx q[2];
rz(-2.731583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.820206) q[1];
sx q[1];
rz(-0.69637978) q[1];
sx q[1];
rz(3.1192944) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24210838) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-2.4702934) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(-1.0971309) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
