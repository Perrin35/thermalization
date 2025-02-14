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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(-0.771653) q[0];
rz(2.9906988) q[1];
sx q[1];
rz(3.4540662) q[1];
sx q[1];
rz(11.087505) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85017289) q[0];
sx q[0];
rz(-1.2787191) q[0];
sx q[0];
rz(0.34698457) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0059228) q[2];
sx q[2];
rz(-0.56081334) q[2];
sx q[2];
rz(0.31080526) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7141908) q[1];
sx q[1];
rz(-1.8334249) q[1];
sx q[1];
rz(2.3841969) q[1];
rz(-pi) q[2];
rz(1.3948049) q[3];
sx q[3];
rz(-2.2424091) q[3];
sx q[3];
rz(2.1064389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(-0.83211952) q[2];
rz(-0.65699792) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(-2.7895797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326139) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(0.12282148) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(0.67272433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8685659) q[0];
sx q[0];
rz(-1.582524) q[0];
sx q[0];
rz(-1.5290029) q[0];
rz(-1.2956728) q[2];
sx q[2];
rz(-1.2689991) q[2];
sx q[2];
rz(1.2040621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0716463) q[1];
sx q[1];
rz(-0.59161797) q[1];
sx q[1];
rz(-1.9272734) q[1];
x q[2];
rz(-3.1150713) q[3];
sx q[3];
rz(-2.0904198) q[3];
sx q[3];
rz(-2.5260203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0279205) q[2];
sx q[2];
rz(-1.5936667) q[2];
sx q[2];
rz(0.38635722) q[2];
rz(-1.1682642) q[3];
sx q[3];
rz(-1.9100274) q[3];
sx q[3];
rz(-1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2439709) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(-3.0834055) q[0];
rz(-0.30481401) q[1];
sx q[1];
rz(-1.1331646) q[1];
sx q[1];
rz(2.286639) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8940272) q[0];
sx q[0];
rz(-1.6391048) q[0];
sx q[0];
rz(1.5724284) q[0];
rz(2.0570175) q[2];
sx q[2];
rz(-0.9409608) q[2];
sx q[2];
rz(-1.6597009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.952576) q[1];
sx q[1];
rz(-1.7754648) q[1];
sx q[1];
rz(2.3701028) q[1];
x q[2];
rz(-1.9819354) q[3];
sx q[3];
rz(-0.59468958) q[3];
sx q[3];
rz(-2.2583972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0893687) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(-1.8734056) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190401) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(-2.4265491) q[0];
rz(-2.7252281) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(0.29410902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.607632) q[0];
sx q[0];
rz(-2.2658969) q[0];
sx q[0];
rz(0.075420336) q[0];
x q[1];
rz(1.8744743) q[2];
sx q[2];
rz(-2.9275002) q[2];
sx q[2];
rz(-1.0372727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8079559) q[1];
sx q[1];
rz(-2.9142434) q[1];
sx q[1];
rz(2.5264986) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77463051) q[3];
sx q[3];
rz(-2.3916158) q[3];
sx q[3];
rz(-2.3924654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2405582) q[2];
sx q[2];
rz(-1.7419387) q[2];
sx q[2];
rz(-0.72875363) q[2];
rz(-0.48480222) q[3];
sx q[3];
rz(-1.0032283) q[3];
sx q[3];
rz(0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6840376) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(-2.8628602) q[0];
rz(-3.0740671) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(2.6356437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7161432) q[0];
sx q[0];
rz(-1.6438515) q[0];
sx q[0];
rz(-2.5849033) q[0];
x q[1];
rz(-1.4691817) q[2];
sx q[2];
rz(-1.3739462) q[2];
sx q[2];
rz(-2.5972899) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5927148) q[1];
sx q[1];
rz(-1.028182) q[1];
sx q[1];
rz(0.9701258) q[1];
x q[2];
rz(-0.26624532) q[3];
sx q[3];
rz(-2.2266386) q[3];
sx q[3];
rz(-0.4474934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2013596) q[2];
sx q[2];
rz(-1.4505922) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(-2.4612259) q[0];
rz(0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(-2.5804677) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8770804) q[0];
sx q[0];
rz(-1.7753557) q[0];
sx q[0];
rz(1.9407985) q[0];
rz(0.20724736) q[2];
sx q[2];
rz(-2.2842513) q[2];
sx q[2];
rz(2.0240473) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5270414) q[1];
sx q[1];
rz(-2.8379662) q[1];
sx q[1];
rz(0.18358453) q[1];
x q[2];
rz(0.64021982) q[3];
sx q[3];
rz(-1.9915765) q[3];
sx q[3];
rz(0.32395485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23201021) q[2];
sx q[2];
rz(-2.1640615) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(1.0163418) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(1.1004755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119754) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(0.49402657) q[0];
rz(-1.5771075) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(-0.090242537) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2170885) q[0];
sx q[0];
rz(-2.6441433) q[0];
sx q[0];
rz(0.50636943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87553067) q[2];
sx q[2];
rz(-1.1644496) q[2];
sx q[2];
rz(-1.6900495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.68914946) q[1];
sx q[1];
rz(-1.9724558) q[1];
sx q[1];
rz(-0.31280495) q[1];
rz(-pi) q[2];
rz(-0.97848864) q[3];
sx q[3];
rz(-0.33708015) q[3];
sx q[3];
rz(-0.85142985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6287441) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(0.64344704) q[2];
rz(2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(1.2233618) q[0];
rz(-2.2471097) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(-1.1462513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2558738) q[0];
sx q[0];
rz(-2.4303655) q[0];
sx q[0];
rz(-2.2125986) q[0];
rz(-0.84051657) q[2];
sx q[2];
rz(-2.5926551) q[2];
sx q[2];
rz(-0.38503371) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8967197) q[1];
sx q[1];
rz(-1.280652) q[1];
sx q[1];
rz(1.3194607) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1560539) q[3];
sx q[3];
rz(-1.7904591) q[3];
sx q[3];
rz(2.1287763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2724096) q[2];
sx q[2];
rz(-2.2792008) q[2];
sx q[2];
rz(1.4746846) q[2];
rz(-0.21371755) q[3];
sx q[3];
rz(-1.7361879) q[3];
sx q[3];
rz(-0.86764446) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992582) q[0];
sx q[0];
rz(-2.5304351) q[0];
sx q[0];
rz(-2.4364731) q[0];
rz(-0.57535386) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(-0.56952482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0234719) q[0];
sx q[0];
rz(-0.11205738) q[0];
sx q[0];
rz(-0.81046354) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8649544) q[2];
sx q[2];
rz(-2.0450971) q[2];
sx q[2];
rz(2.2612342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77806738) q[1];
sx q[1];
rz(-2.229564) q[1];
sx q[1];
rz(-0.79658618) q[1];
x q[2];
rz(2.6526057) q[3];
sx q[3];
rz(-1.8400116) q[3];
sx q[3];
rz(-2.0662291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9465955) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(2.1627964) q[2];
rz(2.7106674) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(2.0144958) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2836714) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(1.6790144) q[0];
rz(2.1025533) q[1];
sx q[1];
rz(-1.7580527) q[1];
sx q[1];
rz(1.9336112) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53762728) q[0];
sx q[0];
rz(-0.4930636) q[0];
sx q[0];
rz(0.8925256) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1080148) q[2];
sx q[2];
rz(-1.0462282) q[2];
sx q[2];
rz(1.7941504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0145899) q[1];
sx q[1];
rz(-2.5095438) q[1];
sx q[1];
rz(1.3227425) q[1];
rz(1.1753756) q[3];
sx q[3];
rz(-1.5081769) q[3];
sx q[3];
rz(0.6400125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98564467) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.3204302) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(2.5613274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(-1.6571922) q[1];
sx q[1];
rz(-1.3927554) q[1];
sx q[1];
rz(1.4295255) q[1];
rz(1.4996573) q[2];
sx q[2];
rz(-1.1924469) q[2];
sx q[2];
rz(3.057657) q[2];
rz(2.8410925) q[3];
sx q[3];
rz(-2.2347652) q[3];
sx q[3];
rz(1.4303257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
