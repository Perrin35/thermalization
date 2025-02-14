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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(11.081063) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(-0.73837003) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90621072) q[0];
sx q[0];
rz(-2.8538508) q[0];
sx q[0];
rz(-2.5830027) q[0];
x q[1];
rz(2.9887669) q[2];
sx q[2];
rz(-1.9597561) q[2];
sx q[2];
rz(2.9940384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.112632) q[1];
sx q[1];
rz(-0.46507177) q[1];
sx q[1];
rz(-0.3151703) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5220736) q[3];
sx q[3];
rz(-1.3196919) q[3];
sx q[3];
rz(0.34003231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7585313) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(0.4134678) q[2];
rz(2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41459945) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(0.10398908) q[0];
rz(-0.14101401) q[1];
sx q[1];
rz(-0.68599373) q[1];
sx q[1];
rz(2.7242421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433421) q[0];
sx q[0];
rz(-0.62055991) q[0];
sx q[0];
rz(1.0360361) q[0];
rz(-pi) q[1];
rz(1.9345788) q[2];
sx q[2];
rz(-1.5822389) q[2];
sx q[2];
rz(1.6506843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1097393) q[1];
sx q[1];
rz(-1.8735587) q[1];
sx q[1];
rz(-0.043655386) q[1];
rz(-2.1363229) q[3];
sx q[3];
rz(-1.1124753) q[3];
sx q[3];
rz(-1.2456196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61791164) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(2.9376612) q[2];
rz(1.5150874) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5390891) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(1.0771982) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(-2.5881252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837778) q[0];
sx q[0];
rz(-0.73448616) q[0];
sx q[0];
rz(0.015588394) q[0];
rz(-1.7856423) q[2];
sx q[2];
rz(-2.2707377) q[2];
sx q[2];
rz(-2.6756949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3543713) q[1];
sx q[1];
rz(-1.0007035) q[1];
sx q[1];
rz(-1.9891504) q[1];
x q[2];
rz(-2.9514246) q[3];
sx q[3];
rz(-1.7432508) q[3];
sx q[3];
rz(-1.2631388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.074097721) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(-1.925776) q[2];
rz(-1.0007693) q[3];
sx q[3];
rz(-1.3832904) q[3];
sx q[3];
rz(1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85553402) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(1.704294) q[0];
rz(-1.3358491) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(-0.45509532) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2360596) q[0];
sx q[0];
rz(-2.5192671) q[0];
sx q[0];
rz(3.0041191) q[0];
x q[1];
rz(1.9432151) q[2];
sx q[2];
rz(-1.3758059) q[2];
sx q[2];
rz(0.41620987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14070732) q[1];
sx q[1];
rz(-1.8155087) q[1];
sx q[1];
rz(0.34414704) q[1];
rz(-pi) q[2];
rz(2.9330071) q[3];
sx q[3];
rz(-1.6634395) q[3];
sx q[3];
rz(-1.1535597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(-2.4882107) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(-0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6298544) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(-2.9298933) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(-0.45417085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1566166) q[0];
sx q[0];
rz(-1.0318705) q[0];
sx q[0];
rz(2.6205553) q[0];
x q[1];
rz(0.37068292) q[2];
sx q[2];
rz(-2.1739568) q[2];
sx q[2];
rz(-2.2070845) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5502624) q[1];
sx q[1];
rz(-2.2095517) q[1];
sx q[1];
rz(0.62229054) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16726475) q[3];
sx q[3];
rz(-1.8297046) q[3];
sx q[3];
rz(-0.46239275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5357369) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(-1.9262675) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(2.6281443) q[0];
rz(2.2250941) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(-1.7880012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83780629) q[0];
sx q[0];
rz(-2.6897061) q[0];
sx q[0];
rz(-0.19048283) q[0];
rz(-pi) q[1];
rz(2.2960179) q[2];
sx q[2];
rz(-2.6161199) q[2];
sx q[2];
rz(2.2487244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60823764) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(0.4877301) q[1];
rz(0.065033536) q[3];
sx q[3];
rz(-1.4819996) q[3];
sx q[3];
rz(2.9119583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6650271) q[2];
sx q[2];
rz(-0.55321425) q[2];
sx q[2];
rz(3.108007) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(2.3869042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(1.9710185) q[0];
rz(-0.19078828) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(2.4986787) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211484) q[0];
sx q[0];
rz(-1.6729676) q[0];
sx q[0];
rz(-1.8381071) q[0];
rz(2.229081) q[2];
sx q[2];
rz(-0.40900074) q[2];
sx q[2];
rz(0.3106948) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1509779) q[1];
sx q[1];
rz(-2.7508368) q[1];
sx q[1];
rz(-0.67024605) q[1];
rz(-1.4516109) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(1.6329671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(-1.1419123) q[2];
rz(3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(-0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340045) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(-1.3080066) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(0.21211472) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40529272) q[0];
sx q[0];
rz(-2.2515902) q[0];
sx q[0];
rz(-2.731057) q[0];
rz(-2.8372357) q[2];
sx q[2];
rz(-1.3469014) q[2];
sx q[2];
rz(-2.9606113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5221716) q[1];
sx q[1];
rz(-2.4345349) q[1];
sx q[1];
rz(2.2823384) q[1];
x q[2];
rz(-2.2905825) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(1.855576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1954631) q[2];
sx q[2];
rz(-1.3529494) q[2];
sx q[2];
rz(-1.8606261) q[2];
rz(1.7383176) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(0.72428552) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487149) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(-0.88248673) q[0];
rz(2.8485883) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(2.015347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348081) q[0];
sx q[0];
rz(-0.92807584) q[0];
sx q[0];
rz(1.1086784) q[0];
x q[1];
rz(2.5416667) q[2];
sx q[2];
rz(-1.5740105) q[2];
sx q[2];
rz(-2.7300138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6997776) q[1];
sx q[1];
rz(-0.80301563) q[1];
sx q[1];
rz(-0.074549874) q[1];
rz(-pi) q[2];
rz(-0.65139095) q[3];
sx q[3];
rz(-0.49534251) q[3];
sx q[3];
rz(1.9598531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61665159) q[2];
sx q[2];
rz(-2.5280819) q[2];
sx q[2];
rz(2.9743527) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.1789221) q[3];
sx q[3];
rz(0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328218) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(-2.3311145) q[0];
rz(-1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(0.097361758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7435849) q[0];
sx q[0];
rz(-0.79774374) q[0];
sx q[0];
rz(1.3378851) q[0];
x q[1];
rz(3.112744) q[2];
sx q[2];
rz(-1.1897773) q[2];
sx q[2];
rz(2.888607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3835427) q[1];
sx q[1];
rz(-1.6191747) q[1];
sx q[1];
rz(0.58396062) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9769659) q[3];
sx q[3];
rz(-2.8517234) q[3];
sx q[3];
rz(-2.4932662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8436766) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(1.8664912) q[2];
rz(0.44025931) q[3];
sx q[3];
rz(-2.2831254) q[3];
sx q[3];
rz(0.27537235) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-0.62200017) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(1.4595819) q[2];
sx q[2];
rz(-1.6485452) q[2];
sx q[2];
rz(-2.9948276) q[2];
rz(0.6435271) q[3];
sx q[3];
rz(-1.4326) q[3];
sx q[3];
rz(2.8958733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
