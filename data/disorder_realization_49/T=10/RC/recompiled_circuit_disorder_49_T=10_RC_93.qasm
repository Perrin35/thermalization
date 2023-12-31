OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5306659) q[0];
sx q[0];
rz(-2.2025684) q[0];
sx q[0];
rz(0.0052069081) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(-1.985328) q[1];
sx q[1];
rz(-1.1896689) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1253163) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(-0.54434158) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9615133) q[2];
sx q[2];
rz(-1.7109949) q[2];
sx q[2];
rz(2.7498498) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7835044) q[1];
sx q[1];
rz(-1.7184098) q[1];
sx q[1];
rz(-1.258177) q[1];
rz(-0.94250836) q[3];
sx q[3];
rz(-1.9058459) q[3];
sx q[3];
rz(-0.059700746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71620119) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(-0.5973967) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.3436915) q[3];
sx q[3];
rz(-1.8610154) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(-0.33357099) q[0];
rz(-1.0936273) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-0.11322583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118606) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(3.0990764) q[0];
x q[1];
rz(3.0683238) q[2];
sx q[2];
rz(-1.5904625) q[2];
sx q[2];
rz(0.71436963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7624224) q[1];
sx q[1];
rz(-1.9813073) q[1];
sx q[1];
rz(-1.1344086) q[1];
x q[2];
rz(2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(-1.295134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.018628) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-0.22228995) q[2];
rz(0.22953454) q[3];
sx q[3];
rz(-2.4042606) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.570913) q[0];
sx q[0];
rz(-0.88071841) q[0];
rz(-2.0942028) q[1];
sx q[1];
rz(-1.2068345) q[1];
sx q[1];
rz(-3.0139794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-2.2420609) q[0];
sx q[0];
rz(0.22247252) q[0];
rz(-pi) q[1];
rz(0.80673809) q[2];
sx q[2];
rz(-0.85959896) q[2];
sx q[2];
rz(-0.15653175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0652005) q[1];
sx q[1];
rz(-0.25307357) q[1];
sx q[1];
rz(1.2364926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25281275) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(-1.52724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7814653) q[2];
sx q[2];
rz(-2.0121622) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(2.7919853) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.67868245) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(2.7691832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2686553) q[0];
sx q[0];
rz(-2.4252709) q[0];
sx q[0];
rz(-1.3852081) q[0];
rz(-pi) q[1];
rz(1.0834951) q[2];
sx q[2];
rz(-0.63939017) q[2];
sx q[2];
rz(-1.0207748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1540893) q[1];
sx q[1];
rz(-1.1945063) q[1];
sx q[1];
rz(0.16102468) q[1];
rz(0.72752556) q[3];
sx q[3];
rz(-1.1480867) q[3];
sx q[3];
rz(-2.5142575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-0.70251846) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(-0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(0.4549543) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882778) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-0.89748853) q[1];
sx q[1];
rz(-0.67684832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880496) q[0];
sx q[0];
rz(-1.6497668) q[0];
sx q[0];
rz(2.700564) q[0];
rz(-pi) q[1];
rz(1.3555688) q[2];
sx q[2];
rz(-2.2903246) q[2];
sx q[2];
rz(-2.1446251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.10830282) q[1];
sx q[1];
rz(-1.8255594) q[1];
sx q[1];
rz(1.5441896) q[1];
x q[2];
rz(2.7480514) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(-1.9963025) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695078) q[0];
sx q[0];
rz(-2.5818995) q[0];
sx q[0];
rz(-2.1493388) q[0];
rz(-0.51322333) q[2];
sx q[2];
rz(-0.66868082) q[2];
sx q[2];
rz(-0.19323397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7416523) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(-1.9386577) q[1];
x q[2];
rz(-0.17721456) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(-0.23055102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(-0.28665001) q[2];
rz(-0.24946985) q[3];
sx q[3];
rz(-1.9727861) q[3];
sx q[3];
rz(-0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(-1.4666784) q[0];
rz(1.02007) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(-1.0010304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0953656) q[0];
sx q[0];
rz(-1.8512748) q[0];
sx q[0];
rz(3.1026955) q[0];
rz(-pi) q[1];
rz(-1.2008576) q[2];
sx q[2];
rz(-1.8897893) q[2];
sx q[2];
rz(-1.621304) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66623464) q[1];
sx q[1];
rz(-1.959414) q[1];
sx q[1];
rz(-2.9734441) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27490297) q[3];
sx q[3];
rz(-0.66920815) q[3];
sx q[3];
rz(1.6717403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(3.0380847) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(2.3039968) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(0.57530595) q[0];
rz(1.9837603) q[2];
sx q[2];
rz(-2.0418842) q[2];
sx q[2];
rz(-2.4177891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0583565) q[1];
sx q[1];
rz(-0.55587686) q[1];
sx q[1];
rz(-0.91169375) q[1];
rz(-0.89160664) q[3];
sx q[3];
rz(-1.3636175) q[3];
sx q[3];
rz(0.87219119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-0.71425444) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927239) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(2.2946987) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(0.27639595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4682639) q[0];
sx q[0];
rz(-1.4110761) q[0];
sx q[0];
rz(0.41202742) q[0];
rz(-pi) q[1];
rz(2.9133965) q[2];
sx q[2];
rz(-1.3820634) q[2];
sx q[2];
rz(-2.2224094) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.285187) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(-1.7941979) q[1];
rz(-1.8180088) q[3];
sx q[3];
rz(-0.97655481) q[3];
sx q[3];
rz(-1.1288527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.85236621) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(3.0140871) q[2];
rz(-0.036711983) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7413095) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(2.1113077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0477662) q[0];
sx q[0];
rz(-2.2802135) q[0];
sx q[0];
rz(1.8295656) q[0];
rz(-pi) q[1];
rz(-1.2376386) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(-1.9106939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-1.1332129) q[1];
sx q[1];
rz(1.9220819) q[1];
rz(0.45566166) q[3];
sx q[3];
rz(-2.4767498) q[3];
sx q[3];
rz(-1.957422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0673922) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(0.34118787) q[2];
rz(0.7406922) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-3.0497131) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0108903) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.40217051) q[2];
sx q[2];
rz(-0.46605863) q[2];
sx q[2];
rz(-3.0083187) q[2];
rz(2.1445198) q[3];
sx q[3];
rz(-0.96521796) q[3];
sx q[3];
rz(-0.60346606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
