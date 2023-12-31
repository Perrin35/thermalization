OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.1975343) q[0];
rz(3.0543047) q[2];
sx q[2];
rz(-0.44863551) q[2];
sx q[2];
rz(-1.0686312) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.672294) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(-1.6507571) q[1];
rz(-pi) q[2];
rz(0.2557405) q[3];
sx q[3];
rz(-1.2107953) q[3];
sx q[3];
rz(-0.4482625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6204651) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(2.492766) q[0];
rz(-pi) q[1];
rz(-2.8380727) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(0.34740651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38824575) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(1.5335598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.84045029) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-0.99951807) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32854983) q[0];
sx q[0];
rz(-1.2721491) q[0];
sx q[0];
rz(-0.15655984) q[0];
x q[1];
rz(-2.2314084) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(-0.78188932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(1.2815777) q[1];
x q[2];
rz(0.079654982) q[3];
sx q[3];
rz(-2.0783391) q[3];
sx q[3];
rz(-1.5505276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(2.8667563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4764458) q[0];
sx q[0];
rz(-0.11867141) q[0];
sx q[0];
rz(-2.2975886) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5722256) q[2];
sx q[2];
rz(-1.187385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3596613) q[1];
sx q[1];
rz(-2.0205824) q[1];
sx q[1];
rz(-0.91314544) q[1];
rz(1.4196017) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(-2.7694626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(-2.143798) q[0];
rz(-0.18355852) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.957513) q[0];
sx q[0];
rz(-0.74680579) q[0];
sx q[0];
rz(2.1225131) q[0];
x q[1];
rz(-1.1797656) q[2];
sx q[2];
rz(-2.070825) q[2];
sx q[2];
rz(-0.21991877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5114054) q[1];
sx q[1];
rz(-1.9747707) q[1];
sx q[1];
rz(-1.9460815) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3715641) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(-1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-2.8009159) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6402123) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(1.120938) q[0];
rz(-0.75603007) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(2.0331969) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.016547116) q[1];
sx q[1];
rz(-2.6066337) q[1];
sx q[1];
rz(2.6782481) q[1];
x q[2];
rz(-1.0814704) q[3];
sx q[3];
rz(-1.193207) q[3];
sx q[3];
rz(-1.9469572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-0.1691549) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15360399) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(0.21501712) q[0];
x q[1];
rz(-2.1452227) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(0.043957274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.834224) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(2.6878396) q[1];
rz(-2.081359) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(-2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6358444) q[0];
sx q[0];
rz(-0.76857476) q[0];
sx q[0];
rz(-1.4710674) q[0];
x q[1];
rz(2.3381091) q[2];
sx q[2];
rz(-2.8065971) q[2];
sx q[2];
rz(-1.4575046) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2625418) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(2.2294728) q[1];
rz(-pi) q[2];
rz(1.7275229) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(-2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-2.267568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4955935) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(1.5587224) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1633515) q[2];
sx q[2];
rz(-1.904084) q[2];
sx q[2];
rz(2.8826706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8585426) q[1];
sx q[1];
rz(-1.0714604) q[1];
sx q[1];
rz(3.0602171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.9627409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261242) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(1.8490851) q[0];
x q[1];
rz(-3.0145698) q[2];
sx q[2];
rz(-1.6633031) q[2];
sx q[2];
rz(-2.9542343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8995754) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(2.1980397) q[1];
x q[2];
rz(-2.6188649) q[3];
sx q[3];
rz(-1.3739112) q[3];
sx q[3];
rz(2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(1.4964676) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
