OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(4.5067956) q[0];
sx q[0];
rz(11.542008) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60351935) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(-pi) q[1];
rz(2.3697882) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-0.31847218) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.484326) q[1];
sx q[1];
rz(-0.53521672) q[1];
sx q[1];
rz(0.80429299) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4952881) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26596507) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-0.47505501) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(2.1038726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514725) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-pi) q[1];
rz(2.6056387) q[2];
sx q[2];
rz(-2.0336656) q[2];
sx q[2];
rz(-3.0218389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84320074) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(0.028224736) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59229895) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-0.084687106) q[2];
rz(0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-0.57317615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37704913) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(-3.0717875) q[0];
rz(-2.3689752) q[2];
sx q[2];
rz(-1.801991) q[2];
sx q[2];
rz(2.2881743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9925476) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(2.4942314) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8543386) q[3];
sx q[3];
rz(-2.5594098) q[3];
sx q[3];
rz(2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(0.36270025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3804647) q[0];
sx q[0];
rz(-1.8109545) q[0];
sx q[0];
rz(2.2855177) q[0];
rz(-1.6022801) q[2];
sx q[2];
rz(-1.5891979) q[2];
sx q[2];
rz(0.035426332) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.6994233) q[1];
sx q[1];
rz(-1.4675092) q[1];
sx q[1];
rz(-2.5073754) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1953771) q[3];
sx q[3];
rz(-1.6295625) q[3];
sx q[3];
rz(-2.0277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(3.1385699) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(3.1242127) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.4594706) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634954) q[0];
sx q[0];
rz(-2.8601544) q[0];
sx q[0];
rz(2.1241758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1826477) q[2];
sx q[2];
rz(-2.358837) q[2];
sx q[2];
rz(0.81629717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0121213) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(0.71838897) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0586117) q[3];
sx q[3];
rz(-1.8564463) q[3];
sx q[3];
rz(0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370699) q[0];
sx q[0];
rz(-1.5504019) q[0];
sx q[0];
rz(1.5674595) q[0];
rz(1.6476829) q[2];
sx q[2];
rz(-1.7711519) q[2];
sx q[2];
rz(-0.69918699) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12339679) q[1];
sx q[1];
rz(-1.0180078) q[1];
sx q[1];
rz(2.4559896) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74216446) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-0.62430635) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(-1.2929582) q[0];
rz(2.675266) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(1.2517267) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.149951) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(2.409694) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1948366) q[3];
sx q[3];
rz(-1.1723926) q[3];
sx q[3];
rz(1.0362253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(0.12891842) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55420586) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(-1.9681853) q[0];
rz(2.2419937) q[2];
sx q[2];
rz(-2.0313782) q[2];
sx q[2];
rz(-0.34924289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.376437) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(0.57904412) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4017879) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(-1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63697469) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46247813) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(1.0366584) q[0];
rz(1.9753014) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(-1.5664958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1180601) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(-0.94783028) q[1];
rz(-1.8072855) q[3];
sx q[3];
rz(-0.83179501) q[3];
sx q[3];
rz(-2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-0.54668033) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71781681) q[0];
sx q[0];
rz(-0.86200889) q[0];
sx q[0];
rz(0.29341673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28995138) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(-1.314756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(0.22565266) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(-1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.5236241) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
rz(2.8125433) q[3];
sx q[3];
rz(-2.0206738) q[3];
sx q[3];
rz(-2.1856367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];