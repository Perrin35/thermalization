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
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(2.6095698) q[1];
sx q[1];
rz(11.397059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7919851) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(-0.7138568) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77180441) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-0.31847218) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.2092023) q[1];
sx q[1];
rz(-1.1671288) q[1];
x q[2];
rz(1.6463046) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(-0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1319565) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(2.1038726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.422721) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(-1.8296815) q[0];
x q[1];
rz(2.3677164) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(0.80640031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4251551) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(1.1643216) q[1];
x q[2];
rz(-0.81289566) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(-1.9302492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-2.5684165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5432376) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(1.7217365) q[0];
rz(0.77261749) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(-2.2881743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.9531986) q[1];
sx q[1];
rz(1.8797727) q[1];
rz(-0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(-0.69609387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(0.20955071) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(0.36270025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.761128) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(2.2855177) q[0];
x q[1];
rz(-1.5393125) q[2];
sx q[2];
rz(-1.5523947) q[2];
sx q[2];
rz(-3.1061663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6994233) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(-2.5073754) q[1];
x q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(-2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5495816) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(-2.9907988) q[0];
rz(-pi) q[1];
rz(0.51909165) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(-1.5965243) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(0.71838897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0829809) q[3];
sx q[3];
rz(-1.8564463) q[3];
sx q[3];
rz(0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(3.0029283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370699) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(-1.5674595) q[0];
rz(-2.7799941) q[2];
sx q[2];
rz(-0.21441678) q[2];
sx q[2];
rz(-0.32985652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87673346) q[1];
sx q[1];
rz(-0.85163341) q[1];
sx q[1];
rz(0.77244669) q[1];
rz(0.45913978) q[3];
sx q[3];
rz(-0.79677478) q[3];
sx q[3];
rz(-1.2129984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-3.0464723) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(2.6050674) q[0];
rz(-0.58553186) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3793959) q[0];
sx q[0];
rz(-1.7805829) q[0];
sx q[0];
rz(-2.4136153) q[0];
x q[1];
rz(2.1937218) q[2];
sx q[2];
rz(-1.1827173) q[2];
sx q[2];
rz(0.58448234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0956456) q[1];
sx q[1];
rz(-1.026517) q[1];
sx q[1];
rz(-0.73927684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6631213) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0607973) q[0];
sx q[0];
rz(-1.9658488) q[0];
sx q[0];
rz(-0.11443826) q[0];
rz(-pi) q[1];
rz(-2.2419937) q[2];
sx q[2];
rz(-2.0313782) q[2];
sx q[2];
rz(0.34924289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0205295) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(1.8706277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.586732) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(0.15077886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794466) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(1.790766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(-1.3895967) q[0];
x q[1];
rz(2.5647854) q[2];
sx q[2];
rz(-1.2262605) q[2];
sx q[2];
rz(-0.21201269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13547922) q[1];
sx q[1];
rz(-1.0891501) q[1];
sx q[1];
rz(-2.3856132) q[1];
rz(-pi) q[2];
rz(-0.25165598) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(-2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(0.54668033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(2.3012327) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4184065) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(1.6413123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4320847) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.6179686) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
rz(2.1605282) q[3];
sx q[3];
rz(-2.5909501) q[3];
sx q[3];
rz(0.29028374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
