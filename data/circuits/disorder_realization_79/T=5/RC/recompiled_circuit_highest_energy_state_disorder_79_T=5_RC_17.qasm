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
rz(1.5647178) q[0];
sx q[0];
rz(1.8160507) q[0];
sx q[0];
rz(9.9845822) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(1.8140732) q[1];
sx q[1];
rz(10.811364) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5965393) q[0];
sx q[0];
rz(-0.60205108) q[0];
sx q[0];
rz(1.1132147) q[0];
rz(-pi) q[1];
rz(1.4580016) q[2];
sx q[2];
rz(-1.8646282) q[2];
sx q[2];
rz(0.71982671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6416642) q[1];
sx q[1];
rz(-1.6545424) q[1];
sx q[1];
rz(-2.9332764) q[1];
rz(1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(2.5293722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6036512) q[2];
sx q[2];
rz(-2.5141022) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(-0.80080992) q[3];
sx q[3];
rz(-1.2468612) q[3];
sx q[3];
rz(0.39113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9073198) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(0.32670879) q[0];
rz(-1.3455343) q[1];
sx q[1];
rz(-1.5252472) q[1];
sx q[1];
rz(-0.84652841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162982) q[0];
sx q[0];
rz(-2.7413053) q[0];
sx q[0];
rz(1.2918703) q[0];
x q[1];
rz(1.2319599) q[2];
sx q[2];
rz(-2.4299271) q[2];
sx q[2];
rz(0.5063048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1607779) q[1];
sx q[1];
rz(-0.52923584) q[1];
sx q[1];
rz(1.5917626) q[1];
x q[2];
rz(1.2760824) q[3];
sx q[3];
rz(-2.3008153) q[3];
sx q[3];
rz(1.7129912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.232406) q[2];
sx q[2];
rz(-1.0973944) q[2];
sx q[2];
rz(2.6573913) q[2];
rz(-1.8353315) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(1.1854712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66121286) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(0.54006201) q[0];
rz(1.1506608) q[1];
sx q[1];
rz(-1.8692503) q[1];
sx q[1];
rz(0.59404341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0399678) q[0];
sx q[0];
rz(-1.7408082) q[0];
sx q[0];
rz(-3.0307458) q[0];
rz(-1.8833022) q[2];
sx q[2];
rz(-0.93949003) q[2];
sx q[2];
rz(2.3986651) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.328209) q[1];
sx q[1];
rz(-1.6594619) q[1];
sx q[1];
rz(0.81221272) q[1];
x q[2];
rz(2.0527127) q[3];
sx q[3];
rz(-1.6301117) q[3];
sx q[3];
rz(0.35463599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8192886) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(0.6907531) q[2];
rz(0.47762075) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(-0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531042) q[0];
sx q[0];
rz(-2.7957343) q[0];
sx q[0];
rz(-2.389287) q[0];
rz(0.90657702) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(-0.22901542) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9045712) q[0];
sx q[0];
rz(-1.8638049) q[0];
sx q[0];
rz(1.8051534) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36831736) q[2];
sx q[2];
rz(-0.57007705) q[2];
sx q[2];
rz(-1.7961479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6743857) q[1];
sx q[1];
rz(-1.3027281) q[1];
sx q[1];
rz(1.8169136) q[1];
x q[2];
rz(2.5563779) q[3];
sx q[3];
rz(-1.8793595) q[3];
sx q[3];
rz(-1.9208966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(-1.0329186) q[2];
rz(2.0756663) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(-1.3037995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8835939) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(-1.3294543) q[0];
rz(2.6138002) q[1];
sx q[1];
rz(-1.8254435) q[1];
sx q[1];
rz(-2.7425308) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2891727) q[0];
sx q[0];
rz(-1.6483083) q[0];
sx q[0];
rz(-1.9231251) q[0];
rz(-1.9827323) q[2];
sx q[2];
rz(-0.88232458) q[2];
sx q[2];
rz(-2.7044123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7150103) q[1];
sx q[1];
rz(-1.8816598) q[1];
sx q[1];
rz(-2.3764627) q[1];
x q[2];
rz(-2.2656029) q[3];
sx q[3];
rz(-1.4238947) q[3];
sx q[3];
rz(1.5796652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.993678) q[2];
sx q[2];
rz(-0.75239158) q[2];
sx q[2];
rz(2.1368775) q[2];
rz(-2.8435454) q[3];
sx q[3];
rz(-1.7883251) q[3];
sx q[3];
rz(-2.0731488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816958) q[0];
sx q[0];
rz(-1.6258465) q[0];
sx q[0];
rz(-1.1915278) q[0];
rz(-0.53939348) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(-0.76621145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8704925) q[0];
sx q[0];
rz(-1.2249761) q[0];
sx q[0];
rz(2.150751) q[0];
x q[1];
rz(0.91583385) q[2];
sx q[2];
rz(-2.1127491) q[2];
sx q[2];
rz(0.65108991) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0141361) q[1];
sx q[1];
rz(-1.3434504) q[1];
sx q[1];
rz(2.3037687) q[1];
rz(3.0449103) q[3];
sx q[3];
rz(-1.9396693) q[3];
sx q[3];
rz(2.3547805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.811565) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(-2.2583029) q[2];
rz(-0.53089321) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(-0.77597165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(2.4921932) q[0];
rz(0.3736639) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(-1.9992453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9050347) q[0];
sx q[0];
rz(-1.243374) q[0];
sx q[0];
rz(0.16752807) q[0];
rz(1.5863922) q[2];
sx q[2];
rz(-1.7806539) q[2];
sx q[2];
rz(-0.02053989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5264633) q[1];
sx q[1];
rz(-0.81306902) q[1];
sx q[1];
rz(-0.43083664) q[1];
x q[2];
rz(1.4799398) q[3];
sx q[3];
rz(-2.6808219) q[3];
sx q[3];
rz(-0.49660578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1335527) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(-2.7827061) q[2];
rz(2.533203) q[3];
sx q[3];
rz(-1.5549992) q[3];
sx q[3];
rz(2.2801094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95978874) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(0.96610075) q[0];
rz(2.6995662) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(0.40837049) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934487) q[0];
sx q[0];
rz(-1.7033332) q[0];
sx q[0];
rz(-3.114633) q[0];
rz(-1.5277065) q[2];
sx q[2];
rz(-1.7183964) q[2];
sx q[2];
rz(-1.6809631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1392532) q[1];
sx q[1];
rz(-2.1531257) q[1];
sx q[1];
rz(2.9703388) q[1];
rz(1.8281874) q[3];
sx q[3];
rz(-2.2857749) q[3];
sx q[3];
rz(1.4090713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24892204) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(1.9425617) q[2];
rz(-1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(3.1395932) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93160981) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(1.1071052) q[0];
rz(0.18868407) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(-1.8170549) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134091) q[0];
sx q[0];
rz(-1.285901) q[0];
sx q[0];
rz(2.6064742) q[0];
rz(-pi) q[1];
rz(2.1679384) q[2];
sx q[2];
rz(-1.5882379) q[2];
sx q[2];
rz(-2.9090925) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81428567) q[1];
sx q[1];
rz(-1.5734356) q[1];
sx q[1];
rz(-2.1111108) q[1];
rz(-2.4408165) q[3];
sx q[3];
rz(-1.1808529) q[3];
sx q[3];
rz(1.571947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3482427) q[2];
sx q[2];
rz(-0.57955727) q[2];
sx q[2];
rz(0.14774518) q[2];
rz(2.31367) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350236) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(-1.6178004) q[0];
rz(-1.891547) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(-2.7203383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68116248) q[0];
sx q[0];
rz(-1.958485) q[0];
sx q[0];
rz(-1.1739385) q[0];
x q[1];
rz(0.6217072) q[2];
sx q[2];
rz(-1.499147) q[2];
sx q[2];
rz(0.10393427) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6727461) q[1];
sx q[1];
rz(-1.4839659) q[1];
sx q[1];
rz(-1.5803017) q[1];
x q[2];
rz(-1.3280844) q[3];
sx q[3];
rz(-1.1303217) q[3];
sx q[3];
rz(2.3005568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0201575) q[2];
sx q[2];
rz(-2.406106) q[2];
sx q[2];
rz(-0.45041931) q[2];
rz(1.8592698) q[3];
sx q[3];
rz(-0.6855135) q[3];
sx q[3];
rz(2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0192075) q[0];
sx q[0];
rz(-1.6078124) q[0];
sx q[0];
rz(-1.5851198) q[0];
rz(-0.71313329) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(2.6789011) q[2];
sx q[2];
rz(-1.0924871) q[2];
sx q[2];
rz(0.15747216) q[2];
rz(0.73591821) q[3];
sx q[3];
rz(-0.10781846) q[3];
sx q[3];
rz(2.929647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
