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
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(-3.1402052) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.815925) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(1.146846) q[0];
rz(-pi) q[1];
rz(-0.44714655) q[2];
sx q[2];
rz(-1.5329754) q[2];
sx q[2];
rz(0.58085261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4692987) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(-1.6507571) q[1];
rz(0.2557405) q[3];
sx q[3];
rz(-1.9307973) q[3];
sx q[3];
rz(-2.6933302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.095123) q[0];
sx q[0];
rz(-0.9233343) q[0];
sx q[0];
rz(-1.6460653) q[0];
rz(-pi) q[1];
rz(-0.30351992) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(-2.7941861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26317877) q[1];
sx q[1];
rz(-0.30305949) q[1];
sx q[1];
rz(-0.11942272) q[1];
x q[2];
rz(1.3482413) q[3];
sx q[3];
rz(-0.9369623) q[3];
sx q[3];
rz(-0.27130493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-2.5088076) q[2];
rz(-1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(0.99951807) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777269) q[0];
sx q[0];
rz(-2.8054872) q[0];
sx q[0];
rz(2.0396114) q[0];
x q[1];
rz(-2.2314084) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(-0.78188932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1516583) q[1];
sx q[1];
rz(-2.5807568) q[1];
sx q[1];
rz(2.6480617) q[1];
rz(2.0796892) q[3];
sx q[3];
rz(-1.5012) q[3];
sx q[3];
rz(-3.0825465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9592322) q[0];
sx q[0];
rz(-1.649547) q[0];
sx q[0];
rz(1.6596646) q[0];
rz(-pi) q[1];
rz(3.1399973) q[2];
sx q[2];
rz(-1.110382) q[2];
sx q[2];
rz(-0.38412016) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8400164) q[1];
sx q[1];
rz(-0.77743545) q[1];
sx q[1];
rz(-2.2393054) q[1];
rz(-1.7219909) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(0.99779469) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.516974) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18407962) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(1.0190796) q[0];
rz(-2.5324608) q[2];
sx q[2];
rz(-0.62437781) q[2];
sx q[2];
rz(-2.6513211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(2.7108971) q[1];
rz(0.59668221) q[3];
sx q[3];
rz(-2.0697069) q[3];
sx q[3];
rz(2.1614206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(0.34067672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6402123) q[0];
sx q[0];
rz(-0.073177241) q[0];
sx q[0];
rz(2.0206547) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31008115) q[2];
sx q[2];
rz(-2.3611464) q[2];
sx q[2];
rz(2.9030637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50960474) q[1];
sx q[1];
rz(-2.0443516) q[1];
sx q[1];
rz(1.8297086) q[1];
x q[2];
rz(-2.0601222) q[3];
sx q[3];
rz(-1.193207) q[3];
sx q[3];
rz(1.9469572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(-2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-0.46494928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60663215) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(-1.9726994) q[0];
rz(-pi) q[1];
rz(0.49352383) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.7871737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3073687) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(-0.4537531) q[1];
x q[2];
rz(2.081359) q[3];
sx q[3];
rz(-0.052882346) q[3];
sx q[3];
rz(-2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-3.0632339) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975472) q[0];
sx q[0];
rz(-0.80701485) q[0];
sx q[0];
rz(0.09597309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9044754) q[2];
sx q[2];
rz(-1.3318921) q[2];
sx q[2];
rz(2.4799926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17532149) q[1];
sx q[1];
rz(-0.92370719) q[1];
sx q[1];
rz(2.9233169) q[1];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-0.16031081) q[3];
sx q[3];
rz(2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3437929) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(-2.2180166) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(2.267568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0747781) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(-0.76807036) q[0];
x q[1];
rz(-2.2888695) q[2];
sx q[2];
rz(-0.52041473) q[2];
sx q[2];
rz(-2.4783217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8928788) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(-2.0715269) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(2.4172879) q[0];
rz(0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26319474) q[0];
sx q[0];
rz(-2.7976755) q[0];
sx q[0];
rz(0.92349903) q[0];
rz(-pi) q[1];
rz(-3.0145698) q[2];
sx q[2];
rz(-1.4782895) q[2];
sx q[2];
rz(-0.18735838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8995754) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(-2.1980397) q[1];
x q[2];
rz(1.7970656) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(-0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-2.010871) q[2];
sx q[2];
rz(-1.5390839) q[2];
sx q[2];
rz(2.5577953) q[2];
rz(-1.7189797) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];