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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697203) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(-1.1138492) q[0];
rz(-pi) q[1];
rz(2.4835303) q[2];
sx q[2];
rz(-2.1076638) q[2];
sx q[2];
rz(-1.8368349) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(1.1671288) q[1];
x q[2];
rz(1.4952881) q[3];
sx q[3];
rz(-1.8823874) q[3];
sx q[3];
rz(-0.15234767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-1.0377201) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(0.56818509) q[0];
rz(-pi) q[1];
rz(-2.3677164) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(2.3351923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71643752) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(1.9772711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(-0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-0.084687106) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111074) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(0.95570046) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(-1.7217365) q[0];
rz(-pi) q[1];
rz(2.8163221) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(-0.48649597) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5830071) q[1];
sx q[1];
rz(-1.9531986) q[1];
sx q[1];
rz(1.8797727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1099986) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(-0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3804647) q[0];
sx q[0];
rz(-1.8109545) q[0];
sx q[0];
rz(0.856075) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0417468) q[2];
sx q[2];
rz(-0.036465557) q[2];
sx q[2];
rz(-2.1349825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1945222) q[1];
sx q[1];
rz(-0.94049373) q[1];
sx q[1];
rz(1.4428201) q[1];
x q[2];
rz(1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(3.127393) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057102324) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(1.3296207) q[0];
rz(-0.95894496) q[2];
sx q[2];
rz(-0.78275567) q[2];
sx q[2];
rz(-0.81629717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7534605) q[1];
sx q[1];
rz(-1.3452483) q[1];
sx q[1];
rz(1.3688341) q[1];
rz(-pi) q[2];
rz(-0.32096433) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7370699) q[0];
sx q[0];
rz(-1.5504019) q[0];
sx q[0];
rz(-1.5674595) q[0];
x q[1];
rz(-2.7799941) q[2];
sx q[2];
rz(-0.21441678) q[2];
sx q[2];
rz(-0.32985652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2648592) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(-0.77244669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1451374) q[3];
sx q[3];
rz(-0.87493757) q[3];
sx q[3];
rz(1.828572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7666703) q[0];
sx q[0];
rz(-2.2793988) q[0];
sx q[0];
rz(-1.8486345) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46632669) q[2];
sx q[2];
rz(-2.1412686) q[2];
sx q[2];
rz(-1.2517267) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(2.2570616) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94675605) q[3];
sx q[3];
rz(-1.9692) q[3];
sx q[3];
rz(1.0362253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873868) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(1.9681853) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2419937) q[2];
sx q[2];
rz(-2.0313782) q[2];
sx q[2];
rz(-0.34924289) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32313777) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(-0.44169359) q[1];
rz(-0.017459004) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(3.0144514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(-1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(-1.3508266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006376) q[0];
sx q[0];
rz(-1.3906286) q[0];
sx q[0];
rz(3.0336477) q[0];
rz(1.9753014) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(-1.5664958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(-2.4904817) q[1];
x q[2];
rz(2.8899367) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(2.3264399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(0.11432153) q[2];
rz(-1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-0.84035994) q[0];
rz(-pi) q[1];
rz(-1.4184065) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(1.5002804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9676799) q[1];
sx q[1];
rz(-2.850107) q[1];
sx q[1];
rz(-0.12853865) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4320847) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(-0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(-2.0426345) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
