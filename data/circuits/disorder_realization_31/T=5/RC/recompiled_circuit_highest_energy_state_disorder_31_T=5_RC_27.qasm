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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(-1.2555726) q[0];
rz(1.7805055) q[1];
sx q[1];
rz(-2.1828916) q[1];
sx q[1];
rz(-0.60671848) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8162168) q[0];
sx q[0];
rz(-2.0097549) q[0];
sx q[0];
rz(-2.3828242) q[0];
rz(1.6902906) q[2];
sx q[2];
rz(-2.8748389) q[2];
sx q[2];
rz(-0.96386749) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42598924) q[1];
sx q[1];
rz(-1.2066325) q[1];
sx q[1];
rz(-0.13873546) q[1];
rz(1.6912724) q[3];
sx q[3];
rz(-1.7306149) q[3];
sx q[3];
rz(-0.32630703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9386193) q[2];
sx q[2];
rz(-1.1569269) q[2];
sx q[2];
rz(2.970676) q[2];
rz(1.1275229) q[3];
sx q[3];
rz(-0.5564965) q[3];
sx q[3];
rz(3.0881622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663986) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(-2.4531181) q[0];
rz(-1.7636048) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(0.14150208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6287891) q[0];
sx q[0];
rz(-2.4806285) q[0];
sx q[0];
rz(2.4557607) q[0];
rz(-1.6703105) q[2];
sx q[2];
rz(-1.9984102) q[2];
sx q[2];
rz(-2.7832246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32530871) q[1];
sx q[1];
rz(-2.6806596) q[1];
sx q[1];
rz(-2.465425) q[1];
x q[2];
rz(1.9555952) q[3];
sx q[3];
rz(-1.6452922) q[3];
sx q[3];
rz(2.7423679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7576393) q[2];
sx q[2];
rz(-0.44718224) q[2];
sx q[2];
rz(-3.1128856) q[2];
rz(2.7590397) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(2.9401275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89556995) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(2.7699455) q[0];
rz(-2.9314575) q[1];
sx q[1];
rz(-0.34966436) q[1];
sx q[1];
rz(1.9336611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9555617) q[0];
sx q[0];
rz(-0.81350858) q[0];
sx q[0];
rz(0.71788089) q[0];
rz(0.25044424) q[2];
sx q[2];
rz(-3.0565002) q[2];
sx q[2];
rz(-1.9988572) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.713088) q[1];
sx q[1];
rz(-2.2853973) q[1];
sx q[1];
rz(2.3355161) q[1];
rz(3.0284027) q[3];
sx q[3];
rz(-0.34750156) q[3];
sx q[3];
rz(-2.6828535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62260425) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(2.5034215) q[2];
rz(0.41247621) q[3];
sx q[3];
rz(-1.6715489) q[3];
sx q[3];
rz(2.0392058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0647122) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(-1.3077211) q[0];
rz(0.67963302) q[1];
sx q[1];
rz(-2.4145587) q[1];
sx q[1];
rz(1.1839428) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1935496) q[0];
sx q[0];
rz(-1.2468766) q[0];
sx q[0];
rz(-0.69635038) q[0];
x q[1];
rz(0.23573622) q[2];
sx q[2];
rz(-1.6833458) q[2];
sx q[2];
rz(1.3723918) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19540994) q[1];
sx q[1];
rz(-1.2210746) q[1];
sx q[1];
rz(1.6433952) q[1];
x q[2];
rz(1.697942) q[3];
sx q[3];
rz(-1.7004439) q[3];
sx q[3];
rz(-1.4722908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.022698553) q[2];
sx q[2];
rz(-0.1476295) q[2];
sx q[2];
rz(2.0856608) q[2];
rz(2.8583156) q[3];
sx q[3];
rz(-1.2669468) q[3];
sx q[3];
rz(-1.383708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088575514) q[0];
sx q[0];
rz(-2.179189) q[0];
sx q[0];
rz(-1.7330633) q[0];
rz(2.9119496) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.1913258) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9339122) q[0];
sx q[0];
rz(-1.4290853) q[0];
sx q[0];
rz(1.974154) q[0];
rz(1.8764087) q[2];
sx q[2];
rz(-1.6907881) q[2];
sx q[2];
rz(-0.12745276) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7745668) q[1];
sx q[1];
rz(-1.8205678) q[1];
sx q[1];
rz(2.7109409) q[1];
rz(-pi) q[2];
rz(-2.2681885) q[3];
sx q[3];
rz(-2.3583989) q[3];
sx q[3];
rz(0.37885538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4952937) q[2];
sx q[2];
rz(-0.36094347) q[2];
sx q[2];
rz(-1.4781282) q[2];
rz(0.16119257) q[3];
sx q[3];
rz(-1.6094004) q[3];
sx q[3];
rz(0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6394871) q[0];
sx q[0];
rz(-2.7171071) q[0];
sx q[0];
rz(0.91833997) q[0];
rz(0.25686747) q[1];
sx q[1];
rz(-0.99595064) q[1];
sx q[1];
rz(-2.0253983) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5115806) q[0];
sx q[0];
rz(-1.0718379) q[0];
sx q[0];
rz(-2.5415238) q[0];
rz(-2.844238) q[2];
sx q[2];
rz(-2.8562244) q[2];
sx q[2];
rz(3.0653846) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5815365) q[1];
sx q[1];
rz(-1.7651084) q[1];
sx q[1];
rz(2.8981853) q[1];
x q[2];
rz(-1.3969054) q[3];
sx q[3];
rz(-1.3309722) q[3];
sx q[3];
rz(-2.2563129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18672289) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(0.40697971) q[2];
rz(-0.61378971) q[3];
sx q[3];
rz(-1.5301306) q[3];
sx q[3];
rz(-0.49155244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75519049) q[0];
sx q[0];
rz(-2.3004005) q[0];
sx q[0];
rz(-2.2156773) q[0];
rz(-0.47261247) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(-0.47017631) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0553207) q[0];
sx q[0];
rz(-2.4013858) q[0];
sx q[0];
rz(1.6946778) q[0];
x q[1];
rz(2.7514303) q[2];
sx q[2];
rz(-2.1794126) q[2];
sx q[2];
rz(2.2450992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31071102) q[1];
sx q[1];
rz(-2.355) q[1];
sx q[1];
rz(1.6221694) q[1];
x q[2];
rz(0.84104611) q[3];
sx q[3];
rz(-1.5435092) q[3];
sx q[3];
rz(-2.5430665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4818695) q[2];
sx q[2];
rz(-1.0708151) q[2];
sx q[2];
rz(0.83079633) q[2];
rz(3.0610436) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(0.95675937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0944716) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(0.15810529) q[0];
rz(-0.72599167) q[1];
sx q[1];
rz(-0.43015614) q[1];
sx q[1];
rz(2.7714738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86240101) q[0];
sx q[0];
rz(-1.6625615) q[0];
sx q[0];
rz(-2.807995) q[0];
x q[1];
rz(2.3838777) q[2];
sx q[2];
rz(-1.9279459) q[2];
sx q[2];
rz(1.0432318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2850879) q[1];
sx q[1];
rz(-0.56724945) q[1];
sx q[1];
rz(-0.7817661) q[1];
x q[2];
rz(0.40970476) q[3];
sx q[3];
rz(-2.0783882) q[3];
sx q[3];
rz(-1.0445802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(0.13713169) q[2];
rz(-1.6454227) q[3];
sx q[3];
rz(-2.4479595) q[3];
sx q[3];
rz(2.6579198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(-0.12635669) q[0];
rz(-2.5670746) q[1];
sx q[1];
rz(-0.9205598) q[1];
sx q[1];
rz(-0.7472907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3355646) q[0];
sx q[0];
rz(-1.3845516) q[0];
sx q[0];
rz(1.9646909) q[0];
rz(-0.70238446) q[2];
sx q[2];
rz(-1.1051854) q[2];
sx q[2];
rz(1.7543751) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8992363) q[1];
sx q[1];
rz(-1.1287424) q[1];
sx q[1];
rz(-0.41190179) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46183698) q[3];
sx q[3];
rz(-0.29134068) q[3];
sx q[3];
rz(-2.9614465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1025461) q[2];
sx q[2];
rz(-2.0206082) q[2];
sx q[2];
rz(-0.54753629) q[2];
rz(-1.3012137) q[3];
sx q[3];
rz(-2.0042714) q[3];
sx q[3];
rz(0.12714061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8319156) q[0];
sx q[0];
rz(-1.7597821) q[0];
sx q[0];
rz(0.58050138) q[0];
rz(-0.25513395) q[1];
sx q[1];
rz(-2.3105123) q[1];
sx q[1];
rz(1.1855804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0190576) q[0];
sx q[0];
rz(-0.92611971) q[0];
sx q[0];
rz(-2.3700847) q[0];
rz(-1.2078778) q[2];
sx q[2];
rz(-1.5722154) q[2];
sx q[2];
rz(-3.0848173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23981389) q[1];
sx q[1];
rz(-2.0396775) q[1];
sx q[1];
rz(-1.6704876) q[1];
x q[2];
rz(1.3316292) q[3];
sx q[3];
rz(-0.73430919) q[3];
sx q[3];
rz(-1.7074761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(2.1346788) q[2];
rz(1.5236731) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.9699008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11378743) q[0];
sx q[0];
rz(-1.2600949) q[0];
sx q[0];
rz(-2.8366198) q[0];
rz(1.2420568) q[1];
sx q[1];
rz(-1.2868953) q[1];
sx q[1];
rz(-2.3931265) q[1];
rz(2.0025675) q[2];
sx q[2];
rz(-1.7285657) q[2];
sx q[2];
rz(0.68539455) q[2];
rz(-2.5438944) q[3];
sx q[3];
rz(-2.2168474) q[3];
sx q[3];
rz(-0.26099152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
