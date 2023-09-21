OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(-2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7919851) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(-0.7138568) q[0];
rz(2.4835303) q[2];
sx q[2];
rz(-2.1076638) q[2];
sx q[2];
rz(1.3047578) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(1.9744639) q[1];
rz(0.23006769) q[3];
sx q[3];
rz(-0.32031968) q[3];
sx q[3];
rz(-2.7473118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(2.1038726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0449045) q[2];
sx q[2];
rz(-2.0453339) q[2];
sx q[2];
rz(1.7102807) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2983919) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(-3.1133679) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95275767) q[3];
sx q[3];
rz(-0.85988322) q[3];
sx q[3];
rz(0.074912138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5432376) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(-1.7217365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32527058) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(2.6550967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9925476) q[1];
sx q[1];
rz(-0.48679513) q[1];
sx q[1];
rz(-0.64736127) q[1];
rz(1.1099986) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(-2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(0.20955071) q[3];
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
rz(-0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(0.55066806) q[0];
rz(-2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(0.36270025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6837316) q[0];
sx q[0];
rz(-2.3944003) q[0];
sx q[0];
rz(1.213221) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(1.6056431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1945222) q[1];
sx q[1];
rz(-0.94049373) q[1];
sx q[1];
rz(1.6987726) q[1];
x q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.592011) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(0.15079389) q[0];
rz(-2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(2.8487157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38813218) q[1];
sx q[1];
rz(-1.3452483) q[1];
sx q[1];
rz(1.7727586) q[1];
x q[2];
rz(-2.8206283) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(-0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16620557) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(-3.1211981) q[0];
x q[1];
rz(-1.6476829) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(-0.69918699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8530635) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(0.89785518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9964553) q[3];
sx q[3];
rz(-0.87493757) q[3];
sx q[3];
rz(-1.828572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.8072051) q[2];
rz(1.1602317) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(2.5172863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7621967) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(2.4136153) q[0];
rz(-pi) q[1];
rz(2.1937218) q[2];
sx q[2];
rz(-1.1827173) q[2];
sx q[2];
rz(-2.5571103) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1075322) q[1];
sx q[1];
rz(-0.95648396) q[1];
sx q[1];
rz(-0.88453102) q[1];
rz(-pi) q[2];
rz(-0.94675605) q[3];
sx q[3];
rz(-1.9692) q[3];
sx q[3];
rz(-1.0362253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(3.110102) q[3];
sx q[3];
rz(-2.4523906) q[3];
sx q[3];
rz(0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55420586) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(-1.9681853) q[0];
x q[1];
rz(0.56477408) q[2];
sx q[2];
rz(-0.97988765) q[2];
sx q[2];
rz(-1.5806944) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32313777) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(-2.6998991) q[1];
rz(-pi) q[2];
rz(-1.5548607) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(-0.15077886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(-2.074266) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(-1.790766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352167) q[0];
sx q[0];
rz(-1.3906286) q[0];
sx q[0];
rz(-0.10794497) q[0];
rz(2.5596041) q[2];
sx q[2];
rz(-2.4798923) q[2];
sx q[2];
rz(2.2616507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2495888) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(2.4904817) q[1];
rz(-2.3885214) q[3];
sx q[3];
rz(-1.3967447) q[3];
sx q[3];
rz(-0.93833246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262912) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-2.9283438) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-2.5949123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-2.3843117) q[0];
sx q[0];
rz(1.8961294) q[0];
rz(-2.8516413) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(1.8268367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17391275) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(-0.12853865) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7095079) q[3];
sx q[3];
rz(-1.7943534) q[3];
sx q[3];
rz(-0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.5236241) q[2];
sx q[2];
rz(-1.5559559) q[2];
sx q[2];
rz(2.422239) q[2];
rz(2.0426345) q[3];
sx q[3];
rz(-1.8660587) q[3];
sx q[3];
rz(2.3793424) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
