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
rz(1.024363) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0697203) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(-1.1138492) q[0];
rz(-pi) q[1];
rz(2.2157482) q[2];
sx q[2];
rz(-1.0174417) q[2];
sx q[2];
rz(0.64253053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(-2.3372997) q[1];
x q[2];
rz(-0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(-1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(-1.0377201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25181928) q[0];
sx q[0];
rz(-0.60476859) q[0];
sx q[0];
rz(-2.7483727) q[0];
rz(-pi) q[1];
rz(-2.3677164) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(2.3351923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71643752) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(1.1643216) q[1];
rz(-pi) q[2];
rz(-0.81289566) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(-1.9302492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(2.7509007) q[1];
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
rz(2.5983551) q[0];
sx q[0];
rz(-0.43524536) q[0];
sx q[0];
rz(-1.7217365) q[0];
rz(0.32527058) q[2];
sx q[2];
rz(-2.3420482) q[2];
sx q[2];
rz(0.48649597) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5830071) q[1];
sx q[1];
rz(-1.9531986) q[1];
sx q[1];
rz(1.26182) q[1];
rz(-pi) q[2];
rz(0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(0.69609387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(1.9283717) q[0];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(-1.5359495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1945222) q[1];
sx q[1];
rz(-0.94049373) q[1];
sx q[1];
rz(1.6987726) q[1];
x q[2];
rz(-2.1953771) q[3];
sx q[3];
rz(-1.6295625) q[3];
sx q[3];
rz(-2.0277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(-2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844903) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(1.811972) q[0];
rz(-pi) q[1];
rz(2.622501) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(-1.5450684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0121213) q[1];
sx q[1];
rz(-0.3016037) q[1];
sx q[1];
rz(2.4232037) q[1];
x q[2];
rz(2.8206283) q[3];
sx q[3];
rz(-2.0372314) q[3];
sx q[3];
rz(-0.62063673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-0.84189502) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(-0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(0.13866436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16620557) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-3.1211981) q[0];
x q[1];
rz(-2.7799941) q[2];
sx q[2];
rz(-0.21441678) q[2];
sx q[2];
rz(2.8117361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87673346) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(-0.77244669) q[1];
rz(-pi) q[2];
rz(-1.1451374) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(-1.828572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038135197) q[0];
sx q[0];
rz(-0.75224829) q[0];
sx q[0];
rz(0.30970807) q[0];
rz(-pi) q[1];
rz(0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(-0.50146539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.045947) q[1];
sx q[1];
rz(-2.1150757) q[1];
sx q[1];
rz(-0.73927684) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55420586) q[0];
sx q[0];
rz(-1.6763858) q[0];
sx q[0];
rz(1.9681853) q[0];
rz(0.56477408) q[2];
sx q[2];
rz(-0.97988765) q[2];
sx q[2];
rz(1.5608982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.376437) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(-0.57904412) q[1];
x q[2];
rz(2.4017879) q[3];
sx q[3];
rz(-1.5590258) q[3];
sx q[3];
rz(-1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-2.1010667) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006376) q[0];
sx q[0];
rz(-1.750964) q[0];
sx q[0];
rz(0.10794497) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1662912) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(1.5664958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0235325) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(-2.1937624) q[1];
rz(1.3343072) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(-0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(2.9283438) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4237758) q[0];
sx q[0];
rz(-2.2795838) q[0];
sx q[0];
rz(-0.29341673) q[0];
x q[1];
rz(-1.4184065) q[2];
sx q[2];
rz(-1.1001462) q[2];
sx q[2];
rz(1.6413123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.039779546) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7095079) q[3];
sx q[3];
rz(-1.7943534) q[3];
sx q[3];
rz(2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(-3.1372916) q[2];
rz(-0.99758482) q[3];
sx q[3];
rz(-0.49013609) q[3];
sx q[3];
rz(2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(2.6976363) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(3.1267358) q[2];
sx q[2];
rz(-1.6179634) q[2];
sx q[2];
rz(0.85214324) q[2];
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
