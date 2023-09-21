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
rz(-2.536474) q[1];
sx q[1];
rz(-2.6095698) q[1];
sx q[1];
rz(-1.1693118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3496075) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(-2.4277359) q[0];
rz(-0.77180441) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(0.31847218) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(-2.3372997) q[1];
x q[2];
rz(-0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(0.39428082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
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
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
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
rz(-1.0377201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71887165) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(1.3119112) q[0];
rz(2.0966881) q[2];
sx q[2];
rz(-2.0453339) q[2];
sx q[2];
rz(1.431312) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4251551) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(1.1643216) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(-2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-0.084687106) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.39069191) q[1];
sx q[1];
rz(-2.5707468) q[1];
sx q[1];
rz(0.57317615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37704913) q[0];
sx q[0];
rz(-1.1408313) q[0];
sx q[0];
rz(-3.0717875) q[0];
x q[1];
rz(-0.77261749) q[2];
sx q[2];
rz(-1.801991) q[2];
sx q[2];
rz(0.85341838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14904505) q[1];
sx q[1];
rz(-2.6547975) q[1];
sx q[1];
rz(-2.4942314) q[1];
rz(-2.0315941) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-2.9320419) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594584) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3804647) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(-0.856075) q[0];
rz(-0.018410725) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(-1.6056431) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73211654) q[1];
sx q[1];
rz(-2.5001657) q[1];
sx q[1];
rz(-0.17318053) q[1];
x q[2];
rz(1.4705212) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(-0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.592011) q[0];
sx q[0];
rz(-1.8093384) q[0];
sx q[0];
rz(0.15079389) q[0];
rz(-pi) q[1];
rz(-2.622501) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(-1.5965243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-2.8399889) q[1];
sx q[1];
rz(2.4232037) q[1];
rz(-pi) q[2];
rz(-0.32096433) q[3];
sx q[3];
rz(-1.1043613) q[3];
sx q[3];
rz(-2.5209559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-0.84189502) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(1.5766597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(-3.0029283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667144) q[0];
sx q[0];
rz(-0.020665558) q[0];
sx q[0];
rz(-0.1621577) q[0];
x q[1];
rz(1.6476829) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(0.69918699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12339679) q[1];
sx q[1];
rz(-2.1235848) q[1];
sx q[1];
rz(-2.4559896) q[1];
rz(-1.1451374) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3665294) q[2];
sx q[2];
rz(-2.0584006) q[2];
sx q[2];
rz(1.8072051) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(1.2929582) q[0];
rz(-pi) q[1];
rz(-0.94787089) q[2];
sx q[2];
rz(-1.9588753) q[2];
sx q[2];
rz(2.5571103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.034060409) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(0.88453102) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1950486) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.265825) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(-3.0126742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55420586) q[0];
sx q[0];
rz(-1.4652068) q[0];
sx q[0];
rz(1.9681853) q[0];
x q[1];
rz(0.56477408) q[2];
sx q[2];
rz(-2.161705) q[2];
sx q[2];
rz(1.5806944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1210632) q[1];
sx q[1];
rz(-1.9951092) q[1];
sx q[1];
rz(-1.8706277) q[1];
x q[2];
rz(-2.4017879) q[3];
sx q[3];
rz(-1.5825669) q[3];
sx q[3];
rz(-1.4307601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352167) q[0];
sx q[0];
rz(-1.750964) q[0];
sx q[0];
rz(-3.0336477) q[0];
rz(-pi) q[1];
rz(1.1662912) q[2];
sx q[2];
rz(-1.0317689) q[2];
sx q[2];
rz(-1.5750969) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0235325) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(-0.94783028) q[1];
rz(-pi) q[2];
rz(1.8072855) q[3];
sx q[3];
rz(-0.83179501) q[3];
sx q[3];
rz(2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(-1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.91530144) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-2.5949123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71781681) q[0];
sx q[0];
rz(-0.86200889) q[0];
sx q[0];
rz(2.8481759) q[0];
rz(-pi) q[1];
rz(-1.4184065) q[2];
sx q[2];
rz(-2.0414464) q[2];
sx q[2];
rz(-1.6413123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(1.5323557) q[1];
x q[2];
rz(0.22565266) q[3];
sx q[3];
rz(-1.7060346) q[3];
sx q[3];
rz(1.3868563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.5590177) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-2.6314541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-1.2658723) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
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
