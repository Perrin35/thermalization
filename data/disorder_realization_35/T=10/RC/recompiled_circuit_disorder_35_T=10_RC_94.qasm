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
rz(-0.53202283) q[1];
sx q[1];
rz(1.1693118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5380733) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(-2.6253683) q[0];
rz(-0.77180441) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(2.8231205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4989657) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(-2.7514003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8291679) q[3];
sx q[3];
rz(-1.6426622) q[3];
sx q[3];
rz(1.7463328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26596507) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(2.8406075) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1319565) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(-1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99012016) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-2.5734076) q[0];
rz(-pi) q[1];
rz(-2.0966881) q[2];
sx q[2];
rz(-2.0453339) q[2];
sx q[2];
rz(1.7102807) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.84320074) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(-3.1133679) q[1];
rz(2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66118801) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(2.5684165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7645435) q[0];
sx q[0];
rz(-2.0007613) q[0];
sx q[0];
rz(0.069805108) q[0];
rz(-pi) q[1];
rz(1.8884044) q[2];
sx q[2];
rz(-0.8237969) q[2];
sx q[2];
rz(0.93712805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(-1.8797727) q[1];
rz(-pi) q[2];
rz(-2.2872541) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(0.69609387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-2.8373575) q[2];
sx q[2];
rz(0.19392459) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(0.20955071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4578611) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(-1.213221) q[0];
rz(-pi) q[1];
rz(0.018410725) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(-1.6056431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.6994233) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(2.5073754) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0691931) q[3];
sx q[3];
rz(-2.1941333) q[3];
sx q[3];
rz(-0.49923957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-2.3390884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(0.014199646) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844903) q[0];
sx q[0];
rz(-1.4243037) q[0];
sx q[0];
rz(-1.3296207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.622501) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(-1.5965243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1369049) q[1];
sx q[1];
rz(-1.7675752) q[1];
sx q[1];
rz(-2.911527) q[1];
rz(-2.0586117) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(2.1250336) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16620557) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-3.1211981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9406592) q[2];
sx q[2];
rz(-1.6461419) q[2];
sx q[2];
rz(-0.88694015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2885292) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(0.89785518) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74216446) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(3.1165251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.8072051) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(-2.5172863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37492232) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(-1.2929582) q[0];
rz(-pi) q[1];
rz(-2.675266) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(1.2517267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0956456) q[1];
sx q[1];
rz(-1.026517) q[1];
sx q[1];
rz(2.4023158) q[1];
x q[2];
rz(0.47847139) q[3];
sx q[3];
rz(-1.0020743) q[3];
sx q[3];
rz(0.80696054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77057225) q[0];
sx q[0];
rz(-0.41045529) q[0];
sx q[0];
rz(1.8380941) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89739563) q[2];
sx q[2];
rz(-0.79332966) q[2];
sx q[2];
rz(-2.430254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7651556) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(0.57904412) q[1];
x q[2];
rz(-1.5548607) q[3];
sx q[3];
rz(-2.3105379) q[3];
sx q[3];
rz(2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(-1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(-1.3508266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58383656) q[0];
sx q[0];
rz(-1.4646052) q[0];
sx q[0];
rz(1.7519959) q[0];
x q[1];
rz(-2.5647854) q[2];
sx q[2];
rz(-1.9153321) q[2];
sx q[2];
rz(-0.21201269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2495888) q[1];
sx q[1];
rz(-2.2715886) q[1];
sx q[1];
rz(2.4904817) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3343072) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(0.21324883) q[0];
rz(-2.7217216) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(-2.5949123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-0.75728098) q[0];
sx q[0];
rz(1.2454633) q[0];
rz(-2.8516413) q[2];
sx q[2];
rz(-2.6486514) q[2];
sx q[2];
rz(1.314756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5200561) q[1];
sx q[1];
rz(-1.6076419) q[1];
sx q[1];
rz(2.8523793) q[1];
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
rz(-pi) q[1];
x q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(3.1372916) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.8757204) q[2];
sx q[2];
rz(-3.0921428) q[2];
sx q[2];
rz(-2.5947239) q[2];
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
