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
rz(-2.3581714) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(2.3697882) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(-2.8231205) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6572666) q[1];
sx q[1];
rz(-0.53521672) q[1];
sx q[1];
rz(0.80429299) q[1];
rz(-0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(-2.7473118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8756276) q[2];
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
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1319565) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(1.3985727) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(2.1038726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99012016) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-2.5734076) q[0];
x q[1];
rz(0.53595397) q[2];
sx q[2];
rz(-1.1079271) q[2];
sx q[2];
rz(0.11975372) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4251551) q[1];
sx q[1];
rz(-1.5448703) q[1];
sx q[1];
rz(-1.1643216) q[1];
rz(-pi) q[2];
rz(-0.59229895) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(-0.37880138) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-0.57317615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1646106) q[0];
sx q[0];
rz(-1.6342388) q[0];
sx q[0];
rz(2.001686) q[0];
rz(-pi) q[1];
rz(-1.2531883) q[2];
sx q[2];
rz(-0.8237969) q[2];
sx q[2];
rz(-2.2044646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(1.26182) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0315941) q[3];
sx q[3];
rz(-1.9402383) q[3];
sx q[3];
rz(-2.8957469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(0.097269416) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(-0.20955071) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(2.6837316) q[0];
sx q[0];
rz(-0.74719238) q[0];
sx q[0];
rz(-1.9283717) q[0];
rz(-pi) q[1];
rz(3.1231819) q[2];
sx q[2];
rz(-1.6022748) q[2];
sx q[2];
rz(1.5359495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.6994233) q[1];
sx q[1];
rz(-1.4675092) q[1];
sx q[1];
rz(2.5073754) q[1];
x q[2];
rz(1.4705212) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4576733) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(-0.65888843) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(3.127393) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.4594706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844903) q[0];
sx q[0];
rz(-1.717289) q[0];
sx q[0];
rz(1.3296207) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51909165) q[2];
sx q[2];
rz(-2.1862098) q[2];
sx q[2];
rz(1.5965243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0046878) q[1];
sx q[1];
rz(-1.3740174) q[1];
sx q[1];
rz(-0.23006567) q[1];
rz(-pi) q[2];
rz(-2.1305389) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(-1.884348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4045227) q[0];
sx q[0];
rz(-1.5911907) q[0];
sx q[0];
rz(1.5674595) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3615985) q[2];
sx q[2];
rz(-2.9271759) q[2];
sx q[2];
rz(-0.32985652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2648592) q[1];
sx q[1];
rz(-0.85163341) q[1];
sx q[1];
rz(0.77244669) q[1];
x q[2];
rz(0.45913978) q[3];
sx q[3];
rz(-0.79677478) q[3];
sx q[3];
rz(1.9285942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(-1.9813609) q[3];
sx q[3];
rz(-1.7778054) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5360864) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(-0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(2.5172863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37492232) q[0];
sx q[0];
rz(-2.2793988) q[0];
sx q[0];
rz(-1.2929582) q[0];
rz(-0.95958556) q[2];
sx q[2];
rz(-2.4215536) q[2];
sx q[2];
rz(-2.6401273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.149951) q[1];
sx q[1];
rz(-0.88643622) q[1];
sx q[1];
rz(-2.409694) q[1];
rz(-pi) q[2];
rz(0.47847139) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(-0.80696054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-3.0016622) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(3.0126742) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873868) q[0];
sx q[0];
rz(-1.6763858) q[0];
sx q[0];
rz(1.1734074) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.244197) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(2.430254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7651556) q[1];
sx q[1];
rz(-2.6273478) q[1];
sx q[1];
rz(-0.57904412) q[1];
rz(-pi) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(3.0144514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(-1.0673267) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794466) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(0.016013913) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(-1.790766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(2.1049343) q[0];
rz(-0.58198858) q[2];
sx q[2];
rz(-0.66170035) q[2];
sx q[2];
rz(-2.2616507) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8920039) q[1];
sx q[1];
rz(-0.87000404) q[1];
sx q[1];
rz(0.65111098) q[1];
rz(-1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(2.6700499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0525557) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-0.54668033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(-2.3012327) q[0];
x q[1];
rz(0.47537739) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(0.00098468653) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6215366) q[1];
sx q[1];
rz(-1.5339508) q[1];
sx q[1];
rz(-2.8523793) q[1];
rz(-pi) q[2];
rz(0.54639001) q[3];
sx q[3];
rz(-2.8791109) q[3];
sx q[3];
rz(-0.7149834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(0.004301087) q[2];
rz(2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(3.1267358) q[2];
sx q[2];
rz(-1.6179634) q[2];
sx q[2];
rz(0.85214324) q[2];
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