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
rz(0.53606755) q[0];
sx q[0];
rz(-1.9789088) q[0];
sx q[0];
rz(-0.43098488) q[0];
rz(-2.7886136) q[1];
sx q[1];
rz(-1.2231491) q[1];
sx q[1];
rz(-0.42981848) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777464) q[0];
sx q[0];
rz(-2.060685) q[0];
sx q[0];
rz(1.4686945) q[0];
rz(0.016021803) q[2];
sx q[2];
rz(-2.0031158) q[2];
sx q[2];
rz(2.1825618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42001549) q[1];
sx q[1];
rz(-2.7215981) q[1];
sx q[1];
rz(-0.53424044) q[1];
rz(0.63685959) q[3];
sx q[3];
rz(-0.59266337) q[3];
sx q[3];
rz(-0.67837472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.063804403) q[2];
sx q[2];
rz(-2.6142075) q[2];
sx q[2];
rz(-2.1303614) q[2];
rz(2.1003335) q[3];
sx q[3];
rz(-0.51308378) q[3];
sx q[3];
rz(0.43963715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0332727) q[0];
sx q[0];
rz(-0.99650538) q[0];
sx q[0];
rz(0.74747768) q[0];
rz(-0.29391089) q[1];
sx q[1];
rz(-1.1312609) q[1];
sx q[1];
rz(-2.8864313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145484) q[0];
sx q[0];
rz(-0.60227312) q[0];
sx q[0];
rz(-1.8667579) q[0];
x q[1];
rz(2.8361144) q[2];
sx q[2];
rz(-1.7490088) q[2];
sx q[2];
rz(-2.5199948) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1595711) q[1];
sx q[1];
rz(-0.82948331) q[1];
sx q[1];
rz(-2.1950566) q[1];
rz(0.18682602) q[3];
sx q[3];
rz(-0.95462026) q[3];
sx q[3];
rz(-2.0957859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1110288) q[2];
sx q[2];
rz(-0.054628987) q[2];
sx q[2];
rz(-0.89152208) q[2];
rz(0.24957481) q[3];
sx q[3];
rz(-0.88281837) q[3];
sx q[3];
rz(-2.7454624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075277) q[0];
sx q[0];
rz(-2.0052818) q[0];
sx q[0];
rz(-0.044128142) q[0];
rz(1.2491501) q[1];
sx q[1];
rz(-2.7154778) q[1];
sx q[1];
rz(2.6557907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6225548) q[0];
sx q[0];
rz(-1.7960894) q[0];
sx q[0];
rz(1.5507429) q[0];
x q[1];
rz(3.0032512) q[2];
sx q[2];
rz(-1.3855583) q[2];
sx q[2];
rz(-2.0092162) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1399849) q[1];
sx q[1];
rz(-2.1473653) q[1];
sx q[1];
rz(-2.216061) q[1];
x q[2];
rz(-0.39704571) q[3];
sx q[3];
rz(-1.9730803) q[3];
sx q[3];
rz(2.8232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8140063) q[2];
sx q[2];
rz(-2.0531211) q[2];
sx q[2];
rz(2.4468454) q[2];
rz(-2.5658549) q[3];
sx q[3];
rz(-1.8703987) q[3];
sx q[3];
rz(-1.7339138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95076743) q[0];
sx q[0];
rz(-0.29470834) q[0];
sx q[0];
rz(-2.8488391) q[0];
rz(0.5689019) q[1];
sx q[1];
rz(-0.82745537) q[1];
sx q[1];
rz(2.2946766) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900565) q[0];
sx q[0];
rz(-2.7256294) q[0];
sx q[0];
rz(1.7638168) q[0];
rz(-0.072242468) q[2];
sx q[2];
rz(-0.90257593) q[2];
sx q[2];
rz(0.39230686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9255815) q[1];
sx q[1];
rz(-0.30864172) q[1];
sx q[1];
rz(2.4587548) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7645232) q[3];
sx q[3];
rz(-1.6370557) q[3];
sx q[3];
rz(2.3364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34951052) q[2];
sx q[2];
rz(-1.8831207) q[2];
sx q[2];
rz(-2.9328031) q[2];
rz(0.71792349) q[3];
sx q[3];
rz(-2.8585298) q[3];
sx q[3];
rz(2.2010402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013041) q[0];
sx q[0];
rz(-2.6191819) q[0];
sx q[0];
rz(-2.8029602) q[0];
rz(-0.70308095) q[1];
sx q[1];
rz(-0.73892006) q[1];
sx q[1];
rz(0.83831659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037443) q[0];
sx q[0];
rz(-0.5068584) q[0];
sx q[0];
rz(2.5985107) q[0];
x q[1];
rz(-1.7846858) q[2];
sx q[2];
rz(-0.75816064) q[2];
sx q[2];
rz(1.2616829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.053239659) q[1];
sx q[1];
rz(-0.87118282) q[1];
sx q[1];
rz(0.043645388) q[1];
x q[2];
rz(1.1456212) q[3];
sx q[3];
rz(-2.1347649) q[3];
sx q[3];
rz(-0.26568809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68359739) q[2];
sx q[2];
rz(-1.0883051) q[2];
sx q[2];
rz(-2.1387157) q[2];
rz(0.14497997) q[3];
sx q[3];
rz(-0.093791157) q[3];
sx q[3];
rz(-3.0785479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.1231287) q[0];
sx q[0];
rz(-0.58582425) q[0];
sx q[0];
rz(0.96328324) q[0];
rz(-2.9604498) q[1];
sx q[1];
rz(-0.89814848) q[1];
sx q[1];
rz(1.9021665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604804) q[0];
sx q[0];
rz(-1.1976997) q[0];
sx q[0];
rz(0.65848559) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6685817) q[2];
sx q[2];
rz(-1.3227229) q[2];
sx q[2];
rz(2.6864227) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9628004) q[1];
sx q[1];
rz(-1.1988678) q[1];
sx q[1];
rz(0.40403251) q[1];
x q[2];
rz(1.7224947) q[3];
sx q[3];
rz(-1.6983508) q[3];
sx q[3];
rz(0.90326819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2565101) q[2];
sx q[2];
rz(-0.91826597) q[2];
sx q[2];
rz(0.93290848) q[2];
rz(1.7289303) q[3];
sx q[3];
rz(-1.4224854) q[3];
sx q[3];
rz(-2.9322114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4466208) q[0];
sx q[0];
rz(-0.3681204) q[0];
sx q[0];
rz(-2.3387261) q[0];
rz(2.39957) q[1];
sx q[1];
rz(-1.4890393) q[1];
sx q[1];
rz(-0.41025695) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3543676) q[0];
sx q[0];
rz(-0.64866186) q[0];
sx q[0];
rz(0.58726585) q[0];
x q[1];
rz(3.018961) q[2];
sx q[2];
rz(-2.3735614) q[2];
sx q[2];
rz(2.1036069) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.81933016) q[1];
sx q[1];
rz(-0.17824379) q[1];
sx q[1];
rz(1.1398386) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5531261) q[3];
sx q[3];
rz(-0.7477254) q[3];
sx q[3];
rz(2.6315053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2731169) q[2];
sx q[2];
rz(-1.5831524) q[2];
sx q[2];
rz(-0.22129076) q[2];
rz(-0.41915974) q[3];
sx q[3];
rz(-2.6594888) q[3];
sx q[3];
rz(2.6656849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23685037) q[0];
sx q[0];
rz(-2.1433647) q[0];
sx q[0];
rz(-1.5297484) q[0];
rz(-1.3329685) q[1];
sx q[1];
rz(-0.95104853) q[1];
sx q[1];
rz(1.6882247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26448956) q[0];
sx q[0];
rz(-0.96142381) q[0];
sx q[0];
rz(-2.2731056) q[0];
rz(1.4508993) q[2];
sx q[2];
rz(-1.3812307) q[2];
sx q[2];
rz(1.8360863) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.087880922) q[1];
sx q[1];
rz(-0.86999991) q[1];
sx q[1];
rz(-1.3739963) q[1];
rz(-pi) q[2];
rz(0.15590053) q[3];
sx q[3];
rz(-1.2955126) q[3];
sx q[3];
rz(0.017214765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2754485) q[2];
sx q[2];
rz(-2.8871705) q[2];
sx q[2];
rz(-3.1308543) q[2];
rz(0.3206611) q[3];
sx q[3];
rz(-1.2977707) q[3];
sx q[3];
rz(-0.58967725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43309942) q[0];
sx q[0];
rz(-0.95516094) q[0];
sx q[0];
rz(0.47472111) q[0];
rz(-0.2933329) q[1];
sx q[1];
rz(-1.1055929) q[1];
sx q[1];
rz(-1.6620103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06874456) q[0];
sx q[0];
rz(-1.6075875) q[0];
sx q[0];
rz(1.6108914) q[0];
x q[1];
rz(-1.8117254) q[2];
sx q[2];
rz(-2.2325667) q[2];
sx q[2];
rz(0.84217254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8822215) q[1];
sx q[1];
rz(-2.6135635) q[1];
sx q[1];
rz(0.47059437) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18085705) q[3];
sx q[3];
rz(-1.81362) q[3];
sx q[3];
rz(2.9202785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.044363) q[2];
sx q[2];
rz(-1.2791415) q[2];
sx q[2];
rz(-0.43107671) q[2];
rz(0.19032446) q[3];
sx q[3];
rz(-2.7908235) q[3];
sx q[3];
rz(2.2569807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97656074) q[0];
sx q[0];
rz(-2.850552) q[0];
sx q[0];
rz(0.2188368) q[0];
rz(-0.95895514) q[1];
sx q[1];
rz(-2.405165) q[1];
sx q[1];
rz(1.3921907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16698027) q[0];
sx q[0];
rz(-1.7638399) q[0];
sx q[0];
rz(-1.5184605) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7213412) q[2];
sx q[2];
rz(-1.7027579) q[2];
sx q[2];
rz(-0.38334639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4569623) q[1];
sx q[1];
rz(-0.70856386) q[1];
sx q[1];
rz(2.674391) q[1];
x q[2];
rz(2.1851319) q[3];
sx q[3];
rz(-1.2117361) q[3];
sx q[3];
rz(-2.7484855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86768156) q[2];
sx q[2];
rz(-0.76649222) q[2];
sx q[2];
rz(1.319818) q[2];
rz(-0.53226081) q[3];
sx q[3];
rz(-1.7936734) q[3];
sx q[3];
rz(1.5925315) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.559448) q[0];
sx q[0];
rz(-1.5912709) q[0];
sx q[0];
rz(-2.7406319) q[0];
rz(2.9136912) q[1];
sx q[1];
rz(-1.0453929) q[1];
sx q[1];
rz(1.4452404) q[1];
rz(-1.8310905) q[2];
sx q[2];
rz(-1.1681265) q[2];
sx q[2];
rz(0.10377965) q[2];
rz(-2.9798672) q[3];
sx q[3];
rz(-0.99746708) q[3];
sx q[3];
rz(0.4494638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
