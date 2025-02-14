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
rz(0.70145488) q[0];
sx q[0];
rz(3.6297807) q[0];
sx q[0];
rz(9.0253137) q[0];
rz(-1.0723298) q[1];
sx q[1];
rz(-2.2430099) q[1];
sx q[1];
rz(2.8314765) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7880023) q[0];
sx q[0];
rz(-1.912477) q[0];
sx q[0];
rz(0.41312878) q[0];
rz(-pi) q[1];
rz(-2.1793876) q[2];
sx q[2];
rz(-0.96675379) q[2];
sx q[2];
rz(2.5317328) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0779841) q[1];
sx q[1];
rz(-1.0922752) q[1];
sx q[1];
rz(-1.4565181) q[1];
rz(3.0324494) q[3];
sx q[3];
rz(-1.9535084) q[3];
sx q[3];
rz(-0.34527895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.255379) q[2];
sx q[2];
rz(-1.9510521) q[2];
sx q[2];
rz(-1.0211241) q[2];
rz(1.6518263) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(-0.79045734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0806231) q[0];
sx q[0];
rz(-1.0401833) q[0];
sx q[0];
rz(-2.8676497) q[0];
rz(0.25320369) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(0.23987548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2027953) q[0];
sx q[0];
rz(-2.1442841) q[0];
sx q[0];
rz(2.5198562) q[0];
rz(-pi) q[1];
rz(-2.5794058) q[2];
sx q[2];
rz(-1.8442684) q[2];
sx q[2];
rz(-2.4386465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4644313) q[1];
sx q[1];
rz(-1.3422692) q[1];
sx q[1];
rz(0.34739574) q[1];
rz(-pi) q[2];
rz(0.20407014) q[3];
sx q[3];
rz(-0.7925472) q[3];
sx q[3];
rz(2.8784405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46951744) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(-1.391927) q[2];
rz(-2.4972656) q[3];
sx q[3];
rz(-1.143012) q[3];
sx q[3];
rz(-0.91359514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.367584) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(-2.9665663) q[0];
rz(-1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(0.79644901) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7541695) q[0];
sx q[0];
rz(-1.8909374) q[0];
sx q[0];
rz(0.27227905) q[0];
rz(0.59993125) q[2];
sx q[2];
rz(-1.33432) q[2];
sx q[2];
rz(-1.8452871) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15115498) q[1];
sx q[1];
rz(-2.5271522) q[1];
sx q[1];
rz(-2.4555456) q[1];
rz(-0.50694179) q[3];
sx q[3];
rz(-1.4846509) q[3];
sx q[3];
rz(-2.994842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1548057) q[2];
sx q[2];
rz(-3.0974168) q[2];
sx q[2];
rz(3.0955443) q[2];
rz(-1.1151399) q[3];
sx q[3];
rz(-1.0067078) q[3];
sx q[3];
rz(-1.0435638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.030877) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(-0.16511551) q[0];
rz(1.1394399) q[1];
sx q[1];
rz(-2.20729) q[1];
sx q[1];
rz(-1.6868235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6580312) q[0];
sx q[0];
rz(-1.665977) q[0];
sx q[0];
rz(-0.24834085) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3619992) q[2];
sx q[2];
rz(-1.1399802) q[2];
sx q[2];
rz(-2.5519759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6858721) q[1];
sx q[1];
rz(-2.2903951) q[1];
sx q[1];
rz(-3.0830129) q[1];
x q[2];
rz(1.9971354) q[3];
sx q[3];
rz(-1.3071663) q[3];
sx q[3];
rz(0.80962777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7016865) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(2.9627724) q[2];
rz(1.5170826) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(-1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521859) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(1.7328316) q[0];
rz(0.57213712) q[1];
sx q[1];
rz(-1.1502384) q[1];
sx q[1];
rz(-0.22452721) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1869769) q[0];
sx q[0];
rz(-2.3659334) q[0];
sx q[0];
rz(-2.5903754) q[0];
x q[1];
rz(-1.2351745) q[2];
sx q[2];
rz(-0.76424176) q[2];
sx q[2];
rz(-2.3847876) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46426168) q[1];
sx q[1];
rz(-0.867093) q[1];
sx q[1];
rz(1.000885) q[1];
rz(-pi) q[2];
rz(1.270431) q[3];
sx q[3];
rz(-1.1704418) q[3];
sx q[3];
rz(2.3107236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81395927) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(-0.28780469) q[2];
rz(-0.52590251) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(-2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.41878685) q[0];
sx q[0];
rz(-0.64544183) q[0];
sx q[0];
rz(3.1268815) q[0];
rz(2.9171433) q[1];
sx q[1];
rz(-1.4597471) q[1];
sx q[1];
rz(0.86123484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5435641) q[0];
sx q[0];
rz(-0.042994067) q[0];
sx q[0];
rz(0.046648101) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9056896) q[2];
sx q[2];
rz(-1.5318724) q[2];
sx q[2];
rz(-0.30920497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3726681) q[1];
sx q[1];
rz(-2.3834991) q[1];
sx q[1];
rz(2.7102986) q[1];
x q[2];
rz(0.92361633) q[3];
sx q[3];
rz(-1.4391393) q[3];
sx q[3];
rz(-3.0513958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8031249) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(3.0306446) q[2];
rz(1.126193) q[3];
sx q[3];
rz(-1.5092827) q[3];
sx q[3];
rz(-0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647389) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(-1.2891084) q[0];
rz(-1.8076757) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(1.8905554) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.23952) q[0];
sx q[0];
rz(-0.7131359) q[0];
sx q[0];
rz(1.9608742) q[0];
rz(-2.6625113) q[2];
sx q[2];
rz(-0.63036171) q[2];
sx q[2];
rz(0.97409526) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2129647) q[1];
sx q[1];
rz(-1.3273094) q[1];
sx q[1];
rz(0.89694556) q[1];
rz(-pi) q[2];
rz(-0.9034702) q[3];
sx q[3];
rz(-0.37411896) q[3];
sx q[3];
rz(-1.4033409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.053607792) q[2];
sx q[2];
rz(-1.3259652) q[2];
sx q[2];
rz(1.717022) q[2];
rz(0.21814957) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(2.654259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80213928) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(-2.6735027) q[0];
rz(0.83174902) q[1];
sx q[1];
rz(-2.2305198) q[1];
sx q[1];
rz(-2.1591878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42280722) q[0];
sx q[0];
rz(-0.095346538) q[0];
sx q[0];
rz(-2.3991292) q[0];
rz(1.8562819) q[2];
sx q[2];
rz(-2.4997556) q[2];
sx q[2];
rz(1.8571719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0259798) q[1];
sx q[1];
rz(-1.8048688) q[1];
sx q[1];
rz(1.267141) q[1];
rz(-1.7782061) q[3];
sx q[3];
rz(-1.6230427) q[3];
sx q[3];
rz(2.047153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9437647) q[2];
sx q[2];
rz(-2.6524537) q[2];
sx q[2];
rz(0.93734199) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26405239) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(-2.383655) q[0];
rz(-2.0856048) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(-0.24076375) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7215181) q[0];
sx q[0];
rz(-1.2596247) q[0];
sx q[0];
rz(-0.56367409) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40177138) q[2];
sx q[2];
rz(-2.2770498) q[2];
sx q[2];
rz(2.1101949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5220753) q[1];
sx q[1];
rz(-0.9065827) q[1];
sx q[1];
rz(2.7665576) q[1];
x q[2];
rz(-1.9338172) q[3];
sx q[3];
rz(-1.8891834) q[3];
sx q[3];
rz(2.8045079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15464887) q[2];
sx q[2];
rz(-1.1661531) q[2];
sx q[2];
rz(-1.1692952) q[2];
rz(-0.5430921) q[3];
sx q[3];
rz(-1.8648059) q[3];
sx q[3];
rz(0.57524601) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6146522) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(2.8618983) q[0];
rz(2.0277297) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(-2.9934771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62559267) q[0];
sx q[0];
rz(-2.2315488) q[0];
sx q[0];
rz(1.8126677) q[0];
rz(-1.7314379) q[2];
sx q[2];
rz(-0.88692188) q[2];
sx q[2];
rz(1.4768254) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0462936) q[1];
sx q[1];
rz(-1.0314086) q[1];
sx q[1];
rz(2.1124243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.801228) q[3];
sx q[3];
rz(-0.90022579) q[3];
sx q[3];
rz(0.55055314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0327882) q[2];
sx q[2];
rz(-0.84454909) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(1.4480715) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(-1.1650803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447822) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(1.0983559) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(0.57571147) q[2];
sx q[2];
rz(-1.4288708) q[2];
sx q[2];
rz(2.5484011) q[2];
rz(1.4357812) q[3];
sx q[3];
rz(-1.5207471) q[3];
sx q[3];
rz(-0.38226939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
