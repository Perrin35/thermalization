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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(-1.0022751) q[1];
sx q[1];
rz(0.60751539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6014221) q[0];
sx q[0];
rz(-2.5937732) q[0];
sx q[0];
rz(1.1146077) q[0];
rz(2.0801699) q[2];
sx q[2];
rz(-1.4685838) q[2];
sx q[2];
rz(-1.4776023) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.059587) q[1];
sx q[1];
rz(-1.4114037) q[1];
sx q[1];
rz(-2.6578147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3729975) q[3];
sx q[3];
rz(-1.1524876) q[3];
sx q[3];
rz(-1.7278863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6174378) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(-1.9614356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2422159) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(-1.7346409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55757421) q[0];
sx q[0];
rz(-1.0913335) q[0];
sx q[0];
rz(2.8842501) q[0];
x q[1];
rz(1.9029877) q[2];
sx q[2];
rz(-1.536035) q[2];
sx q[2];
rz(1.8046276) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21995658) q[1];
sx q[1];
rz(-1.0590648) q[1];
sx q[1];
rz(-0.10351609) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5464393) q[3];
sx q[3];
rz(-0.66616026) q[3];
sx q[3];
rz(-2.7056138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0591639) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.3272237) q[2];
rz(1.0446769) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5752207) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(-0.4253934) q[0];
rz(-1.7644024) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(-1.707071) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4664513) q[0];
sx q[0];
rz(-1.4527307) q[0];
sx q[0];
rz(-2.4982014) q[0];
x q[1];
rz(-0.038952053) q[2];
sx q[2];
rz(-0.70868451) q[2];
sx q[2];
rz(1.5295636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19615281) q[1];
sx q[1];
rz(-1.9827794) q[1];
sx q[1];
rz(-0.71228551) q[1];
rz(-pi) q[2];
rz(1.1528963) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(-1.0583056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7003358) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(-1.8505992) q[2];
rz(1.7839606) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(-1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(-2.3714016) q[0];
rz(0.99984804) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(1.8005449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65085852) q[0];
sx q[0];
rz(-2.5528209) q[0];
sx q[0];
rz(-2.7522699) q[0];
rz(-pi) q[1];
rz(2.6821939) q[2];
sx q[2];
rz(-2.2041577) q[2];
sx q[2];
rz(2.6448696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.136678) q[1];
sx q[1];
rz(-1.305849) q[1];
sx q[1];
rz(-0.30989225) q[1];
rz(-pi) q[2];
x q[2];
rz(0.099319066) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.2885531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(-2.3731903) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(-1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860344) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(2.9146063) q[0];
rz(1.7598033) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(-1.3963799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67302746) q[0];
sx q[0];
rz(-1.919152) q[0];
sx q[0];
rz(-2.7569016) q[0];
x q[1];
rz(-2.7740741) q[2];
sx q[2];
rz(-2.705239) q[2];
sx q[2];
rz(2.9038716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3815617) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(-1.6595592) q[1];
rz(-pi) q[2];
rz(-0.030524039) q[3];
sx q[3];
rz(-2.3168193) q[3];
sx q[3];
rz(-1.33547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607894) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(0.18641557) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064998) q[0];
sx q[0];
rz(-0.71487521) q[0];
sx q[0];
rz(2.2660648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3720296) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(1.1793062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0708145) q[1];
sx q[1];
rz(-2.1970587) q[1];
sx q[1];
rz(2.8423944) q[1];
x q[2];
rz(1.0746434) q[3];
sx q[3];
rz(-2.0130139) q[3];
sx q[3];
rz(2.7520455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67941252) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(1.818044) q[2];
rz(1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(-2.0967261) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(-2.1814836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7564108) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(1.5002285) q[0];
x q[1];
rz(-1.5627925) q[2];
sx q[2];
rz(-1.2697313) q[2];
sx q[2];
rz(2.0640399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3030678) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(1.5938894) q[1];
rz(-pi) q[2];
rz(-2.3784901) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(-2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94379696) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(-1.5956399) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6087795) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(-3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.1688165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10442142) q[0];
sx q[0];
rz(-1.1544495) q[0];
sx q[0];
rz(-2.092931) q[0];
rz(0.76185267) q[2];
sx q[2];
rz(-1.5708062) q[2];
sx q[2];
rz(-0.61864432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(-1.8046677) q[1];
rz(0.49174546) q[3];
sx q[3];
rz(-1.5534658) q[3];
sx q[3];
rz(-0.31760707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.03269) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-2.4857793) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-0.18812215) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(1.6498097) q[0];
rz(2.1022294) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(-2.6731491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.428498) q[0];
sx q[0];
rz(-1.575043) q[0];
sx q[0];
rz(-0.18761401) q[0];
x q[1];
rz(1.2861112) q[2];
sx q[2];
rz(-1.1435025) q[2];
sx q[2];
rz(2.9197249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5891582) q[1];
sx q[1];
rz(-0.8048519) q[1];
sx q[1];
rz(2.9167152) q[1];
x q[2];
rz(-2.8599817) q[3];
sx q[3];
rz(-1.3002031) q[3];
sx q[3];
rz(2.5789346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(-1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(-1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(-2.0637276) q[0];
rz(-2.7742591) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(-1.2841388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76332247) q[0];
sx q[0];
rz(-2.5978932) q[0];
sx q[0];
rz(2.0732353) q[0];
rz(0.91725332) q[2];
sx q[2];
rz(-1.6628569) q[2];
sx q[2];
rz(-1.1600509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4176191) q[1];
sx q[1];
rz(-1.4129169) q[1];
sx q[1];
rz(-1.6579614) q[1];
rz(-pi) q[2];
rz(2.908091) q[3];
sx q[3];
rz(-2.3931008) q[3];
sx q[3];
rz(-0.12602636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1401356) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(-1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975288) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(-0.33959099) q[3];
sx q[3];
rz(-0.97151269) q[3];
sx q[3];
rz(1.1480939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
