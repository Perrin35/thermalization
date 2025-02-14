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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(3.1290913) q[1];
sx q[1];
rz(-1.8379509) q[1];
sx q[1];
rz(1.571203) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16906315) q[0];
sx q[0];
rz(-2.1816945) q[0];
sx q[0];
rz(-1.1283895) q[0];
x q[1];
rz(1.4123795) q[2];
sx q[2];
rz(-2.7374501) q[2];
sx q[2];
rz(0.86106419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.058152288) q[1];
sx q[1];
rz(-1.5136763) q[1];
sx q[1];
rz(-1.9316462) q[1];
rz(-pi) q[2];
rz(1.38405) q[3];
sx q[3];
rz(-0.31436008) q[3];
sx q[3];
rz(-1.4255185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(2.7008936) q[2];
rz(-3.0618482) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(2.0174446) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663064) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(0.16800198) q[0];
rz(-3.1213144) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(-1.5365938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0609866) q[0];
sx q[0];
rz(-1.3362552) q[0];
sx q[0];
rz(1.3479665) q[0];
x q[1];
rz(0.0024944024) q[2];
sx q[2];
rz(-1.5606784) q[2];
sx q[2];
rz(1.5960787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3602996) q[1];
sx q[1];
rz(-1.5754682) q[1];
sx q[1];
rz(0.0010990573) q[1];
rz(-pi) q[2];
rz(2.9320967) q[3];
sx q[3];
rz(-3.0983546) q[3];
sx q[3];
rz(0.87856495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-2.2296843) q[2];
sx q[2];
rz(-1.4076642) q[2];
rz(-2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8318091) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(0.56104863) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(1.8337839) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1501859) q[0];
sx q[0];
rz(-1.5298109) q[0];
sx q[0];
rz(1.2738373) q[0];
rz(-pi) q[1];
rz(-1.5707914) q[2];
sx q[2];
rz(-1.563407) q[2];
sx q[2];
rz(-2.0781197) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0289291) q[1];
sx q[1];
rz(-1.512292) q[1];
sx q[1];
rz(0.99716352) q[1];
x q[2];
rz(-2.2482292) q[3];
sx q[3];
rz(-1.2033389) q[3];
sx q[3];
rz(1.2585672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3448559) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(2.5522088) q[2];
rz(-1.899259) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(-1.7813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347443) q[0];
sx q[0];
rz(-0.51181781) q[0];
sx q[0];
rz(1.7929329) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-0.033500813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.861167) q[0];
sx q[0];
rz(-0.67641947) q[0];
sx q[0];
rz(1.8489714) q[0];
rz(1.5752931) q[2];
sx q[2];
rz(-1.6904313) q[2];
sx q[2];
rz(-0.41882354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4134994) q[1];
sx q[1];
rz(-1.5809605) q[1];
sx q[1];
rz(-1.8397016) q[1];
x q[2];
rz(0.12199282) q[3];
sx q[3];
rz(-1.7339524) q[3];
sx q[3];
rz(1.4515296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1091619) q[2];
sx q[2];
rz(-3.1350632) q[2];
sx q[2];
rz(2.8429441) q[2];
rz(-1.2700891) q[3];
sx q[3];
rz(-0.01615571) q[3];
sx q[3];
rz(-0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.771516) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(-2.694743) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-3.0797854) q[1];
sx q[1];
rz(-1.4131379) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639358) q[0];
sx q[0];
rz(-1.8015588) q[0];
sx q[0];
rz(0.034548977) q[0];
x q[1];
rz(-1.1485708) q[2];
sx q[2];
rz(-1.3729501) q[2];
sx q[2];
rz(0.6126709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2286716) q[1];
sx q[1];
rz(-3.0745709) q[1];
sx q[1];
rz(-2.3093501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90425332) q[3];
sx q[3];
rz(-2.7891141) q[3];
sx q[3];
rz(-2.3942687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3073005) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(-2.5636766) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(-0.59670603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96108288) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(-2.7951796) q[0];
rz(-2.5397884) q[1];
sx q[1];
rz(-1.5608984) q[1];
sx q[1];
rz(0.7535038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6974073) q[0];
sx q[0];
rz(-1.7427398) q[0];
sx q[0];
rz(-0.19449046) q[0];
x q[1];
rz(0.97917231) q[2];
sx q[2];
rz(-0.1340999) q[2];
sx q[2];
rz(0.025452415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9484472) q[1];
sx q[1];
rz(-2.3420534) q[1];
sx q[1];
rz(0.97481291) q[1];
x q[2];
rz(-2.0830882) q[3];
sx q[3];
rz(-2.9550094) q[3];
sx q[3];
rz(-0.032840289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(-1.5241148) q[2];
rz(3.0213455) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8728802) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(-2.9203316) q[0];
rz(-1.6809173) q[1];
sx q[1];
rz(-2.2082081) q[1];
sx q[1];
rz(-3.0642919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54228264) q[0];
sx q[0];
rz(-3.0998383) q[0];
sx q[0];
rz(1.6058654) q[0];
rz(-pi) q[1];
rz(-0.004782025) q[2];
sx q[2];
rz(-1.5796697) q[2];
sx q[2];
rz(-2.9298669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52409808) q[1];
sx q[1];
rz(-1.3957983) q[1];
sx q[1];
rz(1.4982759) q[1];
rz(-pi) q[2];
rz(-1.3568001) q[3];
sx q[3];
rz(-1.6955175) q[3];
sx q[3];
rz(-0.41984841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7920502) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(2.1816317) q[2];
rz(-2.8121484) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(-0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9462117) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(-0.10928133) q[0];
rz(0.37336135) q[1];
sx q[1];
rz(-2.3318718) q[1];
sx q[1];
rz(1.2304617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4407935) q[0];
sx q[0];
rz(-0.95958086) q[0];
sx q[0];
rz(0.81328765) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9443342) q[2];
sx q[2];
rz(-1.3840809) q[2];
sx q[2];
rz(0.017053617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27442481) q[1];
sx q[1];
rz(-1.503957) q[1];
sx q[1];
rz(3.1301366) q[1];
x q[2];
rz(-0.91992141) q[3];
sx q[3];
rz(-1.3215142) q[3];
sx q[3];
rz(1.1703619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(1.8129978) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(-1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0777271) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(-2.5685837) q[0];
rz(-2.833448) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(2.134197) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44237374) q[0];
sx q[0];
rz(-1.464252) q[0];
sx q[0];
rz(-3.1359948) q[0];
x q[1];
rz(3.0976866) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(-0.034417987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7511661) q[1];
sx q[1];
rz(-1.6114283) q[1];
sx q[1];
rz(-1.4875814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6009523) q[3];
sx q[3];
rz(-1.584225) q[3];
sx q[3];
rz(1.7755938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.31743) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(-0.38995788) q[2];
rz(-3.0725078) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(2.8100584) q[3];
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
rz(0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(0.48594117) q[0];
rz(0.87156975) q[1];
sx q[1];
rz(-1.8337245) q[1];
sx q[1];
rz(-1.647324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594306) q[0];
sx q[0];
rz(-2.5100187) q[0];
sx q[0];
rz(-0.53379121) q[0];
rz(-pi) q[1];
rz(-1.6149841) q[2];
sx q[2];
rz(-0.95643294) q[2];
sx q[2];
rz(3.0725266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68077129) q[1];
sx q[1];
rz(-2.6893565) q[1];
sx q[1];
rz(-2.3994563) q[1];
x q[2];
rz(-0.11313862) q[3];
sx q[3];
rz(-1.5290878) q[3];
sx q[3];
rz(2.0737518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(0.031878397) q[2];
rz(-2.3632862) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(-0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185709) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(0.12693916) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(1.610582) q[2];
sx q[2];
rz(-3.0018158) q[2];
sx q[2];
rz(-2.9374585) q[2];
rz(3.0849135) q[3];
sx q[3];
rz(-0.89691848) q[3];
sx q[3];
rz(0.46700029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
