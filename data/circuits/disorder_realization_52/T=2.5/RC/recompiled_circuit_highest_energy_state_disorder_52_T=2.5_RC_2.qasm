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
rz(-1.988451) q[0];
sx q[0];
rz(-2.3260131) q[0];
sx q[0];
rz(-2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(4.9795436) q[1];
sx q[1];
rz(10.995168) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2822097) q[0];
sx q[0];
rz(-2.4042685) q[0];
sx q[0];
rz(2.5928549) q[0];
rz(-1.4123795) q[2];
sx q[2];
rz(-0.40414256) q[2];
sx q[2];
rz(-2.2805285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4788073) q[1];
sx q[1];
rz(-2.7764455) q[1];
sx q[1];
rz(1.7313596) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.38405) q[3];
sx q[3];
rz(-2.8272326) q[3];
sx q[3];
rz(1.7160742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5300753) q[2];
sx q[2];
rz(-0.0081743058) q[2];
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
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663064) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(2.9735907) q[0];
rz(-3.1213144) q[1];
sx q[1];
rz(-2.8333277) q[1];
sx q[1];
rz(-1.6049989) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.080606) q[0];
sx q[0];
rz(-1.3362552) q[0];
sx q[0];
rz(1.3479665) q[0];
x q[1];
rz(-1.3290908) q[2];
sx q[2];
rz(-3.1311718) q[2];
sx q[2];
rz(-1.7872321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5913493) q[1];
sx q[1];
rz(-3.1367932) q[1];
sx q[1];
rz(1.8018434) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20949598) q[3];
sx q[3];
rz(-0.043238021) q[3];
sx q[3];
rz(-0.87856495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(1.4076642) q[2];
rz(-1.0450854) q[3];
sx q[3];
rz(-0.049523517) q[3];
sx q[3];
rz(2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.3097836) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(-0.56104863) q[0];
rz(2.8647515) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(-1.3078088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9914068) q[0];
sx q[0];
rz(-1.5298109) q[0];
sx q[0];
rz(-1.2738373) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.00066999992) q[2];
sx q[2];
rz(-3.1342034) q[2];
sx q[2];
rz(-1.0628029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4958932) q[1];
sx q[1];
rz(-0.99826854) q[1];
sx q[1];
rz(-3.0719724) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2482292) q[3];
sx q[3];
rz(-1.2033389) q[3];
sx q[3];
rz(-1.8830255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3448559) q[2];
sx q[2];
rz(-3.1414746) q[2];
sx q[2];
rz(-0.5893839) q[2];
rz(1.2423337) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(-3.1349365) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-3.1080918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5101227) q[0];
sx q[0];
rz(-2.2167248) q[0];
sx q[0];
rz(-2.9246246) q[0];
x q[1];
rz(0.037390402) q[2];
sx q[2];
rz(-3.0218736) q[2];
sx q[2];
rz(0.38116383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7280933) q[1];
sx q[1];
rz(-1.5809605) q[1];
sx q[1];
rz(-1.3018911) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12199282) q[3];
sx q[3];
rz(-1.4076403) q[3];
sx q[3];
rz(1.6900631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1091619) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(1.2700891) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(3.0884009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700767) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(2.694743) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(1.4131379) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2010445) q[0];
sx q[0];
rz(-1.5371635) q[0];
sx q[0];
rz(-1.3399009) q[0];
rz(-pi) q[1];
rz(1.9930219) q[2];
sx q[2];
rz(-1.3729501) q[2];
sx q[2];
rz(-2.5289218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6525938) q[1];
sx q[1];
rz(-1.5212544) q[1];
sx q[1];
rz(-0.045157305) q[1];
x q[2];
rz(1.8522315) q[3];
sx q[3];
rz(-1.7858943) q[3];
sx q[3];
rz(-0.18751442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83429217) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(0.49868047) q[2];
rz(0.57791609) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(-0.59670603) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96108288) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(-2.7951796) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(-2.3880889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0929209) q[0];
sx q[0];
rz(-1.3792097) q[0];
sx q[0];
rz(1.3956153) q[0];
x q[1];
rz(3.0664938) q[2];
sx q[2];
rz(-1.4595928) q[2];
sx q[2];
rz(2.5203343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9648887) q[1];
sx q[1];
rz(-2.2060925) q[1];
sx q[1];
rz(-0.52365644) q[1];
rz(-pi) q[2];
x q[2];
rz(0.092272225) q[3];
sx q[3];
rz(-1.408395) q[3];
sx q[3];
rz(-2.5889461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-0.0034905958) q[2];
sx q[2];
rz(-1.6174779) q[2];
rz(0.12024719) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(1.4606754) q[1];
sx q[1];
rz(-2.2082081) q[1];
sx q[1];
rz(-3.0642919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.59931) q[0];
sx q[0];
rz(-0.041754384) q[0];
sx q[0];
rz(1.6058654) q[0];
rz(-pi) q[1];
rz(1.5796698) q[2];
sx q[2];
rz(-1.5660145) q[2];
sx q[2];
rz(1.3591131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12880023) q[1];
sx q[1];
rz(-2.9523053) q[1];
sx q[1];
rz(0.38893338) q[1];
x q[2];
rz(0.12760136) q[3];
sx q[3];
rz(-1.3584879) q[3];
sx q[3];
rz(1.177976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7920502) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(0.95996094) q[2];
rz(-2.8121484) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(-0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195381) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(-3.0323113) q[0];
rz(-0.37336135) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.9111309) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7007992) q[0];
sx q[0];
rz(-0.95958086) q[0];
sx q[0];
rz(-0.81328765) q[0];
rz(1.7611165) q[2];
sx q[2];
rz(-1.3770104) q[2];
sx q[2];
rz(1.5908257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8671678) q[1];
sx q[1];
rz(-1.503957) q[1];
sx q[1];
rz(3.1301366) q[1];
x q[2];
rz(-0.30970311) q[3];
sx q[3];
rz(-2.1983357) q[3];
sx q[3];
rz(-2.5554339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(1.8129978) q[2];
rz(-1.396842) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(2.1197135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.063865572) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(2.5685837) q[0];
rz(2.833448) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(1.0073957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38978168) q[0];
sx q[0];
rz(-0.10669076) q[0];
sx q[0];
rz(1.5185028) q[0];
rz(-pi) q[1];
rz(-1.7083659) q[2];
sx q[2];
rz(-1.5273046) q[2];
sx q[2];
rz(1.5303591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5080164) q[1];
sx q[1];
rz(-3.0490082) q[1];
sx q[1];
rz(2.0256945) q[1];
rz(3.1281578) q[3];
sx q[3];
rz(-1.5406431) q[3];
sx q[3];
rz(-0.20520255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.31743) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(0.38995788) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214509) q[0];
sx q[0];
rz(-0.75103253) q[0];
sx q[0];
rz(-2.6556515) q[0];
rz(2.2700229) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(1.4942687) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594306) q[0];
sx q[0];
rz(-2.5100187) q[0];
sx q[0];
rz(-0.53379121) q[0];
x q[1];
rz(-2.5267692) q[2];
sx q[2];
rz(-1.6069001) q[2];
sx q[2];
rz(-1.5272128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4758265) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(-1.8880512) q[1];
x q[2];
rz(1.5288197) q[3];
sx q[3];
rz(-1.4577565) q[3];
sx q[3];
rz(2.6338995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-3.1097143) q[2];
rz(2.3632862) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42302172) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(-3.0146535) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(-1.610582) q[2];
sx q[2];
rz(-0.13977681) q[2];
sx q[2];
rz(0.20413414) q[2];
rz(-1.6416141) q[3];
sx q[3];
rz(-2.4657065) q[3];
sx q[3];
rz(0.55767725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
