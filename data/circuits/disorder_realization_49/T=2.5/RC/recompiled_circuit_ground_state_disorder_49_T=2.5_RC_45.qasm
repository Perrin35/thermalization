OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0030274) q[0];
sx q[0];
rz(4.0199036) q[0];
sx q[0];
rz(10.255788) q[0];
rz(2.5660958) q[1];
sx q[1];
rz(-0.73605186) q[1];
sx q[1];
rz(-0.73986685) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66798009) q[0];
sx q[0];
rz(-1.6428609) q[0];
sx q[0];
rz(2.9524809) q[0];
rz(-1.3492435) q[2];
sx q[2];
rz(-2.5622517) q[2];
sx q[2];
rz(0.18150961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4783096) q[1];
sx q[1];
rz(-1.4245207) q[1];
sx q[1];
rz(-2.0520794) q[1];
rz(0.44828592) q[3];
sx q[3];
rz(-0.42196754) q[3];
sx q[3];
rz(1.2039876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9326707) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(-0.33898655) q[2];
rz(0.50421667) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(1.1241815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(0.50952953) q[0];
rz(-1.5024028) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(2.1972919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4128542) q[0];
sx q[0];
rz(-1.7704238) q[0];
sx q[0];
rz(2.9908604) q[0];
rz(-pi) q[1];
rz(1.239907) q[2];
sx q[2];
rz(-2.9941786) q[2];
sx q[2];
rz(2.7380163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2067926) q[1];
sx q[1];
rz(-1.0369318) q[1];
sx q[1];
rz(2.4660048) q[1];
rz(-pi) q[2];
rz(-1.1491165) q[3];
sx q[3];
rz(-0.72315589) q[3];
sx q[3];
rz(0.92484091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3462191) q[2];
sx q[2];
rz(-1.5920762) q[2];
sx q[2];
rz(0.53595558) q[2];
rz(2.0761944) q[3];
sx q[3];
rz(-2.8841116) q[3];
sx q[3];
rz(-0.65142256) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(-0.73597062) q[0];
rz(-2.60587) q[1];
sx q[1];
rz(-1.842272) q[1];
sx q[1];
rz(-2.2419498) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9743703) q[0];
sx q[0];
rz(-1.456454) q[0];
sx q[0];
rz(0.13916441) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2399142) q[2];
sx q[2];
rz(-2.5106259) q[2];
sx q[2];
rz(-3.0598989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50541836) q[1];
sx q[1];
rz(-2.6405425) q[1];
sx q[1];
rz(1.4062642) q[1];
rz(-pi) q[2];
rz(2.8923762) q[3];
sx q[3];
rz(-1.2050516) q[3];
sx q[3];
rz(0.62376991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54124093) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(0.80292732) q[2];
rz(0.7488572) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.51711851) q[0];
sx q[0];
rz(-1.4117389) q[0];
sx q[0];
rz(0.13036048) q[0];
rz(-1.8761926) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(-0.1098384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30194908) q[0];
sx q[0];
rz(-0.57900864) q[0];
sx q[0];
rz(0.92424519) q[0];
rz(-pi) q[1];
rz(-2.7214987) q[2];
sx q[2];
rz(-2.5330628) q[2];
sx q[2];
rz(2.9836754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.500587) q[1];
sx q[1];
rz(-2.2595123) q[1];
sx q[1];
rz(-0.50488774) q[1];
rz(-pi) q[2];
rz(0.62105824) q[3];
sx q[3];
rz(-2.406236) q[3];
sx q[3];
rz(1.7572973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5023697) q[2];
sx q[2];
rz(-1.1797735) q[2];
sx q[2];
rz(0.42312527) q[2];
rz(1.0387748) q[3];
sx q[3];
rz(-0.3813425) q[3];
sx q[3];
rz(1.0264621) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410685) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-0.50222841) q[1];
sx q[1];
rz(-2.8531029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5768439) q[0];
sx q[0];
rz(-1.1001407) q[0];
sx q[0];
rz(0.0708357) q[0];
x q[1];
rz(1.1568858) q[2];
sx q[2];
rz(-0.26662808) q[2];
sx q[2];
rz(-1.8265343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.57196778) q[1];
sx q[1];
rz(-1.5020084) q[1];
sx q[1];
rz(1.6602762) q[1];
rz(0.27933018) q[3];
sx q[3];
rz(-2.175176) q[3];
sx q[3];
rz(1.6681748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5491817) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(-0.16331095) q[2];
rz(1.1641938) q[3];
sx q[3];
rz(-0.67622286) q[3];
sx q[3];
rz(-1.7827079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613681) q[0];
sx q[0];
rz(-0.017266406) q[0];
sx q[0];
rz(-1.226271) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-0.93003479) q[1];
sx q[1];
rz(0.61940449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7295655) q[0];
sx q[0];
rz(-0.71456996) q[0];
sx q[0];
rz(-0.78972915) q[0];
x q[1];
rz(1.7707497) q[2];
sx q[2];
rz(-1.8189578) q[2];
sx q[2];
rz(-1.8979817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0521026) q[1];
sx q[1];
rz(-1.5137496) q[1];
sx q[1];
rz(1.4228348) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9483637) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.61249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9161393) q[2];
sx q[2];
rz(-1.4832387) q[2];
sx q[2];
rz(-2.7062866) q[2];
rz(0.12886038) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(2.784957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72047609) q[0];
sx q[0];
rz(-2.5832376) q[0];
sx q[0];
rz(-2.5893353) q[0];
rz(2.5241959) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(1.5026106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4308145) q[0];
sx q[0];
rz(-2.7511247) q[0];
sx q[0];
rz(-1.4560519) q[0];
rz(-2.0755956) q[2];
sx q[2];
rz(-1.4420954) q[2];
sx q[2];
rz(-1.7247049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0459111) q[1];
sx q[1];
rz(-2.0740182) q[1];
sx q[1];
rz(2.9614425) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076528744) q[3];
sx q[3];
rz(-0.26004836) q[3];
sx q[3];
rz(-1.6116774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5611394) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(-1.0882161) q[2];
rz(3.1214118) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(-1.5599686) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35767558) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(2.1141323) q[0];
rz(2.0383535) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(1.9267513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2526557) q[0];
sx q[0];
rz(-1.2644469) q[0];
sx q[0];
rz(-3.0333748) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80963366) q[2];
sx q[2];
rz(-0.61695951) q[2];
sx q[2];
rz(0.074658565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.852927) q[1];
sx q[1];
rz(-0.92354362) q[1];
sx q[1];
rz(0.44026466) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5276466) q[3];
sx q[3];
rz(-1.1505732) q[3];
sx q[3];
rz(2.9752034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2601629) q[2];
sx q[2];
rz(-1.9378928) q[2];
sx q[2];
rz(-1.0677968) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-0.46025899) q[3];
sx q[3];
rz(2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(-2.886014) q[0];
rz(-0.74343395) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(1.5240634) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8403064) q[0];
sx q[0];
rz(-0.69536007) q[0];
sx q[0];
rz(-0.82454234) q[0];
rz(-0.41478283) q[2];
sx q[2];
rz(-2.4786502) q[2];
sx q[2];
rz(-0.65185968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8614304) q[1];
sx q[1];
rz(-1.4373684) q[1];
sx q[1];
rz(0.89806865) q[1];
rz(2.8761707) q[3];
sx q[3];
rz(-0.94688225) q[3];
sx q[3];
rz(1.8056295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0301547) q[2];
sx q[2];
rz(-1.7058426) q[2];
sx q[2];
rz(-1.5343522) q[2];
rz(2.0700908) q[3];
sx q[3];
rz(-1.6000308) q[3];
sx q[3];
rz(-1.1736386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(0.80825835) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.6917276) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6351371) q[0];
sx q[0];
rz(-2.8525557) q[0];
sx q[0];
rz(1.6063547) q[0];
rz(-2.6053455) q[2];
sx q[2];
rz(-1.9384346) q[2];
sx q[2];
rz(0.73208955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7781642) q[1];
sx q[1];
rz(-1.320667) q[1];
sx q[1];
rz(0.56285897) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4042481) q[3];
sx q[3];
rz(-0.72190815) q[3];
sx q[3];
rz(-1.192821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(-0.9355363) q[2];
rz(-1.4714636) q[3];
sx q[3];
rz(-1.6664489) q[3];
sx q[3];
rz(-2.0475625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6846631) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(-2.6928071) q[1];
sx q[1];
rz(-1.9955336) q[1];
sx q[1];
rz(1.21036) q[1];
rz(1.2032897) q[2];
sx q[2];
rz(-1.5264282) q[2];
sx q[2];
rz(-2.1849968) q[2];
rz(2.5815677) q[3];
sx q[3];
rz(-1.2784132) q[3];
sx q[3];
rz(-1.9412435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
