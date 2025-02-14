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
rz(0.74066585) q[0];
sx q[0];
rz(-1.292115) q[0];
sx q[0];
rz(2.3722755) q[0];
rz(-1.5462592) q[1];
sx q[1];
rz(-0.21472628) q[1];
sx q[1];
rz(-2.5200342) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21827605) q[0];
sx q[0];
rz(-2.0686604) q[0];
sx q[0];
rz(2.455869) q[0];
x q[1];
rz(1.1530959) q[2];
sx q[2];
rz(-1.5868895) q[2];
sx q[2];
rz(-2.180445) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5468419) q[1];
sx q[1];
rz(-1.8545281) q[1];
sx q[1];
rz(-2.3111514) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4906472) q[3];
sx q[3];
rz(-1.7043796) q[3];
sx q[3];
rz(-2.7770673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8605211) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(0.93112913) q[2];
rz(2.9562601) q[3];
sx q[3];
rz(-2.1087746) q[3];
sx q[3];
rz(-1.4144271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.0952048) q[0];
sx q[0];
rz(-2.8075908) q[0];
sx q[0];
rz(-2.1686676) q[0];
rz(-2.3880549) q[1];
sx q[1];
rz(-0.85644186) q[1];
sx q[1];
rz(-0.10261745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0282183) q[0];
sx q[0];
rz(-1.4732142) q[0];
sx q[0];
rz(-0.085208864) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4717536) q[2];
sx q[2];
rz(-1.4375039) q[2];
sx q[2];
rz(-1.006056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2778138) q[1];
sx q[1];
rz(-1.4881388) q[1];
sx q[1];
rz(-1.7277) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1903193) q[3];
sx q[3];
rz(-0.35934908) q[3];
sx q[3];
rz(2.7294266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40833452) q[2];
sx q[2];
rz(-1.7034986) q[2];
sx q[2];
rz(-2.5737838) q[2];
rz(-2.7693977) q[3];
sx q[3];
rz(-1.3293543) q[3];
sx q[3];
rz(1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55783015) q[0];
sx q[0];
rz(-0.74302858) q[0];
sx q[0];
rz(1.8999735) q[0];
rz(0.81635967) q[1];
sx q[1];
rz(-2.7528449) q[1];
sx q[1];
rz(0.98966086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802487) q[0];
sx q[0];
rz(-1.3575866) q[0];
sx q[0];
rz(0.38589392) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7317361) q[2];
sx q[2];
rz(-1.66203) q[2];
sx q[2];
rz(-0.021878069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.951713) q[1];
sx q[1];
rz(-2.6939309) q[1];
sx q[1];
rz(0.42231222) q[1];
rz(-pi) q[2];
rz(-1.3157522) q[3];
sx q[3];
rz(-2.7560184) q[3];
sx q[3];
rz(-2.5346867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.376754) q[2];
sx q[2];
rz(-1.9003442) q[2];
sx q[2];
rz(1.6468916) q[2];
rz(0.74690789) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(2.1493105) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.295739) q[0];
sx q[0];
rz(-0.60765147) q[0];
sx q[0];
rz(0.94957748) q[0];
rz(-1.9751366) q[1];
sx q[1];
rz(-1.2587073) q[1];
sx q[1];
rz(-0.20225254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40914772) q[0];
sx q[0];
rz(-1.5209001) q[0];
sx q[0];
rz(1.6014535) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4030898) q[2];
sx q[2];
rz(-2.2532092) q[2];
sx q[2];
rz(-0.81425999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65860393) q[1];
sx q[1];
rz(-2.3301396) q[1];
sx q[1];
rz(-0.8354018) q[1];
x q[2];
rz(-2.3712986) q[3];
sx q[3];
rz(-2.3912734) q[3];
sx q[3];
rz(1.6495153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4919081) q[2];
sx q[2];
rz(-1.4338355) q[2];
sx q[2];
rz(1.3362308) q[2];
rz(-2.6797471) q[3];
sx q[3];
rz(-2.0838085) q[3];
sx q[3];
rz(-2.113078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.860054) q[0];
sx q[0];
rz(-3.0335732) q[0];
sx q[0];
rz(-0.39529133) q[0];
rz(-0.95834243) q[1];
sx q[1];
rz(-1.7610565) q[1];
sx q[1];
rz(0.39452943) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4636587) q[0];
sx q[0];
rz(-1.5531085) q[0];
sx q[0];
rz(-1.836548) q[0];
rz(2.5960693) q[2];
sx q[2];
rz(-0.50456556) q[2];
sx q[2];
rz(-1.2009753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1198152) q[1];
sx q[1];
rz(-1.2477326) q[1];
sx q[1];
rz(-2.3141239) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8452931) q[3];
sx q[3];
rz(-0.63929108) q[3];
sx q[3];
rz(-0.019867181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40256527) q[2];
sx q[2];
rz(-1.5274916) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(1.3022276) q[3];
sx q[3];
rz(-3.0890833) q[3];
sx q[3];
rz(-0.70560613) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5810982) q[0];
sx q[0];
rz(-2.5769233) q[0];
sx q[0];
rz(2.1288921) q[0];
rz(-0.46214354) q[1];
sx q[1];
rz(-1.4116762) q[1];
sx q[1];
rz(-3.0949458) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.843294) q[0];
sx q[0];
rz(-3.0180535) q[0];
sx q[0];
rz(-1.4840829) q[0];
rz(-1.2744802) q[2];
sx q[2];
rz(-1.8159397) q[2];
sx q[2];
rz(-2.5398538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.093864249) q[1];
sx q[1];
rz(-1.0055491) q[1];
sx q[1];
rz(-0.39611343) q[1];
x q[2];
rz(1.9479475) q[3];
sx q[3];
rz(-1.5483678) q[3];
sx q[3];
rz(-0.42322657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3095653) q[2];
sx q[2];
rz(-0.088531606) q[2];
sx q[2];
rz(-0.95477611) q[2];
rz(-0.66983062) q[3];
sx q[3];
rz(-1.1845183) q[3];
sx q[3];
rz(-2.7772016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1466115) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(-1.1501508) q[0];
rz(3.0448044) q[1];
sx q[1];
rz(-2.0014747) q[1];
sx q[1];
rz(1.4801625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.262275) q[0];
sx q[0];
rz(-1.9497383) q[0];
sx q[0];
rz(-2.3728601) q[0];
x q[1];
rz(1.3315013) q[2];
sx q[2];
rz(-1.8836438) q[2];
sx q[2];
rz(0.32101908) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5642359) q[1];
sx q[1];
rz(-0.32315394) q[1];
sx q[1];
rz(1.1250457) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0050762) q[3];
sx q[3];
rz(-1.7107114) q[3];
sx q[3];
rz(-0.29014465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73411554) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(2.5131098) q[2];
rz(-0.34144044) q[3];
sx q[3];
rz(-0.98613685) q[3];
sx q[3];
rz(-2.3054874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2899365) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(-0.063902721) q[0];
rz(-2.5996161) q[1];
sx q[1];
rz(-1.130645) q[1];
sx q[1];
rz(-0.37905395) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0444894) q[0];
sx q[0];
rz(-1.3987204) q[0];
sx q[0];
rz(1.340036) q[0];
rz(-pi) q[1];
rz(-3.0115773) q[2];
sx q[2];
rz(-0.81323821) q[2];
sx q[2];
rz(2.7521135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8547092) q[1];
sx q[1];
rz(-1.5552551) q[1];
sx q[1];
rz(-2.5980224) q[1];
rz(-0.42551252) q[3];
sx q[3];
rz(-2.1149016) q[3];
sx q[3];
rz(-1.5489006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9717676) q[2];
sx q[2];
rz(-1.1173893) q[2];
sx q[2];
rz(2.8525412) q[2];
rz(1.0158094) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(-0.22739205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.10580258) q[0];
sx q[0];
rz(-2.7547024) q[0];
sx q[0];
rz(-2.8111358) q[0];
rz(0.48078787) q[1];
sx q[1];
rz(-1.582522) q[1];
sx q[1];
rz(2.6343583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.071166) q[0];
sx q[0];
rz(-0.88001635) q[0];
sx q[0];
rz(0.69657495) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4749996) q[2];
sx q[2];
rz(-2.5799516) q[2];
sx q[2];
rz(3.0009342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54785905) q[1];
sx q[1];
rz(-2.6630109) q[1];
sx q[1];
rz(-1.0793769) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6303924) q[3];
sx q[3];
rz(-0.97138849) q[3];
sx q[3];
rz(0.58356482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5656188) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-2.7536075) q[2];
rz(-3.0787789) q[3];
sx q[3];
rz(-2.4020436) q[3];
sx q[3];
rz(1.1886103) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0509309) q[0];
sx q[0];
rz(-0.020314038) q[0];
sx q[0];
rz(0.055572979) q[0];
rz(0.22124258) q[1];
sx q[1];
rz(-2.1317) q[1];
sx q[1];
rz(-1.6411068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.198198) q[0];
sx q[0];
rz(-2.072081) q[0];
sx q[0];
rz(2.8265165) q[0];
rz(-0.13672853) q[2];
sx q[2];
rz(-1.6199582) q[2];
sx q[2];
rz(1.8135742) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21162578) q[1];
sx q[1];
rz(-1.1911508) q[1];
sx q[1];
rz(-0.34348653) q[1];
rz(-pi) q[2];
rz(0.9110578) q[3];
sx q[3];
rz(-2.3039936) q[3];
sx q[3];
rz(-0.78105951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9556433) q[2];
sx q[2];
rz(-2.9538437) q[2];
sx q[2];
rz(-2.032568) q[2];
rz(-0.016131314) q[3];
sx q[3];
rz(-1.4343836) q[3];
sx q[3];
rz(2.6764638) q[3];
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
rz(-pi/2) q[3];
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
rz(-1.7598509) q[0];
sx q[0];
rz(-1.009059) q[0];
sx q[0];
rz(-0.41643634) q[0];
rz(-2.3460559) q[1];
sx q[1];
rz(-1.3858613) q[1];
sx q[1];
rz(-0.081079986) q[1];
rz(1.9151081) q[2];
sx q[2];
rz(-1.7034265) q[2];
sx q[2];
rz(-1.8272057) q[2];
rz(2.1063188) q[3];
sx q[3];
rz(-2.2761619) q[3];
sx q[3];
rz(-0.63009562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
