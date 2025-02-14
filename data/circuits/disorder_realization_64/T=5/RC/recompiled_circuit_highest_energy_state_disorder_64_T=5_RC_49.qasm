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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(-0.91941961) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(3.0923831) q[1];
sx q[1];
rz(10.13848) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714234) q[0];
sx q[0];
rz(-1.4235625) q[0];
sx q[0];
rz(1.8482313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.117341) q[2];
sx q[2];
rz(-0.14278097) q[2];
sx q[2];
rz(2.873024) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9956869) q[1];
sx q[1];
rz(-0.45098454) q[1];
sx q[1];
rz(-3.1383951) q[1];
rz(-2.0421175) q[3];
sx q[3];
rz(-1.5402392) q[3];
sx q[3];
rz(1.2191284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.087622341) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-2.7719882) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(0.37231529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0358148) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(-1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(-0.47164741) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99804634) q[0];
sx q[0];
rz(-0.61206383) q[0];
sx q[0];
rz(1.7630811) q[0];
rz(-0.94309965) q[2];
sx q[2];
rz(-2.1617956) q[2];
sx q[2];
rz(-2.2636388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14395606) q[1];
sx q[1];
rz(-2.156971) q[1];
sx q[1];
rz(-2.2464941) q[1];
x q[2];
rz(-2.2120958) q[3];
sx q[3];
rz(-0.62493443) q[3];
sx q[3];
rz(0.24111023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(2.0527077) q[2];
rz(-2.0770843) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(0.3256807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(0.017070008) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(2.5753218) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44646548) q[0];
sx q[0];
rz(-0.74120616) q[0];
sx q[0];
rz(-1.5873853) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0300358) q[2];
sx q[2];
rz(-0.93387253) q[2];
sx q[2];
rz(0.43535296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10668392) q[1];
sx q[1];
rz(-1.3843754) q[1];
sx q[1];
rz(1.4843462) q[1];
rz(1.6376466) q[3];
sx q[3];
rz(-1.8354356) q[3];
sx q[3];
rz(-2.8422249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8038586) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(3.1324978) q[2];
rz(3.005262) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71753865) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(2.9754382) q[0];
rz(-1.1135788) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(0.16214935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067569392) q[0];
sx q[0];
rz(-0.20299202) q[0];
sx q[0];
rz(2.0503886) q[0];
x q[1];
rz(-0.040302868) q[2];
sx q[2];
rz(-2.5922814) q[2];
sx q[2];
rz(-1.2072762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2723586) q[1];
sx q[1];
rz(-1.2125535) q[1];
sx q[1];
rz(-0.76613249) q[1];
rz(-pi) q[2];
rz(1.9982073) q[3];
sx q[3];
rz(-0.95164883) q[3];
sx q[3];
rz(2.038508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98264155) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(-2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.270179) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(-0.80108368) q[0];
rz(-0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(-2.1606826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50612569) q[0];
sx q[0];
rz(-0.46403316) q[0];
sx q[0];
rz(2.6543174) q[0];
rz(2.7655914) q[2];
sx q[2];
rz(-2.0472976) q[2];
sx q[2];
rz(-2.8636572) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0438761) q[1];
sx q[1];
rz(-1.6642769) q[1];
sx q[1];
rz(3.082117) q[1];
rz(0.24215908) q[3];
sx q[3];
rz(-2.1091614) q[3];
sx q[3];
rz(-0.70555609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.058978) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(-0.8832461) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.9969295) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(0.99217478) q[0];
rz(-2.3400173) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(-1.4206402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0376037) q[0];
sx q[0];
rz(-2.1891356) q[0];
sx q[0];
rz(1.9644587) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0384485) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(-0.20648512) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0040505) q[1];
sx q[1];
rz(-2.1909449) q[1];
sx q[1];
rz(-1.1637079) q[1];
rz(1.4735953) q[3];
sx q[3];
rz(-1.5335113) q[3];
sx q[3];
rz(2.4241222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27986032) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2456197) q[0];
sx q[0];
rz(-0.27181044) q[0];
sx q[0];
rz(-2.8017092) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(-1.9702912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1780258) q[0];
sx q[0];
rz(-1.8676571) q[0];
sx q[0];
rz(1.5328477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7843855) q[2];
sx q[2];
rz(-1.3554975) q[2];
sx q[2];
rz(0.70206308) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.724546) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.538365) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62452353) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(-2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0168212) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(-0.67145124) q[2];
rz(-2.6050383) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-2.0898537) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(-0.89624727) q[0];
rz(-0.25262901) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3327961) q[0];
sx q[0];
rz(-2.0967297) q[0];
sx q[0];
rz(1.512708) q[0];
rz(-pi) q[1];
rz(-1.7017548) q[2];
sx q[2];
rz(-0.90990657) q[2];
sx q[2];
rz(-2.5051528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72880367) q[1];
sx q[1];
rz(-1.1096242) q[1];
sx q[1];
rz(-2.9739266) q[1];
x q[2];
rz(-0.75783055) q[3];
sx q[3];
rz(-1.5327454) q[3];
sx q[3];
rz(-1.1421957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62319055) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(-1.5625578) q[2];
rz(-1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(-0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(-2.7516464) q[0];
rz(0.82603106) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(1.0521851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9814126) q[0];
sx q[0];
rz(-2.3773068) q[0];
sx q[0];
rz(2.0212964) q[0];
rz(-2.7446943) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(-1.4805178) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6893951) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(-2.2188975) q[1];
rz(-pi) q[2];
rz(-2.1901699) q[3];
sx q[3];
rz(-2.4149899) q[3];
sx q[3];
rz(-0.28765905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0612001) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(0.52919394) q[2];
rz(-0.75096327) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(-2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879409) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.1645114) q[0];
rz(-1.8544082) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(-0.50416344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8541478) q[0];
sx q[0];
rz(-1.7896255) q[0];
sx q[0];
rz(-1.8356505) q[0];
x q[1];
rz(1.5515366) q[2];
sx q[2];
rz(-1.1941205) q[2];
sx q[2];
rz(-2.6835359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87889844) q[1];
sx q[1];
rz(-0.32890186) q[1];
sx q[1];
rz(2.1961588) q[1];
x q[2];
rz(-0.21728094) q[3];
sx q[3];
rz(-2.2437527) q[3];
sx q[3];
rz(2.6126044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9265499) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(-3.1070993) q[2];
rz(-1.6132332) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(-2.8300986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.9813949) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(1.002671) q[1];
sx q[1];
rz(-1.4613338) q[1];
sx q[1];
rz(2.7636539) q[1];
rz(-0.38242292) q[2];
sx q[2];
rz(-2.6527846) q[2];
sx q[2];
rz(-2.6058886) q[2];
rz(-2.814449) q[3];
sx q[3];
rz(-1.5113304) q[3];
sx q[3];
rz(-2.3979901) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
