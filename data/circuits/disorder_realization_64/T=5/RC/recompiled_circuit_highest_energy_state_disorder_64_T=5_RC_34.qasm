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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714234) q[0];
sx q[0];
rz(-1.4235625) q[0];
sx q[0];
rz(-1.2933613) q[0];
rz(2.0242516) q[2];
sx q[2];
rz(-0.14278097) q[2];
sx q[2];
rz(2.873024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9956869) q[1];
sx q[1];
rz(-0.45098454) q[1];
sx q[1];
rz(3.1383951) q[1];
rz(-pi) q[2];
rz(-1.0994751) q[3];
sx q[3];
rz(-1.5402392) q[3];
sx q[3];
rz(-1.2191284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.087622341) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(-2.7719882) q[2];
rz(1.9965648) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10577781) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(0.57304397) q[0];
rz(-1.0967968) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(0.47164741) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9100138) q[0];
sx q[0];
rz(-2.1699804) q[0];
sx q[0];
rz(-3.0082361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4493982) q[2];
sx q[2];
rz(-2.080215) q[2];
sx q[2];
rz(0.30864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9976366) q[1];
sx q[1];
rz(-2.156971) q[1];
sx q[1];
rz(2.2464941) q[1];
rz(-0.40741323) q[3];
sx q[3];
rz(-1.0828567) q[3];
sx q[3];
rz(2.1564302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7219217) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(2.0527077) q[2];
rz(2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(0.3256807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0786781) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(-0.017070008) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(0.56627083) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46895263) q[0];
sx q[0];
rz(-0.8297161) q[0];
sx q[0];
rz(-3.1264105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60795075) q[2];
sx q[2];
rz(-2.3312116) q[2];
sx q[2];
rz(-1.9160401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6614187) q[1];
sx q[1];
rz(-1.4858477) q[1];
sx q[1];
rz(2.9544891) q[1];
rz(-2.8763883) q[3];
sx q[3];
rz(-1.6353161) q[3];
sx q[3];
rz(1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.337734) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(3.1324978) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.9280619) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71753865) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(-0.16615443) q[0];
rz(-1.1135788) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(0.16214935) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067569392) q[0];
sx q[0];
rz(-0.20299202) q[0];
sx q[0];
rz(2.0503886) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5926431) q[2];
sx q[2];
rz(-1.5918343) q[2];
sx q[2];
rz(0.32914135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6238193) q[1];
sx q[1];
rz(-2.2775067) q[1];
sx q[1];
rz(1.0916187) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.477263) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(2.9323102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98264155) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(-0.52784935) q[2];
rz(0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.270179) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(0.98091006) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.635467) q[0];
sx q[0];
rz(-0.46403316) q[0];
sx q[0];
rz(-2.6543174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0641685) q[2];
sx q[2];
rz(-1.9031823) q[2];
sx q[2];
rz(1.6695963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0438761) q[1];
sx q[1];
rz(-1.6642769) q[1];
sx q[1];
rz(-3.082117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.122287) q[3];
sx q[3];
rz(-1.3634342) q[3];
sx q[3];
rz(-2.1503748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(0.8832461) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.14466318) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-0.99217478) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(1.7209524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29692263) q[0];
sx q[0];
rz(-1.2529182) q[0];
sx q[0];
rz(2.4852089) q[0];
x q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(-2.9351075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7758549) q[1];
sx q[1];
rz(-2.4147662) q[1];
sx q[1];
rz(0.50623843) q[1];
rz(-pi) q[2];
rz(1.6679974) q[3];
sx q[3];
rz(-1.6080813) q[3];
sx q[3];
rz(-0.71747045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27986032) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(-3.1308657) q[3];
sx q[3];
rz(-0.74672943) q[3];
sx q[3];
rz(-2.2906394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.895973) q[0];
sx q[0];
rz(-0.27181044) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(1.3019568) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(-1.9702912) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1780258) q[0];
sx q[0];
rz(-1.2739355) q[0];
sx q[0];
rz(1.6087449) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22014648) q[2];
sx q[2];
rz(-1.7793806) q[2];
sx q[2];
rz(2.2265547) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3893632) q[1];
sx q[1];
rz(-3.1023355) q[1];
sx q[1];
rz(-2.1696349) q[1];
rz(-pi) q[2];
rz(2.9006935) q[3];
sx q[3];
rz(-2.5031075) q[3];
sx q[3];
rz(1.0257667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0168212) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(2.4701414) q[2];
rz(2.6050383) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3143828) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(0.89624727) q[0];
rz(2.8889636) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79118585) q[0];
sx q[0];
rz(-1.5205654) q[0];
sx q[0];
rz(-2.6149261) q[0];
rz(-1.4398378) q[2];
sx q[2];
rz(-2.2316861) q[2];
sx q[2];
rz(-2.5051528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0493361) q[1];
sx q[1];
rz(-0.4886371) q[1];
sx q[1];
rz(1.2467924) q[1];
x q[2];
rz(1.518431) q[3];
sx q[3];
rz(-0.81365055) q[3];
sx q[3];
rz(-0.46458515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62319055) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(-1.5790348) q[2];
rz(1.7744428) q[3];
sx q[3];
rz(-2.0953777) q[3];
sx q[3];
rz(2.3794543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(0.38994625) q[0];
rz(-2.3155616) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(-2.0894076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5717333) q[0];
sx q[0];
rz(-2.2433407) q[0];
sx q[0];
rz(-0.39542266) q[0];
x q[1];
rz(2.7414514) q[2];
sx q[2];
rz(-2.7145128) q[2];
sx q[2];
rz(-2.6838322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4521976) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(-0.92269519) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1972606) q[3];
sx q[3];
rz(-1.1748702) q[3];
sx q[3];
rz(-2.348071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0612001) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(0.52919394) q[2];
rz(2.3906294) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.1645114) q[0];
rz(-1.8544082) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(-2.6374292) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39162779) q[0];
sx q[0];
rz(-0.34191439) q[0];
sx q[0];
rz(0.86655273) q[0];
rz(1.590056) q[2];
sx q[2];
rz(-1.9474721) q[2];
sx q[2];
rz(-2.6835359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.092445651) q[1];
sx q[1];
rz(-1.3805674) q[1];
sx q[1];
rz(-1.3008429) q[1];
rz(-2.2553798) q[3];
sx q[3];
rz(-1.4014114) q[3];
sx q[3];
rz(-2.2365295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21504271) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(-0.034493383) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(2.8300986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9813949) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(-2.1389217) q[1];
sx q[1];
rz(-1.4613338) q[1];
sx q[1];
rz(2.7636539) q[1];
rz(-2.6832081) q[2];
sx q[2];
rz(-1.3946563) q[2];
sx q[2];
rz(-1.3763225) q[2];
rz(2.814449) q[3];
sx q[3];
rz(-1.6302623) q[3];
sx q[3];
rz(0.74360256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
