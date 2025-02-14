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
rz(-0.79688537) q[0];
sx q[0];
rz(4.8719811) q[0];
sx q[0];
rz(13.48579) q[0];
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
rz(-2.8714234) q[0];
sx q[0];
rz(-1.7180301) q[0];
sx q[0];
rz(-1.2933613) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0242516) q[2];
sx q[2];
rz(-2.9988117) q[2];
sx q[2];
rz(0.26856865) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9956869) q[1];
sx q[1];
rz(-2.6906081) q[1];
sx q[1];
rz(0.0031975071) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0994751) q[3];
sx q[3];
rz(-1.5402392) q[3];
sx q[3];
rz(-1.9224642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.087622341) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-2.7719882) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10577781) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(2.5685487) q[0];
rz(-1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(2.6699452) q[1];
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
x q[1];
rz(-0.69219442) q[2];
sx q[2];
rz(-1.0613777) q[2];
sx q[2];
rz(0.30864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9976366) q[1];
sx q[1];
rz(-0.9846217) q[1];
sx q[1];
rz(-0.8950986) q[1];
rz(-1.0466659) q[3];
sx q[3];
rz(-1.2132304) q[3];
sx q[3];
rz(0.78525728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4196709) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(2.0527077) q[2];
rz(-1.0645083) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(2.815912) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(0.017070008) q[0];
rz(-2.4315289) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(-0.56627083) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6951272) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(1.5873853) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5336419) q[2];
sx q[2];
rz(-0.81038108) q[2];
sx q[2];
rz(1.2255526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8108845) q[1];
sx q[1];
rz(-0.20527923) q[1];
sx q[1];
rz(0.42930557) q[1];
x q[2];
rz(2.8763883) q[3];
sx q[3];
rz(-1.5062766) q[3];
sx q[3];
rz(1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.337734) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(-0.0090948661) q[2];
rz(0.13633063) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424054) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(2.0280139) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(2.9794433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1095162) q[0];
sx q[0];
rz(-1.6639532) q[0];
sx q[0];
rz(1.3901802) q[0];
rz(-2.5926431) q[2];
sx q[2];
rz(-1.5918343) q[2];
sx q[2];
rz(2.8124513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5177734) q[1];
sx q[1];
rz(-0.86408593) q[1];
sx q[1];
rz(-2.049974) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6643296) q[3];
sx q[3];
rz(-1.9151315) q[3];
sx q[3];
rz(-2.9323102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1589511) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8714137) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(2.357645) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(2.1606826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62215319) q[0];
sx q[0];
rz(-1.781917) q[0];
sx q[0];
rz(-2.7252174) q[0];
rz(-2.1891315) q[2];
sx q[2];
rz(-0.59788579) q[2];
sx q[2];
rz(2.7088504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.09771654) q[1];
sx q[1];
rz(-1.6642769) q[1];
sx q[1];
rz(-3.082117) q[1];
rz(2.122287) q[3];
sx q[3];
rz(-1.7781584) q[3];
sx q[3];
rz(2.1503748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(0.8832461) q[2];
rz(-0.20032459) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(1.0909572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9969295) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(-2.1494179) q[0];
rz(-0.80157533) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(1.7209524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039889) q[0];
sx q[0];
rz(-0.9524571) q[0];
sx q[0];
rz(-1.9644587) q[0];
rz(2.1031441) q[2];
sx q[2];
rz(-1.4803998) q[2];
sx q[2];
rz(-0.20648512) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9538624) q[1];
sx q[1];
rz(-1.2427325) q[1];
sx q[1];
rz(2.4806054) q[1];
x q[2];
rz(1.4735953) q[3];
sx q[3];
rz(-1.6080813) q[3];
sx q[3];
rz(-2.4241222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27986032) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(0.53708616) q[2];
rz(-0.01072695) q[3];
sx q[3];
rz(-0.74672943) q[3];
sx q[3];
rz(2.2906394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(-1.1713015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9635668) q[0];
sx q[0];
rz(-1.8676571) q[0];
sx q[0];
rz(1.6087449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3572071) q[2];
sx q[2];
rz(-1.7860951) q[2];
sx q[2];
rz(-2.4395296) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4170467) q[1];
sx q[1];
rz(-1.5486716) q[1];
sx q[1];
rz(1.538365) q[1];
rz(-pi) q[2];
rz(-2.5170691) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0168212) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-2.4701414) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(2.0898537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(-2.2453454) q[0];
rz(-0.25262901) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(1.0505229) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3327961) q[0];
sx q[0];
rz(-2.0967297) q[0];
sx q[0];
rz(-1.6288847) q[0];
x q[1];
rz(-2.4765313) q[2];
sx q[2];
rz(-1.4675234) q[2];
sx q[2];
rz(2.2879083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3747762) q[1];
sx q[1];
rz(-1.4207867) q[1];
sx q[1];
rz(1.1039724) q[1];
rz(2.3837621) q[3];
sx q[3];
rz(-1.6088472) q[3];
sx q[3];
rz(-1.999397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5184021) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(-1.5625578) q[2];
rz(-1.3671499) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(-2.7516464) q[0];
rz(2.3155616) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(1.0521851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74653502) q[0];
sx q[0];
rz(-1.2647226) q[0];
sx q[0];
rz(0.85880441) q[0];
x q[1];
rz(0.40014123) q[2];
sx q[2];
rz(-2.7145128) q[2];
sx q[2];
rz(2.6838322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4521976) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(0.92269519) q[1];
rz(-0.94433208) q[3];
sx q[3];
rz(-1.9667224) q[3];
sx q[3];
rz(-0.79352165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(-2.6123987) q[2];
rz(-2.3906294) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879409) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(1.1645114) q[0];
rz(-1.8544082) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(2.6374292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8541478) q[0];
sx q[0];
rz(-1.3519671) q[0];
sx q[0];
rz(1.3059421) q[0];
x q[1];
rz(0.048647886) q[2];
sx q[2];
rz(-2.7644483) q[2];
sx q[2];
rz(-0.51037558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87889844) q[1];
sx q[1];
rz(-0.32890186) q[1];
sx q[1];
rz(0.94543381) q[1];
x q[2];
rz(-0.88621288) q[3];
sx q[3];
rz(-1.7401812) q[3];
sx q[3];
rz(-2.2365295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21504271) q[2];
sx q[2];
rz(-1.5861009) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9813949) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(2.1389217) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(1.3748693) q[2];
sx q[2];
rz(-1.1200323) q[2];
sx q[2];
rz(-3.0333698) q[2];
rz(0.18319753) q[3];
sx q[3];
rz(-0.33231322) q[3];
sx q[3];
rz(-0.65381369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
