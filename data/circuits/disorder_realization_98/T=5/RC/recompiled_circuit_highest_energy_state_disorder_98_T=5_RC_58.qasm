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
rz(2.4636318) q[0];
sx q[0];
rz(-2.8673661) q[0];
sx q[0];
rz(1.2908614) q[0];
rz(-3.3554606) q[1];
sx q[1];
rz(5.9670347) q[1];
sx q[1];
rz(11.066758) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052657) q[0];
sx q[0];
rz(-1.8417497) q[0];
sx q[0];
rz(1.8150369) q[0];
rz(-pi) q[1];
rz(0.77969061) q[2];
sx q[2];
rz(-1.3609972) q[2];
sx q[2];
rz(-0.22127929) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5668402) q[1];
sx q[1];
rz(-2.2623203) q[1];
sx q[1];
rz(-1.994446) q[1];
rz(2.1898299) q[3];
sx q[3];
rz(-0.72975273) q[3];
sx q[3];
rz(0.24977926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1715673) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(0.901326) q[2];
rz(-0.46402913) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-0.06981167) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0483911) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.5346079) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8044455) q[0];
sx q[0];
rz(-2.5408816) q[0];
sx q[0];
rz(-1.5128193) q[0];
rz(1.2708475) q[2];
sx q[2];
rz(-2.9510164) q[2];
sx q[2];
rz(1.9924763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6834641) q[1];
sx q[1];
rz(-2.5741815) q[1];
sx q[1];
rz(0.50000425) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8022763) q[3];
sx q[3];
rz(-2.5893988) q[3];
sx q[3];
rz(1.4331762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.719912) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-2.9812532) q[2];
rz(0.73873377) q[3];
sx q[3];
rz(-2.2952047) q[3];
sx q[3];
rz(1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(2.6318188) q[0];
rz(1.7104507) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4228446) q[0];
sx q[0];
rz(-1.5507409) q[0];
sx q[0];
rz(2.123318) q[0];
x q[1];
rz(-2.4116729) q[2];
sx q[2];
rz(-1.0906719) q[2];
sx q[2];
rz(-1.7034704) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.989403) q[1];
sx q[1];
rz(-2.2446989) q[1];
sx q[1];
rz(2.6722355) q[1];
x q[2];
rz(-2.7535149) q[3];
sx q[3];
rz(-0.92736926) q[3];
sx q[3];
rz(2.3449947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9868682) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(-1.9179087) q[2];
rz(1.0736505) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(-2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735485) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(0.37937382) q[0];
rz(-0.1768449) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-2.3462229) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75126266) q[0];
sx q[0];
rz(-0.38026938) q[0];
sx q[0];
rz(1.8230536) q[0];
rz(0.81669871) q[2];
sx q[2];
rz(-1.7672605) q[2];
sx q[2];
rz(-2.8981371) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2242608) q[1];
sx q[1];
rz(-0.59228173) q[1];
sx q[1];
rz(-2.0439953) q[1];
rz(-1.5136817) q[3];
sx q[3];
rz(-1.6962255) q[3];
sx q[3];
rz(0.40298395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4312326) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(1.360652) q[2];
rz(-1.0294754) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.3220538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5929476) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.7318447) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.627305) q[0];
sx q[0];
rz(-2.3804733) q[0];
sx q[0];
rz(-0.29157761) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35946741) q[2];
sx q[2];
rz(-1.4233982) q[2];
sx q[2];
rz(-2.0161674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4673938) q[1];
sx q[1];
rz(-1.129303) q[1];
sx q[1];
rz(-0.91798325) q[1];
rz(-pi) q[2];
rz(-1.665776) q[3];
sx q[3];
rz(-2.5174369) q[3];
sx q[3];
rz(0.74172663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8283525) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(-2.7102615) q[2];
rz(-0.055015419) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(-0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3738275) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(-1.025169) q[0];
rz(0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-0.90389171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4898281) q[0];
sx q[0];
rz(-1.8790885) q[0];
sx q[0];
rz(-0.90109189) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3329817) q[2];
sx q[2];
rz(-2.3249194) q[2];
sx q[2];
rz(-2.6220235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76662582) q[1];
sx q[1];
rz(-1.300525) q[1];
sx q[1];
rz(-0.75281669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7359832) q[3];
sx q[3];
rz(-1.4017229) q[3];
sx q[3];
rz(1.0195093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30770939) q[2];
sx q[2];
rz(-1.6857952) q[2];
sx q[2];
rz(-2.9252388) q[2];
rz(2.2458535) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(1.5007796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(-2.2898477) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18139938) q[0];
sx q[0];
rz(-1.6087247) q[0];
sx q[0];
rz(2.1071803) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5642407) q[2];
sx q[2];
rz(-1.5923579) q[2];
sx q[2];
rz(-0.93761629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5761583) q[1];
sx q[1];
rz(-1.2093739) q[1];
sx q[1];
rz(2.4225077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64166358) q[3];
sx q[3];
rz(-1.9681276) q[3];
sx q[3];
rz(1.404055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8375887) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(-2.6336929) q[2];
rz(-1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(-0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4429338) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(-0.0096631924) q[0];
rz(-1.5254321) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-0.38472167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5670861) q[0];
sx q[0];
rz(-0.28038803) q[0];
sx q[0];
rz(1.0100421) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9392249) q[2];
sx q[2];
rz(-1.191539) q[2];
sx q[2];
rz(-2.8059562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5009969) q[1];
sx q[1];
rz(-1.2944376) q[1];
sx q[1];
rz(-2.9261057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1423903) q[3];
sx q[3];
rz(-2.2561142) q[3];
sx q[3];
rz(0.1911605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.559451) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(3.1257296) q[2];
rz(-1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(-0.62172833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7568307) q[0];
sx q[0];
rz(-0.66487304) q[0];
sx q[0];
rz(2.050198) q[0];
rz(0.81820828) q[1];
sx q[1];
rz(-1.3312157) q[1];
sx q[1];
rz(-0.6689201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2125137) q[0];
sx q[0];
rz(-0.90796472) q[0];
sx q[0];
rz(0.79848358) q[0];
rz(1.0878272) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(-0.20394606) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8525772) q[1];
sx q[1];
rz(-2.1271879) q[1];
sx q[1];
rz(0.08361967) q[1];
x q[2];
rz(0.98812859) q[3];
sx q[3];
rz(-1.291953) q[3];
sx q[3];
rz(-0.13259618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6341256) q[2];
sx q[2];
rz(-2.6125245) q[2];
sx q[2];
rz(-0.96366209) q[2];
rz(-1.730718) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(1.3655519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249579) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(-2.813521) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011400819) q[0];
sx q[0];
rz(-1.6284124) q[0];
sx q[0];
rz(-1.685623) q[0];
rz(-pi) q[1];
rz(-2.5727082) q[2];
sx q[2];
rz(-2.9279104) q[2];
sx q[2];
rz(1.793022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0340424) q[1];
sx q[1];
rz(-0.75354105) q[1];
sx q[1];
rz(0.57489245) q[1];
rz(-pi) q[2];
rz(-1.7441196) q[3];
sx q[3];
rz(-1.4651235) q[3];
sx q[3];
rz(-0.25143501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6803117) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(-2.7247562) q[2];
rz(-2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1998491) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(-0.70882123) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(-1.3016635) q[2];
sx q[2];
rz(-0.65432815) q[2];
sx q[2];
rz(2.339044) q[2];
rz(2.57614) q[3];
sx q[3];
rz(-0.91309091) q[3];
sx q[3];
rz(2.0500195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
