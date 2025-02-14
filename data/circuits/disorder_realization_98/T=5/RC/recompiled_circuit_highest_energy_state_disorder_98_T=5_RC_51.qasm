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
rz(-1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0555178) q[0];
sx q[0];
rz(-0.3627622) q[0];
sx q[0];
rz(2.4253009) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2798645) q[2];
sx q[2];
rz(-0.81255823) q[2];
sx q[2];
rz(1.1464553) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95924458) q[1];
sx q[1];
rz(-0.79234353) q[1];
sx q[1];
rz(2.6807129) q[1];
x q[2];
rz(-0.95176272) q[3];
sx q[3];
rz(-2.4118399) q[3];
sx q[3];
rz(-0.24977926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1715673) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(-2.2402666) q[2];
rz(-2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0932015) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(1.2777591) q[0];
rz(-0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.6069848) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33714715) q[0];
sx q[0];
rz(-0.60071105) q[0];
sx q[0];
rz(1.5128193) q[0];
rz(1.8707451) q[2];
sx q[2];
rz(-2.9510164) q[2];
sx q[2];
rz(1.1491164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1086191) q[1];
sx q[1];
rz(-1.0796282) q[1];
sx q[1];
rz(-1.2742548) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8022763) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(1.4331762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4216807) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(2.9812532) q[2];
rz(0.73873377) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(-1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(-2.6318188) q[0];
rz(1.431142) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0261542) q[0];
sx q[0];
rz(-0.55284772) q[0];
sx q[0];
rz(1.5325969) q[0];
x q[1];
rz(-0.72991972) q[2];
sx q[2];
rz(-2.0509208) q[2];
sx q[2];
rz(1.4381222) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.4697527) q[1];
sx q[1];
rz(-2.3418509) q[1];
sx q[1];
rz(-1.0554764) q[1];
x q[2];
rz(0.38807773) q[3];
sx q[3];
rz(-2.2142234) q[3];
sx q[3];
rz(-2.3449947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15472445) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(2.0679421) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(0.79536974) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.11957) q[0];
sx q[0];
rz(-1.2031462) q[0];
sx q[0];
rz(-0.099438503) q[0];
x q[1];
rz(-1.8537117) q[2];
sx q[2];
rz(-2.3672315) q[2];
sx q[2];
rz(2.0191569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6713555) q[1];
sx q[1];
rz(-2.0908326) q[1];
sx q[1];
rz(-2.8440471) q[1];
rz(-pi) q[2];
rz(3.0159608) q[3];
sx q[3];
rz(-1.5141308) q[3];
sx q[3];
rz(1.1606596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7103601) q[2];
sx q[2];
rz(-1.3453588) q[2];
sx q[2];
rz(1.360652) q[2];
rz(-2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(-1.3220538) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(1.5486451) q[0];
rz(-0.40052888) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(-1.7318447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.627305) q[0];
sx q[0];
rz(-0.76111932) q[0];
sx q[0];
rz(0.29157761) q[0];
rz(-pi) q[1];
rz(1.7280987) q[2];
sx q[2];
rz(-1.2154007) q[2];
sx q[2];
rz(0.50050625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21239195) q[1];
sx q[1];
rz(-2.1521795) q[1];
sx q[1];
rz(-2.6049032) q[1];
rz(-pi) q[2];
rz(1.665776) q[3];
sx q[3];
rz(-0.62415571) q[3];
sx q[3];
rz(0.74172663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3132402) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(-0.43133119) q[2];
rz(3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(2.5855248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7677652) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(2.1164236) q[0];
rz(0.45285666) q[1];
sx q[1];
rz(-2.5621474) q[1];
sx q[1];
rz(0.90389171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1548123) q[0];
sx q[0];
rz(-0.93789369) q[0];
sx q[0];
rz(2.7557719) q[0];
rz(-pi) q[1];
rz(-0.80861096) q[2];
sx q[2];
rz(-2.3249194) q[2];
sx q[2];
rz(2.6220235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3749668) q[1];
sx q[1];
rz(-1.8410676) q[1];
sx q[1];
rz(-0.75281669) q[1];
rz(-pi) q[2];
rz(-0.76670209) q[3];
sx q[3];
rz(-2.9057716) q[3];
sx q[3];
rz(-0.23877777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8338833) q[2];
sx q[2];
rz(-1.6857952) q[2];
sx q[2];
rz(-0.21635381) q[2];
rz(-2.2458535) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(-0.69881451) q[0];
rz(-1.9578594) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(0.85174495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18139938) q[0];
sx q[0];
rz(-1.532868) q[0];
sx q[0];
rz(1.0344124) q[0];
rz(-1.5965271) q[2];
sx q[2];
rz(-2.1479969) q[2];
sx q[2];
rz(0.647223) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5761583) q[1];
sx q[1];
rz(-1.9322188) q[1];
sx q[1];
rz(2.4225077) q[1];
rz(2.5301039) q[3];
sx q[3];
rz(-2.4019255) q[3];
sx q[3];
rz(0.31114331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.304004) q[2];
sx q[2];
rz(-1.8303266) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4429338) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(-0.0096631924) q[0];
rz(-1.5254321) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(0.38472167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54668173) q[0];
sx q[0];
rz(-1.4230886) q[0];
sx q[0];
rz(1.8099996) q[0];
x q[1];
rz(-2.9392249) q[2];
sx q[2];
rz(-1.191539) q[2];
sx q[2];
rz(-0.33563644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8249121) q[1];
sx q[1];
rz(-2.7928565) q[1];
sx q[1];
rz(2.2168) q[1];
rz(-pi) q[2];
rz(-2.3703897) q[3];
sx q[3];
rz(-2.0029541) q[3];
sx q[3];
rz(-1.766253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5821417) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(-3.1257296) q[2];
rz(1.2265685) q[3];
sx q[3];
rz(-2.3358986) q[3];
sx q[3];
rz(-2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3847619) q[0];
sx q[0];
rz(-0.66487304) q[0];
sx q[0];
rz(1.0913947) q[0];
rz(-2.3233844) q[1];
sx q[1];
rz(-1.3312157) q[1];
sx q[1];
rz(-0.6689201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2050987) q[0];
sx q[0];
rz(-0.97081796) q[0];
sx q[0];
rz(-0.72941248) q[0];
x q[1];
rz(1.0878272) q[2];
sx q[2];
rz(-0.7730128) q[2];
sx q[2];
rz(0.20394606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.131625) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(1.4373006) q[1];
rz(-pi) q[2];
rz(-0.33031611) q[3];
sx q[3];
rz(-2.1282176) q[3];
sx q[3];
rz(1.6176318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5074671) q[2];
sx q[2];
rz(-2.6125245) q[2];
sx q[2];
rz(-0.96366209) q[2];
rz(1.4108747) q[3];
sx q[3];
rz(-1.7129292) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016634781) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(-1.4547263) q[0];
rz(-1.1085054) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(2.813521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011400819) q[0];
sx q[0];
rz(-1.5131803) q[0];
sx q[0];
rz(-1.4559697) q[0];
x q[1];
rz(-0.180822) q[2];
sx q[2];
rz(-1.6852813) q[2];
sx q[2];
rz(-0.78071981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83392622) q[1];
sx q[1];
rz(-0.95912479) q[1];
sx q[1];
rz(2.0425379) q[1];
rz(-pi) q[2];
rz(2.1222018) q[3];
sx q[3];
rz(-0.20272045) q[3];
sx q[3];
rz(-0.77714506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46128094) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(0.4168365) q[2];
rz(2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(-0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9417435) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(2.4327714) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(1.3016635) q[2];
sx q[2];
rz(-2.4872645) q[2];
sx q[2];
rz(-0.80254868) q[2];
rz(-0.56545267) q[3];
sx q[3];
rz(-0.91309091) q[3];
sx q[3];
rz(2.0500195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
