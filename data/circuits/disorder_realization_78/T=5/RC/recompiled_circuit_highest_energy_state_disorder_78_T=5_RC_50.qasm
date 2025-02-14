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
rz(-1.6260835) q[0];
sx q[0];
rz(-0.27684394) q[0];
sx q[0];
rz(-3.0247363) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(0.49122223) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304162) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(-1.6947794) q[0];
rz(-0.48727076) q[2];
sx q[2];
rz(-1.3830162) q[2];
sx q[2];
rz(-1.0010546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95211103) q[1];
sx q[1];
rz(-0.64095381) q[1];
sx q[1];
rz(2.9753995) q[1];
rz(-pi) q[2];
rz(1.9141044) q[3];
sx q[3];
rz(-1.5945537) q[3];
sx q[3];
rz(-2.6809965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3789856) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-2.911705) q[2];
rz(-0.057279438) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(3.1087061) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(0.65509534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65194535) q[0];
sx q[0];
rz(-0.86344749) q[0];
sx q[0];
rz(0.26443692) q[0];
rz(1.7704324) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(2.9889929) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078121878) q[1];
sx q[1];
rz(-2.1912716) q[1];
sx q[1];
rz(-0.93777754) q[1];
x q[2];
rz(-1.412185) q[3];
sx q[3];
rz(-0.74071032) q[3];
sx q[3];
rz(2.9155758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(-0.61400145) q[0];
rz(-1.7738495) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(0.49555379) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7921831) q[0];
sx q[0];
rz(-1.111545) q[0];
sx q[0];
rz(-0.71190603) q[0];
rz(-pi) q[1];
rz(-0.88286384) q[2];
sx q[2];
rz(-1.9830215) q[2];
sx q[2];
rz(-2.8799873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5835201) q[1];
sx q[1];
rz(-2.8360543) q[1];
sx q[1];
rz(0.90684857) q[1];
x q[2];
rz(1.9248149) q[3];
sx q[3];
rz(-1.1292158) q[3];
sx q[3];
rz(1.3819249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(2.5490254) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(-2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-0.52111202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65893764) q[0];
sx q[0];
rz(-1.5636496) q[0];
sx q[0];
rz(0.004645017) q[0];
rz(-1.5278339) q[2];
sx q[2];
rz(-0.60138541) q[2];
sx q[2];
rz(-3.0850981) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8940034) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(-2.4127059) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24346015) q[3];
sx q[3];
rz(-0.32103466) q[3];
sx q[3];
rz(0.76317549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(2.4288154) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.6596863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7707959) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(-3.1053012) q[0];
rz(-0.52641422) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(0.2821736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68809915) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(-2.009925) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93499489) q[2];
sx q[2];
rz(-2.3381669) q[2];
sx q[2];
rz(-0.73551169) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5314922) q[1];
sx q[1];
rz(-2.2740915) q[1];
sx q[1];
rz(1.478757) q[1];
rz(-pi) q[2];
rz(-0.80970851) q[3];
sx q[3];
rz(-0.96862176) q[3];
sx q[3];
rz(-0.47145999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7669749) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(0.96595079) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.1740603) q[0];
rz(2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-3.1016268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46972457) q[0];
sx q[0];
rz(-1.7789808) q[0];
sx q[0];
rz(-1.2428268) q[0];
x q[1];
rz(-2.2117679) q[2];
sx q[2];
rz(-1.3017968) q[2];
sx q[2];
rz(-0.19942927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3921622) q[1];
sx q[1];
rz(-2.431173) q[1];
sx q[1];
rz(0.53954069) q[1];
rz(1.3371991) q[3];
sx q[3];
rz(-0.75832483) q[3];
sx q[3];
rz(0.17994623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2018955) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(0.8197909) q[2];
rz(2.36006) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(-1.7614822) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(-1.3167205) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(-2.7046611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1986982) q[0];
sx q[0];
rz(-1.534053) q[0];
sx q[0];
rz(1.7612639) q[0];
rz(-pi) q[1];
rz(2.390449) q[2];
sx q[2];
rz(-0.1607543) q[2];
sx q[2];
rz(-3.1243665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97368846) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(-2.3547947) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42610069) q[3];
sx q[3];
rz(-2.0668732) q[3];
sx q[3];
rz(-2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(-1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(-1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(0.18184161) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-2.9827319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358064) q[0];
sx q[0];
rz(-1.5022359) q[0];
sx q[0];
rz(0.44802702) q[0];
rz(-pi) q[1];
rz(1.8552637) q[2];
sx q[2];
rz(-1.163223) q[2];
sx q[2];
rz(1.4409232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5103858) q[1];
sx q[1];
rz(-2.800731) q[1];
sx q[1];
rz(-2.4426961) q[1];
rz(-pi) q[2];
rz(-1.7233061) q[3];
sx q[3];
rz(-2.1195863) q[3];
sx q[3];
rz(1.9085409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2813256) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(-0.9935317) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(-0.55307585) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4489514) q[0];
sx q[0];
rz(-1.7990489) q[0];
sx q[0];
rz(-0.3897764) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0680066) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(-0.70895665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.304405) q[1];
sx q[1];
rz(-2.4585729) q[1];
sx q[1];
rz(-2.7262615) q[1];
x q[2];
rz(2.210238) q[3];
sx q[3];
rz(-1.7903312) q[3];
sx q[3];
rz(-0.94506782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9370344) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(0.23703144) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700664) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(0.1189098) q[0];
rz(-0.43926829) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(0.063442245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7164511) q[0];
sx q[0];
rz(-2.2277895) q[0];
sx q[0];
rz(0.69752661) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62427743) q[2];
sx q[2];
rz(-1.7221525) q[2];
sx q[2];
rz(2.4454012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6193004) q[1];
sx q[1];
rz(-1.5792773) q[1];
sx q[1];
rz(-1.2441166) q[1];
x q[2];
rz(2.2584553) q[3];
sx q[3];
rz(-2.4534907) q[3];
sx q[3];
rz(2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1305609) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.4981184) q[2];
rz(-2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7508004) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(0.87286585) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(0.93449705) q[2];
sx q[2];
rz(-2.1515982) q[2];
sx q[2];
rz(1.2610255) q[2];
rz(2.5641702) q[3];
sx q[3];
rz(-1.631177) q[3];
sx q[3];
rz(-2.9384818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
