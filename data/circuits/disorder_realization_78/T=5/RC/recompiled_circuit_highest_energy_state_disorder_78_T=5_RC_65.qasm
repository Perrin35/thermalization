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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304162) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(-1.6947794) q[0];
rz(-pi) q[1];
rz(0.38552706) q[2];
sx q[2];
rz(-2.6221199) q[2];
sx q[2];
rz(0.2310209) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74581214) q[1];
sx q[1];
rz(-0.94008259) q[1];
sx q[1];
rz(-1.4480026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5003203) q[3];
sx q[3];
rz(-2.797496) q[3];
sx q[3];
rz(1.1765574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3789856) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(2.911705) q[2];
rz(0.057279438) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-0.21695319) q[0];
rz(-3.1087061) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(0.65509534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.396916) q[0];
sx q[0];
rz(-1.370805) q[0];
sx q[0];
rz(2.295664) q[0];
x q[1];
rz(2.6397471) q[2];
sx q[2];
rz(-1.7463533) q[2];
sx q[2];
rz(-1.6278536) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82258525) q[1];
sx q[1];
rz(-2.2862541) q[1];
sx q[1];
rz(0.69147488) q[1];
x q[2];
rz(-0.83637303) q[3];
sx q[3];
rz(-1.4640088) q[3];
sx q[3];
rz(1.2272782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(-1.2201747) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(-2.5275912) q[0];
rz(-1.7738495) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(-2.6460389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3494095) q[0];
sx q[0];
rz(-2.0300477) q[0];
sx q[0];
rz(2.4296866) q[0];
rz(2.6265385) q[2];
sx q[2];
rz(-2.1916789) q[2];
sx q[2];
rz(-0.99109288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55807251) q[1];
sx q[1];
rz(-2.8360543) q[1];
sx q[1];
rz(-2.2347441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2167778) q[3];
sx q[3];
rz(-2.0123768) q[3];
sx q[3];
rz(-1.3819249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3765748) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-2.5490254) q[2];
rz(0.23558922) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(2.6872915) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6465103) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(0.5799154) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(-2.6204806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65893764) q[0];
sx q[0];
rz(-1.5636496) q[0];
sx q[0];
rz(-0.004645017) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6137587) q[2];
sx q[2];
rz(-0.60138541) q[2];
sx q[2];
rz(-0.056494519) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8940034) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(2.4127059) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(-2.8233675) q[2];
rz(2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.6596863) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(-3.1053012) q[0];
rz(-2.6151784) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-2.859419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4534935) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(-2.009925) q[0];
rz(-pi) q[1];
rz(-2.2065978) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(2.406081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6101004) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(1.478757) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75936486) q[3];
sx q[3];
rz(-2.1753484) q[3];
sx q[3];
rz(-2.5370363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7669749) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(-0.10147632) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(-0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58979708) q[0];
sx q[0];
rz(-0.17687251) q[0];
sx q[0];
rz(-1.1740603) q[0];
rz(-0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.171282) q[0];
sx q[0];
rz(-1.891428) q[0];
sx q[0];
rz(2.9220504) q[0];
x q[1];
rz(-2.2117679) q[2];
sx q[2];
rz(-1.3017968) q[2];
sx q[2];
rz(2.9421634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.894095) q[1];
sx q[1];
rz(-1.9124417) q[1];
sx q[1];
rz(2.5057807) q[1];
rz(-pi) q[2];
rz(-1.8043936) q[3];
sx q[3];
rz(-0.75832483) q[3];
sx q[3];
rz(-2.9616464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(-2.3218018) q[2];
rz(-0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509197) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(1.3167205) q[0];
rz(-1.3680178) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(2.7046611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1986982) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(-1.7612639) q[0];
rz(-0.11796911) q[2];
sx q[2];
rz(-1.461339) q[2];
sx q[2];
rz(-0.80889672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1679042) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(2.3547947) q[1];
x q[2];
rz(2.223001) q[3];
sx q[3];
rz(-0.64213412) q[3];
sx q[3];
rz(1.0561424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(0.13460049) q[2];
rz(1.918119) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80158919) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(-0.34135154) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-2.9827319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6065935) q[0];
sx q[0];
rz(-0.45289055) q[0];
sx q[0];
rz(0.15720982) q[0];
rz(-pi) q[1];
rz(-1.8552637) q[2];
sx q[2];
rz(-1.9783696) q[2];
sx q[2];
rz(1.4409232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73028683) q[1];
sx q[1];
rz(-1.7875693) q[1];
sx q[1];
rz(-2.8764399) q[1];
rz(-1.4182866) q[3];
sx q[3];
rz(-2.1195863) q[3];
sx q[3];
rz(1.2330518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(2.2567828) q[2];
rz(1.7224711) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(0.9935317) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(2.5885168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7599604) q[0];
sx q[0];
rz(-2.6928718) q[0];
sx q[0];
rz(-2.5928709) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6378013) q[2];
sx q[2];
rz(-0.83344668) q[2];
sx q[2];
rz(2.532117) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40402624) q[1];
sx q[1];
rz(-1.8282923) q[1];
sx q[1];
rz(0.63996686) q[1];
rz(-pi) q[2];
rz(-0.27121765) q[3];
sx q[3];
rz(-0.94908774) q[3];
sx q[3];
rz(-0.78628899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20455827) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(2.9045612) q[2];
rz(-1.6009181) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(-1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700664) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(0.43926829) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(0.063442245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251415) q[0];
sx q[0];
rz(-0.91380318) q[0];
sx q[0];
rz(-0.69752661) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.756606) q[2];
sx q[2];
rz(-0.95473052) q[2];
sx q[2];
rz(-0.98279976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6193004) q[1];
sx q[1];
rz(-1.5792773) q[1];
sx q[1];
rz(-1.897476) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48096913) q[3];
sx q[3];
rz(-2.0837373) q[3];
sx q[3];
rz(-2.6875567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0110317) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(-1.4981184) q[2];
rz(-2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(2.2070956) q[2];
sx q[2];
rz(-0.9899944) q[2];
sx q[2];
rz(-1.8805671) q[2];
rz(3.031293) q[3];
sx q[3];
rz(-2.5613789) q[3];
sx q[3];
rz(-1.2753244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
