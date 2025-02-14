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
rz(1.0951618) q[0];
sx q[0];
rz(-2.996063) q[0];
sx q[0];
rz(-0.50330436) q[0];
rz(1.1777999) q[1];
sx q[1];
rz(-1.5947394) q[1];
sx q[1];
rz(1.4454747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88004011) q[0];
sx q[0];
rz(-1.6140811) q[0];
sx q[0];
rz(2.7425984) q[0];
x q[1];
rz(2.4769449) q[2];
sx q[2];
rz(-0.87021135) q[2];
sx q[2];
rz(-2.075891) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7912476) q[1];
sx q[1];
rz(-0.71202946) q[1];
sx q[1];
rz(-2.6256034) q[1];
rz(-pi) q[2];
rz(0.30685356) q[3];
sx q[3];
rz(-0.05159353) q[3];
sx q[3];
rz(1.8966228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6426927) q[2];
sx q[2];
rz(-1.2314726) q[2];
sx q[2];
rz(1.5042245) q[2];
rz(-0.73689342) q[3];
sx q[3];
rz(-1.5259909) q[3];
sx q[3];
rz(2.1604497) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40642834) q[0];
sx q[0];
rz(-2.6728215) q[0];
sx q[0];
rz(2.6306187) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.7402486) q[1];
sx q[1];
rz(-2.8151292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0827507) q[0];
sx q[0];
rz(-1.8371757) q[0];
sx q[0];
rz(-0.27348117) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3148221) q[2];
sx q[2];
rz(-0.69967383) q[2];
sx q[2];
rz(-2.889461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67124704) q[1];
sx q[1];
rz(-1.874028) q[1];
sx q[1];
rz(0.95945759) q[1];
rz(-1.1835008) q[3];
sx q[3];
rz(-1.8346458) q[3];
sx q[3];
rz(1.6843759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.1702801) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(-1.2962606) q[2];
rz(2.7235459) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(-0.70438284) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0728077) q[0];
sx q[0];
rz(-2.7588221) q[0];
sx q[0];
rz(2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-0.03104041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.870793) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(-0.55678456) q[0];
rz(-1.8607742) q[2];
sx q[2];
rz(-0.6781247) q[2];
sx q[2];
rz(1.4416308) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87796383) q[1];
sx q[1];
rz(-0.97743703) q[1];
sx q[1];
rz(0.29747648) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5865745) q[3];
sx q[3];
rz(-1.5567829) q[3];
sx q[3];
rz(-2.009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62640181) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(0.82702965) q[2];
rz(0.51860297) q[3];
sx q[3];
rz(-2.0293472) q[3];
sx q[3];
rz(-1.6270858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9943635) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.9299141) q[0];
rz(2.0934824) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.3406219) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82523924) q[0];
sx q[0];
rz(-1.2636501) q[0];
sx q[0];
rz(-0.65748837) q[0];
x q[1];
rz(1.7197837) q[2];
sx q[2];
rz(-1.3366615) q[2];
sx q[2];
rz(2.9643167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87340277) q[1];
sx q[1];
rz(-0.65618578) q[1];
sx q[1];
rz(2.6471958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8302835) q[3];
sx q[3];
rz(-1.177586) q[3];
sx q[3];
rz(-1.1953199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94053215) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57946) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(-2.539047) q[0];
rz(2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(-2.3604732) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3042219) q[0];
sx q[0];
rz(-1.8580372) q[0];
sx q[0];
rz(-2.9648613) q[0];
rz(0.65004543) q[2];
sx q[2];
rz(-0.87832762) q[2];
sx q[2];
rz(1.6672857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8153861) q[1];
sx q[1];
rz(-1.4844196) q[1];
sx q[1];
rz(2.2934329) q[1];
rz(-pi) q[2];
rz(-2.9424465) q[3];
sx q[3];
rz(-0.73561397) q[3];
sx q[3];
rz(-0.74935952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1802804) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(2.9126634) q[2];
rz(-1.3198352) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.9127386) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(-1.9629021) q[0];
rz(-1.2294058) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(-1.2961402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4018009) q[0];
sx q[0];
rz(-1.290088) q[0];
sx q[0];
rz(0.14174353) q[0];
rz(-pi) q[1];
rz(-0.73298323) q[2];
sx q[2];
rz(-2.146581) q[2];
sx q[2];
rz(2.1819262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.069433) q[1];
sx q[1];
rz(-1.0899223) q[1];
sx q[1];
rz(-2.062501) q[1];
rz(-pi) q[2];
rz(0.18978203) q[3];
sx q[3];
rz(-0.49473195) q[3];
sx q[3];
rz(-3.1056946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.487315) q[2];
sx q[2];
rz(-2.2389905) q[2];
sx q[2];
rz(2.8793092) q[2];
rz(0.30019635) q[3];
sx q[3];
rz(-1.9191977) q[3];
sx q[3];
rz(0.20849553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009785) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(-1.8390919) q[0];
rz(0.66849661) q[1];
sx q[1];
rz(-1.786247) q[1];
sx q[1];
rz(2.1065333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4651606) q[0];
sx q[0];
rz(-0.4598099) q[0];
sx q[0];
rz(-2.0681429) q[0];
x q[1];
rz(0.18685734) q[2];
sx q[2];
rz(-1.9344988) q[2];
sx q[2];
rz(1.1101369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8866402) q[1];
sx q[1];
rz(-1.6463841) q[1];
sx q[1];
rz(1.2264538) q[1];
rz(2.5104654) q[3];
sx q[3];
rz(-2.2911117) q[3];
sx q[3];
rz(-0.18319727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0237026) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(-1.6560076) q[2];
rz(-0.28247908) q[3];
sx q[3];
rz(-1.7829203) q[3];
sx q[3];
rz(0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(-1.2046332) q[1];
sx q[1];
rz(-2.5614673) q[1];
sx q[1];
rz(2.5801632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1875136) q[0];
sx q[0];
rz(-0.8773548) q[0];
sx q[0];
rz(-2.2827143) q[0];
rz(-pi) q[1];
rz(-1.2677293) q[2];
sx q[2];
rz(-0.91721877) q[2];
sx q[2];
rz(-0.67959626) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9778487) q[1];
sx q[1];
rz(-1.1625966) q[1];
sx q[1];
rz(0.331649) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2116978) q[3];
sx q[3];
rz(-0.42503438) q[3];
sx q[3];
rz(0.89034058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5795035) q[2];
sx q[2];
rz(-1.1056489) q[2];
sx q[2];
rz(-1.7378463) q[2];
rz(-1.3459407) q[3];
sx q[3];
rz(-0.74508777) q[3];
sx q[3];
rz(2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3812934) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(-0.90989939) q[0];
rz(-1.1970041) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(2.2534456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040374856) q[0];
sx q[0];
rz(-1.4618317) q[0];
sx q[0];
rz(3.0473112) q[0];
rz(-pi) q[1];
rz(1.6309225) q[2];
sx q[2];
rz(-0.91183582) q[2];
sx q[2];
rz(-1.2252285) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3067813) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(-2.4624834) q[1];
rz(-pi) q[2];
x q[2];
rz(0.03753438) q[3];
sx q[3];
rz(-1.8662631) q[3];
sx q[3];
rz(0.42643828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1278648) q[2];
sx q[2];
rz(-1.4415386) q[2];
sx q[2];
rz(2.0590651) q[2];
rz(2.9108289) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(-0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3103264) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(0.58151522) q[0];
rz(0.61559081) q[1];
sx q[1];
rz(-0.94771996) q[1];
sx q[1];
rz(1.8448578) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0727507) q[0];
sx q[0];
rz(-1.5889935) q[0];
sx q[0];
rz(0.010404603) q[0];
rz(-1.5611951) q[2];
sx q[2];
rz(-1.2266955) q[2];
sx q[2];
rz(-2.8167958) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3677926) q[1];
sx q[1];
rz(-1.8270711) q[1];
sx q[1];
rz(-1.132375) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1704135) q[3];
sx q[3];
rz(-2.2030911) q[3];
sx q[3];
rz(0.99611547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41946188) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(0.44931832) q[2];
rz(-2.8214473) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(1.2464574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.37353361) q[0];
sx q[0];
rz(-0.67714416) q[0];
sx q[0];
rz(-1.9376391) q[0];
rz(1.7569348) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(0.25945406) q[2];
sx q[2];
rz(-2.1726407) q[2];
sx q[2];
rz(3.0214027) q[2];
rz(2.2438335) q[3];
sx q[3];
rz(-0.57716341) q[3];
sx q[3];
rz(-2.9152277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
