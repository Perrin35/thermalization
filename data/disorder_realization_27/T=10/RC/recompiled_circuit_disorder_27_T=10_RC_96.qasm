OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.1383837) q[0];
sx q[0];
rz(-2.9870343) q[0];
sx q[0];
rz(2.4490693) q[0];
rz(-1.2094296) q[1];
sx q[1];
rz(-1.8930607) q[1];
sx q[1];
rz(-1.7564397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801075) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(2.0018342) q[0];
x q[1];
rz(0.42983774) q[2];
sx q[2];
rz(-2.5463383) q[2];
sx q[2];
rz(1.0560448) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7117118) q[1];
sx q[1];
rz(-1.8036588) q[1];
sx q[1];
rz(-1.7377322) q[1];
rz(-pi) q[2];
rz(-0.65269835) q[3];
sx q[3];
rz(-1.9737118) q[3];
sx q[3];
rz(-2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7339864) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.9143547) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7176712) q[0];
sx q[0];
rz(-2.9917891) q[0];
sx q[0];
rz(-1.040209) q[0];
x q[1];
rz(-0.51867698) q[2];
sx q[2];
rz(-1.9762602) q[2];
sx q[2];
rz(-0.60207089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7397346) q[1];
sx q[1];
rz(-1.8622073) q[1];
sx q[1];
rz(-1.3009562) q[1];
rz(-pi) q[2];
x q[2];
rz(2.482588) q[3];
sx q[3];
rz(-0.25203029) q[3];
sx q[3];
rz(-2.8641831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.041302117) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(-2.7764017) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(0.89865249) q[0];
rz(-0.99575106) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(-2.8083037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3086739) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(-2.4426016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0731508) q[2];
sx q[2];
rz(-0.2012673) q[2];
sx q[2];
rz(-1.4301436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7283199) q[1];
sx q[1];
rz(-1.4188758) q[1];
sx q[1];
rz(-2.1876213) q[1];
rz(1.0201449) q[3];
sx q[3];
rz(-1.2216976) q[3];
sx q[3];
rz(-2.246644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68625346) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(0.68471318) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-2.6338449) q[1];
sx q[1];
rz(1.205014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3048153) q[0];
sx q[0];
rz(-0.69201058) q[0];
sx q[0];
rz(-1.6230323) q[0];
rz(-2.2150546) q[2];
sx q[2];
rz(-1.7346003) q[2];
sx q[2];
rz(2.0870199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.64441427) q[1];
sx q[1];
rz(-0.89343151) q[1];
sx q[1];
rz(-2.7888984) q[1];
rz(-1.6550001) q[3];
sx q[3];
rz(-0.70662543) q[3];
sx q[3];
rz(-0.018761793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8923607) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-2.0223117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7147987) q[0];
sx q[0];
rz(-1.5835276) q[0];
sx q[0];
rz(-1.2554332) q[0];
rz(-pi) q[1];
rz(1.985717) q[2];
sx q[2];
rz(-1.6332111) q[2];
sx q[2];
rz(-1.2635363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.811378) q[1];
sx q[1];
rz(-2.0159617) q[1];
sx q[1];
rz(-3.1203169) q[1];
rz(-pi) q[2];
rz(1.5186148) q[3];
sx q[3];
rz(-1.6515886) q[3];
sx q[3];
rz(0.83524708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30620265) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(2.561835) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.6019843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.471305) q[0];
sx q[0];
rz(-1.7738713) q[0];
sx q[0];
rz(-2.2905486) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6164262) q[2];
sx q[2];
rz(-1.7012193) q[2];
sx q[2];
rz(1.4554731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0582038) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(2.8298488) q[1];
rz(-pi) q[2];
rz(-2.2555389) q[3];
sx q[3];
rz(-1.3430809) q[3];
sx q[3];
rz(2.7453604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6948029) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(2.457298) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(-0.51876846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25012384) q[0];
sx q[0];
rz(-0.84202535) q[0];
sx q[0];
rz(0.17983371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22612818) q[2];
sx q[2];
rz(-0.78352189) q[2];
sx q[2];
rz(2.4353611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7873951) q[1];
sx q[1];
rz(-2.0320315) q[1];
sx q[1];
rz(-0.72871491) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026168907) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(0.27163423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.96034399) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(-2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(-2.3972437) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9104011) q[0];
sx q[0];
rz(-1.7185128) q[0];
sx q[0];
rz(1.6258679) q[0];
x q[1];
rz(1.3049576) q[2];
sx q[2];
rz(-0.53005866) q[2];
sx q[2];
rz(0.30345464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8213615) q[1];
sx q[1];
rz(-1.3470955) q[1];
sx q[1];
rz(2.6279468) q[1];
x q[2];
rz(1.5472502) q[3];
sx q[3];
rz(-1.1163201) q[3];
sx q[3];
rz(2.8794895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-0.60950935) q[2];
rz(-0.65731796) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(1.6375861) q[0];
rz(-1.9001182) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(-0.77493587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79561728) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(-2.6443291) q[0];
rz(-pi) q[1];
rz(2.9853285) q[2];
sx q[2];
rz(-1.1322349) q[2];
sx q[2];
rz(1.4807448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(2.7369569) q[1];
x q[2];
rz(1.3516515) q[3];
sx q[3];
rz(-1.7100167) q[3];
sx q[3];
rz(-0.62300357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(-2.4841323) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(1.0459895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25046529) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(1.6842708) q[0];
rz(-pi) q[1];
rz(3.0984512) q[2];
sx q[2];
rz(-2.5501745) q[2];
sx q[2];
rz(0.45515781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18098772) q[1];
sx q[1];
rz(-2.1121896) q[1];
sx q[1];
rz(0.42591806) q[1];
rz(-3.0797327) q[3];
sx q[3];
rz(-0.40611551) q[3];
sx q[3];
rz(1.5554242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7252698) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(-2.0521169) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(2.7325148) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(-3.0388721) q[3];
sx q[3];
rz(-0.552388) q[3];
sx q[3];
rz(-1.296476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
