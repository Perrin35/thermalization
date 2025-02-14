OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9893875) q[0];
sx q[0];
rz(-2.9289991) q[0];
sx q[0];
rz(2.4074182) q[0];
rz(1.2165767) q[1];
sx q[1];
rz(-2.7956378) q[1];
sx q[1];
rz(0.024756519) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5261794) q[0];
sx q[0];
rz(-3.0657112) q[0];
sx q[0];
rz(-2.5350201) q[0];
x q[1];
rz(-2.0169746) q[2];
sx q[2];
rz(-0.25814787) q[2];
sx q[2];
rz(-0.1218957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9563958) q[1];
sx q[1];
rz(-1.6231771) q[1];
sx q[1];
rz(-2.8422794) q[1];
rz(-3.0856783) q[3];
sx q[3];
rz(-1.7828807) q[3];
sx q[3];
rz(0.66842118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2572702) q[2];
sx q[2];
rz(-1.6486282) q[2];
sx q[2];
rz(-2.8327668) q[2];
rz(-0.8257927) q[3];
sx q[3];
rz(-2.1335996) q[3];
sx q[3];
rz(0.76396137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17077133) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(-1.0276851) q[0];
rz(2.3254501) q[1];
sx q[1];
rz(-2.0232537) q[1];
sx q[1];
rz(2.3291086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1493133) q[0];
sx q[0];
rz(-1.9198787) q[0];
sx q[0];
rz(0.16916738) q[0];
rz(2.7682253) q[2];
sx q[2];
rz(-1.5575711) q[2];
sx q[2];
rz(2.2227299) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8390439) q[1];
sx q[1];
rz(-1.7086281) q[1];
sx q[1];
rz(-1.0562357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3489209) q[3];
sx q[3];
rz(-1.9066208) q[3];
sx q[3];
rz(-0.15062697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57427788) q[2];
sx q[2];
rz(-1.652176) q[2];
sx q[2];
rz(2.2226649) q[2];
rz(-2.610176) q[3];
sx q[3];
rz(-2.0960977) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70319217) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(-1.137314) q[0];
rz(-1.3905585) q[1];
sx q[1];
rz(-2.5847692) q[1];
sx q[1];
rz(2.4086319) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175388) q[0];
sx q[0];
rz(-1.8630872) q[0];
sx q[0];
rz(2.7089416) q[0];
rz(-pi) q[1];
rz(-2.175019) q[2];
sx q[2];
rz(-2.0546277) q[2];
sx q[2];
rz(-2.1343729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.809124) q[1];
sx q[1];
rz(-2.4753503) q[1];
sx q[1];
rz(0.82102832) q[1];
rz(-pi) q[2];
rz(-3.1350144) q[3];
sx q[3];
rz(-1.1220782) q[3];
sx q[3];
rz(2.6980163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0040032337) q[2];
sx q[2];
rz(-1.7294451) q[2];
sx q[2];
rz(1.0388733) q[2];
rz(-2.1393356) q[3];
sx q[3];
rz(-1.8398617) q[3];
sx q[3];
rz(0.29016289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.20659474) q[0];
sx q[0];
rz(-1.0863786) q[0];
sx q[0];
rz(2.2344053) q[0];
rz(-2.9838003) q[1];
sx q[1];
rz(-2.1267499) q[1];
sx q[1];
rz(1.2812322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32203755) q[0];
sx q[0];
rz(-2.0896308) q[0];
sx q[0];
rz(2.4253393) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2864986) q[2];
sx q[2];
rz(-1.2629328) q[2];
sx q[2];
rz(1.5474873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7042027) q[1];
sx q[1];
rz(-1.1466195) q[1];
sx q[1];
rz(2.4553039) q[1];
rz(0.63381845) q[3];
sx q[3];
rz(-2.1396643) q[3];
sx q[3];
rz(-2.947629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24718757) q[2];
sx q[2];
rz(-0.93623585) q[2];
sx q[2];
rz(-1.2897162) q[2];
rz(1.3442518) q[3];
sx q[3];
rz(-1.4569837) q[3];
sx q[3];
rz(2.4414506) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46071389) q[0];
sx q[0];
rz(-0.49837708) q[0];
sx q[0];
rz(1.0182678) q[0];
rz(-2.0385888) q[1];
sx q[1];
rz(-2.4296727) q[1];
sx q[1];
rz(1.8011372) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0809652) q[0];
sx q[0];
rz(-0.61765352) q[0];
sx q[0];
rz(-0.65184848) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6685772) q[2];
sx q[2];
rz(-1.9488397) q[2];
sx q[2];
rz(-0.60286544) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.05100194) q[1];
sx q[1];
rz(-0.18185073) q[1];
sx q[1];
rz(-0.79504063) q[1];
rz(-pi) q[2];
rz(-0.057647905) q[3];
sx q[3];
rz(-1.0312005) q[3];
sx q[3];
rz(0.12069139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(1.0020024) q[2];
rz(2.4501948) q[3];
sx q[3];
rz(-1.1804429) q[3];
sx q[3];
rz(1.5207759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6028676) q[0];
sx q[0];
rz(-1.069101) q[0];
sx q[0];
rz(-0.62527239) q[0];
rz(-1.193115) q[1];
sx q[1];
rz(-2.4517877) q[1];
sx q[1];
rz(-0.56732059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35463542) q[0];
sx q[0];
rz(-2.3428095) q[0];
sx q[0];
rz(1.8171993) q[0];
rz(-pi) q[1];
rz(-2.3151933) q[2];
sx q[2];
rz(-1.5652839) q[2];
sx q[2];
rz(-2.7720087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4299791) q[1];
sx q[1];
rz(-2.1890854) q[1];
sx q[1];
rz(3.1097163) q[1];
rz(-pi) q[2];
rz(0.17622275) q[3];
sx q[3];
rz(-1.3941951) q[3];
sx q[3];
rz(3.0295102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.304473) q[2];
sx q[2];
rz(-1.2834872) q[2];
sx q[2];
rz(1.5817969) q[2];
rz(0.61257735) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(-2.9858203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(1.1268536) q[0];
rz(1.8719748) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(-0.52106214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57018897) q[0];
sx q[0];
rz(-0.91556433) q[0];
sx q[0];
rz(-0.44555026) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2418141) q[2];
sx q[2];
rz(-1.6018724) q[2];
sx q[2];
rz(0.46361332) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5875257) q[1];
sx q[1];
rz(-0.8272285) q[1];
sx q[1];
rz(0.83283333) q[1];
rz(-0.43989681) q[3];
sx q[3];
rz(-2.5554113) q[3];
sx q[3];
rz(0.90914721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0236987) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(1.914631) q[2];
rz(1.6795233) q[3];
sx q[3];
rz(-1.6242124) q[3];
sx q[3];
rz(2.7000361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284477) q[0];
sx q[0];
rz(-1.0297091) q[0];
sx q[0];
rz(2.5944769) q[0];
rz(0.088134915) q[1];
sx q[1];
rz(-1.793975) q[1];
sx q[1];
rz(-0.42253447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20132658) q[0];
sx q[0];
rz(-0.28183386) q[0];
sx q[0];
rz(-1.6495709) q[0];
x q[1];
rz(-1.2014548) q[2];
sx q[2];
rz(-1.0779194) q[2];
sx q[2];
rz(-2.3702459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6531038) q[1];
sx q[1];
rz(-3.0526027) q[1];
sx q[1];
rz(-2.2527534) q[1];
x q[2];
rz(2.2219031) q[3];
sx q[3];
rz(-0.22907478) q[3];
sx q[3];
rz(-0.13815115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2253458) q[2];
sx q[2];
rz(-1.5044745) q[2];
sx q[2];
rz(-2.8543191) q[2];
rz(2.5987127) q[3];
sx q[3];
rz(-2.3714239) q[3];
sx q[3];
rz(0.058102593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3995689) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(-2.3642819) q[0];
rz(0.9264535) q[1];
sx q[1];
rz(-1.1421721) q[1];
sx q[1];
rz(-1.3949589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.004121) q[0];
sx q[0];
rz(-2.6290659) q[0];
sx q[0];
rz(0.88051535) q[0];
rz(-1.3205166) q[2];
sx q[2];
rz(-2.1553073) q[2];
sx q[2];
rz(0.78148851) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.3753739) q[1];
sx q[1];
rz(-1.9002943) q[1];
sx q[1];
rz(-1.0887515) q[1];
rz(-0.80628245) q[3];
sx q[3];
rz(-0.7976992) q[3];
sx q[3];
rz(2.9258941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7159783) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(1.0905637) q[2];
rz(-2.1770832) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(-1.9264268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5665117) q[0];
sx q[0];
rz(-0.33841857) q[0];
sx q[0];
rz(-1.6315208) q[0];
rz(3.1072726) q[1];
sx q[1];
rz(-1.8020554) q[1];
sx q[1];
rz(2.7095749) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1572947) q[0];
sx q[0];
rz(-2.108886) q[0];
sx q[0];
rz(-3.0731766) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3506579) q[2];
sx q[2];
rz(-2.7825522) q[2];
sx q[2];
rz(-0.15147789) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68709438) q[1];
sx q[1];
rz(-1.7155572) q[1];
sx q[1];
rz(-1.471038) q[1];
rz(-0.26507399) q[3];
sx q[3];
rz(-1.005583) q[3];
sx q[3];
rz(0.044818002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67046514) q[2];
sx q[2];
rz(-1.1521143) q[2];
sx q[2];
rz(1.0738037) q[2];
rz(0.18887575) q[3];
sx q[3];
rz(-2.5329068) q[3];
sx q[3];
rz(-2.7591738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0976681) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-0.64622579) q[2];
sx q[2];
rz(-0.85474174) q[2];
sx q[2];
rz(-1.3893736) q[2];
rz(-1.5021642) q[3];
sx q[3];
rz(-2.0292239) q[3];
sx q[3];
rz(-2.6440764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
