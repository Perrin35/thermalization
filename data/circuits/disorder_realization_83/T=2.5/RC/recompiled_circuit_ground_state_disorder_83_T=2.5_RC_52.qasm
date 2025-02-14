OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72252005) q[0];
sx q[0];
rz(2.4604586) q[0];
sx q[0];
rz(13.606986) q[0];
rz(-0.71086565) q[1];
sx q[1];
rz(-0.65029538) q[1];
sx q[1];
rz(-1.8898036) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8914362) q[0];
sx q[0];
rz(-0.52816872) q[0];
sx q[0];
rz(-2.3528966) q[0];
rz(-2.5898143) q[2];
sx q[2];
rz(-1.1349196) q[2];
sx q[2];
rz(0.18407719) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7573479) q[1];
sx q[1];
rz(-0.66435821) q[1];
sx q[1];
rz(0.039907736) q[1];
x q[2];
rz(-0.85504882) q[3];
sx q[3];
rz(-2.2576393) q[3];
sx q[3];
rz(-0.035720197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9616482) q[2];
sx q[2];
rz(-1.0893931) q[2];
sx q[2];
rz(-1.3570448) q[2];
rz(-2.0876136) q[3];
sx q[3];
rz(-1.6097924) q[3];
sx q[3];
rz(2.0747144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.2422975) q[0];
sx q[0];
rz(-0.59277788) q[0];
sx q[0];
rz(1.584126) q[0];
rz(-1.8862995) q[1];
sx q[1];
rz(-0.86304945) q[1];
sx q[1];
rz(3.0994298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52707878) q[0];
sx q[0];
rz(-1.9111484) q[0];
sx q[0];
rz(-0.69410364) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0980646) q[2];
sx q[2];
rz(-1.0587436) q[2];
sx q[2];
rz(1.5498424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.592432) q[1];
sx q[1];
rz(-1.7322408) q[1];
sx q[1];
rz(-1.8668951) q[1];
rz(-pi) q[2];
rz(1.0125431) q[3];
sx q[3];
rz(-1.7413845) q[3];
sx q[3];
rz(2.0743897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92370799) q[2];
sx q[2];
rz(-1.2486685) q[2];
sx q[2];
rz(-0.62151796) q[2];
rz(0.38226852) q[3];
sx q[3];
rz(-1.1416124) q[3];
sx q[3];
rz(-0.83466616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37654787) q[0];
sx q[0];
rz(-1.5060197) q[0];
sx q[0];
rz(1.3733093) q[0];
rz(0.37257591) q[1];
sx q[1];
rz(-2.617045) q[1];
sx q[1];
rz(3.0212044) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62154462) q[0];
sx q[0];
rz(-2.1839474) q[0];
sx q[0];
rz(0.91213062) q[0];
x q[1];
rz(-1.7104501) q[2];
sx q[2];
rz(-2.8595028) q[2];
sx q[2];
rz(2.8733746) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38594151) q[1];
sx q[1];
rz(-2.3214968) q[1];
sx q[1];
rz(2.7044317) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1266842) q[3];
sx q[3];
rz(-0.57775195) q[3];
sx q[3];
rz(1.8803949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7239712) q[2];
sx q[2];
rz(-1.3105023) q[2];
sx q[2];
rz(1.2869147) q[2];
rz(0.78479615) q[3];
sx q[3];
rz(-2.0311821) q[3];
sx q[3];
rz(-0.77814656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6374636) q[0];
sx q[0];
rz(-2.6478719) q[0];
sx q[0];
rz(2.8959287) q[0];
rz(3.1371112) q[1];
sx q[1];
rz(-0.5537529) q[1];
sx q[1];
rz(-0.063668879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144062) q[0];
sx q[0];
rz(-0.92824358) q[0];
sx q[0];
rz(1.8802934) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9826239) q[2];
sx q[2];
rz(-2.5971488) q[2];
sx q[2];
rz(3.009764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6323327) q[1];
sx q[1];
rz(-0.89069388) q[1];
sx q[1];
rz(-2.9010309) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40157179) q[3];
sx q[3];
rz(-2.662559) q[3];
sx q[3];
rz(-3.1216321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8513546) q[2];
sx q[2];
rz(-1.3685702) q[2];
sx q[2];
rz(-1.2151388) q[2];
rz(-0.75567192) q[3];
sx q[3];
rz(-0.98545051) q[3];
sx q[3];
rz(0.22172609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3168685) q[0];
sx q[0];
rz(-1.2687954) q[0];
sx q[0];
rz(-0.49279898) q[0];
rz(-2.1700676) q[1];
sx q[1];
rz(-1.3243472) q[1];
sx q[1];
rz(2.9528101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014009091) q[0];
sx q[0];
rz(-0.32453093) q[0];
sx q[0];
rz(0.55248673) q[0];
rz(-pi) q[1];
rz(0.014585693) q[2];
sx q[2];
rz(-1.5684394) q[2];
sx q[2];
rz(-2.6429526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9663736) q[1];
sx q[1];
rz(-1.0639251) q[1];
sx q[1];
rz(-3.0543249) q[1];
rz(-pi) q[2];
rz(-2.4972992) q[3];
sx q[3];
rz(-1.7731993) q[3];
sx q[3];
rz(-2.6930893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8652953) q[2];
sx q[2];
rz(-1.6388288) q[2];
sx q[2];
rz(0.067528188) q[2];
rz(-1.4589795) q[3];
sx q[3];
rz(-2.429481) q[3];
sx q[3];
rz(-2.0290831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825298) q[0];
sx q[0];
rz(-1.9628061) q[0];
sx q[0];
rz(-1.4558526) q[0];
rz(0.23446941) q[1];
sx q[1];
rz(-0.5916943) q[1];
sx q[1];
rz(0.04105982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0309248) q[0];
sx q[0];
rz(-2.2718856) q[0];
sx q[0];
rz(-1.5337471) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1424871) q[2];
sx q[2];
rz(-1.6950144) q[2];
sx q[2];
rz(1.6060843) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8248661) q[1];
sx q[1];
rz(-1.1371964) q[1];
sx q[1];
rz(1.0528013) q[1];
rz(-pi) q[2];
rz(-1.4047296) q[3];
sx q[3];
rz(-0.73544805) q[3];
sx q[3];
rz(-2.0388132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2264651) q[2];
sx q[2];
rz(-0.76639405) q[2];
sx q[2];
rz(2.8506632) q[2];
rz(2.9229524) q[3];
sx q[3];
rz(-0.69059697) q[3];
sx q[3];
rz(-2.2712928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640294) q[0];
sx q[0];
rz(-0.95528269) q[0];
sx q[0];
rz(0.56354228) q[0];
rz(-2.3130985) q[1];
sx q[1];
rz(-0.68873134) q[1];
sx q[1];
rz(0.73961893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002454438) q[0];
sx q[0];
rz(-2.1147303) q[0];
sx q[0];
rz(0.55940658) q[0];
x q[1];
rz(-1.2752259) q[2];
sx q[2];
rz(-0.53662125) q[2];
sx q[2];
rz(1.581096) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.0060192903) q[1];
sx q[1];
rz(-1.3050627) q[1];
sx q[1];
rz(-1.0510848) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4007934) q[3];
sx q[3];
rz(-2.1925266) q[3];
sx q[3];
rz(-0.19843693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9699817) q[2];
sx q[2];
rz(-0.7526528) q[2];
sx q[2];
rz(2.8438711) q[2];
rz(-0.026084829) q[3];
sx q[3];
rz(-1.1931158) q[3];
sx q[3];
rz(-2.397876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.17224269) q[0];
sx q[0];
rz(-1.5650711) q[0];
sx q[0];
rz(1.7283537) q[0];
rz(-0.50624943) q[1];
sx q[1];
rz(-1.0289611) q[1];
sx q[1];
rz(1.1594353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6613061) q[0];
sx q[0];
rz(-1.4692069) q[0];
sx q[0];
rz(0.30509588) q[0];
rz(-pi) q[1];
rz(1.5250823) q[2];
sx q[2];
rz(-0.9897784) q[2];
sx q[2];
rz(2.8948262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22190801) q[1];
sx q[1];
rz(-1.8076928) q[1];
sx q[1];
rz(-2.6568031) q[1];
rz(3.0049802) q[3];
sx q[3];
rz(-2.3834565) q[3];
sx q[3];
rz(-1.6389169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1737698) q[2];
sx q[2];
rz(-2.2161667) q[2];
sx q[2];
rz(1.7769163) q[2];
rz(-0.81554282) q[3];
sx q[3];
rz(-1.2929076) q[3];
sx q[3];
rz(1.7447785) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57104617) q[0];
sx q[0];
rz(-2.1111574) q[0];
sx q[0];
rz(-1.1676769) q[0];
rz(-0.22444621) q[1];
sx q[1];
rz(-2.3851604) q[1];
sx q[1];
rz(-2.6969553) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8380008) q[0];
sx q[0];
rz(-1.5193257) q[0];
sx q[0];
rz(-0.055368467) q[0];
rz(-pi) q[1];
rz(2.2811057) q[2];
sx q[2];
rz(-1.0376467) q[2];
sx q[2];
rz(1.8955961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2336572) q[1];
sx q[1];
rz(-0.3668712) q[1];
sx q[1];
rz(-1.1790102) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4280714) q[3];
sx q[3];
rz(-1.1253353) q[3];
sx q[3];
rz(-0.023825432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4148407) q[2];
sx q[2];
rz(-1.3123935) q[2];
sx q[2];
rz(-1.4078338) q[2];
rz(-1.6057181) q[3];
sx q[3];
rz(-1.1955669) q[3];
sx q[3];
rz(2.2685952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6294412) q[0];
sx q[0];
rz(-0.81708556) q[0];
sx q[0];
rz(2.0508811) q[0];
rz(1.0461944) q[1];
sx q[1];
rz(-0.9895784) q[1];
sx q[1];
rz(-0.88669056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42536286) q[0];
sx q[0];
rz(-1.9337448) q[0];
sx q[0];
rz(-1.1858398) q[0];
rz(2.1616814) q[2];
sx q[2];
rz(-1.3279898) q[2];
sx q[2];
rz(-1.7942015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1613933) q[1];
sx q[1];
rz(-0.43159136) q[1];
sx q[1];
rz(3.1100247) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90417741) q[3];
sx q[3];
rz(-1.2632367) q[3];
sx q[3];
rz(-2.9136861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0118759) q[2];
sx q[2];
rz(-2.1265714) q[2];
sx q[2];
rz(-1.2474308) q[2];
rz(0.85793197) q[3];
sx q[3];
rz(-1.4638289) q[3];
sx q[3];
rz(-1.5425382) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7958551) q[0];
sx q[0];
rz(-1.7248187) q[0];
sx q[0];
rz(3.0378573) q[0];
rz(-0.70032447) q[1];
sx q[1];
rz(-1.2979957) q[1];
sx q[1];
rz(-1.8630984) q[1];
rz(1.4486031) q[2];
sx q[2];
rz(-1.3753825) q[2];
sx q[2];
rz(-1.5254979) q[2];
rz(0.35459749) q[3];
sx q[3];
rz(-1.1413367) q[3];
sx q[3];
rz(-1.8324112) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
