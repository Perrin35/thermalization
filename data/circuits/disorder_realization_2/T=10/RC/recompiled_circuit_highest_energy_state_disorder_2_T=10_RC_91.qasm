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
rz(-0.18345565) q[0];
sx q[0];
rz(2.3975211) q[0];
sx q[0];
rz(10.476713) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(0.57300353) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0664803) q[0];
sx q[0];
rz(-1.5085992) q[0];
sx q[0];
rz(2.1960427) q[0];
rz(-pi) q[1];
rz(2.6296205) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(-1.3775502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1243368) q[1];
sx q[1];
rz(-1.5005497) q[1];
sx q[1];
rz(2.6439366) q[1];
rz(-pi) q[2];
rz(-2.4450847) q[3];
sx q[3];
rz(-1.4918431) q[3];
sx q[3];
rz(0.95206184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(1.0563043) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.2317155) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(-0.43014446) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8970222) q[0];
sx q[0];
rz(-1.8720227) q[0];
sx q[0];
rz(-1.3006849) q[0];
rz(-pi) q[1];
rz(1.6280074) q[2];
sx q[2];
rz(-1.621641) q[2];
sx q[2];
rz(-1.0982996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3900657) q[1];
sx q[1];
rz(-1.2428778) q[1];
sx q[1];
rz(-2.4310914) q[1];
rz(-pi) q[2];
rz(1.7875769) q[3];
sx q[3];
rz(-1.0415029) q[3];
sx q[3];
rz(2.2640219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0251856) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(0.44542584) q[2];
rz(-0.40924254) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(-2.4267922) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540368) q[0];
sx q[0];
rz(-1.4041413) q[0];
sx q[0];
rz(-1.0221857) q[0];
rz(0.94560854) q[2];
sx q[2];
rz(-2.473712) q[2];
sx q[2];
rz(2.5651074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2043874) q[1];
sx q[1];
rz(-1.722285) q[1];
sx q[1];
rz(-3.1207496) q[1];
rz(-pi) q[2];
rz(0.9483665) q[3];
sx q[3];
rz(-2.7223058) q[3];
sx q[3];
rz(-0.47800999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0032234) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(0.98201069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0729436) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(0.15596998) q[0];
rz(2.7929557) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(-1.9085931) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6207244) q[0];
sx q[0];
rz(-1.8308795) q[0];
sx q[0];
rz(0.12518945) q[0];
x q[1];
rz(-2.9301823) q[2];
sx q[2];
rz(-0.93358913) q[2];
sx q[2];
rz(1.5906972) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8089227) q[1];
sx q[1];
rz(-0.70007174) q[1];
sx q[1];
rz(1.1247404) q[1];
rz(1.5560914) q[3];
sx q[3];
rz(-1.3193519) q[3];
sx q[3];
rz(-0.78212839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6984581) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(1.5114463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(2.1512234) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.9885709) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(-2.8797454) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76806289) q[0];
sx q[0];
rz(-0.090001194) q[0];
sx q[0];
rz(1.9530208) q[0];
rz(2.8072281) q[2];
sx q[2];
rz(-1.0131256) q[2];
sx q[2];
rz(-1.1209436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5763229) q[1];
sx q[1];
rz(-0.63332958) q[1];
sx q[1];
rz(-2.759658) q[1];
rz(-pi) q[2];
rz(-1.366794) q[3];
sx q[3];
rz(-2.3424934) q[3];
sx q[3];
rz(-2.4287756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61949817) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(-0.74553472) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0010506823) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-0.35181272) q[0];
rz(-0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(-0.17131677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.139411) q[0];
sx q[0];
rz(-1.3277963) q[0];
sx q[0];
rz(0.67968207) q[0];
rz(-pi) q[1];
rz(-1.8753042) q[2];
sx q[2];
rz(-2.300548) q[2];
sx q[2];
rz(-2.8685121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.379516) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(-1.5980491) q[1];
rz(-0.76254179) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(1.9866379) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.6963541) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(0.0048986991) q[0];
rz(2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-2.4536224) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60058355) q[0];
sx q[0];
rz(-2.2760411) q[0];
sx q[0];
rz(1.3441104) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46779386) q[2];
sx q[2];
rz(-0.67343283) q[2];
sx q[2];
rz(1.3040257) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0258642) q[1];
sx q[1];
rz(-2.0021571) q[1];
sx q[1];
rz(-2.906745) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79598456) q[3];
sx q[3];
rz(-0.86382691) q[3];
sx q[3];
rz(-0.10315264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(2.2116057) q[2];
rz(1.2215349) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(0.090106877) q[0];
rz(-1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(-2.7065014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85283576) q[0];
sx q[0];
rz(-0.52271862) q[0];
sx q[0];
rz(0.29490348) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7516364) q[2];
sx q[2];
rz(-2.2308084) q[2];
sx q[2];
rz(-1.3446128) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1798517) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(-0.95203103) q[1];
rz(-pi) q[2];
rz(-1.5699205) q[3];
sx q[3];
rz(-1.9951576) q[3];
sx q[3];
rz(2.3800338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(0.62620658) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-0.99072376) q[3];
sx q[3];
rz(1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597252) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(-3.0873121) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(1.8064226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8322382) q[0];
sx q[0];
rz(-1.0296427) q[0];
sx q[0];
rz(0.65017976) q[0];
x q[1];
rz(0.35690709) q[2];
sx q[2];
rz(-2.1015328) q[2];
sx q[2];
rz(2.4414947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3980435) q[1];
sx q[1];
rz(-2.365592) q[1];
sx q[1];
rz(0.77568027) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0306251) q[3];
sx q[3];
rz(-2.9879192) q[3];
sx q[3];
rz(0.75017649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52428952) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(-2.3824298) q[0];
rz(-1.0836541) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(2.0738475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5006951) q[0];
sx q[0];
rz(-1.0984549) q[0];
sx q[0];
rz(0.91820902) q[0];
rz(-1.8482089) q[2];
sx q[2];
rz(-1.8545863) q[2];
sx q[2];
rz(1.0671771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0800277) q[1];
sx q[1];
rz(-2.1127709) q[1];
sx q[1];
rz(-0.47596495) q[1];
rz(0.33371933) q[3];
sx q[3];
rz(-2.2401143) q[3];
sx q[3];
rz(-3.1143318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4296253) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982211) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(-1.632985) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(2.9859424) q[3];
sx q[3];
rz(-1.8425103) q[3];
sx q[3];
rz(1.507148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
