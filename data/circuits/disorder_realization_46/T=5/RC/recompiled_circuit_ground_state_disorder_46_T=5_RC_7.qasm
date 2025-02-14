OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98193869) q[0];
sx q[0];
rz(-1.5768134) q[0];
sx q[0];
rz(-1.2793581) q[0];
rz(0.98969412) q[1];
sx q[1];
rz(5.0862105) q[1];
sx q[1];
rz(12.401019) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757614) q[0];
sx q[0];
rz(-2.1935757) q[0];
sx q[0];
rz(0.39259194) q[0];
rz(-1.1520391) q[2];
sx q[2];
rz(-1.9809696) q[2];
sx q[2];
rz(-0.23677408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27439427) q[1];
sx q[1];
rz(-1.9748678) q[1];
sx q[1];
rz(-1.0342802) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8824485) q[3];
sx q[3];
rz(-1.8168279) q[3];
sx q[3];
rz(-0.450799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87236658) q[2];
sx q[2];
rz(-1.9908315) q[2];
sx q[2];
rz(1.1671789) q[2];
rz(-2.732318) q[3];
sx q[3];
rz(-2.8661178) q[3];
sx q[3];
rz(-0.95612139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-0.15412155) q[0];
sx q[0];
rz(-1.7896205) q[0];
rz(1.4714454) q[1];
sx q[1];
rz(-2.2189326) q[1];
sx q[1];
rz(2.2460489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038747) q[0];
sx q[0];
rz(-0.26709891) q[0];
sx q[0];
rz(-3.1338723) q[0];
x q[1];
rz(-0.47426736) q[2];
sx q[2];
rz(-1.6401498) q[2];
sx q[2];
rz(-2.4276395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1481759) q[1];
sx q[1];
rz(-1.9636781) q[1];
sx q[1];
rz(2.3200289) q[1];
x q[2];
rz(-1.1510405) q[3];
sx q[3];
rz(-1.7872056) q[3];
sx q[3];
rz(-1.2866502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46513778) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(-1.4466393) q[2];
rz(-2.1220574) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(-2.8442966) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119706) q[0];
sx q[0];
rz(-1.977704) q[0];
sx q[0];
rz(-3.1373366) q[0];
rz(2.440522) q[1];
sx q[1];
rz(-2.6316167) q[1];
sx q[1];
rz(3.1140936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.741318) q[0];
sx q[0];
rz(-1.4158465) q[0];
sx q[0];
rz(-1.9727895) q[0];
rz(2.296059) q[2];
sx q[2];
rz(-1.4730244) q[2];
sx q[2];
rz(-2.1883983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3361275) q[1];
sx q[1];
rz(-1.256811) q[1];
sx q[1];
rz(-2.9538395) q[1];
rz(3.1276567) q[3];
sx q[3];
rz(-1.2346141) q[3];
sx q[3];
rz(-0.61830904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6444401) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(0.39786097) q[2];
rz(-1.0128939) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(0.40867543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560646) q[0];
sx q[0];
rz(-1.8714454) q[0];
sx q[0];
rz(0.84554607) q[0];
rz(0.0088508765) q[1];
sx q[1];
rz(-2.2933941) q[1];
sx q[1];
rz(1.1525851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801013) q[0];
sx q[0];
rz(-1.9753755) q[0];
sx q[0];
rz(2.220551) q[0];
rz(-pi) q[1];
rz(-2.0699224) q[2];
sx q[2];
rz(-0.50412699) q[2];
sx q[2];
rz(-1.1386516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2680451) q[1];
sx q[1];
rz(-1.598804) q[1];
sx q[1];
rz(0.59240474) q[1];
rz(-2.939393) q[3];
sx q[3];
rz(-0.31980896) q[3];
sx q[3];
rz(-0.6347189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74959603) q[2];
sx q[2];
rz(-1.0276724) q[2];
sx q[2];
rz(-1.9409625) q[2];
rz(0.47731733) q[3];
sx q[3];
rz(-2.6715607) q[3];
sx q[3];
rz(2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93477997) q[0];
sx q[0];
rz(-2.2581357) q[0];
sx q[0];
rz(2.2373037) q[0];
rz(-0.9710871) q[1];
sx q[1];
rz(-1.9458408) q[1];
sx q[1];
rz(-2.8579874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0440031) q[0];
sx q[0];
rz(-1.6416802) q[0];
sx q[0];
rz(-1.5868091) q[0];
rz(2.3432557) q[2];
sx q[2];
rz(-2.5588648) q[2];
sx q[2];
rz(-0.59205627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3963822) q[1];
sx q[1];
rz(-1.3644427) q[1];
sx q[1];
rz(-2.2554805) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6568703) q[3];
sx q[3];
rz(-2.1987763) q[3];
sx q[3];
rz(2.064765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7859555) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(2.9034485) q[2];
rz(-2.1613878) q[3];
sx q[3];
rz(-1.9832059) q[3];
sx q[3];
rz(2.3908206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439483) q[0];
sx q[0];
rz(-1.5357966) q[0];
sx q[0];
rz(-0.407298) q[0];
rz(2.7382355) q[1];
sx q[1];
rz(-2.401001) q[1];
sx q[1];
rz(-0.86722803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2326538) q[0];
sx q[0];
rz(-2.5728658) q[0];
sx q[0];
rz(1.6816116) q[0];
rz(-pi) q[1];
rz(-1.3620141) q[2];
sx q[2];
rz(-2.1174413) q[2];
sx q[2];
rz(-2.7289313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85485172) q[1];
sx q[1];
rz(-0.84567243) q[1];
sx q[1];
rz(-1.9560567) q[1];
rz(-1.2739185) q[3];
sx q[3];
rz(-1.6921077) q[3];
sx q[3];
rz(-2.6685126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0200218) q[2];
sx q[2];
rz(-0.82841221) q[2];
sx q[2];
rz(1.7106445) q[2];
rz(0.53267789) q[3];
sx q[3];
rz(-1.0360274) q[3];
sx q[3];
rz(-1.7238341) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.983404) q[0];
sx q[0];
rz(-2.991365) q[0];
sx q[0];
rz(-1.0839373) q[0];
rz(-0.051941959) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(-3.1165677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1442102) q[0];
sx q[0];
rz(-0.48244993) q[0];
sx q[0];
rz(-0.061003322) q[0];
rz(2.04966) q[2];
sx q[2];
rz(-0.81170481) q[2];
sx q[2];
rz(2.1213437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4551291) q[1];
sx q[1];
rz(-1.741998) q[1];
sx q[1];
rz(-1.135395) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42633485) q[3];
sx q[3];
rz(-1.4407743) q[3];
sx q[3];
rz(-3.1257747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0673125) q[2];
sx q[2];
rz(-1.4558027) q[2];
sx q[2];
rz(2.32453) q[2];
rz(0.95064154) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0364712) q[0];
sx q[0];
rz(-0.28577411) q[0];
sx q[0];
rz(-0.60086077) q[0];
rz(0.034424456) q[1];
sx q[1];
rz(-1.8079146) q[1];
sx q[1];
rz(1.6011802) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9118155) q[0];
sx q[0];
rz(-0.86965771) q[0];
sx q[0];
rz(2.3100353) q[0];
x q[1];
rz(0.73141269) q[2];
sx q[2];
rz(-2.3225975) q[2];
sx q[2];
rz(-2.7186944) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34791495) q[1];
sx q[1];
rz(-1.131212) q[1];
sx q[1];
rz(-1.5949834) q[1];
x q[2];
rz(0.65798379) q[3];
sx q[3];
rz(-0.65449981) q[3];
sx q[3];
rz(-2.2006048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2711266) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(1.4938483) q[2];
rz(-2.3624524) q[3];
sx q[3];
rz(-1.4804163) q[3];
sx q[3];
rz(2.486865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036024) q[0];
sx q[0];
rz(-2.0605189) q[0];
sx q[0];
rz(-2.3533452) q[0];
rz(1.2429271) q[1];
sx q[1];
rz(-1.582076) q[1];
sx q[1];
rz(-2.8020614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718455) q[0];
sx q[0];
rz(-1.7632428) q[0];
sx q[0];
rz(-3.0234973) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8571157) q[2];
sx q[2];
rz(-1.0272046) q[2];
sx q[2];
rz(-1.7329777) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.371468) q[1];
sx q[1];
rz(-1.1684019) q[1];
sx q[1];
rz(-1.2312908) q[1];
rz(-pi) q[2];
rz(-0.65776627) q[3];
sx q[3];
rz(-2.2493746) q[3];
sx q[3];
rz(-0.15746169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.01123151) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.2569024) q[2];
rz(2.0508749) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(2.2244661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.286769) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(2.0538034) q[0];
rz(0.65025672) q[1];
sx q[1];
rz(-1.6842664) q[1];
sx q[1];
rz(0.021171721) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41474671) q[0];
sx q[0];
rz(-2.0170037) q[0];
sx q[0];
rz(-1.3394974) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7672054) q[2];
sx q[2];
rz(-1.1078664) q[2];
sx q[2];
rz(-0.48505515) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52439201) q[1];
sx q[1];
rz(-0.78689303) q[1];
sx q[1];
rz(1.2715497) q[1];
rz(-2.1032501) q[3];
sx q[3];
rz(-1.5335011) q[3];
sx q[3];
rz(1.2398014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0040969) q[2];
sx q[2];
rz(-1.7986412) q[2];
sx q[2];
rz(-0.98319483) q[2];
rz(-0.99299562) q[3];
sx q[3];
rz(-2.0296622) q[3];
sx q[3];
rz(-0.23428169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642333) q[0];
sx q[0];
rz(-2.3491884) q[0];
sx q[0];
rz(1.0607251) q[0];
rz(2.0296774) q[1];
sx q[1];
rz(-2.8363375) q[1];
sx q[1];
rz(-3.0189966) q[1];
rz(-2.0531102) q[2];
sx q[2];
rz(-1.7822722) q[2];
sx q[2];
rz(1.7581802) q[2];
rz(-0.52295104) q[3];
sx q[3];
rz(-1.5802438) q[3];
sx q[3];
rz(2.6680744) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
