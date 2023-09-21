OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2250741) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(0.9401456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2034982) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(-3.0245568) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(0.64330805) q[1];
rz(-pi) q[2];
rz(0.75818054) q[3];
sx q[3];
rz(-1.4361793) q[3];
sx q[3];
rz(1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(-2.2629471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4449578) q[0];
sx q[0];
rz(-0.96224552) q[0];
sx q[0];
rz(0.47136013) q[0];
rz(-pi) q[1];
rz(2.5231045) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-2.6478812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6138184) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(-2.7223177) q[1];
rz(-1.2611748) q[3];
sx q[3];
rz(-1.505758) q[3];
sx q[3];
rz(1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(1.1211959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7911581) q[0];
sx q[0];
rz(-2.3179623) q[0];
sx q[0];
rz(-0.50823786) q[0];
rz(-pi) q[1];
rz(1.9636743) q[2];
sx q[2];
rz(-1.5117206) q[2];
sx q[2];
rz(-1.7977561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0908302) q[1];
sx q[1];
rz(-1.8659667) q[1];
sx q[1];
rz(-0.34543085) q[1];
rz(-1.4533914) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(0.11638164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375835) q[0];
sx q[0];
rz(-1.596367) q[0];
sx q[0];
rz(0.63347647) q[0];
rz(0.94159796) q[2];
sx q[2];
rz(-0.59855748) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0013106) q[1];
sx q[1];
rz(-1.1661068) q[1];
sx q[1];
rz(2.6005122) q[1];
x q[2];
rz(3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(-1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.1486357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0557077) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(-1.6968615) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45102851) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(0.23652467) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8685638) q[1];
sx q[1];
rz(-1.7595353) q[1];
sx q[1];
rz(0.36004685) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7025181) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(-0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.59153581) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059522) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(-0.12775001) q[0];
rz(2.3328853) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(1.8546113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1604662) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(-1.6756945) q[1];
rz(-pi) q[2];
rz(-2.0993125) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(-0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(2.8708141) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(2.0628827) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966184) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-0.78290126) q[0];
sx q[0];
rz(2.3379675) q[0];
x q[1];
rz(-0.48970512) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(-0.44039886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83963517) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(1.1595999) q[1];
x q[2];
rz(0.32254036) q[3];
sx q[3];
rz(-0.99223677) q[3];
sx q[3];
rz(-1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85212612) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(0.37129398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4079354) q[2];
sx q[2];
rz(-1.3961627) q[2];
sx q[2];
rz(2.3101431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.029433) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(0.68316858) q[1];
rz(-pi) q[2];
rz(-2.8052748) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(-1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(2.5382036) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98709479) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(-0.82018606) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7681098) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(-1.2158074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77482624) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(1.9701387) q[1];
x q[2];
rz(2.845876) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(0.28702345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-2.4323145) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-0.15670776) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126497) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(-1.2055231) q[0];
x q[1];
rz(-2.8814949) q[2];
sx q[2];
rz(-1.1220699) q[2];
sx q[2];
rz(-1.5425494) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1117489) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-0.44058056) q[1];
rz(-pi) q[2];
rz(0.79348989) q[3];
sx q[3];
rz(-1.7975382) q[3];
sx q[3];
rz(1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.043914) q[2];
sx q[2];
rz(-1.1259148) q[2];
sx q[2];
rz(-1.7150707) q[2];
rz(-2.5614212) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
