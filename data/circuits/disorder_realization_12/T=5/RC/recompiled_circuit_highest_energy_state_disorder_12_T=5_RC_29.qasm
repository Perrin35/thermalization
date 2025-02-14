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
rz(-0.15859088) q[0];
sx q[0];
rz(-0.51704419) q[0];
sx q[0];
rz(1.3102732) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(-2.1769843) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83822891) q[0];
sx q[0];
rz(-1.705929) q[0];
sx q[0];
rz(-0.63146308) q[0];
x q[1];
rz(2.7959441) q[2];
sx q[2];
rz(-1.7549577) q[2];
sx q[2];
rz(-3.057827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8339616) q[1];
sx q[1];
rz(-1.8677068) q[1];
sx q[1];
rz(-2.3351257) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4258641) q[3];
sx q[3];
rz(-0.67620899) q[3];
sx q[3];
rz(-1.165426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1356807) q[2];
sx q[2];
rz(-0.19401208) q[2];
sx q[2];
rz(1.676959) q[2];
rz(1.3095193) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(0.32002282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784073) q[0];
sx q[0];
rz(-1.139737) q[0];
sx q[0];
rz(1.7742668) q[0];
rz(-1.4890081) q[1];
sx q[1];
rz(-0.79843489) q[1];
sx q[1];
rz(2.8489825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83787888) q[0];
sx q[0];
rz(-2.7927164) q[0];
sx q[0];
rz(0.29033355) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66537036) q[2];
sx q[2];
rz(-0.20466802) q[2];
sx q[2];
rz(2.7432347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83963385) q[1];
sx q[1];
rz(-1.7758766) q[1];
sx q[1];
rz(-0.2409711) q[1];
x q[2];
rz(1.6368565) q[3];
sx q[3];
rz(-1.8365897) q[3];
sx q[3];
rz(3.0817666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9802398) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(0.18388595) q[2];
rz(-2.6889177) q[3];
sx q[3];
rz(-1.6190745) q[3];
sx q[3];
rz(3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1995131) q[0];
sx q[0];
rz(-2.4266854) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(-1.0792271) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(-1.315518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41322177) q[0];
sx q[0];
rz(-2.411619) q[0];
sx q[0];
rz(0.088800207) q[0];
rz(-pi) q[1];
rz(2.6805515) q[2];
sx q[2];
rz(-0.88141841) q[2];
sx q[2];
rz(-0.79511666) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9971738) q[1];
sx q[1];
rz(-2.8447629) q[1];
sx q[1];
rz(1.2522015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30267164) q[3];
sx q[3];
rz(-1.4883452) q[3];
sx q[3];
rz(2.3668805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2074073) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(-2.5918813) q[2];
rz(3.1304729) q[3];
sx q[3];
rz(-0.79922262) q[3];
sx q[3];
rz(1.3219249) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84871197) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(2.8380561) q[0];
rz(0.15829463) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(2.1319481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5925668) q[0];
sx q[0];
rz(-1.1404622) q[0];
sx q[0];
rz(-0.31220147) q[0];
rz(-pi) q[1];
rz(-2.539297) q[2];
sx q[2];
rz(-1.8297075) q[2];
sx q[2];
rz(-2.7342755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19349603) q[1];
sx q[1];
rz(-1.0973012) q[1];
sx q[1];
rz(1.3934474) q[1];
x q[2];
rz(2.6437838) q[3];
sx q[3];
rz(-1.5890997) q[3];
sx q[3];
rz(0.25204424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2046854) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(-2.8724907) q[2];
rz(-1.6875632) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(-1.5627741) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27075416) q[0];
sx q[0];
rz(-2.0814867) q[0];
sx q[0];
rz(-2.7794072) q[0];
rz(-0.65131342) q[1];
sx q[1];
rz(-1.4910699) q[1];
sx q[1];
rz(-1.6168894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633879) q[0];
sx q[0];
rz(-1.4414865) q[0];
sx q[0];
rz(0.61152258) q[0];
rz(0.94793041) q[2];
sx q[2];
rz(-1.6375223) q[2];
sx q[2];
rz(-0.78363505) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73495364) q[1];
sx q[1];
rz(-1.2124774) q[1];
sx q[1];
rz(1.0124769) q[1];
x q[2];
rz(2.7157341) q[3];
sx q[3];
rz(-0.79090624) q[3];
sx q[3];
rz(-0.87825852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0321956) q[2];
sx q[2];
rz(-0.95462644) q[2];
sx q[2];
rz(1.3738013) q[2];
rz(2.7023756) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(-0.60905987) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1010901) q[0];
sx q[0];
rz(-2.4093565) q[0];
sx q[0];
rz(-2.8771583) q[0];
rz(2.4624372) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(-0.98495475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59439182) q[0];
sx q[0];
rz(-1.8972881) q[0];
sx q[0];
rz(0.31983967) q[0];
rz(-0.28552766) q[2];
sx q[2];
rz(-1.1368903) q[2];
sx q[2];
rz(1.7800231) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15036525) q[1];
sx q[1];
rz(-2.5179407) q[1];
sx q[1];
rz(-0.6536478) q[1];
x q[2];
rz(1.1285176) q[3];
sx q[3];
rz(-1.3849713) q[3];
sx q[3];
rz(3.1215661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(0.43231535) q[2];
rz(-2.1281706) q[3];
sx q[3];
rz(-1.6067959) q[3];
sx q[3];
rz(-0.54949808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6559615) q[0];
sx q[0];
rz(-2.7483342) q[0];
sx q[0];
rz(2.0320758) q[0];
rz(-0.79311496) q[1];
sx q[1];
rz(-1.7879281) q[1];
sx q[1];
rz(-1.5951593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6731747) q[0];
sx q[0];
rz(-2.1270848) q[0];
sx q[0];
rz(-0.89222096) q[0];
rz(-pi) q[1];
rz(-1.3605453) q[2];
sx q[2];
rz(-0.86019197) q[2];
sx q[2];
rz(-2.9637314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3137957) q[1];
sx q[1];
rz(-0.89875752) q[1];
sx q[1];
rz(-0.25700493) q[1];
x q[2];
rz(3.0199354) q[3];
sx q[3];
rz(-2.4215048) q[3];
sx q[3];
rz(2.4984604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3194797) q[2];
sx q[2];
rz(-2.959368) q[2];
sx q[2];
rz(-1.7721843) q[2];
rz(-0.93044126) q[3];
sx q[3];
rz(-1.3481827) q[3];
sx q[3];
rz(1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.1996138) q[0];
sx q[0];
rz(-1.1811341) q[0];
sx q[0];
rz(2.8379295) q[0];
rz(0.97967255) q[1];
sx q[1];
rz(-2.517608) q[1];
sx q[1];
rz(0.40506515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2361778) q[0];
sx q[0];
rz(-1.9472031) q[0];
sx q[0];
rz(-1.7787399) q[0];
x q[1];
rz(1.151565) q[2];
sx q[2];
rz(-1.7287325) q[2];
sx q[2];
rz(-1.7309703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0593589) q[1];
sx q[1];
rz(-0.59525604) q[1];
sx q[1];
rz(2.8221376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1478506) q[3];
sx q[3];
rz(-2.6788524) q[3];
sx q[3];
rz(3.1266318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8651809) q[2];
sx q[2];
rz(-2.1750735) q[2];
sx q[2];
rz(2.442404) q[2];
rz(2.3326323) q[3];
sx q[3];
rz(-1.8374846) q[3];
sx q[3];
rz(2.8185524) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9998099) q[0];
sx q[0];
rz(-1.9843822) q[0];
sx q[0];
rz(-0.13741563) q[0];
rz(0.049086463) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(-0.5079937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415599) q[0];
sx q[0];
rz(-1.5440436) q[0];
sx q[0];
rz(0.49515611) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069641308) q[2];
sx q[2];
rz(-0.17270522) q[2];
sx q[2];
rz(2.791419) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6235828) q[1];
sx q[1];
rz(-1.8933663) q[1];
sx q[1];
rz(-2.0913893) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2035962) q[3];
sx q[3];
rz(-0.15060234) q[3];
sx q[3];
rz(-1.6737398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5892443) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(-2.153896) q[2];
rz(0.47932953) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(-2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4722897) q[0];
sx q[0];
rz(-2.9030114) q[0];
sx q[0];
rz(-0.14078374) q[0];
rz(-0.36551481) q[1];
sx q[1];
rz(-0.82217685) q[1];
sx q[1];
rz(-2.4929094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539584) q[0];
sx q[0];
rz(-1.2231069) q[0];
sx q[0];
rz(-0.68490048) q[0];
rz(-pi) q[1];
rz(-3.0200483) q[2];
sx q[2];
rz(-0.76833188) q[2];
sx q[2];
rz(-2.2939081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71350559) q[1];
sx q[1];
rz(-0.97555842) q[1];
sx q[1];
rz(1.5028605) q[1];
rz(-0.62913904) q[3];
sx q[3];
rz(-0.39689343) q[3];
sx q[3];
rz(-1.6861738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4623798) q[2];
sx q[2];
rz(-2.1361735) q[2];
sx q[2];
rz(0.072754808) q[2];
rz(-0.66271979) q[3];
sx q[3];
rz(-1.3136256) q[3];
sx q[3];
rz(-0.16547468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872021) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(1.5351334) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(-2.6680505) q[2];
sx q[2];
rz(-2.8859856) q[2];
sx q[2];
rz(-2.4398266) q[2];
rz(0.88964626) q[3];
sx q[3];
rz(-1.2500661) q[3];
sx q[3];
rz(1.4483857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
