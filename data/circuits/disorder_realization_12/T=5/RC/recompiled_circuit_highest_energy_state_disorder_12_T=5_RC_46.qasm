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
rz(2.9830018) q[0];
sx q[0];
rz(-2.6245485) q[0];
sx q[0];
rz(-1.3102732) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(1.9603536) q[1];
sx q[1];
rz(11.601762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2268596) q[0];
sx q[0];
rz(-0.64382416) q[0];
sx q[0];
rz(-2.9152246) q[0];
rz(-pi) q[1];
rz(-0.50268051) q[2];
sx q[2];
rz(-0.38990228) q[2];
sx q[2];
rz(-2.1250059) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.010299461) q[1];
sx q[1];
rz(-2.2939766) q[1];
sx q[1];
rz(2.7406969) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2418666) q[3];
sx q[3];
rz(-1.6613071) q[3];
sx q[3];
rz(-0.29202005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0059119314) q[2];
sx q[2];
rz(-2.9475806) q[2];
sx q[2];
rz(-1.676959) q[2];
rz(-1.3095193) q[3];
sx q[3];
rz(-0.94683164) q[3];
sx q[3];
rz(2.8215698) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86318535) q[0];
sx q[0];
rz(-2.0018556) q[0];
sx q[0];
rz(-1.3673258) q[0];
rz(-1.4890081) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(-2.8489825) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6824043) q[0];
sx q[0];
rz(-1.4727797) q[0];
sx q[0];
rz(-2.8062264) q[0];
rz(-pi) q[1];
rz(-2.9797249) q[2];
sx q[2];
rz(-1.6965995) q[2];
sx q[2];
rz(-0.51728546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78116592) q[1];
sx q[1];
rz(-1.3349717) q[1];
sx q[1];
rz(1.7818013) q[1];
rz(-2.8752453) q[3];
sx q[3];
rz(-1.6345336) q[3];
sx q[3];
rz(1.6132465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1613529) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(2.9577067) q[2];
rz(-0.45267496) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-0.12990738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1995131) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(2.0623656) q[1];
sx q[1];
rz(-0.77003038) q[1];
sx q[1];
rz(-1.315518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7283709) q[0];
sx q[0];
rz(-2.411619) q[0];
sx q[0];
rz(0.088800207) q[0];
x q[1];
rz(-0.82683021) q[2];
sx q[2];
rz(-1.9212124) q[2];
sx q[2];
rz(-1.0817127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1444188) q[1];
sx q[1];
rz(-0.29682973) q[1];
sx q[1];
rz(1.8893912) q[1];
rz(0.27045111) q[3];
sx q[3];
rz(-0.31336774) q[3];
sx q[3];
rz(-0.53821401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2074073) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(-0.54971131) q[2];
rz(0.011119757) q[3];
sx q[3];
rz(-0.79922262) q[3];
sx q[3];
rz(1.8196677) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2928807) q[0];
sx q[0];
rz(-2.1826545) q[0];
sx q[0];
rz(-2.8380561) q[0];
rz(-0.15829463) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(-2.1319481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9859814) q[0];
sx q[0];
rz(-1.2878875) q[0];
sx q[0];
rz(-1.1213746) q[0];
x q[1];
rz(0.60229566) q[2];
sx q[2];
rz(-1.8297075) q[2];
sx q[2];
rz(0.4073172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6827439) q[1];
sx q[1];
rz(-1.7284596) q[1];
sx q[1];
rz(-0.47994061) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1032755) q[3];
sx q[3];
rz(-2.6434757) q[3];
sx q[3];
rz(1.2850873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2046854) q[2];
sx q[2];
rz(-0.82131177) q[2];
sx q[2];
rz(0.26910195) q[2];
rz(1.6875632) q[3];
sx q[3];
rz(-1.5917835) q[3];
sx q[3];
rz(1.5788186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27075416) q[0];
sx q[0];
rz(-1.0601059) q[0];
sx q[0];
rz(-0.36218542) q[0];
rz(0.65131342) q[1];
sx q[1];
rz(-1.4910699) q[1];
sx q[1];
rz(1.6168894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633879) q[0];
sx q[0];
rz(-1.7001061) q[0];
sx q[0];
rz(2.5300701) q[0];
rz(-pi) q[1];
rz(-0.082090898) q[2];
sx q[2];
rz(-0.94952784) q[2];
sx q[2];
rz(-2.4022849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.32423577) q[1];
sx q[1];
rz(-2.488616) q[1];
sx q[1];
rz(0.95545902) q[1];
x q[2];
rz(-2.397419) q[3];
sx q[3];
rz(-1.2726882) q[3];
sx q[3];
rz(-1.0013195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0321956) q[2];
sx q[2];
rz(-2.1869662) q[2];
sx q[2];
rz(-1.7677914) q[2];
rz(0.43921709) q[3];
sx q[3];
rz(-2.4924811) q[3];
sx q[3];
rz(-0.60905987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1010901) q[0];
sx q[0];
rz(-0.73223615) q[0];
sx q[0];
rz(2.8771583) q[0];
rz(-0.6791555) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(-0.98495475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0593582) q[0];
sx q[0];
rz(-1.8731887) q[0];
sx q[0];
rz(-1.2281657) q[0];
rz(-0.28552766) q[2];
sx q[2];
rz(-1.1368903) q[2];
sx q[2];
rz(-1.3615695) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90675844) q[1];
sx q[1];
rz(-1.088716) q[1];
sx q[1];
rz(-1.9831897) q[1];
rz(1.1285176) q[3];
sx q[3];
rz(-1.3849713) q[3];
sx q[3];
rz(3.1215661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(-0.43231535) q[2];
rz(-2.1281706) q[3];
sx q[3];
rz(-1.5347967) q[3];
sx q[3];
rz(-2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48563114) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(1.1095169) q[0];
rz(-0.79311496) q[1];
sx q[1];
rz(-1.3536645) q[1];
sx q[1];
rz(-1.5464334) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4594134) q[0];
sx q[0];
rz(-0.84852444) q[0];
sx q[0];
rz(-2.3514777) q[0];
x q[1];
rz(-0.23792365) q[2];
sx q[2];
rz(-0.73582651) q[2];
sx q[2];
rz(3.0032681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90518752) q[1];
sx q[1];
rz(-1.3705472) q[1];
sx q[1];
rz(2.2591641) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4647382) q[3];
sx q[3];
rz(-0.85717382) q[3];
sx q[3];
rz(2.3372363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3194797) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(1.7721843) q[2];
rz(-0.93044126) q[3];
sx q[3];
rz(-1.7934099) q[3];
sx q[3];
rz(-1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996138) q[0];
sx q[0];
rz(-1.1811341) q[0];
sx q[0];
rz(0.30366316) q[0];
rz(-0.97967255) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(0.40506515) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7295718) q[0];
sx q[0];
rz(-1.7639909) q[0];
sx q[0];
rz(2.7576819) q[0];
rz(-pi) q[1];
x q[1];
rz(1.151565) q[2];
sx q[2];
rz(-1.7287325) q[2];
sx q[2];
rz(-1.7309703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0822337) q[1];
sx q[1];
rz(-2.5463366) q[1];
sx q[1];
rz(0.31945503) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9977244) q[3];
sx q[3];
rz(-1.7550623) q[3];
sx q[3];
rz(-1.9686521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8651809) q[2];
sx q[2];
rz(-2.1750735) q[2];
sx q[2];
rz(2.442404) q[2];
rz(-0.80896038) q[3];
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
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.14178273) q[0];
sx q[0];
rz(-1.1572105) q[0];
sx q[0];
rz(-3.004177) q[0];
rz(-0.049086463) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(0.5079937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72599365) q[0];
sx q[0];
rz(-1.5440436) q[0];
sx q[0];
rz(0.49515611) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.069641308) q[2];
sx q[2];
rz(-0.17270522) q[2];
sx q[2];
rz(0.35017362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6235828) q[1];
sx q[1];
rz(-1.8933663) q[1];
sx q[1];
rz(-2.0913893) q[1];
rz(-pi) q[2];
rz(-1.6925595) q[3];
sx q[3];
rz(-1.4819488) q[3];
sx q[3];
rz(0.73032398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5892443) q[2];
sx q[2];
rz(-1.2184703) q[2];
sx q[2];
rz(-2.153896) q[2];
rz(2.6622631) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(2.2492669) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.669303) q[0];
sx q[0];
rz(-0.23858128) q[0];
sx q[0];
rz(-3.0008089) q[0];
rz(0.36551481) q[1];
sx q[1];
rz(-2.3194158) q[1];
sx q[1];
rz(0.64868322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1295531) q[0];
sx q[0];
rz(-2.3863992) q[0];
sx q[0];
rz(2.6213403) q[0];
rz(-1.4541523) q[2];
sx q[2];
rz(-0.80956993) q[2];
sx q[2];
rz(-0.67942373) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4280871) q[1];
sx q[1];
rz(-0.97555842) q[1];
sx q[1];
rz(1.6387322) q[1];
rz(0.32673959) q[3];
sx q[3];
rz(-1.800273) q[3];
sx q[3];
rz(0.47577259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6792128) q[2];
sx q[2];
rz(-1.0054192) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.4521815) q[2];
sx q[2];
rz(-1.3438423) q[2];
sx q[2];
rz(-1.9527506) q[2];
rz(2.0562511) q[3];
sx q[3];
rz(-2.3997636) q[3];
sx q[3];
rz(2.6482481) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
