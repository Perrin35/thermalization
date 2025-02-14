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
rz(0.17231365) q[0];
sx q[0];
rz(3.3495164) q[0];
sx q[0];
rz(10.715545) q[0];
rz(-2.9895904) q[1];
sx q[1];
rz(-0.64270371) q[1];
sx q[1];
rz(2.7132577) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9865532) q[0];
sx q[0];
rz(-1.4615834) q[0];
sx q[0];
rz(-0.10515736) q[0];
rz(2.4691888) q[2];
sx q[2];
rz(-1.3383631) q[2];
sx q[2];
rz(-1.642579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2178035) q[1];
sx q[1];
rz(-0.83309595) q[1];
sx q[1];
rz(1.6290725) q[1];
rz(-pi) q[2];
rz(-2.9228404) q[3];
sx q[3];
rz(-1.3748589) q[3];
sx q[3];
rz(0.71813717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20994818) q[2];
sx q[2];
rz(-0.95982176) q[2];
sx q[2];
rz(-1.9317365) q[2];
rz(-3.0620388) q[3];
sx q[3];
rz(-2.7456561) q[3];
sx q[3];
rz(-0.52566665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.88456589) q[0];
sx q[0];
rz(-0.50978065) q[0];
sx q[0];
rz(2.7575745) q[0];
rz(-0.89029038) q[1];
sx q[1];
rz(-1.9046116) q[1];
sx q[1];
rz(0.30337897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10643364) q[0];
sx q[0];
rz(-1.037853) q[0];
sx q[0];
rz(-0.22179166) q[0];
x q[1];
rz(0.25945865) q[2];
sx q[2];
rz(-1.0676117) q[2];
sx q[2];
rz(-0.34215701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1002101) q[1];
sx q[1];
rz(-3.0315371) q[1];
sx q[1];
rz(-2.1724114) q[1];
rz(1.1214662) q[3];
sx q[3];
rz(-2.1746965) q[3];
sx q[3];
rz(0.7627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7552135) q[2];
sx q[2];
rz(-2.0342125) q[2];
sx q[2];
rz(2.4066822) q[2];
rz(0.84878659) q[3];
sx q[3];
rz(-0.74257344) q[3];
sx q[3];
rz(-2.7809704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3291149) q[0];
sx q[0];
rz(-2.766093) q[0];
sx q[0];
rz(0.078711674) q[0];
rz(0.82376897) q[1];
sx q[1];
rz(-1.0853826) q[1];
sx q[1];
rz(-1.8471921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1549653) q[0];
sx q[0];
rz(-1.6324348) q[0];
sx q[0];
rz(0.060527965) q[0];
rz(2.1903388) q[2];
sx q[2];
rz(-1.143537) q[2];
sx q[2];
rz(-0.45742971) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6523167) q[1];
sx q[1];
rz(-1.4536263) q[1];
sx q[1];
rz(-2.8588309) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8494959) q[3];
sx q[3];
rz(-0.72360814) q[3];
sx q[3];
rz(-1.1147182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2417629) q[2];
sx q[2];
rz(-0.37909847) q[2];
sx q[2];
rz(-2.3251593) q[2];
rz(1.624931) q[3];
sx q[3];
rz(-2.1818325) q[3];
sx q[3];
rz(2.4412156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4888332) q[0];
sx q[0];
rz(-2.3214898) q[0];
sx q[0];
rz(-0.56831992) q[0];
rz(-1.142451) q[1];
sx q[1];
rz(-1.7894952) q[1];
sx q[1];
rz(1.6894587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84684337) q[0];
sx q[0];
rz(-0.94943014) q[0];
sx q[0];
rz(-1.2970379) q[0];
x q[1];
rz(2.6645489) q[2];
sx q[2];
rz(-2.9397268) q[2];
sx q[2];
rz(-2.9143726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4751079) q[1];
sx q[1];
rz(-2.1748073) q[1];
sx q[1];
rz(2.758065) q[1];
x q[2];
rz(1.0978974) q[3];
sx q[3];
rz(-2.7644025) q[3];
sx q[3];
rz(-2.8998714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51231724) q[2];
sx q[2];
rz(-2.6111111) q[2];
sx q[2];
rz(0.076889195) q[2];
rz(0.39938375) q[3];
sx q[3];
rz(-0.9136343) q[3];
sx q[3];
rz(-2.5975749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15453108) q[0];
sx q[0];
rz(-2.7847544) q[0];
sx q[0];
rz(-0.29990184) q[0];
rz(1.9918282) q[1];
sx q[1];
rz(-1.0118142) q[1];
sx q[1];
rz(-2.6531175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804084) q[0];
sx q[0];
rz(-1.6035491) q[0];
sx q[0];
rz(1.7561654) q[0];
rz(1.064642) q[2];
sx q[2];
rz(-1.0380259) q[2];
sx q[2];
rz(-1.2974933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3255182) q[1];
sx q[1];
rz(-1.5518922) q[1];
sx q[1];
rz(2.5579648) q[1];
rz(-pi) q[2];
rz(0.3533479) q[3];
sx q[3];
rz(-0.67253695) q[3];
sx q[3];
rz(-1.3143244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38146314) q[2];
sx q[2];
rz(-3.0035786) q[2];
sx q[2];
rz(-1.4786973) q[2];
rz(2.0315157) q[3];
sx q[3];
rz(-0.9980945) q[3];
sx q[3];
rz(-2.4983675) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396624) q[0];
sx q[0];
rz(-1.9597541) q[0];
sx q[0];
rz(2.8231743) q[0];
rz(-0.034612522) q[1];
sx q[1];
rz(-2.5739659) q[1];
sx q[1];
rz(-0.45077032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50758963) q[0];
sx q[0];
rz(-1.5566155) q[0];
sx q[0];
rz(-2.0812278) q[0];
x q[1];
rz(-2.4813406) q[2];
sx q[2];
rz(-2.1532032) q[2];
sx q[2];
rz(-0.16124111) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8776226) q[1];
sx q[1];
rz(-1.1671986) q[1];
sx q[1];
rz(0.91124673) q[1];
x q[2];
rz(0.81574635) q[3];
sx q[3];
rz(-1.295608) q[3];
sx q[3];
rz(-1.311917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0938809) q[2];
sx q[2];
rz(-2.9517951) q[2];
sx q[2];
rz(0.06165687) q[2];
rz(3.0430074) q[3];
sx q[3];
rz(-2.3969789) q[3];
sx q[3];
rz(-2.4390167) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75673574) q[0];
sx q[0];
rz(-2.0282133) q[0];
sx q[0];
rz(-0.039948832) q[0];
rz(2.3710251) q[1];
sx q[1];
rz(-0.71208411) q[1];
sx q[1];
rz(-2.9133453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9271999) q[0];
sx q[0];
rz(-1.6480371) q[0];
sx q[0];
rz(1.4820251) q[0];
rz(-pi) q[1];
rz(1.575861) q[2];
sx q[2];
rz(-0.3061848) q[2];
sx q[2];
rz(-1.926946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7668263) q[1];
sx q[1];
rz(-0.09775459) q[1];
sx q[1];
rz(-0.575211) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2049115) q[3];
sx q[3];
rz(-1.2290658) q[3];
sx q[3];
rz(-0.069119819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.063529) q[2];
sx q[2];
rz(-2.0827796) q[2];
sx q[2];
rz(2.883319) q[2];
rz(0.20049788) q[3];
sx q[3];
rz(-0.79810464) q[3];
sx q[3];
rz(-2.958278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16266009) q[0];
sx q[0];
rz(-1.3717835) q[0];
sx q[0];
rz(-0.3072511) q[0];
rz(1.2112674) q[1];
sx q[1];
rz(-2.6478421) q[1];
sx q[1];
rz(-2.5603851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085614) q[0];
sx q[0];
rz(-1.1662467) q[0];
sx q[0];
rz(3.1076215) q[0];
rz(-0.18980726) q[2];
sx q[2];
rz(-0.71868374) q[2];
sx q[2];
rz(0.41658066) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4631066) q[1];
sx q[1];
rz(-2.3358279) q[1];
sx q[1];
rz(-1.8496978) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56116207) q[3];
sx q[3];
rz(-0.61120874) q[3];
sx q[3];
rz(1.2964013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6792949) q[2];
sx q[2];
rz(-1.1335979) q[2];
sx q[2];
rz(-0.031166859) q[2];
rz(-2.9160685) q[3];
sx q[3];
rz(-1.7799107) q[3];
sx q[3];
rz(2.1112736) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088260055) q[0];
sx q[0];
rz(-3.0971165) q[0];
sx q[0];
rz(0.70892507) q[0];
rz(0.21253474) q[1];
sx q[1];
rz(-2.1268851) q[1];
sx q[1];
rz(2.2841891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014985059) q[0];
sx q[0];
rz(-1.4410271) q[0];
sx q[0];
rz(-0.0051104498) q[0];
rz(2.6231758) q[2];
sx q[2];
rz(-0.5974434) q[2];
sx q[2];
rz(-3.0068676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45891201) q[1];
sx q[1];
rz(-2.4574728) q[1];
sx q[1];
rz(-2.4980809) q[1];
x q[2];
rz(1.7836501) q[3];
sx q[3];
rz(-2.2567333) q[3];
sx q[3];
rz(-2.6712772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4400441) q[2];
sx q[2];
rz(-0.35228071) q[2];
sx q[2];
rz(0.80553833) q[2];
rz(0.37197477) q[3];
sx q[3];
rz(-1.493908) q[3];
sx q[3];
rz(-0.39723799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-3.1086248) q[0];
sx q[0];
rz(-1.2297577) q[0];
sx q[0];
rz(0.97880542) q[0];
rz(2.4188614) q[1];
sx q[1];
rz(-1.1284072) q[1];
sx q[1];
rz(-2.5236948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2513951) q[0];
sx q[0];
rz(-0.82444901) q[0];
sx q[0];
rz(-1.192481) q[0];
rz(-1.0523791) q[2];
sx q[2];
rz(-2.1627277) q[2];
sx q[2];
rz(0.059378864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1714461) q[1];
sx q[1];
rz(-1.3916124) q[1];
sx q[1];
rz(-2.9703364) q[1];
rz(-1.4862107) q[3];
sx q[3];
rz(-1.654155) q[3];
sx q[3];
rz(2.3918834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2021947) q[2];
sx q[2];
rz(-1.3556182) q[2];
sx q[2];
rz(3.1414269) q[2];
rz(-2.557142) q[3];
sx q[3];
rz(-2.1355459) q[3];
sx q[3];
rz(-2.5463026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7547739) q[0];
sx q[0];
rz(-1.5703572) q[0];
sx q[0];
rz(1.5686709) q[0];
rz(1.8008925) q[1];
sx q[1];
rz(-1.1005713) q[1];
sx q[1];
rz(1.5250199) q[1];
rz(2.2095815) q[2];
sx q[2];
rz(-1.1975653) q[2];
sx q[2];
rz(-1.4690659) q[2];
rz(0.67047337) q[3];
sx q[3];
rz(-2.5994876) q[3];
sx q[3];
rz(2.2710298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
