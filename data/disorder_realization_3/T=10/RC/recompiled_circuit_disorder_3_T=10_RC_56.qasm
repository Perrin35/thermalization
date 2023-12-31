OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(0.41419849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.480455) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(2.7829091) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0010927) q[2];
sx q[2];
rz(-1.4909407) q[2];
sx q[2];
rz(-2.9542838) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1297076) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(0.070406291) q[1];
rz(-pi) q[2];
rz(1.6730509) q[3];
sx q[3];
rz(-1.8342606) q[3];
sx q[3];
rz(1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-0.53952113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9233421) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(1.1211066) q[0];
x q[1];
rz(-2.4194854) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(-2.9177641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19903781) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(1.1953137) q[1];
rz(-pi) q[2];
rz(0.94445618) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-0.32523793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494333) q[0];
sx q[0];
rz(-1.306635) q[0];
sx q[0];
rz(0.63412068) q[0];
rz(-pi) q[1];
rz(-0.82661144) q[2];
sx q[2];
rz(-2.5275143) q[2];
sx q[2];
rz(2.5342864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33830723) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(2.0558946) q[1];
rz(-pi) q[2];
rz(2.9647397) q[3];
sx q[3];
rz(-1.6079418) q[3];
sx q[3];
rz(2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(1.3556708) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.2333966) q[3];
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
rz(1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-0.25948778) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-2.4096699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5607802) q[0];
sx q[0];
rz(-0.73045759) q[0];
sx q[0];
rz(0.74303445) q[0];
rz(-pi) q[1];
rz(-0.71728431) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(0.35536534) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(0.64756067) q[1];
x q[2];
rz(2.6287574) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(-1.4506838) q[0];
rz(-pi) q[1];
rz(-1.7902137) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(2.9079633) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.080575374) q[1];
sx q[1];
rz(-1.6643545) q[1];
sx q[1];
rz(0.94228014) q[1];
rz(2.5161414) q[3];
sx q[3];
rz(-0.52895412) q[3];
sx q[3];
rz(1.7367425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-3.0715122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698343) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(-1.397875) q[0];
rz(-1.0191392) q[2];
sx q[2];
rz(-0.86691228) q[2];
sx q[2];
rz(-0.64955074) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97092123) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(-1.7621653) q[1];
rz(1.1280941) q[3];
sx q[3];
rz(-3.0590995) q[3];
sx q[3];
rz(-3.1174297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(-0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746349) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(0.5704244) q[0];
rz(2.7107312) q[2];
sx q[2];
rz(-2.2164946) q[2];
sx q[2];
rz(-0.25260392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4850033) q[1];
sx q[1];
rz(-0.90066972) q[1];
sx q[1];
rz(-0.90799241) q[1];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(0.34208959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4456711) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2973328) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(2.8914333) q[0];
rz(1.0566696) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(2.4664997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(-2.3676141) q[1];
rz(3.0495166) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(-0.74842194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(0.41697821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.5331393) q[0];
sx q[0];
rz(-2.2249939) q[0];
rz(-0.68848227) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(0.85306963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6777842) q[1];
sx q[1];
rz(-2.6208372) q[1];
sx q[1];
rz(1.2893454) q[1];
rz(-pi) q[2];
rz(3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.5243994) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(-1.1907207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0199213) q[2];
sx q[2];
rz(-1.7977062) q[2];
sx q[2];
rz(-0.34044468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(0.91099693) q[1];
x q[2];
rz(-2.2262276) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.1673499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(2.1255169) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-2.6323742) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
