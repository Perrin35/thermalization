OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.769387) q[0];
sx q[0];
rz(-1.1507478) q[0];
sx q[0];
rz(-0.37680349) q[0];
rz(-pi) q[1];
rz(2.8777962) q[2];
sx q[2];
rz(-0.88636878) q[2];
sx q[2];
rz(-0.85927187) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3914324) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(-1.0990012) q[1];
rz(-pi) q[2];
rz(-1.5815758) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(-2.74995) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(-0.02123775) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613974) q[0];
sx q[0];
rz(-1.178) q[0];
sx q[0];
rz(-2.0307721) q[0];
rz(-0.86512489) q[2];
sx q[2];
rz(-1.428831) q[2];
sx q[2];
rz(2.1324468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.9351164) q[1];
sx q[1];
rz(-1.4644535) q[1];
rz(-pi) q[2];
rz(0.89800091) q[3];
sx q[3];
rz(-1.7963396) q[3];
sx q[3];
rz(-2.0190092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(0.35456625) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.59042674) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(-0.61808008) q[0];
rz(3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(1.9504257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-2.7911721) q[0];
sx q[0];
rz(-2.0273897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8580655) q[2];
sx q[2];
rz(-1.8054188) q[2];
sx q[2];
rz(-1.7262176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1515854) q[1];
sx q[1];
rz(-1.7813492) q[1];
sx q[1];
rz(-1.7819808) q[1];
rz(-pi) q[2];
rz(1.1690087) q[3];
sx q[3];
rz(-1.0649293) q[3];
sx q[3];
rz(1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(-1.7640242) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-2.2185982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3337293) q[0];
sx q[0];
rz(-0.45818746) q[0];
sx q[0];
rz(-2.3966167) q[0];
rz(-2.8795305) q[2];
sx q[2];
rz(-3.0960961) q[2];
sx q[2];
rz(-0.87693518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2642306) q[1];
sx q[1];
rz(-1.1738395) q[1];
sx q[1];
rz(-3.0287659) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5444078) q[3];
sx q[3];
rz(-2.4370586) q[3];
sx q[3];
rz(1.964331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(0.50951177) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.3607508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960984) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(1.7646837) q[0];
rz(-pi) q[1];
rz(2.4945716) q[2];
sx q[2];
rz(-0.87304742) q[2];
sx q[2];
rz(-0.24017142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9779224) q[1];
sx q[1];
rz(-2.0241892) q[1];
sx q[1];
rz(-1.8015253) q[1];
rz(0.97134437) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(-1.5980699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3551066) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(3.0267267) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-0.27522603) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(0.37822322) q[0];
rz(-0.49595828) q[2];
sx q[2];
rz(-0.99018103) q[2];
sx q[2];
rz(2.1883287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59906193) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(-2.6229726) q[1];
rz(-1.4554548) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(1.8231525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39367166) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(-2.6829524) q[2];
rz(-0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28073072) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(2.4218959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170383) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(-0.65529234) q[0];
rz(-pi) q[1];
rz(-2.3709488) q[2];
sx q[2];
rz(-2.1241803) q[2];
sx q[2];
rz(-1.1469974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91832671) q[1];
sx q[1];
rz(-1.0263138) q[1];
sx q[1];
rz(1.5189927) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3004488) q[3];
sx q[3];
rz(-0.90793228) q[3];
sx q[3];
rz(-2.8519252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.4950745) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197501) q[0];
sx q[0];
rz(-1.2312206) q[0];
sx q[0];
rz(1.7975438) q[0];
rz(-0.69627701) q[2];
sx q[2];
rz(-1.7211203) q[2];
sx q[2];
rz(-2.6339649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5852768) q[1];
sx q[1];
rz(-1.3116515) q[1];
sx q[1];
rz(-1.4633333) q[1];
x q[2];
rz(-0.90663119) q[3];
sx q[3];
rz(-2.2615221) q[3];
sx q[3];
rz(-0.074113473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(3.126826) q[2];
rz(-1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(-1.7285255) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.7810129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1454865) q[0];
sx q[0];
rz(-1.6084451) q[0];
sx q[0];
rz(-2.1426175) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0948823) q[2];
sx q[2];
rz(-0.94418664) q[2];
sx q[2];
rz(1.4057297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66879241) q[1];
sx q[1];
rz(-1.2392514) q[1];
sx q[1];
rz(-1.214295) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9901048) q[3];
sx q[3];
rz(-1.4290223) q[3];
sx q[3];
rz(-1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(2.2488135) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-0.75004309) q[3];
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
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(-3.0950586) q[0];
rz(-0.90905601) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68144804) q[0];
sx q[0];
rz(-1.8143192) q[0];
sx q[0];
rz(0.081865099) q[0];
x q[1];
rz(0.30883046) q[2];
sx q[2];
rz(-1.575982) q[2];
sx q[2];
rz(-1.7538824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9626999) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(0.015951338) q[1];
x q[2];
rz(0.7474483) q[3];
sx q[3];
rz(-2.2564853) q[3];
sx q[3];
rz(-1.3626584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-0.026467888) q[2];
rz(1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(2.9251255) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-2.7770476) q[2];
sx q[2];
rz(-0.74082965) q[2];
sx q[2];
rz(-0.34168591) q[2];
rz(-3.1090267) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
