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
rz(-0.42521617) q[0];
sx q[0];
rz(-1.8073616) q[0];
sx q[0];
rz(-0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(2.9549197) q[1];
sx q[1];
rz(8.5001707) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036669) q[0];
sx q[0];
rz(-0.80090085) q[0];
sx q[0];
rz(-0.4842224) q[0];
x q[1];
rz(1.0546646) q[2];
sx q[2];
rz(-1.7164111) q[2];
sx q[2];
rz(-0.21805412) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9052995) q[1];
sx q[1];
rz(-0.84262203) q[1];
sx q[1];
rz(3.0889126) q[1];
rz(-3.0136631) q[3];
sx q[3];
rz(-1.5797714) q[3];
sx q[3];
rz(2.8913446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0164464) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(-2.9347349) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(2.8894506) q[0];
rz(0.02154669) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-0.85539877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148672) q[0];
sx q[0];
rz(-1.5692737) q[0];
sx q[0];
rz(1.5978088) q[0];
rz(-0.55412519) q[2];
sx q[2];
rz(-1.9337855) q[2];
sx q[2];
rz(-2.3100694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.336074) q[1];
sx q[1];
rz(-1.0401498) q[1];
sx q[1];
rz(-2.8855399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68287373) q[3];
sx q[3];
rz(-1.3078017) q[3];
sx q[3];
rz(-2.9337286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97027332) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(2.8805736) q[2];
rz(3.1018992) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(-2.4470827) q[0];
rz(1.1272686) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(-2.1515501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5607669) q[0];
sx q[0];
rz(-0.35835727) q[0];
sx q[0];
rz(0.23507765) q[0];
rz(-1.2743397) q[2];
sx q[2];
rz(-2.5387089) q[2];
sx q[2];
rz(0.68494132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5824008) q[1];
sx q[1];
rz(-1.9291996) q[1];
sx q[1];
rz(2.6646975) q[1];
rz(2.2740836) q[3];
sx q[3];
rz(-1.1425619) q[3];
sx q[3];
rz(-0.060286418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29828829) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.8753768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0855584) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(2.339824) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(2.3559949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39488897) q[0];
sx q[0];
rz(-2.3869042) q[0];
sx q[0];
rz(-2.7695023) q[0];
x q[1];
rz(-1.7217721) q[2];
sx q[2];
rz(-2.2822126) q[2];
sx q[2];
rz(0.47455088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6735437) q[1];
sx q[1];
rz(-2.1077413) q[1];
sx q[1];
rz(-1.931744) q[1];
rz(-0.81964747) q[3];
sx q[3];
rz(-1.420212) q[3];
sx q[3];
rz(2.8898847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12513146) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(-1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9652902) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(-1.0705795) q[0];
rz(-2.3640682) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(1.0106962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010653822) q[0];
sx q[0];
rz(-1.21216) q[0];
sx q[0];
rz(0.50827311) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0762482) q[2];
sx q[2];
rz(-1.7453219) q[2];
sx q[2];
rz(-0.068152817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3123828) q[1];
sx q[1];
rz(-2.0264033) q[1];
sx q[1];
rz(-2.4036951) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56212496) q[3];
sx q[3];
rz(-1.486612) q[3];
sx q[3];
rz(-2.4474194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(-1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(-1.4122081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509336) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(2.3405128) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(-0.4020234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659365) q[0];
sx q[0];
rz(-2.2992327) q[0];
sx q[0];
rz(-0.45439646) q[0];
x q[1];
rz(1.5779675) q[2];
sx q[2];
rz(-2.8507887) q[2];
sx q[2];
rz(2.4343642) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1014598) q[1];
sx q[1];
rz(-1.8454843) q[1];
sx q[1];
rz(1.0708315) q[1];
x q[2];
rz(0.16120637) q[3];
sx q[3];
rz(-2.4771002) q[3];
sx q[3];
rz(2.4624069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5037527) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(-0.77581882) q[2];
rz(1.4555629) q[3];
sx q[3];
rz(-1.173923) q[3];
sx q[3];
rz(1.0078526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0966454) q[0];
sx q[0];
rz(-2.2438887) q[0];
sx q[0];
rz(2.4626379) q[0];
rz(-1.2848162) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(1.3791893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5011936) q[0];
sx q[0];
rz(-1.3974579) q[0];
sx q[0];
rz(2.3480183) q[0];
rz(-pi) q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.2922568) q[2];
sx q[2];
rz(-2.4065774) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0028982) q[1];
sx q[1];
rz(-2.0636673) q[1];
sx q[1];
rz(2.0763796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3016765) q[3];
sx q[3];
rz(-2.9342071) q[3];
sx q[3];
rz(0.086750448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3369559) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(1.2389368) q[2];
rz(-1.3321446) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(3.0031257) q[0];
rz(-2.9578517) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(0.26652452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375232) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(-2.0384602) q[0];
rz(-pi) q[1];
rz(2.1378072) q[2];
sx q[2];
rz(-1.8006174) q[2];
sx q[2];
rz(2.9858495) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6147242) q[1];
sx q[1];
rz(-0.42335948) q[1];
sx q[1];
rz(-1.7080659) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1317461) q[3];
sx q[3];
rz(-0.72692633) q[3];
sx q[3];
rz(1.3971412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0900241) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-0.23452342) q[2];
rz(-1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33673564) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(-2.3418703) q[0];
rz(-2.5526478) q[1];
sx q[1];
rz(-0.84732333) q[1];
sx q[1];
rz(-0.2074997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5867347) q[0];
sx q[0];
rz(-1.4401146) q[0];
sx q[0];
rz(-2.7452637) q[0];
rz(-0.3090119) q[2];
sx q[2];
rz(-2.2591619) q[2];
sx q[2];
rz(0.38028827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7991528) q[1];
sx q[1];
rz(-2.229433) q[1];
sx q[1];
rz(-2.6299146) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1010246) q[3];
sx q[3];
rz(-1.54514) q[3];
sx q[3];
rz(-0.86650833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80692446) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(-0.36726382) q[2];
rz(-1.4372829) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87840286) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(-0.81018418) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(2.5883163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5636175) q[0];
sx q[0];
rz(-1.6262691) q[0];
sx q[0];
rz(-0.93240057) q[0];
rz(-pi) q[1];
x q[1];
rz(2.192334) q[2];
sx q[2];
rz(-2.261603) q[2];
sx q[2];
rz(2.031428) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77299632) q[1];
sx q[1];
rz(-2.8008411) q[1];
sx q[1];
rz(0.98253886) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7891879) q[3];
sx q[3];
rz(-0.81768546) q[3];
sx q[3];
rz(-1.6116097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0193923) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(-0.23966399) q[2];
rz(-1.8137118) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(0.46752587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079530579) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(-0.65555864) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(2.4782933) q[2];
sx q[2];
rz(-2.6313783) q[2];
sx q[2];
rz(0.22102697) q[2];
rz(-0.58047337) q[3];
sx q[3];
rz(-1.7844641) q[3];
sx q[3];
rz(-2.6337558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
