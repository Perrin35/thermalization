OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(-2.2401016) q[0];
rz(1.3540406) q[1];
sx q[1];
rz(-1.5259589) q[1];
sx q[1];
rz(1.807133) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5164099) q[0];
sx q[0];
rz(-1.3526655) q[0];
sx q[0];
rz(-1.5152009) q[0];
rz(-pi) q[1];
rz(1.5152211) q[2];
sx q[2];
rz(-1.7211421) q[2];
sx q[2];
rz(-1.0483345) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5049274) q[1];
sx q[1];
rz(-1.0475155) q[1];
sx q[1];
rz(0.06097554) q[1];
rz(-1.5235352) q[3];
sx q[3];
rz(-1.3287373) q[3];
sx q[3];
rz(-0.74064613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91288519) q[2];
sx q[2];
rz(-0.097147377) q[2];
sx q[2];
rz(2.482282) q[2];
rz(-2.3653476) q[3];
sx q[3];
rz(-3.1228437) q[3];
sx q[3];
rz(-2.4565878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2289675) q[0];
sx q[0];
rz(-1.9321059) q[0];
sx q[0];
rz(-1.9483161) q[0];
rz(0.044366447) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(-0.22656974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9539174) q[0];
sx q[0];
rz(-1.5724475) q[0];
sx q[0];
rz(-1.5673914) q[0];
x q[1];
rz(1.5537276) q[2];
sx q[2];
rz(-1.9498132) q[2];
sx q[2];
rz(1.5660777) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0262278) q[1];
sx q[1];
rz(-1.4997002) q[1];
sx q[1];
rz(2.709862) q[1];
x q[2];
rz(-1.3188521) q[3];
sx q[3];
rz(-0.74352194) q[3];
sx q[3];
rz(1.3321259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6988397) q[2];
sx q[2];
rz(-1.5696462) q[2];
sx q[2];
rz(1.5296096) q[2];
rz(-2.2085341) q[3];
sx q[3];
rz(-1.6589087) q[3];
sx q[3];
rz(0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642692) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(-2.9413057) q[0];
rz(3.1411723) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(-0.013484152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4448252) q[0];
sx q[0];
rz(-0.42967859) q[0];
sx q[0];
rz(1.4151156) q[0];
rz(-pi) q[1];
rz(0.066670074) q[2];
sx q[2];
rz(-1.5278897) q[2];
sx q[2];
rz(1.5868452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2993907) q[1];
sx q[1];
rz(-2.9889571) q[1];
sx q[1];
rz(1.0959036) q[1];
x q[2];
rz(-1.5192658) q[3];
sx q[3];
rz(-1.7331707) q[3];
sx q[3];
rz(-2.3414617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1991594) q[2];
sx q[2];
rz(-1.5853256) q[2];
sx q[2];
rz(1.5419434) q[2];
rz(2.13983) q[3];
sx q[3];
rz(-2.9250513) q[3];
sx q[3];
rz(-0.17222968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7672985) q[0];
sx q[0];
rz(-0.16574398) q[0];
sx q[0];
rz(-0.34749183) q[0];
rz(2.5047498) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(-1.1905131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74032846) q[0];
sx q[0];
rz(-1.51121) q[0];
sx q[0];
rz(-1.4387095) q[0];
rz(0.0077771386) q[2];
sx q[2];
rz(-1.5019122) q[2];
sx q[2];
rz(-1.5780256) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5305883) q[1];
sx q[1];
rz(-1.3827033) q[1];
sx q[1];
rz(0.76716073) q[1];
rz(-pi) q[2];
rz(-2.7772831) q[3];
sx q[3];
rz(-1.7765883) q[3];
sx q[3];
rz(2.8575773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5570598) q[2];
sx q[2];
rz(-3.1019042) q[2];
sx q[2];
rz(1.7812799) q[2];
rz(-1.6654061) q[3];
sx q[3];
rz(-1.5740266) q[3];
sx q[3];
rz(2.7358957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0880459) q[0];
sx q[0];
rz(-0.73943728) q[0];
sx q[0];
rz(-0.0060225688) q[0];
rz(1.7209523) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(-0.086070148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0417418) q[0];
sx q[0];
rz(-1.6779416) q[0];
sx q[0];
rz(-1.5539377) q[0];
rz(3.086523) q[2];
sx q[2];
rz(-1.5447642) q[2];
sx q[2];
rz(0.30773417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2246612) q[1];
sx q[1];
rz(-1.5337394) q[1];
sx q[1];
rz(0.091450973) q[1];
x q[2];
rz(2.213201) q[3];
sx q[3];
rz(-0.07559055) q[3];
sx q[3];
rz(-2.9380611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1083131) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(1.3839728) q[2];
rz(0.39517394) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(1.1672195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5351717) q[0];
sx q[0];
rz(-0.21927729) q[0];
sx q[0];
rz(-1.0004591) q[0];
rz(-0.74042997) q[1];
sx q[1];
rz(-0.43559566) q[1];
sx q[1];
rz(2.7620517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6839121) q[0];
sx q[0];
rz(-0.098794071) q[0];
sx q[0];
rz(0.26215078) q[0];
rz(-0.47490317) q[2];
sx q[2];
rz(-1.3273718) q[2];
sx q[2];
rz(-2.7096675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15495519) q[1];
sx q[1];
rz(-1.0670245) q[1];
sx q[1];
rz(0.030435199) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0256512) q[3];
sx q[3];
rz(-1.1244457) q[3];
sx q[3];
rz(1.4614814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1257816) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(-0.077032653) q[2];
rz(-3.1214664) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(-0.7974112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064442) q[0];
sx q[0];
rz(-0.012520944) q[0];
sx q[0];
rz(-1.7092108) q[0];
rz(-0.2969946) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(-0.061554734) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650353) q[0];
sx q[0];
rz(-1.1662959) q[0];
sx q[0];
rz(0.96488953) q[0];
rz(-0.02587895) q[2];
sx q[2];
rz(-2.9270083) q[2];
sx q[2];
rz(1.1550922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24352001) q[1];
sx q[1];
rz(-2.8419333) q[1];
sx q[1];
rz(0.23345848) q[1];
rz(1.6084303) q[3];
sx q[3];
rz(-2.0052344) q[3];
sx q[3];
rz(-1.5805336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84294549) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(1.418815) q[2];
rz(-1.336054) q[3];
sx q[3];
rz(-3.1137443) q[3];
sx q[3];
rz(1.7347887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.0272738) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(1.5007716) q[0];
rz(-2.0188792) q[1];
sx q[1];
rz(-2.8148459) q[1];
sx q[1];
rz(1.8480802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1294848) q[0];
sx q[0];
rz(-0.92381682) q[0];
sx q[0];
rz(-2.7478474) q[0];
x q[1];
rz(0.69399909) q[2];
sx q[2];
rz(-1.8376164) q[2];
sx q[2];
rz(2.9129183) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39513341) q[1];
sx q[1];
rz(-0.32276216) q[1];
sx q[1];
rz(-2.5844896) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84728665) q[3];
sx q[3];
rz(-2.1186817) q[3];
sx q[3];
rz(-0.63036608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5844172) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(1.3593675) q[2];
rz(-0.44698295) q[3];
sx q[3];
rz(-3.0990661) q[3];
sx q[3];
rz(0.16714787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3249224) q[0];
sx q[0];
rz(-1.8055547) q[0];
sx q[0];
rz(2.0409806) q[0];
rz(-1.0758411) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(-0.67695391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3269126) q[0];
sx q[0];
rz(-1.956358) q[0];
sx q[0];
rz(-1.7837693) q[0];
rz(-pi) q[1];
rz(0.81440429) q[2];
sx q[2];
rz(-0.16781092) q[2];
sx q[2];
rz(-1.7240768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1558324) q[1];
sx q[1];
rz(-1.5703771) q[1];
sx q[1];
rz(9.930519e-05) q[1];
rz(-pi) q[2];
rz(0.30218924) q[3];
sx q[3];
rz(-2.6679278) q[3];
sx q[3];
rz(-0.74407265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0000275) q[2];
sx q[2];
rz(-3.1391149) q[2];
sx q[2];
rz(-2.4139717) q[2];
rz(-1.8992807) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(1.1718933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2321371) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(-2.452028) q[0];
rz(-1.4762956) q[1];
sx q[1];
rz(-2.8910525) q[1];
sx q[1];
rz(0.15588674) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4613381) q[0];
sx q[0];
rz(-0.91947407) q[0];
sx q[0];
rz(-1.7375577) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8901398) q[2];
sx q[2];
rz(-1.4551468) q[2];
sx q[2];
rz(2.419099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5248651) q[1];
sx q[1];
rz(-1.5666125) q[1];
sx q[1];
rz(1.5692488) q[1];
rz(1.0267657) q[3];
sx q[3];
rz(-2.9599905) q[3];
sx q[3];
rz(-1.3917519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9547687) q[2];
sx q[2];
rz(-0.12207741) q[2];
sx q[2];
rz(0.99511498) q[2];
rz(-0.22594813) q[3];
sx q[3];
rz(-0.048361691) q[3];
sx q[3];
rz(-2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.7786998) q[0];
sx q[0];
rz(-2.1669372) q[0];
sx q[0];
rz(-1.7397276) q[0];
rz(-1.7097991) q[1];
sx q[1];
rz(-1.2964389) q[1];
sx q[1];
rz(-2.5242205) q[1];
rz(-2.974343) q[2];
sx q[2];
rz(-2.4495301) q[2];
sx q[2];
rz(-2.4514103) q[2];
rz(1.2479242) q[3];
sx q[3];
rz(-0.40798305) q[3];
sx q[3];
rz(1.8508607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
