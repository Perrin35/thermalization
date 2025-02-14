OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(-2.5929911) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34615883) q[0];
sx q[0];
rz(-1.2761269) q[0];
sx q[0];
rz(-0.3064038) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59738417) q[2];
sx q[2];
rz(-1.5533012) q[2];
sx q[2];
rz(1.04523) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59034172) q[1];
sx q[1];
rz(-1.1635374) q[1];
sx q[1];
rz(-2.1339416) q[1];
rz(-1.2554507) q[3];
sx q[3];
rz(-1.3419125) q[3];
sx q[3];
rz(-1.6381519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7301664) q[2];
sx q[2];
rz(-2.7226518) q[2];
sx q[2];
rz(-2.1298998) q[2];
rz(2.2569979) q[3];
sx q[3];
rz(-1.1583637) q[3];
sx q[3];
rz(-1.2456892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0737792) q[0];
sx q[0];
rz(-2.1429006) q[0];
sx q[0];
rz(-2.3538537) q[0];
rz(0.9785606) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(1.8914793) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7536016) q[0];
sx q[0];
rz(-1.4630254) q[0];
sx q[0];
rz(2.0600256) q[0];
rz(1.1680295) q[2];
sx q[2];
rz(-1.313123) q[2];
sx q[2];
rz(-2.7855025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3652894) q[1];
sx q[1];
rz(-2.3765916) q[1];
sx q[1];
rz(1.5359775) q[1];
x q[2];
rz(2.9216209) q[3];
sx q[3];
rz(-1.7860054) q[3];
sx q[3];
rz(-0.53088354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1193739) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(-0.12164965) q[2];
rz(-2.4308448) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(2.5392883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0372593) q[0];
sx q[0];
rz(-0.72838825) q[0];
sx q[0];
rz(-0.65993586) q[0];
rz(-0.51689369) q[1];
sx q[1];
rz(-2.4362502) q[1];
sx q[1];
rz(-2.0491811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67319523) q[0];
sx q[0];
rz(-1.930113) q[0];
sx q[0];
rz(-0.27984377) q[0];
x q[1];
rz(2.8404929) q[2];
sx q[2];
rz(-0.89576713) q[2];
sx q[2];
rz(1.1464034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1221083) q[1];
sx q[1];
rz(-1.0398163) q[1];
sx q[1];
rz(2.2506258) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64869439) q[3];
sx q[3];
rz(-1.7770542) q[3];
sx q[3];
rz(-2.1163396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2163781) q[2];
sx q[2];
rz(-0.34319147) q[2];
sx q[2];
rz(-1.421831) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.657594) q[3];
sx q[3];
rz(1.4181731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5928818) q[0];
sx q[0];
rz(-0.084267862) q[0];
sx q[0];
rz(-1.8362057) q[0];
rz(-0.0072172324) q[1];
sx q[1];
rz(-0.22166285) q[1];
sx q[1];
rz(0.73297393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77196593) q[0];
sx q[0];
rz(-1.4384171) q[0];
sx q[0];
rz(0.11970971) q[0];
rz(-pi) q[1];
rz(-0.29455955) q[2];
sx q[2];
rz(-1.6798875) q[2];
sx q[2];
rz(-0.2738758) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8165255) q[1];
sx q[1];
rz(-2.1754334) q[1];
sx q[1];
rz(0.52765347) q[1];
rz(-1.6812612) q[3];
sx q[3];
rz(-0.85064954) q[3];
sx q[3];
rz(-3.0283749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8492154) q[2];
sx q[2];
rz(-1.7376309) q[2];
sx q[2];
rz(2.2461069) q[2];
rz(0.53564566) q[3];
sx q[3];
rz(-1.5786542) q[3];
sx q[3];
rz(0.061804684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8931005) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(-0.50791159) q[0];
rz(1.4503362) q[1];
sx q[1];
rz(-1.0117057) q[1];
sx q[1];
rz(-1.3571665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0471418) q[0];
sx q[0];
rz(-3.019382) q[0];
sx q[0];
rz(-0.41441865) q[0];
rz(2.6772887) q[2];
sx q[2];
rz(-1.1912734) q[2];
sx q[2];
rz(0.10019856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30631599) q[1];
sx q[1];
rz(-0.40566722) q[1];
sx q[1];
rz(1.8464645) q[1];
rz(-pi) q[2];
x q[2];
rz(2.761854) q[3];
sx q[3];
rz(-2.1354699) q[3];
sx q[3];
rz(-2.5329451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7148529) q[2];
sx q[2];
rz(-1.667495) q[2];
sx q[2];
rz(2.0392141) q[2];
rz(-1.1357931) q[3];
sx q[3];
rz(-1.9250684) q[3];
sx q[3];
rz(-0.77146012) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256846) q[0];
sx q[0];
rz(-0.9698292) q[0];
sx q[0];
rz(2.2127175) q[0];
rz(2.7595787) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(1.4487723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083915868) q[0];
sx q[0];
rz(-1.8734212) q[0];
sx q[0];
rz(0.50047366) q[0];
x q[1];
rz(-1.5097627) q[2];
sx q[2];
rz(-2.054913) q[2];
sx q[2];
rz(-2.8060437) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12940059) q[1];
sx q[1];
rz(-1.4621648) q[1];
sx q[1];
rz(1.1715749) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17254707) q[3];
sx q[3];
rz(-1.2013519) q[3];
sx q[3];
rz(1.9011232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53943071) q[2];
sx q[2];
rz(-2.7950037) q[2];
sx q[2];
rz(0.76751417) q[2];
rz(1.8481988) q[3];
sx q[3];
rz(-2.4496205) q[3];
sx q[3];
rz(-0.20492157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(2.8906004) q[0];
rz(-2.9639066) q[1];
sx q[1];
rz(-1.4439986) q[1];
sx q[1];
rz(-0.65690717) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11372317) q[0];
sx q[0];
rz(-1.8060469) q[0];
sx q[0];
rz(0.21224888) q[0];
rz(-1.7288293) q[2];
sx q[2];
rz(-0.92301805) q[2];
sx q[2];
rz(0.72731804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4106573) q[1];
sx q[1];
rz(-1.9563795) q[1];
sx q[1];
rz(-1.247184) q[1];
rz(-0.55258379) q[3];
sx q[3];
rz(-2.4173418) q[3];
sx q[3];
rz(0.16952215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3333007) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(0.91326886) q[2];
rz(0.67780668) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(1.0031797) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.102757) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(-2.5586149) q[0];
rz(-1.673117) q[1];
sx q[1];
rz(-0.99248326) q[1];
sx q[1];
rz(2.800422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2058294) q[0];
sx q[0];
rz(-1.7700164) q[0];
sx q[0];
rz(-0.043313839) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24218817) q[2];
sx q[2];
rz(-1.2953087) q[2];
sx q[2];
rz(-1.3812068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.063733) q[1];
sx q[1];
rz(-1.2021515) q[1];
sx q[1];
rz(0.39409448) q[1];
x q[2];
rz(-1.0811133) q[3];
sx q[3];
rz(-0.73990179) q[3];
sx q[3];
rz(-2.5748753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94656452) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(0.80671802) q[2];
rz(-2.6840456) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(0.88361067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44613999) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(0.17644185) q[0];
rz(-2.8314619) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(-0.93607059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061358) q[0];
sx q[0];
rz(-1.8671037) q[0];
sx q[0];
rz(-1.4151876) q[0];
rz(-pi) q[1];
rz(-2.5199982) q[2];
sx q[2];
rz(-2.4450069) q[2];
sx q[2];
rz(-2.7736349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99451274) q[1];
sx q[1];
rz(-1.4777061) q[1];
sx q[1];
rz(-1.2017815) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8027169) q[3];
sx q[3];
rz(-0.48503625) q[3];
sx q[3];
rz(-2.1887961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.04756847) q[2];
sx q[2];
rz(-1.3022283) q[2];
sx q[2];
rz(-3.1330718) q[2];
rz(2.408037) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(3.0146397) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200639) q[0];
sx q[0];
rz(-0.41291741) q[0];
sx q[0];
rz(2.371696) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-2.3513992) q[1];
sx q[1];
rz(1.92164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1320187) q[0];
sx q[0];
rz(-1.2694799) q[0];
sx q[0];
rz(-1.8150041) q[0];
rz(-1.8656857) q[2];
sx q[2];
rz(-0.98433009) q[2];
sx q[2];
rz(0.49582729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97970672) q[1];
sx q[1];
rz(-2.5570013) q[1];
sx q[1];
rz(-2.4120055) q[1];
rz(1.4195326) q[3];
sx q[3];
rz(-1.6409887) q[3];
sx q[3];
rz(2.9541525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28025815) q[2];
sx q[2];
rz(-2.4836149) q[2];
sx q[2];
rz(0.78557837) q[2];
rz(2.806459) q[3];
sx q[3];
rz(-2.1527055) q[3];
sx q[3];
rz(-2.3194763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.32421865) q[0];
sx q[0];
rz(-2.2786409) q[0];
sx q[0];
rz(1.8550158) q[0];
rz(0.68589504) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(2.1519312) q[2];
sx q[2];
rz(-1.5662532) q[2];
sx q[2];
rz(-3.1305542) q[2];
rz(-1.2780581) q[3];
sx q[3];
rz(-1.6061898) q[3];
sx q[3];
rz(1.3502179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
