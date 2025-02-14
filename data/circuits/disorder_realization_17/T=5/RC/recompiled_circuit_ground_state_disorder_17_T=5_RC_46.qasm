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
rz(-0.53417438) q[0];
sx q[0];
rz(-2.0913273) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(-2.5929911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596244) q[0];
sx q[0];
rz(-0.42185703) q[0];
sx q[0];
rz(0.78849383) q[0];
x q[1];
rz(-1.5919551) q[2];
sx q[2];
rz(-2.1680764) q[2];
sx q[2];
rz(-0.51366546) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5512509) q[1];
sx q[1];
rz(-1.9780553) q[1];
sx q[1];
rz(2.1339416) q[1];
rz(-pi) q[2];
rz(-2.9012783) q[3];
sx q[3];
rz(-1.263947) q[3];
sx q[3];
rz(0.14123973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7301664) q[2];
sx q[2];
rz(-2.7226518) q[2];
sx q[2];
rz(-1.0116928) q[2];
rz(0.88459477) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(1.8959034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678134) q[0];
sx q[0];
rz(-0.99869204) q[0];
sx q[0];
rz(-2.3538537) q[0];
rz(-0.9785606) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(1.2501134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.387991) q[0];
sx q[0];
rz(-1.6785673) q[0];
sx q[0];
rz(2.0600256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1627126) q[2];
sx q[2];
rz(-2.667281) q[2];
sx q[2];
rz(1.3878551) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82455222) q[1];
sx q[1];
rz(-2.3352156) q[1];
sx q[1];
rz(-3.1081852) q[1];
rz(-pi) q[2];
rz(-1.7911508) q[3];
sx q[3];
rz(-1.3559794) q[3];
sx q[3];
rz(-2.1493916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1193739) q[2];
sx q[2];
rz(-1.677607) q[2];
sx q[2];
rz(-0.12164965) q[2];
rz(0.71074784) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(-0.60230437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1043333) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(2.4816568) q[0];
rz(0.51689369) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(1.0924115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79688841) q[0];
sx q[0];
rz(-1.3092586) q[0];
sx q[0];
rz(1.198223) q[0];
x q[1];
rz(-2.2683892) q[2];
sx q[2];
rz(-1.8044458) q[2];
sx q[2];
rz(2.9088809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9398492) q[1];
sx q[1];
rz(-2.1437316) q[1];
sx q[1];
rz(0.64676379) q[1];
rz(2.4928983) q[3];
sx q[3];
rz(-1.3645384) q[3];
sx q[3];
rz(-1.0252531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92521459) q[2];
sx q[2];
rz(-2.7984012) q[2];
sx q[2];
rz(1.421831) q[2];
rz(2.6575994) q[3];
sx q[3];
rz(-1.4839987) q[3];
sx q[3];
rz(1.4181731) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5928818) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(1.3053869) q[0];
rz(-3.1343754) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(0.73297393) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6306046) q[0];
sx q[0];
rz(-0.17824358) q[0];
sx q[0];
rz(-2.3018738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4568366) q[2];
sx q[2];
rz(-1.8635529) q[2];
sx q[2];
rz(1.2639015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8165255) q[1];
sx q[1];
rz(-2.1754334) q[1];
sx q[1];
rz(2.6139392) q[1];
x q[2];
rz(0.72317883) q[3];
sx q[3];
rz(-1.4878325) q[3];
sx q[3];
rz(-1.3845598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8492154) q[2];
sx q[2];
rz(-1.4039618) q[2];
sx q[2];
rz(-0.89548573) q[2];
rz(0.53564566) q[3];
sx q[3];
rz(-1.5786542) q[3];
sx q[3];
rz(-3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8931005) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(-2.6336811) q[0];
rz(-1.4503362) q[1];
sx q[1];
rz(-2.129887) q[1];
sx q[1];
rz(-1.3571665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535808) q[0];
sx q[0];
rz(-1.5216899) q[0];
sx q[0];
rz(-0.11195575) q[0];
rz(-pi) q[1];
rz(2.4139348) q[2];
sx q[2];
rz(-2.5508159) q[2];
sx q[2];
rz(-1.0340921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1313981) q[1];
sx q[1];
rz(-1.6784188) q[1];
sx q[1];
rz(-1.1788998) q[1];
x q[2];
rz(-2.761854) q[3];
sx q[3];
rz(-2.1354699) q[3];
sx q[3];
rz(2.5329451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42673972) q[2];
sx q[2];
rz(-1.4740976) q[2];
sx q[2];
rz(1.1023785) q[2];
rz(-2.0057996) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(-0.77146012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.8847467) q[0];
sx q[0];
rz(-2.1717635) q[0];
sx q[0];
rz(-0.92887512) q[0];
rz(2.7595787) q[1];
sx q[1];
rz(-0.99645749) q[1];
sx q[1];
rz(1.6928203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9858157) q[0];
sx q[0];
rz(-0.57812968) q[0];
sx q[0];
rz(-2.5647519) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11547757) q[2];
sx q[2];
rz(-0.48764569) q[2];
sx q[2];
rz(0.46610006) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4871028) q[1];
sx q[1];
rz(-1.9675323) q[1];
sx q[1];
rz(-3.023772) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9881416) q[3];
sx q[3];
rz(-2.7355237) q[3];
sx q[3];
rz(-1.690133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53943071) q[2];
sx q[2];
rz(-0.346589) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(1.8481988) q[3];
sx q[3];
rz(-2.4496205) q[3];
sx q[3];
rz(-0.20492157) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(-2.8906004) q[0];
rz(-2.9639066) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(-2.4846855) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278695) q[0];
sx q[0];
rz(-1.8060469) q[0];
sx q[0];
rz(-2.9293438) q[0];
rz(1.7288293) q[2];
sx q[2];
rz(-0.92301805) q[2];
sx q[2];
rz(2.4142746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4272139) q[1];
sx q[1];
rz(-1.2717112) q[1];
sx q[1];
rz(-2.7370791) q[1];
rz(-pi) q[2];
rz(-0.64539306) q[3];
sx q[3];
rz(-1.2155967) q[3];
sx q[3];
rz(1.8339637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3333007) q[2];
sx q[2];
rz(-2.5623645) q[2];
sx q[2];
rz(2.2283238) q[2];
rz(2.463786) q[3];
sx q[3];
rz(-1.5014239) q[3];
sx q[3];
rz(-2.138413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.102757) q[0];
sx q[0];
rz(-1.9982194) q[0];
sx q[0];
rz(-2.5586149) q[0];
rz(-1.673117) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(0.34117064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7852029) q[0];
sx q[0];
rz(-1.613253) q[0];
sx q[0];
rz(-1.7701985) q[0];
x q[1];
rz(-0.24218817) q[2];
sx q[2];
rz(-1.846284) q[2];
sx q[2];
rz(-1.7603859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64168186) q[1];
sx q[1];
rz(-1.2044831) q[1];
sx q[1];
rz(1.1745499) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0811133) q[3];
sx q[3];
rz(-2.4016909) q[3];
sx q[3];
rz(0.56671732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94656452) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(2.3348746) q[2];
rz(0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(-2.257982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44613999) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(-2.9651508) q[0];
rz(0.31013075) q[1];
sx q[1];
rz(-1.4896432) q[1];
sx q[1];
rz(-2.2055221) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1131206) q[0];
sx q[0];
rz(-0.33361379) q[0];
sx q[0];
rz(2.6718475) q[0];
rz(0.62159448) q[2];
sx q[2];
rz(-0.69658579) q[2];
sx q[2];
rz(2.7736349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99451274) q[1];
sx q[1];
rz(-1.6638866) q[1];
sx q[1];
rz(1.2017815) q[1];
rz(-pi) q[2];
rz(-0.12054875) q[3];
sx q[3];
rz(-1.0998187) q[3];
sx q[3];
rz(0.69191832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.04756847) q[2];
sx q[2];
rz(-1.8393643) q[2];
sx q[2];
rz(3.1330718) q[2];
rz(-0.73355567) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(-0.12695299) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215288) q[0];
sx q[0];
rz(-2.7286752) q[0];
sx q[0];
rz(0.76989663) q[0];
rz(-0.13433111) q[1];
sx q[1];
rz(-0.7901935) q[1];
sx q[1];
rz(1.92164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1320187) q[0];
sx q[0];
rz(-1.2694799) q[0];
sx q[0];
rz(1.8150041) q[0];
rz(-pi) q[1];
rz(-2.5346181) q[2];
sx q[2];
rz(-1.3262981) q[2];
sx q[2];
rz(2.2331657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97970672) q[1];
sx q[1];
rz(-0.58459133) q[1];
sx q[1];
rz(2.4120055) q[1];
rz(-pi) q[2];
rz(0.071000428) q[3];
sx q[3];
rz(-1.7216847) q[3];
sx q[3];
rz(-1.7689266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28025815) q[2];
sx q[2];
rz(-2.4836149) q[2];
sx q[2];
rz(0.78557837) q[2];
rz(2.806459) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(2.3194763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817374) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(0.68589504) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(-1.5790719) q[2];
sx q[2];
rz(-2.5604421) q[2];
sx q[2];
rz(1.5887518) q[2];
rz(3.104628) q[3];
sx q[3];
rz(-1.8633458) q[3];
sx q[3];
rz(2.9316791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
