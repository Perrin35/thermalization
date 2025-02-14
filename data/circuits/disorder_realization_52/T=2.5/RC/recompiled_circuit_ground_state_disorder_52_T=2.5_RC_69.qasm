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
rz(-1.787552) q[1];
sx q[1];
rz(-1.6156337) q[1];
sx q[1];
rz(-1.807133) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37346807) q[0];
sx q[0];
rz(-2.9165977) q[0];
sx q[0];
rz(-0.24562545) q[0];
rz(-pi) q[1];
rz(0.15057474) q[2];
sx q[2];
rz(-1.5158487) q[2];
sx q[2];
rz(2.6274632) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5150979) q[1];
sx q[1];
rz(-0.52649311) q[1];
sx q[1];
rz(-1.676031) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6180575) q[3];
sx q[3];
rz(-1.3287373) q[3];
sx q[3];
rz(-2.4009465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2287075) q[2];
sx q[2];
rz(-3.0444453) q[2];
sx q[2];
rz(2.482282) q[2];
rz(-0.77624503) q[3];
sx q[3];
rz(-0.018748911) q[3];
sx q[3];
rz(-2.4565878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9126251) q[0];
sx q[0];
rz(-1.2094867) q[0];
sx q[0];
rz(1.1932766) q[0];
rz(-3.0972262) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(-0.22656974) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3069585) q[0];
sx q[0];
rz(-3.1378085) q[0];
sx q[0];
rz(2.0223122) q[0];
x q[1];
rz(-2.7625257) q[2];
sx q[2];
rz(-1.5866536) q[2];
sx q[2];
rz(-0.011034688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5331796) q[1];
sx q[1];
rz(-0.4371818) q[1];
sx q[1];
rz(2.9730148) q[1];
x q[2];
rz(-2.2983589) q[3];
sx q[3];
rz(-1.7403462) q[3];
sx q[3];
rz(2.7156626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.442753) q[2];
sx q[2];
rz(-1.5719465) q[2];
sx q[2];
rz(1.6119831) q[2];
rz(0.93305856) q[3];
sx q[3];
rz(-1.482684) q[3];
sx q[3];
rz(2.9034767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17732349) q[0];
sx q[0];
rz(-3.1260335) q[0];
sx q[0];
rz(0.200287) q[0];
rz(0.00042032584) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(-3.1281085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015756135) q[0];
sx q[0];
rz(-1.5061597) q[0];
sx q[0];
rz(-1.1457074) q[0];
rz(-2.5691367) q[2];
sx q[2];
rz(-3.0623263) q[2];
sx q[2];
rz(2.5545189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.779055) q[1];
sx q[1];
rz(-1.7064306) q[1];
sx q[1];
rz(-3.0713697) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30463574) q[3];
sx q[3];
rz(-0.17028642) q[3];
sx q[3];
rz(-2.0326322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1991594) q[2];
sx q[2];
rz(-1.556267) q[2];
sx q[2];
rz(-1.5996492) q[2];
rz(2.13983) q[3];
sx q[3];
rz(-2.9250513) q[3];
sx q[3];
rz(-0.17222968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7672985) q[0];
sx q[0];
rz(-2.9758487) q[0];
sx q[0];
rz(-0.34749183) q[0];
rz(2.5047498) q[1];
sx q[1];
rz(-0.0056191365) q[1];
sx q[1];
rz(1.1905131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83837968) q[0];
sx q[0];
rz(-1.4389453) q[0];
sx q[0];
rz(0.060108713) q[0];
rz(1.4585481) q[2];
sx q[2];
rz(-0.069321037) q[2];
sx q[2];
rz(-1.6905418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0034297) q[1];
sx q[1];
rz(-0.82050475) q[1];
sx q[1];
rz(1.8292887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52996323) q[3];
sx q[3];
rz(-0.41614446) q[3];
sx q[3];
rz(0.79485369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5570598) q[2];
sx q[2];
rz(-0.0396885) q[2];
sx q[2];
rz(-1.3603127) q[2];
rz(1.4761866) q[3];
sx q[3];
rz(-1.5740266) q[3];
sx q[3];
rz(-0.40569693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0535468) q[0];
sx q[0];
rz(-2.4021554) q[0];
sx q[0];
rz(-3.1355701) q[0];
rz(-1.7209523) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(-3.0555225) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0417418) q[0];
sx q[0];
rz(-1.4636511) q[0];
sx q[0];
rz(-1.5876549) q[0];
rz(-1.5968679) q[2];
sx q[2];
rz(-1.5157454) q[2];
sx q[2];
rz(1.2616273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4115178) q[1];
sx q[1];
rz(-0.098654276) q[1];
sx q[1];
rz(0.38563557) q[1];
rz(0.045343355) q[3];
sx q[3];
rz(-1.6312977) q[3];
sx q[3];
rz(2.7013472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.033279557) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(-1.7576199) q[2];
rz(0.39517394) q[3];
sx q[3];
rz(-0.049001781) q[3];
sx q[3];
rz(1.9743732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60642099) q[0];
sx q[0];
rz(-0.21927729) q[0];
sx q[0];
rz(-2.1411335) q[0];
rz(0.74042997) q[1];
sx q[1];
rz(-2.705997) q[1];
sx q[1];
rz(-0.37954095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.767547) q[0];
sx q[0];
rz(-1.5963608) q[0];
sx q[0];
rz(-3.0461531) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6666895) q[2];
sx q[2];
rz(-1.3273718) q[2];
sx q[2];
rz(2.7096675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21793903) q[1];
sx q[1];
rz(-0.50461136) q[1];
sx q[1];
rz(-1.6259471) q[1];
x q[2];
rz(0.11594144) q[3];
sx q[3];
rz(-2.017147) q[3];
sx q[3];
rz(-1.6801113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.015811054) q[2];
sx q[2];
rz(-0.28818473) q[2];
sx q[2];
rz(-0.077032653) q[2];
rz(0.020126255) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(-0.7974112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035148419) q[0];
sx q[0];
rz(-3.1290717) q[0];
sx q[0];
rz(-1.7092108) q[0];
rz(0.2969946) q[1];
sx q[1];
rz(-2.98731) q[1];
sx q[1];
rz(-0.061554734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135573) q[0];
sx q[0];
rz(-1.0197864) q[0];
sx q[0];
rz(-2.6614266) q[0];
rz(-pi) q[1];
rz(1.5764357) q[2];
sx q[2];
rz(-1.3562849) q[2];
sx q[2];
rz(-1.1286061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8980726) q[1];
sx q[1];
rz(-0.29965934) q[1];
sx q[1];
rz(0.23345848) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.2986472) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(1.418815) q[2];
rz(-1.336054) q[3];
sx q[3];
rz(-0.027848363) q[3];
sx q[3];
rz(-1.7347887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0272738) q[0];
sx q[0];
rz(-0.71684664) q[0];
sx q[0];
rz(-1.640821) q[0];
rz(1.1227135) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.2935125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8282765) q[0];
sx q[0];
rz(-1.2596845) q[0];
sx q[0];
rz(-0.88514741) q[0];
rz(2.7377364) q[2];
sx q[2];
rz(-2.4060898) q[2];
sx q[2];
rz(2.1063358) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.499524) q[1];
sx q[1];
rz(-1.4022938) q[1];
sx q[1];
rz(-0.27660714) q[1];
rz(-pi) q[2];
x q[2];
rz(2.294306) q[3];
sx q[3];
rz(-1.022911) q[3];
sx q[3];
rz(-2.5112266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5844172) q[2];
sx q[2];
rz(-2.7807005) q[2];
sx q[2];
rz(-1.3593675) q[2];
rz(-2.6946097) q[3];
sx q[3];
rz(-3.0990661) q[3];
sx q[3];
rz(-0.16714787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3249224) q[0];
sx q[0];
rz(-1.8055547) q[0];
sx q[0];
rz(-1.100612) q[0];
rz(-2.0657516) q[1];
sx q[1];
rz(-0.64972076) q[1];
sx q[1];
rz(2.4646387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16273424) q[0];
sx q[0];
rz(-1.7679201) q[0];
sx q[0];
rz(0.39359025) q[0];
rz(-pi) q[1];
rz(-1.4482037) q[2];
sx q[2];
rz(-1.685678) q[2];
sx q[2];
rz(0.5960532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7266287) q[1];
sx q[1];
rz(-1.570697) q[1];
sx q[1];
rz(-1.5712156) q[1];
x q[2];
rz(-2.8394034) q[3];
sx q[3];
rz(-0.47366484) q[3];
sx q[3];
rz(0.74407265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1415652) q[2];
sx q[2];
rz(-3.1391149) q[2];
sx q[2];
rz(0.72762093) q[2];
rz(1.8992807) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(-1.1718933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2321371) q[0];
sx q[0];
rz(-1.0219034) q[0];
sx q[0];
rz(2.452028) q[0];
rz(-1.4762956) q[1];
sx q[1];
rz(-2.8910525) q[1];
sx q[1];
rz(0.15588674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4613381) q[0];
sx q[0];
rz(-2.2221186) q[0];
sx q[0];
rz(-1.7375577) q[0];
rz(-1.4514267) q[2];
sx q[2];
rz(-1.3210591) q[2];
sx q[2];
rz(-2.2636556) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6167276) q[1];
sx q[1];
rz(-1.5749802) q[1];
sx q[1];
rz(1.5692488) q[1];
x q[2];
rz(1.7266375) q[3];
sx q[3];
rz(-1.6644125) q[3];
sx q[3];
rz(-2.7839212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9547687) q[2];
sx q[2];
rz(-3.0195152) q[2];
sx q[2];
rz(2.1464777) q[2];
rz(0.22594813) q[3];
sx q[3];
rz(-3.093231) q[3];
sx q[3];
rz(-2.3154955) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7786998) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(1.4317935) q[1];
sx q[1];
rz(-1.2964389) q[1];
sx q[1];
rz(-2.5242205) q[1];
rz(2.974343) q[2];
sx q[2];
rz(-0.69206253) q[2];
sx q[2];
rz(0.69018232) q[2];
rz(1.9598087) q[3];
sx q[3];
rz(-1.6970194) q[3];
sx q[3];
rz(0.57803911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
