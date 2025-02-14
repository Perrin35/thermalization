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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37346807) q[0];
sx q[0];
rz(-0.22499496) q[0];
sx q[0];
rz(-0.24562545) q[0];
rz(1.6263715) q[2];
sx q[2];
rz(-1.4204506) q[2];
sx q[2];
rz(2.0932582) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5049274) q[1];
sx q[1];
rz(-2.0940771) q[1];
sx q[1];
rz(3.0806171) q[1];
rz(-pi) q[2];
rz(2.8992735) q[3];
sx q[3];
rz(-1.6166787) q[3];
sx q[3];
rz(2.322779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91288519) q[2];
sx q[2];
rz(-3.0444453) q[2];
sx q[2];
rz(2.482282) q[2];
rz(-0.77624503) q[3];
sx q[3];
rz(-3.1228437) q[3];
sx q[3];
rz(2.4565878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9126251) q[0];
sx q[0];
rz(-1.2094867) q[0];
sx q[0];
rz(1.9483161) q[0];
rz(-3.0972262) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(-0.22656974) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18767522) q[0];
sx q[0];
rz(-1.5724475) q[0];
sx q[0];
rz(-1.5742013) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5878651) q[2];
sx q[2];
rz(-1.1917795) q[2];
sx q[2];
rz(-1.5660777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4227155) q[1];
sx q[1];
rz(-1.1402292) q[1];
sx q[1];
rz(1.4925455) q[1];
x q[2];
rz(1.3188521) q[3];
sx q[3];
rz(-2.3980707) q[3];
sx q[3];
rz(-1.8094667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.442753) q[2];
sx q[2];
rz(-1.5719465) q[2];
sx q[2];
rz(1.5296096) q[2];
rz(-2.2085341) q[3];
sx q[3];
rz(-1.482684) q[3];
sx q[3];
rz(2.9034767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17732349) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(0.200287) q[0];
rz(-3.1411723) q[1];
sx q[1];
rz(-2.2047408) q[1];
sx q[1];
rz(3.1281085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1258365) q[0];
sx q[0];
rz(-1.635433) q[0];
sx q[0];
rz(1.9958853) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5277943) q[2];
sx q[2];
rz(-1.637405) q[2];
sx q[2];
rz(-0.013184908) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2993907) q[1];
sx q[1];
rz(-0.1526356) q[1];
sx q[1];
rz(-1.0959036) q[1];
rz(1.5192658) q[3];
sx q[3];
rz(-1.7331707) q[3];
sx q[3];
rz(-0.80013093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1991594) q[2];
sx q[2];
rz(-1.5853256) q[2];
sx q[2];
rz(-1.5419434) q[2];
rz(1.0017627) q[3];
sx q[3];
rz(-0.21654138) q[3];
sx q[3];
rz(-0.17222968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7672985) q[0];
sx q[0];
rz(-2.9758487) q[0];
sx q[0];
rz(-0.34749183) q[0];
rz(-2.5047498) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(-1.9510795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74032846) q[0];
sx q[0];
rz(-1.51121) q[0];
sx q[0];
rz(1.4387095) q[0];
x q[1];
rz(-1.6396825) q[2];
sx q[2];
rz(-1.578555) q[2];
sx q[2];
rz(0.0077645609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0034297) q[1];
sx q[1];
rz(-0.82050475) q[1];
sx q[1];
rz(1.8292887) q[1];
rz(1.7905964) q[3];
sx q[3];
rz(-1.2145208) q[3];
sx q[3];
rz(-1.7770467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5845329) q[2];
sx q[2];
rz(-0.0396885) q[2];
sx q[2];
rz(1.7812799) q[2];
rz(1.4761866) q[3];
sx q[3];
rz(-1.5675661) q[3];
sx q[3];
rz(0.40569693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0535468) q[0];
sx q[0];
rz(-2.4021554) q[0];
sx q[0];
rz(-0.0060225688) q[0];
rz(-1.7209523) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(0.086070148) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6688441) q[0];
sx q[0];
rz(-1.5875582) q[0];
sx q[0];
rz(-3.0344323) q[0];
rz(2.6997296) q[2];
sx q[2];
rz(-3.0806858) q[2];
sx q[2];
rz(-1.4373844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9169315) q[1];
sx q[1];
rz(-1.6078533) q[1];
sx q[1];
rz(-0.091450973) q[1];
x q[2];
rz(-2.213201) q[3];
sx q[3];
rz(-0.07559055) q[3];
sx q[3];
rz(2.9380611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.033279557) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(-1.7576199) q[2];
rz(-0.39517394) q[3];
sx q[3];
rz(-0.049001781) q[3];
sx q[3];
rz(-1.9743732) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60642099) q[0];
sx q[0];
rz(-0.21927729) q[0];
sx q[0];
rz(-2.1411335) q[0];
rz(-0.74042997) q[1];
sx q[1];
rz(-0.43559566) q[1];
sx q[1];
rz(2.7620517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6839121) q[0];
sx q[0];
rz(-3.0427986) q[0];
sx q[0];
rz(-2.8794419) q[0];
rz(1.8431115) q[2];
sx q[2];
rz(-2.0305995) q[2];
sx q[2];
rz(-2.1260171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21793903) q[1];
sx q[1];
rz(-0.50461136) q[1];
sx q[1];
rz(1.6259471) q[1];
rz(-3.0256512) q[3];
sx q[3];
rz(-1.1244457) q[3];
sx q[3];
rz(1.6801113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.015811054) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(-3.06456) q[2];
rz(-0.020126255) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(-2.3441815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1064442) q[0];
sx q[0];
rz(-3.1290717) q[0];
sx q[0];
rz(1.4323819) q[0];
rz(-2.8445981) q[1];
sx q[1];
rz(-2.98731) q[1];
sx q[1];
rz(3.0800379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81082329) q[0];
sx q[0];
rz(-2.4274438) q[0];
sx q[0];
rz(-0.9261976) q[0];
x q[1];
rz(-0.21451471) q[2];
sx q[2];
rz(-1.5652862) q[2];
sx q[2];
rz(0.44098976) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48744943) q[1];
sx q[1];
rz(-1.2795078) q[1];
sx q[1];
rz(-1.4994411) q[1];
x q[2];
rz(-0.43470862) q[3];
sx q[3];
rz(-1.5366597) q[3];
sx q[3];
rz(3.1354836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2986472) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(-1.418815) q[2];
rz(-1.8055387) q[3];
sx q[3];
rz(-0.027848363) q[3];
sx q[3];
rz(-1.4068039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0272738) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(1.640821) q[0];
rz(2.0188792) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.8480802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155535) q[0];
sx q[0];
rz(-2.3991971) q[0];
sx q[0];
rz(-2.0406988) q[0];
x q[1];
rz(2.7377364) q[2];
sx q[2];
rz(-2.4060898) q[2];
sx q[2];
rz(2.1063358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64206865) q[1];
sx q[1];
rz(-1.7392989) q[1];
sx q[1];
rz(-2.8649855) q[1];
rz(-2.3154852) q[3];
sx q[3];
rz(-2.2650121) q[3];
sx q[3];
rz(2.7339767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5844172) q[2];
sx q[2];
rz(-2.7807005) q[2];
sx q[2];
rz(1.3593675) q[2];
rz(-2.6946097) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(-2.9744448) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2928325) q[0];
sx q[0];
rz(-0.4378795) q[0];
sx q[0];
rz(2.6614702) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6933889) q[2];
sx q[2];
rz(-1.685678) q[2];
sx q[2];
rz(-0.5960532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0648411) q[1];
sx q[1];
rz(-3.1411618) q[1];
sx q[1];
rz(-1.3382124) q[1];
rz(-pi) q[2];
rz(-0.4550981) q[3];
sx q[3];
rz(-1.4346204) q[3];
sx q[3];
rz(1.0973339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0000275) q[2];
sx q[2];
rz(-3.1391149) q[2];
sx q[2];
rz(2.4139717) q[2];
rz(1.242312) q[3];
sx q[3];
rz(-3.1051903) q[3];
sx q[3];
rz(1.9696994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2321371) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(2.452028) q[0];
rz(-1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(-2.9857059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0077654) q[0];
sx q[0];
rz(-1.7031914) q[0];
sx q[0];
rz(-2.4835229) q[0];
x q[1];
rz(-1.6901659) q[2];
sx q[2];
rz(-1.8205336) q[2];
sx q[2];
rz(-2.2636556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.879133) q[1];
sx q[1];
rz(-3.1371318) q[1];
sx q[1];
rz(-0.35426472) q[1];
rz(-pi) q[2];
rz(1.0267657) q[3];
sx q[3];
rz(-0.18160219) q[3];
sx q[3];
rz(-1.7498407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9547687) q[2];
sx q[2];
rz(-3.0195152) q[2];
sx q[2];
rz(0.99511498) q[2];
rz(-2.9156445) q[3];
sx q[3];
rz(-3.093231) q[3];
sx q[3];
rz(0.82609716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7786998) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(1.7097991) q[1];
sx q[1];
rz(-1.8451537) q[1];
sx q[1];
rz(0.61737212) q[1];
rz(0.16724965) q[2];
sx q[2];
rz(-2.4495301) q[2];
sx q[2];
rz(-2.4514103) q[2];
rz(-1.9598087) q[3];
sx q[3];
rz(-1.4445732) q[3];
sx q[3];
rz(-2.5635535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
