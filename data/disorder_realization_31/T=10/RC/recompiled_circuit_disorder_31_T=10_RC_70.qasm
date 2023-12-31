OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9487171) q[0];
sx q[0];
rz(-2.2964381) q[0];
sx q[0];
rz(-1.2851508) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1137142) q[2];
sx q[2];
rz(-0.61402245) q[2];
sx q[2];
rz(-0.30464722) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9021437) q[1];
sx q[1];
rz(-1.7001171) q[1];
sx q[1];
rz(2.7620478) q[1];
x q[2];
rz(-2.9620142) q[3];
sx q[3];
rz(-2.539733) q[3];
sx q[3];
rz(-1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(-2.4029762) q[0];
rz(-pi) q[1];
rz(-0.39312675) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.4002422) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7696015) q[1];
sx q[1];
rz(-0.76906119) q[1];
sx q[1];
rz(0.12188697) q[1];
rz(-1.3080018) q[3];
sx q[3];
rz(-1.3920708) q[3];
sx q[3];
rz(-2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-2.2581805) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.6261684) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0059144817) q[0];
sx q[0];
rz(-2.1164829) q[0];
sx q[0];
rz(0.90555993) q[0];
rz(-0.91471471) q[2];
sx q[2];
rz(-1.9064184) q[2];
sx q[2];
rz(-2.6732973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25635438) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(0.79856915) q[1];
rz(-pi) q[2];
rz(3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(2.8312347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(-1.7128574) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(-3.0435261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1102317) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.6237753) q[1];
rz(-pi) q[2];
x q[2];
rz(0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(-1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244959) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(-0.32075551) q[0];
rz(-0.69663163) q[2];
sx q[2];
rz(-1.4259035) q[2];
sx q[2];
rz(-2.2640413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1340027) q[1];
sx q[1];
rz(-0.78467272) q[1];
sx q[1];
rz(0.44962928) q[1];
rz(-0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(-2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.822086) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(-0.23131891) q[0];
rz(-pi) q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(-2.2350603) q[1];
rz(-pi) q[2];
rz(1.3098573) q[3];
sx q[3];
rz(-0.3762227) q[3];
sx q[3];
rz(2.6747243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-0.033989865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(2.3539761) q[0];
x q[1];
rz(2.6007973) q[2];
sx q[2];
rz(-1.0058837) q[2];
sx q[2];
rz(-2.5164547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9719203) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(-2.2488942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4230698) q[3];
sx q[3];
rz(-2.1805694) q[3];
sx q[3];
rz(-2.8019398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96955339) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(-1.1374377) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82463512) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(0.075721272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8135012) q[1];
sx q[1];
rz(-2.5655167) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94533841) q[3];
sx q[3];
rz(-1.006554) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(-0.30977419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99684925) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(-1.4247308) q[0];
x q[1];
rz(0.0094527761) q[2];
sx q[2];
rz(-1.3938892) q[2];
sx q[2];
rz(-2.0131907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5619547) q[1];
sx q[1];
rz(-1.6288174) q[1];
sx q[1];
rz(-1.3647563) q[1];
x q[2];
rz(1.314332) q[3];
sx q[3];
rz(-0.98955742) q[3];
sx q[3];
rz(2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.6419798) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4831055) q[0];
sx q[0];
rz(-2.8327836) q[0];
sx q[0];
rz(-1.9103861) q[0];
rz(3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(0.27809696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.75913945) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(0.23666246) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1390926) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(0.97027422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(1.0735738) q[2];
sx q[2];
rz(-1.1560658) q[2];
sx q[2];
rz(1.8552468) q[2];
rz(-0.55340135) q[3];
sx q[3];
rz(-0.80080606) q[3];
sx q[3];
rz(-1.3047119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
