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
rz(-0.99011079) q[0];
sx q[0];
rz(-2.4408542) q[0];
sx q[0];
rz(-2.3837958) q[0];
rz(-2.4611729) q[1];
sx q[1];
rz(-1.0935723) q[1];
sx q[1];
rz(3.0419066) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133649) q[0];
sx q[0];
rz(-1.037957) q[0];
sx q[0];
rz(-1.9060978) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9366936) q[2];
sx q[2];
rz(-2.120595) q[2];
sx q[2];
rz(-1.1680574) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1336109) q[1];
sx q[1];
rz(-0.33002285) q[1];
sx q[1];
rz(2.2225478) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7277579) q[3];
sx q[3];
rz(-0.36245868) q[3];
sx q[3];
rz(-2.3316104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0199355) q[2];
sx q[2];
rz(-1.1172373) q[2];
sx q[2];
rz(-0.026148671) q[2];
rz(1.3610241) q[3];
sx q[3];
rz(-2.0676421) q[3];
sx q[3];
rz(1.1092383) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015198135) q[0];
sx q[0];
rz(-0.69271883) q[0];
sx q[0];
rz(-1.3808845) q[0];
rz(-2.8610002) q[1];
sx q[1];
rz(-2.3255489) q[1];
sx q[1];
rz(-0.38995829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1564962) q[0];
sx q[0];
rz(-2.3802983) q[0];
sx q[0];
rz(3.0846315) q[0];
x q[1];
rz(0.89050129) q[2];
sx q[2];
rz(-1.599031) q[2];
sx q[2];
rz(0.75301127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9401318) q[1];
sx q[1];
rz(-1.6319449) q[1];
sx q[1];
rz(0.1810302) q[1];
rz(-pi) q[2];
rz(1.5481351) q[3];
sx q[3];
rz(-1.0575599) q[3];
sx q[3];
rz(3.0340241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0754764) q[2];
sx q[2];
rz(-2.7710997) q[2];
sx q[2];
rz(2.9541435) q[2];
rz(2.1665393) q[3];
sx q[3];
rz(-1.8911898) q[3];
sx q[3];
rz(1.8072849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27075574) q[0];
sx q[0];
rz(-2.2360531) q[0];
sx q[0];
rz(-1.723473) q[0];
rz(-2.1899147) q[1];
sx q[1];
rz(-2.8680809) q[1];
sx q[1];
rz(-2.4804514) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9959629) q[0];
sx q[0];
rz(-1.055777) q[0];
sx q[0];
rz(-1.4670232) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9533526) q[2];
sx q[2];
rz(-0.91209665) q[2];
sx q[2];
rz(3.0223568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74562711) q[1];
sx q[1];
rz(-0.14603171) q[1];
sx q[1];
rz(1.3457937) q[1];
rz(-2.5824359) q[3];
sx q[3];
rz(-1.4803018) q[3];
sx q[3];
rz(1.4714946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4989138) q[2];
sx q[2];
rz(-1.4475409) q[2];
sx q[2];
rz(2.101669) q[2];
rz(-2.6652279) q[3];
sx q[3];
rz(-0.56291181) q[3];
sx q[3];
rz(-2.3001455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845402) q[0];
sx q[0];
rz(-0.81890506) q[0];
sx q[0];
rz(2.1475329) q[0];
rz(-2.0301863) q[1];
sx q[1];
rz(-1.2329654) q[1];
sx q[1];
rz(2.9434189) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4497714) q[0];
sx q[0];
rz(-1.9110962) q[0];
sx q[0];
rz(-1.4744076) q[0];
rz(-pi) q[1];
rz(0.54951602) q[2];
sx q[2];
rz(-1.9870111) q[2];
sx q[2];
rz(0.93577318) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6804643) q[1];
sx q[1];
rz(-1.8682419) q[1];
sx q[1];
rz(-0.72338413) q[1];
rz(-1.0425838) q[3];
sx q[3];
rz(-1.5798395) q[3];
sx q[3];
rz(-2.8419055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17890945) q[2];
sx q[2];
rz(-1.6390025) q[2];
sx q[2];
rz(-2.7965014) q[2];
rz(2.1033449) q[3];
sx q[3];
rz(-2.0227261) q[3];
sx q[3];
rz(3.069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8615231) q[0];
sx q[0];
rz(-0.11174209) q[0];
sx q[0];
rz(0.34568632) q[0];
rz(-3.0142036) q[1];
sx q[1];
rz(-2.7368339) q[1];
sx q[1];
rz(1.5397127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9939748) q[0];
sx q[0];
rz(-0.26122005) q[0];
sx q[0];
rz(1.0077976) q[0];
rz(2.0284838) q[2];
sx q[2];
rz(-2.1808592) q[2];
sx q[2];
rz(-1.318536) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3469763) q[1];
sx q[1];
rz(-1.8759751) q[1];
sx q[1];
rz(0.48726122) q[1];
x q[2];
rz(-2.5681189) q[3];
sx q[3];
rz(-1.4961362) q[3];
sx q[3];
rz(-0.71686137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4442399) q[2];
sx q[2];
rz(-1.3109861) q[2];
sx q[2];
rz(-2.7250807) q[2];
rz(3.0537649) q[3];
sx q[3];
rz(-2.646324) q[3];
sx q[3];
rz(0.4270359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2764565) q[0];
sx q[0];
rz(-0.75053954) q[0];
sx q[0];
rz(-2.1476792) q[0];
rz(-1.7991426) q[1];
sx q[1];
rz(-1.7667222) q[1];
sx q[1];
rz(-1.3329175) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7291746) q[0];
sx q[0];
rz(-1.6977457) q[0];
sx q[0];
rz(0.35611313) q[0];
x q[1];
rz(-2.7287675) q[2];
sx q[2];
rz(-1.7805837) q[2];
sx q[2];
rz(-1.0722463) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2675587) q[1];
sx q[1];
rz(-1.7090952) q[1];
sx q[1];
rz(-0.036768032) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52458472) q[3];
sx q[3];
rz(-0.3915873) q[3];
sx q[3];
rz(3.0767431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6028613) q[2];
sx q[2];
rz(-1.8659325) q[2];
sx q[2];
rz(-2.4219647) q[2];
rz(-1.0031649) q[3];
sx q[3];
rz(-1.4041308) q[3];
sx q[3];
rz(-0.42705944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3908983) q[0];
sx q[0];
rz(-2.7321132) q[0];
sx q[0];
rz(2.3636708) q[0];
rz(1.3748417) q[1];
sx q[1];
rz(-2.1781616) q[1];
sx q[1];
rz(2.8533997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2552528) q[0];
sx q[0];
rz(-1.985441) q[0];
sx q[0];
rz(-0.1840846) q[0];
rz(-pi) q[1];
rz(-0.32518816) q[2];
sx q[2];
rz(-3.1085494) q[2];
sx q[2];
rz(0.059800241) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58449897) q[1];
sx q[1];
rz(-0.58773041) q[1];
sx q[1];
rz(2.6705301) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1918357) q[3];
sx q[3];
rz(-1.4498447) q[3];
sx q[3];
rz(-0.17571501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2253183) q[2];
sx q[2];
rz(-2.306873) q[2];
sx q[2];
rz(1.5073331) q[2];
rz(2.7507014) q[3];
sx q[3];
rz(-1.2619799) q[3];
sx q[3];
rz(-1.619537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36562076) q[0];
sx q[0];
rz(-0.93637744) q[0];
sx q[0];
rz(0.17882624) q[0];
rz(1.3911744) q[1];
sx q[1];
rz(-1.4427253) q[1];
sx q[1];
rz(-2.9475382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5882307) q[0];
sx q[0];
rz(-1.0798619) q[0];
sx q[0];
rz(0.42598283) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3909245) q[2];
sx q[2];
rz(-1.0714415) q[2];
sx q[2];
rz(2.4908092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68264254) q[1];
sx q[1];
rz(-1.7006665) q[1];
sx q[1];
rz(0.81241205) q[1];
rz(-pi) q[2];
rz(-1.3369183) q[3];
sx q[3];
rz(-0.66576427) q[3];
sx q[3];
rz(0.0364937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8630744) q[2];
sx q[2];
rz(-0.14894177) q[2];
sx q[2];
rz(1.8982559) q[2];
rz(0.30531034) q[3];
sx q[3];
rz(-1.6122626) q[3];
sx q[3];
rz(0.39795157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0968001) q[0];
sx q[0];
rz(-3.0718006) q[0];
sx q[0];
rz(-3.0332562) q[0];
rz(-2.4414252) q[1];
sx q[1];
rz(-1.533875) q[1];
sx q[1];
rz(0.38665032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17849191) q[0];
sx q[0];
rz(-1.0722536) q[0];
sx q[0];
rz(0.83591423) q[0];
x q[1];
rz(-2.9896624) q[2];
sx q[2];
rz(-1.4174089) q[2];
sx q[2];
rz(1.5832324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79124444) q[1];
sx q[1];
rz(-1.2859259) q[1];
sx q[1];
rz(2.4288252) q[1];
x q[2];
rz(-0.4824494) q[3];
sx q[3];
rz(-2.1157935) q[3];
sx q[3];
rz(1.4712131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.391523) q[2];
sx q[2];
rz(-1.1015026) q[2];
sx q[2];
rz(-2.2186685) q[2];
rz(2.7643438) q[3];
sx q[3];
rz(-2.173285) q[3];
sx q[3];
rz(1.3737465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7467576) q[0];
sx q[0];
rz(-2.6834798) q[0];
sx q[0];
rz(-2.3858892) q[0];
rz(2.34756) q[1];
sx q[1];
rz(-0.67818063) q[1];
sx q[1];
rz(1.1620577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1235706) q[0];
sx q[0];
rz(-1.7777052) q[0];
sx q[0];
rz(-1.9777954) q[0];
x q[1];
rz(-0.5616582) q[2];
sx q[2];
rz(-2.409081) q[2];
sx q[2];
rz(1.2309625) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2311467) q[1];
sx q[1];
rz(-0.9462463) q[1];
sx q[1];
rz(-2.7471119) q[1];
x q[2];
rz(-2.1697767) q[3];
sx q[3];
rz(-0.47656588) q[3];
sx q[3];
rz(-3.0314881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0843167) q[2];
sx q[2];
rz(-0.82547775) q[2];
sx q[2];
rz(-3.0153073) q[2];
rz(2.5430191) q[3];
sx q[3];
rz(-1.6152265) q[3];
sx q[3];
rz(-1.1461145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79692688) q[0];
sx q[0];
rz(-1.4753337) q[0];
sx q[0];
rz(-1.4519806) q[0];
rz(1.7563734) q[1];
sx q[1];
rz(-1.1744371) q[1];
sx q[1];
rz(-0.074443346) q[1];
rz(-0.66567265) q[2];
sx q[2];
rz(-0.72279983) q[2];
sx q[2];
rz(1.3574251) q[2];
rz(2.606288) q[3];
sx q[3];
rz(-2.6953831) q[3];
sx q[3];
rz(2.7627177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
