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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(1.992835) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(1.6869071) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8073976) q[0];
sx q[0];
rz(-0.79803033) q[0];
sx q[0];
rz(-1.0502771) q[0];
rz(-pi) q[1];
rz(-3.0266255) q[2];
sx q[2];
rz(-0.76412725) q[2];
sx q[2];
rz(2.4121491) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4920884) q[1];
sx q[1];
rz(-1.7028812) q[1];
sx q[1];
rz(2.3852939) q[1];
rz(0.65806453) q[3];
sx q[3];
rz(-1.4941553) q[3];
sx q[3];
rz(1.6057526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9601606) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(-1.9350447) q[2];
rz(-1.1613965) q[3];
sx q[3];
rz(-0.56777081) q[3];
sx q[3];
rz(-1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-1.125536) q[0];
sx q[0];
rz(-0.97208446) q[0];
rz(0.26845911) q[1];
sx q[1];
rz(-1.4413036) q[1];
sx q[1];
rz(0.45480248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80538023) q[0];
sx q[0];
rz(-2.8180362) q[0];
sx q[0];
rz(-1.1624883) q[0];
rz(-1.201638) q[2];
sx q[2];
rz(-2.7281986) q[2];
sx q[2];
rz(-1.3184551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3784967) q[1];
sx q[1];
rz(-0.46751577) q[1];
sx q[1];
rz(0.11390953) q[1];
x q[2];
rz(2.5760994) q[3];
sx q[3];
rz(-1.0988937) q[3];
sx q[3];
rz(1.3513235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1099757) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(-0.99962437) q[2];
rz(-0.085974606) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(3.1302248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.40083945) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(2.9174347) q[0];
rz(1.594918) q[1];
sx q[1];
rz(-1.5432576) q[1];
sx q[1];
rz(-1.8348144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6431893) q[0];
sx q[0];
rz(-2.1191594) q[0];
sx q[0];
rz(-1.7140618) q[0];
x q[1];
rz(2.3036912) q[2];
sx q[2];
rz(-2.1616677) q[2];
sx q[2];
rz(1.5082347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0780454) q[1];
sx q[1];
rz(-1.987146) q[1];
sx q[1];
rz(-2.0267118) q[1];
rz(-0.66931386) q[3];
sx q[3];
rz(-1.770854) q[3];
sx q[3];
rz(-1.1049096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.33009067) q[2];
sx q[2];
rz(-1.6197438) q[2];
sx q[2];
rz(-1.0203934) q[2];
rz(-1.329782) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(1.6167195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.735585) q[0];
sx q[0];
rz(-2.0841053) q[0];
sx q[0];
rz(-1.965858) q[0];
rz(-2.8658087) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(2.5036459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9383134) q[0];
sx q[0];
rz(-1.5601451) q[0];
sx q[0];
rz(3.108884) q[0];
x q[1];
rz(-1.3512072) q[2];
sx q[2];
rz(-0.96115005) q[2];
sx q[2];
rz(2.3860506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.040311) q[1];
sx q[1];
rz(-0.71523803) q[1];
sx q[1];
rz(0.50846993) q[1];
rz(-pi) q[2];
rz(0.27880554) q[3];
sx q[3];
rz(-0.33167095) q[3];
sx q[3];
rz(1.7808033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5059169) q[2];
sx q[2];
rz(-1.5636779) q[2];
sx q[2];
rz(1.6944616) q[2];
rz(2.9315089) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90784812) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(-2.0727169) q[0];
rz(1.8416038) q[1];
sx q[1];
rz(-1.735382) q[1];
sx q[1];
rz(-1.9357505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763472) q[0];
sx q[0];
rz(-1.6505213) q[0];
sx q[0];
rz(-1.1673156) q[0];
x q[1];
rz(-1.5583493) q[2];
sx q[2];
rz(-0.65285359) q[2];
sx q[2];
rz(1.659153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1002378) q[1];
sx q[1];
rz(-2.2861028) q[1];
sx q[1];
rz(0.13731401) q[1];
x q[2];
rz(-2.8338594) q[3];
sx q[3];
rz(-0.50129563) q[3];
sx q[3];
rz(-0.0020108027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.892889) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(-0.068566337) q[2];
rz(2.2019703) q[3];
sx q[3];
rz(-1.4087804) q[3];
sx q[3];
rz(-1.2287963) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3786316) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(-0.5214386) q[0];
rz(1.6261082) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(-0.62561402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8541295) q[0];
sx q[0];
rz(-2.1105327) q[0];
sx q[0];
rz(3.1053393) q[0];
rz(-pi) q[1];
rz(-2.8270673) q[2];
sx q[2];
rz(-0.39681602) q[2];
sx q[2];
rz(-2.0443969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7661356) q[1];
sx q[1];
rz(-1.7682255) q[1];
sx q[1];
rz(1.1079815) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7196521) q[3];
sx q[3];
rz(-1.6263537) q[3];
sx q[3];
rz(-2.0618604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0079415) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(-3.0360743) q[2];
rz(-2.1775235) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(-2.156637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0725726) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(-1.7286812) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(0.55317318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39554292) q[0];
sx q[0];
rz(-0.93727124) q[0];
sx q[0];
rz(-0.49479158) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14153962) q[2];
sx q[2];
rz(-1.6602548) q[2];
sx q[2];
rz(-0.78023655) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87482086) q[1];
sx q[1];
rz(-1.2874585) q[1];
sx q[1];
rz(1.6127519) q[1];
rz(2.9384261) q[3];
sx q[3];
rz(-1.3152514) q[3];
sx q[3];
rz(2.0269964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-0.41213948) q[2];
rz(-2.4814217) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039345) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(2.9397021) q[0];
rz(-2.0837325) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(-2.8813664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5903918) q[0];
sx q[0];
rz(-2.0856212) q[0];
sx q[0];
rz(2.0082672) q[0];
rz(2.260602) q[2];
sx q[2];
rz(-0.63369232) q[2];
sx q[2];
rz(-1.4680167) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8152129) q[1];
sx q[1];
rz(-1.9163487) q[1];
sx q[1];
rz(-1.2735081) q[1];
rz(1.4020355) q[3];
sx q[3];
rz(-2.2243554) q[3];
sx q[3];
rz(0.26168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95320931) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(2.5592213) q[2];
rz(-0.44238704) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(0.83627397) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22170947) q[0];
sx q[0];
rz(-1.6136805) q[0];
sx q[0];
rz(-1.7145994) q[0];
rz(1.3555948) q[1];
sx q[1];
rz(-1.4878863) q[1];
sx q[1];
rz(-1.7048763) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6528931) q[0];
sx q[0];
rz(-2.821229) q[0];
sx q[0];
rz(-1.5071177) q[0];
rz(-pi) q[1];
rz(1.270601) q[2];
sx q[2];
rz(-2.0328902) q[2];
sx q[2];
rz(-1.0754881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7683803) q[1];
sx q[1];
rz(-1.5564959) q[1];
sx q[1];
rz(-1.3296024) q[1];
rz(-1.837265) q[3];
sx q[3];
rz(-2.5992166) q[3];
sx q[3];
rz(-1.8625384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23472486) q[2];
sx q[2];
rz(-1.8610672) q[2];
sx q[2];
rz(-1.6024164) q[2];
rz(-2.5336044) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.7749087) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(-0.82345024) q[0];
rz(2.2460294) q[1];
sx q[1];
rz(-1.3500373) q[1];
sx q[1];
rz(-0.38111883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32981634) q[0];
sx q[0];
rz(-2.1069706) q[0];
sx q[0];
rz(-2.0945805) q[0];
rz(-pi) q[1];
rz(0.37883436) q[2];
sx q[2];
rz(-0.4288579) q[2];
sx q[2];
rz(-0.90012401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19958147) q[1];
sx q[1];
rz(-1.8560272) q[1];
sx q[1];
rz(-0.99883325) q[1];
rz(-1.0166083) q[3];
sx q[3];
rz(-2.1234406) q[3];
sx q[3];
rz(1.2070398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33619189) q[2];
sx q[2];
rz(-2.7298584) q[2];
sx q[2];
rz(2.1545048) q[2];
rz(1.5385212) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(-0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9170452) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(-0.70980258) q[1];
sx q[1];
rz(-0.55955049) q[1];
sx q[1];
rz(1.4581663) q[1];
rz(0.3644402) q[2];
sx q[2];
rz(-1.6779998) q[2];
sx q[2];
rz(-2.6371434) q[2];
rz(-0.74123989) q[3];
sx q[3];
rz(-2.0251353) q[3];
sx q[3];
rz(0.97788772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
