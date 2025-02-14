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
rz(0.69831508) q[0];
sx q[0];
rz(2.7628216) q[0];
sx q[0];
rz(8.2674352) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9357932) q[0];
sx q[0];
rz(-2.4895634) q[0];
sx q[0];
rz(1.0656641) q[0];
x q[1];
rz(-2.4073462) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(0.89370382) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4545645) q[1];
sx q[1];
rz(-2.6978081) q[1];
sx q[1];
rz(-0.084974809) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0105114) q[3];
sx q[3];
rz(-1.5961732) q[3];
sx q[3];
rz(2.199774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-2.989952) q[2];
rz(2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66459429) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(-0.81749302) q[0];
rz(-0.11257653) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(2.2419825) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(1.8194356) q[0];
x q[1];
rz(-0.95006534) q[2];
sx q[2];
rz(-2.143444) q[2];
sx q[2];
rz(0.060163035) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.591685) q[1];
sx q[1];
rz(-0.75061676) q[1];
sx q[1];
rz(-3.0590712) q[1];
rz(-2.1856861) q[3];
sx q[3];
rz(-2.3645176) q[3];
sx q[3];
rz(1.2082548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67109913) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.4772011) q[3];
sx q[3];
rz(3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0055493) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(1.9480202) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.4793414) q[1];
sx q[1];
rz(-1.6568291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1910716) q[0];
sx q[0];
rz(-1.5245887) q[0];
sx q[0];
rz(-3.1165857) q[0];
x q[1];
rz(-0.12642352) q[2];
sx q[2];
rz(-1.7806541) q[2];
sx q[2];
rz(-2.1239547) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79968184) q[1];
sx q[1];
rz(-1.4355625) q[1];
sx q[1];
rz(-0.34457259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2079575) q[3];
sx q[3];
rz(-1.6479744) q[3];
sx q[3];
rz(-0.046260351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9803311) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(2.4397591) q[2];
rz(0.069132239) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(2.5162146) q[0];
rz(2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(-1.8249003) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024361) q[0];
sx q[0];
rz(-2.5768752) q[0];
sx q[0];
rz(-1.8403649) q[0];
x q[1];
rz(1.6676297) q[2];
sx q[2];
rz(-1.2487354) q[2];
sx q[2];
rz(-1.2351954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.229828) q[1];
sx q[1];
rz(-0.87512866) q[1];
sx q[1];
rz(-1.8063481) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98424498) q[3];
sx q[3];
rz(-2.2797425) q[3];
sx q[3];
rz(-1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(3.1114846) q[2];
rz(-2.7249469) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(-0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3629214) q[0];
sx q[0];
rz(-1.2782949) q[0];
sx q[0];
rz(1.9238506) q[0];
rz(2.9171004) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(2.7405558) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7571791) q[0];
sx q[0];
rz(-1.9862439) q[0];
sx q[0];
rz(0.40796221) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6729092) q[2];
sx q[2];
rz(-0.71063738) q[2];
sx q[2];
rz(2.8993894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2485473) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(-0.23706146) q[1];
rz(-pi) q[2];
rz(1.4698074) q[3];
sx q[3];
rz(-1.7324395) q[3];
sx q[3];
rz(-2.9532331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2991221) q[2];
sx q[2];
rz(-1.2229908) q[2];
sx q[2];
rz(-0.50745884) q[2];
rz(0.92042813) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949263) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(-1.0194417) q[0];
rz(1.1722209) q[1];
sx q[1];
rz(-1.9688789) q[1];
sx q[1];
rz(1.2219465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98602277) q[0];
sx q[0];
rz(-1.9180074) q[0];
sx q[0];
rz(-2.272764) q[0];
x q[1];
rz(-1.6331633) q[2];
sx q[2];
rz(-2.6259239) q[2];
sx q[2];
rz(-3.1413648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3675569) q[1];
sx q[1];
rz(-1.8453215) q[1];
sx q[1];
rz(2.5749899) q[1];
rz(-pi) q[2];
rz(-1.1282519) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(-0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4814066) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-1.0931724) q[2];
rz(0.53705755) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(-0.66948906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401232) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(1.4599266) q[0];
rz(1.5096674) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-3.0622838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58188841) q[0];
sx q[0];
rz(-1.4900582) q[0];
sx q[0];
rz(-1.333192) q[0];
rz(-pi) q[1];
rz(-2.5517188) q[2];
sx q[2];
rz(-2.8980719) q[2];
sx q[2];
rz(-2.4517454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39187688) q[1];
sx q[1];
rz(-1.8444711) q[1];
sx q[1];
rz(-0.43507149) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9856624) q[3];
sx q[3];
rz(-2.4178388) q[3];
sx q[3];
rz(-2.2035905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1377533) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(2.7286781) q[2];
rz(1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(0.41954654) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9846648) q[0];
sx q[0];
rz(-0.65021896) q[0];
sx q[0];
rz(-0.28537634) q[0];
rz(3.0889619) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(0.1836798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.015332) q[0];
sx q[0];
rz(-1.53526) q[0];
sx q[0];
rz(-1.4592341) q[0];
rz(-0.48641522) q[2];
sx q[2];
rz(-2.8210495) q[2];
sx q[2];
rz(0.17580168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7990098) q[1];
sx q[1];
rz(-1.4945806) q[1];
sx q[1];
rz(-0.59554312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37723549) q[3];
sx q[3];
rz(-2.4580015) q[3];
sx q[3];
rz(2.205276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8949184) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(-1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(-0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0378549) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(-1.6023585) q[0];
rz(2.1922951) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.4303713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968144) q[0];
sx q[0];
rz(-2.0955293) q[0];
sx q[0];
rz(2.120976) q[0];
rz(1.4865506) q[2];
sx q[2];
rz(-2.4110594) q[2];
sx q[2];
rz(-2.4075395) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9897223) q[1];
sx q[1];
rz(-2.8443326) q[1];
sx q[1];
rz(0.10135915) q[1];
x q[2];
rz(-2.8590081) q[3];
sx q[3];
rz(-1.9842098) q[3];
sx q[3];
rz(-1.9081209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.8370266) q[2];
sx q[2];
rz(-0.12651786) q[2];
rz(-0.33454076) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3009406) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(2.6244923) q[0];
rz(-0.05489796) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(0.089769207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9909279) q[0];
sx q[0];
rz(-2.1164843) q[0];
sx q[0];
rz(0.47713466) q[0];
rz(-pi) q[1];
rz(0.062762063) q[2];
sx q[2];
rz(-1.0107991) q[2];
sx q[2];
rz(-2.7095209) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5489784) q[1];
sx q[1];
rz(-1.4121118) q[1];
sx q[1];
rz(0.25513809) q[1];
rz(0.2281727) q[3];
sx q[3];
rz(-1.9252649) q[3];
sx q[3];
rz(1.4972403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(-2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.1995734) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(0.10689312) q[2];
sx q[2];
rz(-2.2724367) q[2];
sx q[2];
rz(2.5913749) q[2];
rz(0.94115067) q[3];
sx q[3];
rz(-1.0761906) q[3];
sx q[3];
rz(-3.0039345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
