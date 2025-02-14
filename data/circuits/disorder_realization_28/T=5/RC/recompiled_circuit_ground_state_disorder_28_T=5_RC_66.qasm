OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(1.6190785) q[0];
sx q[0];
rz(9.4879307) q[0];
rz(1.8237279) q[1];
sx q[1];
rz(-0.39931077) q[1];
sx q[1];
rz(2.3656486) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467593) q[0];
sx q[0];
rz(-1.5251119) q[0];
sx q[0];
rz(-1.5976357) q[0];
rz(0.49709289) q[2];
sx q[2];
rz(-2.1462206) q[2];
sx q[2];
rz(2.6445553) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5835442) q[1];
sx q[1];
rz(-2.2698102) q[1];
sx q[1];
rz(-0.2322766) q[1];
rz(-0.50252359) q[3];
sx q[3];
rz(-1.3602358) q[3];
sx q[3];
rz(-0.37744194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6905288) q[2];
sx q[2];
rz(-1.3961926) q[2];
sx q[2];
rz(-0.53831354) q[2];
rz(0.29317835) q[3];
sx q[3];
rz(-1.5978975) q[3];
sx q[3];
rz(1.258491) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341768) q[0];
sx q[0];
rz(-2.2947312) q[0];
sx q[0];
rz(-2.9127981) q[0];
rz(-0.070411988) q[1];
sx q[1];
rz(-1.2733302) q[1];
sx q[1];
rz(-1.061903) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5492229) q[0];
sx q[0];
rz(-1.3064991) q[0];
sx q[0];
rz(0.53472526) q[0];
x q[1];
rz(-2.6648971) q[2];
sx q[2];
rz(-1.9856493) q[2];
sx q[2];
rz(-1.7384993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8613597) q[1];
sx q[1];
rz(-1.8953084) q[1];
sx q[1];
rz(1.2636416) q[1];
rz(-pi) q[2];
rz(-2.7135647) q[3];
sx q[3];
rz(-1.3770388) q[3];
sx q[3];
rz(2.0462115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9262907) q[2];
sx q[2];
rz(-2.6965202) q[2];
sx q[2];
rz(0.11623795) q[2];
rz(-2.2262573) q[3];
sx q[3];
rz(-1.6719619) q[3];
sx q[3];
rz(0.62469971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.4794469) q[0];
sx q[0];
rz(-1.5818469) q[0];
sx q[0];
rz(-0.33552718) q[0];
rz(1.4115931) q[1];
sx q[1];
rz(-2.2970707) q[1];
sx q[1];
rz(-1.1988877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035164) q[0];
sx q[0];
rz(-1.8545574) q[0];
sx q[0];
rz(1.0182583) q[0];
rz(0.54036083) q[2];
sx q[2];
rz(-0.66695222) q[2];
sx q[2];
rz(1.4948542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80497031) q[1];
sx q[1];
rz(-1.687007) q[1];
sx q[1];
rz(2.3776965) q[1];
rz(-0.24925225) q[3];
sx q[3];
rz(-2.0468759) q[3];
sx q[3];
rz(-2.2152546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18232839) q[2];
sx q[2];
rz(-2.3458643) q[2];
sx q[2];
rz(-1.1129334) q[2];
rz(0.5111323) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(-2.6326877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6540601) q[0];
sx q[0];
rz(-0.075981058) q[0];
sx q[0];
rz(0.37387601) q[0];
rz(1.182425) q[1];
sx q[1];
rz(-1.1610718) q[1];
sx q[1];
rz(2.9808796) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99237011) q[0];
sx q[0];
rz(-2.9806031) q[0];
sx q[0];
rz(0.54246728) q[0];
x q[1];
rz(-0.66443759) q[2];
sx q[2];
rz(-2.2291846) q[2];
sx q[2];
rz(0.13641549) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.070878167) q[1];
sx q[1];
rz(-1.1859815) q[1];
sx q[1];
rz(1.0139731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8237958) q[3];
sx q[3];
rz(-1.8565912) q[3];
sx q[3];
rz(-0.60493776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60761991) q[2];
sx q[2];
rz(-2.0325568) q[2];
sx q[2];
rz(-0.71471659) q[2];
rz(-2.886582) q[3];
sx q[3];
rz(-1.8111572) q[3];
sx q[3];
rz(2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.052499972) q[0];
sx q[0];
rz(-2.8369501) q[0];
sx q[0];
rz(-2.3783045) q[0];
rz(2.8870562) q[1];
sx q[1];
rz(-1.7691879) q[1];
sx q[1];
rz(-1.6932142) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453706) q[0];
sx q[0];
rz(-1.9217446) q[0];
sx q[0];
rz(0.36349067) q[0];
rz(-pi) q[1];
rz(1.7380452) q[2];
sx q[2];
rz(-2.3838701) q[2];
sx q[2];
rz(-3.0191985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6519439) q[1];
sx q[1];
rz(-1.7484416) q[1];
sx q[1];
rz(0.62477115) q[1];
rz(-0.45381399) q[3];
sx q[3];
rz(-2.3771913) q[3];
sx q[3];
rz(-1.3168471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.280507) q[2];
sx q[2];
rz(-1.0887159) q[2];
sx q[2];
rz(2.1121934) q[2];
rz(-1.6086802) q[3];
sx q[3];
rz(-0.89919388) q[3];
sx q[3];
rz(-2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3935881) q[0];
sx q[0];
rz(-2.2279976) q[0];
sx q[0];
rz(-0.21011259) q[0];
rz(0.70136079) q[1];
sx q[1];
rz(-1.3664093) q[1];
sx q[1];
rz(2.0393541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85066807) q[0];
sx q[0];
rz(-1.0449261) q[0];
sx q[0];
rz(-0.53914244) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1386078) q[2];
sx q[2];
rz(-0.66221234) q[2];
sx q[2];
rz(-2.4110297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.452207) q[1];
sx q[1];
rz(-1.8836792) q[1];
sx q[1];
rz(2.41928) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54164782) q[3];
sx q[3];
rz(-1.4337501) q[3];
sx q[3];
rz(2.8769719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8474951) q[2];
sx q[2];
rz(-1.5832486) q[2];
sx q[2];
rz(1.6451277) q[2];
rz(-1.3930813) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(0.67162544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3944655) q[0];
sx q[0];
rz(-1.8165996) q[0];
sx q[0];
rz(-0.4766683) q[0];
rz(-1.7851625) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(2.2043601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8097157) q[0];
sx q[0];
rz(-0.88167718) q[0];
sx q[0];
rz(3.0948113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5086447) q[2];
sx q[2];
rz(-2.2878433) q[2];
sx q[2];
rz(0.99822068) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37717078) q[1];
sx q[1];
rz(-0.57767361) q[1];
sx q[1];
rz(0.75862268) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20003518) q[3];
sx q[3];
rz(-0.45645255) q[3];
sx q[3];
rz(2.4394622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6095907) q[2];
sx q[2];
rz(-1.8014182) q[2];
sx q[2];
rz(2.0659921) q[2];
rz(-2.7514451) q[3];
sx q[3];
rz(-1.3107927) q[3];
sx q[3];
rz(2.3384317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6049062) q[0];
sx q[0];
rz(-1.0746047) q[0];
sx q[0];
rz(-2.5965776) q[0];
rz(-2.1489428) q[1];
sx q[1];
rz(-0.98572171) q[1];
sx q[1];
rz(-1.4535905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3869868) q[0];
sx q[0];
rz(-2.0029066) q[0];
sx q[0];
rz(-2.291392) q[0];
x q[1];
rz(-1.4382866) q[2];
sx q[2];
rz(-1.6269267) q[2];
sx q[2];
rz(-2.0126337) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0664898) q[1];
sx q[1];
rz(-2.3716092) q[1];
sx q[1];
rz(-1.408073) q[1];
x q[2];
rz(2.4798315) q[3];
sx q[3];
rz(-3.0654717) q[3];
sx q[3];
rz(-0.95835987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40889007) q[2];
sx q[2];
rz(-1.7444892) q[2];
sx q[2];
rz(-2.8453541) q[2];
rz(0.57684165) q[3];
sx q[3];
rz(-0.39671612) q[3];
sx q[3];
rz(1.2805987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532362) q[0];
sx q[0];
rz(-0.78859538) q[0];
sx q[0];
rz(2.9175135) q[0];
rz(-0.60399404) q[1];
sx q[1];
rz(-2.150034) q[1];
sx q[1];
rz(1.4858861) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1224269) q[0];
sx q[0];
rz(-2.6976524) q[0];
sx q[0];
rz(-2.5371518) q[0];
rz(-pi) q[1];
rz(1.3732128) q[2];
sx q[2];
rz(-1.7181529) q[2];
sx q[2];
rz(1.7201529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0476802) q[1];
sx q[1];
rz(-0.96598066) q[1];
sx q[1];
rz(-1.259786) q[1];
rz(-pi) q[2];
rz(2.6194686) q[3];
sx q[3];
rz(-2.587954) q[3];
sx q[3];
rz(-1.6185135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6622582) q[2];
sx q[2];
rz(-1.0666288) q[2];
sx q[2];
rz(0.38702854) q[2];
rz(2.8582063) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(2.7391105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0249483) q[0];
sx q[0];
rz(-2.8706757) q[0];
sx q[0];
rz(-2.0181632) q[0];
rz(-1.6472389) q[1];
sx q[1];
rz(-1.4861743) q[1];
sx q[1];
rz(2.2656238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7520053) q[0];
sx q[0];
rz(-1.8266062) q[0];
sx q[0];
rz(2.385456) q[0];
x q[1];
rz(0.45166679) q[2];
sx q[2];
rz(-2.0658138) q[2];
sx q[2];
rz(1.5507914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3800602) q[1];
sx q[1];
rz(-2.426154) q[1];
sx q[1];
rz(1.9081294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3920636) q[3];
sx q[3];
rz(-1.9944085) q[3];
sx q[3];
rz(2.7152747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7312077) q[2];
sx q[2];
rz(-1.7816252) q[2];
sx q[2];
rz(3.0044921) q[2];
rz(-1.7588663) q[3];
sx q[3];
rz(-2.2189249) q[3];
sx q[3];
rz(-1.9130798) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41505861) q[0];
sx q[0];
rz(-2.6714323) q[0];
sx q[0];
rz(2.5191125) q[0];
rz(2.7525735) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(2.0034267) q[2];
sx q[2];
rz(-2.0180704) q[2];
sx q[2];
rz(1.13712) q[2];
rz(-0.20122726) q[3];
sx q[3];
rz(-1.22898) q[3];
sx q[3];
rz(-2.1635319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
