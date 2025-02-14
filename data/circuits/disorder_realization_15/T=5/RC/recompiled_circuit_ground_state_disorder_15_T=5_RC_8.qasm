OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48159596) q[0];
sx q[0];
rz(-2.7240318) q[0];
sx q[0];
rz(0.041393809) q[0];
rz(0.59250915) q[1];
sx q[1];
rz(-0.23696391) q[1];
sx q[1];
rz(2.2466329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0314908) q[0];
sx q[0];
rz(-0.43927017) q[0];
sx q[0];
rz(1.2554368) q[0];
x q[1];
rz(2.824895) q[2];
sx q[2];
rz(-2.0513655) q[2];
sx q[2];
rz(2.4807535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8348011) q[1];
sx q[1];
rz(-0.39265206) q[1];
sx q[1];
rz(2.3638335) q[1];
rz(-pi) q[2];
rz(1.6144606) q[3];
sx q[3];
rz(-1.4619251) q[3];
sx q[3];
rz(3.0565232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71901739) q[2];
sx q[2];
rz(-0.97192478) q[2];
sx q[2];
rz(-2.1026588) q[2];
rz(-0.59764189) q[3];
sx q[3];
rz(-0.20331764) q[3];
sx q[3];
rz(1.982127) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0036302) q[0];
sx q[0];
rz(-2.2901386) q[0];
sx q[0];
rz(0.33388579) q[0];
rz(2.4489898) q[1];
sx q[1];
rz(-1.2665766) q[1];
sx q[1];
rz(-1.0612706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740299) q[0];
sx q[0];
rz(-2.2252796) q[0];
sx q[0];
rz(0.45463134) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7374836) q[2];
sx q[2];
rz(-1.8066908) q[2];
sx q[2];
rz(-1.0395694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0530653) q[1];
sx q[1];
rz(-2.663041) q[1];
sx q[1];
rz(0.38691945) q[1];
x q[2];
rz(-1.7739822) q[3];
sx q[3];
rz(-1.7082761) q[3];
sx q[3];
rz(0.41609496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4873203) q[2];
sx q[2];
rz(-1.4620179) q[2];
sx q[2];
rz(-2.9576603) q[2];
rz(-3.1390624) q[3];
sx q[3];
rz(-2.3380184) q[3];
sx q[3];
rz(2.7185503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84394395) q[0];
sx q[0];
rz(-0.14744814) q[0];
sx q[0];
rz(-2.3609128) q[0];
rz(0.88790226) q[1];
sx q[1];
rz(-2.0474515) q[1];
sx q[1];
rz(-2.323774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6827777) q[0];
sx q[0];
rz(-1.11894) q[0];
sx q[0];
rz(-2.8917612) q[0];
x q[1];
rz(0.79105241) q[2];
sx q[2];
rz(-2.1169691) q[2];
sx q[2];
rz(-2.5763022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43843514) q[1];
sx q[1];
rz(-0.58227986) q[1];
sx q[1];
rz(0.0011187355) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0637487) q[3];
sx q[3];
rz(-1.972282) q[3];
sx q[3];
rz(2.6428938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7094949) q[2];
sx q[2];
rz(-3.0165387) q[2];
sx q[2];
rz(-2.1091667) q[2];
rz(-0.89140511) q[3];
sx q[3];
rz(-0.67474198) q[3];
sx q[3];
rz(-0.83511043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29435232) q[0];
sx q[0];
rz(-0.27701858) q[0];
sx q[0];
rz(-0.54798049) q[0];
rz(-3.0849988) q[1];
sx q[1];
rz(-1.2013925) q[1];
sx q[1];
rz(3.1170381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5176651) q[0];
sx q[0];
rz(-2.6193287) q[0];
sx q[0];
rz(0.89512093) q[0];
rz(-pi) q[1];
rz(2.265931) q[2];
sx q[2];
rz(-0.49816874) q[2];
sx q[2];
rz(-2.9745575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3520842) q[1];
sx q[1];
rz(-0.86474907) q[1];
sx q[1];
rz(0.45070453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4931498) q[3];
sx q[3];
rz(-2.7973865) q[3];
sx q[3];
rz(1.5045117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5121439) q[2];
sx q[2];
rz(-2.5089743) q[2];
sx q[2];
rz(2.4150685) q[2];
rz(-1.4350545) q[3];
sx q[3];
rz(-1.8972998) q[3];
sx q[3];
rz(-1.4783036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5479945) q[0];
sx q[0];
rz(-1.3837805) q[0];
sx q[0];
rz(1.3872248) q[0];
rz(-0.45135003) q[1];
sx q[1];
rz(-0.74140048) q[1];
sx q[1];
rz(-0.50419921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614642) q[0];
sx q[0];
rz(-2.9315492) q[0];
sx q[0];
rz(-2.5145636) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4455657) q[2];
sx q[2];
rz(-0.88567222) q[2];
sx q[2];
rz(2.7643577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57297774) q[1];
sx q[1];
rz(-2.0842127) q[1];
sx q[1];
rz(-2.364679) q[1];
rz(-pi) q[2];
rz(0.80662722) q[3];
sx q[3];
rz(-0.99116814) q[3];
sx q[3];
rz(-1.5999853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7826358) q[2];
sx q[2];
rz(-2.2613596) q[2];
sx q[2];
rz(3.1202988) q[2];
rz(-3.0535789) q[3];
sx q[3];
rz(-2.8787677) q[3];
sx q[3];
rz(2.9737441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7511895) q[0];
sx q[0];
rz(-2.5683371) q[0];
sx q[0];
rz(-2.7943352) q[0];
rz(-1.685453) q[1];
sx q[1];
rz(-2.8442454) q[1];
sx q[1];
rz(0.28486326) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96831095) q[0];
sx q[0];
rz(-1.8704458) q[0];
sx q[0];
rz(-0.30072995) q[0];
rz(-3.0820644) q[2];
sx q[2];
rz(-0.7452508) q[2];
sx q[2];
rz(1.1203114) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59249944) q[1];
sx q[1];
rz(-0.91895267) q[1];
sx q[1];
rz(2.8054906) q[1];
rz(-pi) q[2];
rz(1.9273734) q[3];
sx q[3];
rz(-1.8468401) q[3];
sx q[3];
rz(0.86892933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.138729) q[2];
sx q[2];
rz(-1.8314654) q[2];
sx q[2];
rz(-1.0339453) q[2];
rz(0.73505861) q[3];
sx q[3];
rz(-1.6481954) q[3];
sx q[3];
rz(2.471931) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0064148684) q[0];
sx q[0];
rz(-2.4786351) q[0];
sx q[0];
rz(-0.906115) q[0];
rz(-2.9754029) q[1];
sx q[1];
rz(-2.6527185) q[1];
sx q[1];
rz(0.30950549) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342732) q[0];
sx q[0];
rz(-1.7953402) q[0];
sx q[0];
rz(-0.17404795) q[0];
x q[1];
rz(1.3679753) q[2];
sx q[2];
rz(-2.3536567) q[2];
sx q[2];
rz(2.1757464) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5823203) q[1];
sx q[1];
rz(-0.67495912) q[1];
sx q[1];
rz(-1.2260502) q[1];
rz(-1.0528492) q[3];
sx q[3];
rz(-0.098909698) q[3];
sx q[3];
rz(2.3581751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6918148) q[2];
sx q[2];
rz(-0.66902995) q[2];
sx q[2];
rz(1.6548033) q[2];
rz(2.6438223) q[3];
sx q[3];
rz(-0.83511746) q[3];
sx q[3];
rz(2.914079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787702) q[0];
sx q[0];
rz(-2.4202122) q[0];
sx q[0];
rz(-0.26807868) q[0];
rz(-2.2396741) q[1];
sx q[1];
rz(-2.0733158) q[1];
sx q[1];
rz(1.3021775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807515) q[0];
sx q[0];
rz(-2.5227418) q[0];
sx q[0];
rz(0.28736823) q[0];
rz(0.44659932) q[2];
sx q[2];
rz(-1.5726358) q[2];
sx q[2];
rz(0.50825495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5561192) q[1];
sx q[1];
rz(-2.1550309) q[1];
sx q[1];
rz(-2.9727436) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6792278) q[3];
sx q[3];
rz(-2.83395) q[3];
sx q[3];
rz(0.82266146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34182081) q[2];
sx q[2];
rz(-1.2398961) q[2];
sx q[2];
rz(2.5566901) q[2];
rz(1.2029485) q[3];
sx q[3];
rz(-1.7479618) q[3];
sx q[3];
rz(-1.5929619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840266) q[0];
sx q[0];
rz(-2.6507222) q[0];
sx q[0];
rz(0.67114818) q[0];
rz(2.8533543) q[1];
sx q[1];
rz(-2.3858374) q[1];
sx q[1];
rz(-0.63405687) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8347522) q[0];
sx q[0];
rz(-0.04148395) q[0];
sx q[0];
rz(1.1783429) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8618635) q[2];
sx q[2];
rz(-0.17269793) q[2];
sx q[2];
rz(-0.23898838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61810869) q[1];
sx q[1];
rz(-1.5963703) q[1];
sx q[1];
rz(2.8167538) q[1];
rz(-pi) q[2];
rz(-2.8170812) q[3];
sx q[3];
rz(-2.4430118) q[3];
sx q[3];
rz(-1.6968264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59342283) q[2];
sx q[2];
rz(-1.5735184) q[2];
sx q[2];
rz(-2.8578952) q[2];
rz(0.13188322) q[3];
sx q[3];
rz(-2.7851084) q[3];
sx q[3];
rz(-0.62294817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9054966) q[0];
sx q[0];
rz(-0.14362366) q[0];
sx q[0];
rz(-0.96963257) q[0];
rz(-2.6456004) q[1];
sx q[1];
rz(-0.6876567) q[1];
sx q[1];
rz(-2.2775441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0522033) q[0];
sx q[0];
rz(-1.5004352) q[0];
sx q[0];
rz(-0.86300993) q[0];
rz(-pi) q[1];
rz(-2.0798565) q[2];
sx q[2];
rz(-2.6531918) q[2];
sx q[2];
rz(1.0906206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8567791) q[1];
sx q[1];
rz(-0.80409324) q[1];
sx q[1];
rz(2.2914679) q[1];
x q[2];
rz(0.91409036) q[3];
sx q[3];
rz(-1.2063081) q[3];
sx q[3];
rz(1.9250223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25829092) q[2];
sx q[2];
rz(-0.51670462) q[2];
sx q[2];
rz(-2.1041763) q[2];
rz(-2.9344905) q[3];
sx q[3];
rz(-0.89829069) q[3];
sx q[3];
rz(-0.5894388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265738) q[0];
sx q[0];
rz(-1.6828231) q[0];
sx q[0];
rz(-1.9376301) q[0];
rz(0.013982458) q[1];
sx q[1];
rz(-1.436186) q[1];
sx q[1];
rz(2.0367429) q[1];
rz(-0.034910708) q[2];
sx q[2];
rz(-1.8014805) q[2];
sx q[2];
rz(-1.0911566) q[2];
rz(-1.1963853) q[3];
sx q[3];
rz(-0.76210124) q[3];
sx q[3];
rz(-2.6740554) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
