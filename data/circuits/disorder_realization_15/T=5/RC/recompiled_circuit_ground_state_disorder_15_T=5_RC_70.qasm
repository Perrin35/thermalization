OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6599967) q[0];
sx q[0];
rz(-0.41756088) q[0];
sx q[0];
rz(-0.041393809) q[0];
rz(-2.5490835) q[1];
sx q[1];
rz(-2.9046287) q[1];
sx q[1];
rz(-2.2466329) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3774619) q[0];
sx q[0];
rz(-1.9870166) q[0];
sx q[0];
rz(-0.14472117) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.069015) q[2];
sx q[2];
rz(-1.291005) q[2];
sx q[2];
rz(-0.75959709) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1393361) q[1];
sx q[1];
rz(-1.8426241) q[1];
sx q[1];
rz(0.28693954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37994639) q[3];
sx q[3];
rz(-0.11726876) q[3];
sx q[3];
rz(-0.46739551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71901739) q[2];
sx q[2];
rz(-2.1696679) q[2];
sx q[2];
rz(-1.0389339) q[2];
rz(-0.59764189) q[3];
sx q[3];
rz(-2.938275) q[3];
sx q[3];
rz(-1.982127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.0036302) q[0];
sx q[0];
rz(-2.2901386) q[0];
sx q[0];
rz(0.33388579) q[0];
rz(2.4489898) q[1];
sx q[1];
rz(-1.2665766) q[1];
sx q[1];
rz(2.0803221) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4505517) q[0];
sx q[0];
rz(-0.77734935) q[0];
sx q[0];
rz(-1.0509899) q[0];
x q[1];
rz(-2.5928934) q[2];
sx q[2];
rz(-0.46462545) q[2];
sx q[2];
rz(1.0313874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1352219) q[1];
sx q[1];
rz(-1.7454444) q[1];
sx q[1];
rz(-2.6937301) q[1];
rz(-pi) q[2];
rz(-1.3676104) q[3];
sx q[3];
rz(-1.4333165) q[3];
sx q[3];
rz(-2.7254977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4873203) q[2];
sx q[2];
rz(-1.6795748) q[2];
sx q[2];
rz(0.1839323) q[2];
rz(-3.1390624) q[3];
sx q[3];
rz(-2.3380184) q[3];
sx q[3];
rz(-0.42304236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2976487) q[0];
sx q[0];
rz(-0.14744814) q[0];
sx q[0];
rz(-0.78067988) q[0];
rz(0.88790226) q[1];
sx q[1];
rz(-1.0941411) q[1];
sx q[1];
rz(-0.81781864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22293689) q[0];
sx q[0];
rz(-1.7950995) q[0];
sx q[0];
rz(-1.1064311) q[0];
x q[1];
rz(-2.3505402) q[2];
sx q[2];
rz(-1.0246236) q[2];
sx q[2];
rz(-0.56529048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7031575) q[1];
sx q[1];
rz(-0.58227986) q[1];
sx q[1];
rz(-3.1404739) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6894916) q[3];
sx q[3];
rz(-2.0341877) q[3];
sx q[3];
rz(-2.2832561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7094949) q[2];
sx q[2];
rz(-0.12505394) q[2];
sx q[2];
rz(-2.1091667) q[2];
rz(2.2501875) q[3];
sx q[3];
rz(-0.67474198) q[3];
sx q[3];
rz(-0.83511043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.29435232) q[0];
sx q[0];
rz(-0.27701858) q[0];
sx q[0];
rz(2.5936122) q[0];
rz(-0.056593865) q[1];
sx q[1];
rz(-1.2013925) q[1];
sx q[1];
rz(0.024554575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62392759) q[0];
sx q[0];
rz(-0.52226394) q[0];
sx q[0];
rz(0.89512093) q[0];
rz(0.33523021) q[2];
sx q[2];
rz(-1.1950699) q[2];
sx q[2];
rz(0.59240985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6645421) q[1];
sx q[1];
rz(-1.9086491) q[1];
sx q[1];
rz(-0.8124688) q[1];
x q[2];
rz(-1.6484429) q[3];
sx q[3];
rz(-0.34420612) q[3];
sx q[3];
rz(1.6370809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5121439) q[2];
sx q[2];
rz(-2.5089743) q[2];
sx q[2];
rz(2.4150685) q[2];
rz(-1.7065382) q[3];
sx q[3];
rz(-1.2442929) q[3];
sx q[3];
rz(1.6632891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5479945) q[0];
sx q[0];
rz(-1.3837805) q[0];
sx q[0];
rz(1.3872248) q[0];
rz(-2.6902426) q[1];
sx q[1];
rz(-2.4001922) q[1];
sx q[1];
rz(-0.50419921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.90550051) q[2];
sx q[2];
rz(-0.93449392) q[2];
sx q[2];
rz(1.8417412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6804769) q[1];
sx q[1];
rz(-2.2409984) q[1];
sx q[1];
rz(0.67733353) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81309115) q[3];
sx q[3];
rz(-2.219354) q[3];
sx q[3];
rz(0.5483707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7826358) q[2];
sx q[2];
rz(-2.2613596) q[2];
sx q[2];
rz(3.1202988) q[2];
rz(3.0535789) q[3];
sx q[3];
rz(-2.8787677) q[3];
sx q[3];
rz(-2.9737441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.39040318) q[0];
sx q[0];
rz(-2.5683371) q[0];
sx q[0];
rz(-0.34725749) q[0];
rz(1.685453) q[1];
sx q[1];
rz(-0.29734722) q[1];
sx q[1];
rz(0.28486326) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96831095) q[0];
sx q[0];
rz(-1.2711469) q[0];
sx q[0];
rz(2.8408627) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6256394) q[2];
sx q[2];
rz(-2.3144139) q[2];
sx q[2];
rz(1.0393927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5490932) q[1];
sx q[1];
rz(-0.91895267) q[1];
sx q[1];
rz(0.33610208) q[1];
rz(1.2142193) q[3];
sx q[3];
rz(-1.8468401) q[3];
sx q[3];
rz(2.2726633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0028637) q[2];
sx q[2];
rz(-1.3101273) q[2];
sx q[2];
rz(-1.0339453) q[2];
rz(-0.73505861) q[3];
sx q[3];
rz(-1.6481954) q[3];
sx q[3];
rz(-2.471931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1351778) q[0];
sx q[0];
rz(-0.66295755) q[0];
sx q[0];
rz(0.906115) q[0];
rz(2.9754029) q[1];
sx q[1];
rz(-2.6527185) q[1];
sx q[1];
rz(-0.30950549) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85649895) q[0];
sx q[0];
rz(-1.740432) q[0];
sx q[0];
rz(-1.3429252) q[0];
x q[1];
rz(1.3679753) q[2];
sx q[2];
rz(-2.3536567) q[2];
sx q[2];
rz(2.1757464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.013405) q[1];
sx q[1];
rz(-2.1995029) q[1];
sx q[1];
rz(2.8774269) q[1];
rz(-0.0490908) q[3];
sx q[3];
rz(-1.6566983) q[3];
sx q[3];
rz(1.3034749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6918148) q[2];
sx q[2];
rz(-0.66902995) q[2];
sx q[2];
rz(1.4867894) q[2];
rz(-2.6438223) q[3];
sx q[3];
rz(-0.83511746) q[3];
sx q[3];
rz(0.2275137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.9437156) q[0];
sx q[0];
rz(-2.4202122) q[0];
sx q[0];
rz(2.873514) q[0];
rz(-2.2396741) q[1];
sx q[1];
rz(-2.0733158) q[1];
sx q[1];
rz(1.3021775) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9807515) q[0];
sx q[0];
rz(-0.61885082) q[0];
sx q[0];
rz(2.8542244) q[0];
rz(0.44659932) q[2];
sx q[2];
rz(-1.5726358) q[2];
sx q[2];
rz(-2.6333377) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2563602) q[1];
sx q[1];
rz(-2.5361885) q[1];
sx q[1];
rz(1.3219236) q[1];
rz(-1.2648495) q[3];
sx q[3];
rz(-1.6035723) q[3];
sx q[3];
rz(-2.2900801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7997718) q[2];
sx q[2];
rz(-1.2398961) q[2];
sx q[2];
rz(2.5566901) q[2];
rz(1.9386442) q[3];
sx q[3];
rz(-1.3936309) q[3];
sx q[3];
rz(-1.5929619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575661) q[0];
sx q[0];
rz(-2.6507222) q[0];
sx q[0];
rz(2.4704445) q[0];
rz(-0.28823832) q[1];
sx q[1];
rz(-2.3858374) q[1];
sx q[1];
rz(2.5075358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3068405) q[0];
sx q[0];
rz(-0.04148395) q[0];
sx q[0];
rz(-1.1783429) q[0];
rz(-1.2797292) q[2];
sx q[2];
rz(-2.9688947) q[2];
sx q[2];
rz(0.23898838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1975164) q[1];
sx q[1];
rz(-1.895525) q[1];
sx q[1];
rz(1.5977809) q[1];
rz(0.32451144) q[3];
sx q[3];
rz(-0.69858089) q[3];
sx q[3];
rz(-1.4447663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59342283) q[2];
sx q[2];
rz(-1.5680743) q[2];
sx q[2];
rz(2.8578952) q[2];
rz(3.0097094) q[3];
sx q[3];
rz(-2.7851084) q[3];
sx q[3];
rz(-2.5186445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236096) q[0];
sx q[0];
rz(-0.14362366) q[0];
sx q[0];
rz(-2.1719601) q[0];
rz(0.49599221) q[1];
sx q[1];
rz(-2.453936) q[1];
sx q[1];
rz(-0.86404854) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5629146) q[0];
sx q[0];
rz(-2.2764674) q[0];
sx q[0];
rz(-3.0490997) q[0];
x q[1];
rz(-2.0052018) q[2];
sx q[2];
rz(-1.8015125) q[2];
sx q[2];
rz(2.2034881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2848136) q[1];
sx q[1];
rz(-2.3374994) q[1];
sx q[1];
rz(-2.2914679) q[1];
rz(-pi) q[2];
x q[2];
rz(0.448927) q[3];
sx q[3];
rz(-2.1777946) q[3];
sx q[3];
rz(-0.62238979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8833017) q[2];
sx q[2];
rz(-0.51670462) q[2];
sx q[2];
rz(-1.0374163) q[2];
rz(-0.20710219) q[3];
sx q[3];
rz(-2.243302) q[3];
sx q[3];
rz(-0.5894388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.0265738) q[0];
sx q[0];
rz(-1.4587695) q[0];
sx q[0];
rz(1.2039626) q[0];
rz(-0.013982458) q[1];
sx q[1];
rz(-1.7054066) q[1];
sx q[1];
rz(-1.1048497) q[1];
rz(-0.034910708) q[2];
sx q[2];
rz(-1.8014805) q[2];
sx q[2];
rz(-1.0911566) q[2];
rz(1.1963853) q[3];
sx q[3];
rz(-2.3794914) q[3];
sx q[3];
rz(0.46753721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
