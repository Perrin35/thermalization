OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(4.6392415) q[0];
sx q[0];
rz(8.5949329) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(2.2652664) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5397545) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(1.3215617) q[0];
rz(-2.1562188) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(-1.0759682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4727002) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(3.1256691) q[1];
x q[2];
rz(-1.7238486) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(-3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(2.412964) q[2];
rz(0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(0.87444011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53423184) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(-0.94632728) q[0];
rz(-1.8057683) q[2];
sx q[2];
rz(-0.58832303) q[2];
sx q[2];
rz(0.36662835) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.892131) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-2.4888121) q[1];
rz(1.3784834) q[3];
sx q[3];
rz(-0.27413878) q[3];
sx q[3];
rz(-2.6499234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4008537) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-0.40288231) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(2.7456465) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1186737) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(2.2454717) q[0];
x q[1];
rz(-0.01404889) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(-0.77169466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.695895) q[1];
sx q[1];
rz(-0.88765111) q[1];
sx q[1];
rz(1.5584857) q[1];
x q[2];
rz(0.25554244) q[3];
sx q[3];
rz(-1.5481755) q[3];
sx q[3];
rz(-0.25657755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-0.23342361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6596286) q[0];
sx q[0];
rz(-0.11419645) q[0];
sx q[0];
rz(2.1301079) q[0];
rz(-pi) q[1];
rz(-1.3350305) q[2];
sx q[2];
rz(-0.94199099) q[2];
sx q[2];
rz(-2.6218888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(0.48308785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78419533) q[3];
sx q[3];
rz(-1.8521706) q[3];
sx q[3];
rz(-0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9984263) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(1.6869998) q[0];
x q[1];
rz(0.20247395) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(1.7970049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4796175) q[1];
sx q[1];
rz(-1.8912589) q[1];
sx q[1];
rz(1.5451876) q[1];
rz(1.8499591) q[3];
sx q[3];
rz(-1.2984635) q[3];
sx q[3];
rz(-2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(3.1150505) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(2.81566) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1074166) q[0];
sx q[0];
rz(-1.0506442) q[0];
sx q[0];
rz(-2.0053894) q[0];
rz(-2.6961156) q[2];
sx q[2];
rz(-1.2947086) q[2];
sx q[2];
rz(-2.7518227) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(-2.8403776) q[1];
rz(-1.5876706) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(2.8942787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5439593) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-0.11925764) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79524604) q[0];
sx q[0];
rz(-0.7659142) q[0];
sx q[0];
rz(0.38608293) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97738816) q[2];
sx q[2];
rz(-2.9466629) q[2];
sx q[2];
rz(-2.5755663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5057999) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(2.990681) q[1];
rz(0.70763208) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(1.5054782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5078686) q[0];
sx q[0];
rz(-1.3390203) q[0];
sx q[0];
rz(-0.50110441) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9912668) q[2];
sx q[2];
rz(-1.5577003) q[2];
sx q[2];
rz(3.1377369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.209621) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(-0.45100905) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2742911) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(-1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.2148946) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58763114) q[0];
sx q[0];
rz(-0.56092867) q[0];
sx q[0];
rz(1.8857303) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43478195) q[2];
sx q[2];
rz(-1.5909305) q[2];
sx q[2];
rz(0.15205631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62577265) q[1];
sx q[1];
rz(-1.9966218) q[1];
sx q[1];
rz(1.6243402) q[1];
rz(-pi) q[2];
rz(-0.92000658) q[3];
sx q[3];
rz(-1.917106) q[3];
sx q[3];
rz(-1.8370093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(2.2686968) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.1970253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.409317) q[0];
sx q[0];
rz(-2.0079552) q[0];
sx q[0];
rz(0.93426312) q[0];
x q[1];
rz(-2.0335774) q[2];
sx q[2];
rz(-0.26460755) q[2];
sx q[2];
rz(0.77361425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5572284) q[1];
sx q[1];
rz(-0.96712501) q[1];
sx q[1];
rz(-0.61324688) q[1];
rz(-pi) q[2];
rz(-1.9479806) q[3];
sx q[3];
rz(-1.0918655) q[3];
sx q[3];
rz(-2.9951028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(0.52195436) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-1.21571) q[2];
sx q[2];
rz(-1.2181031) q[2];
sx q[2];
rz(1.9545771) q[2];
rz(0.79896169) q[3];
sx q[3];
rz(-2.3286455) q[3];
sx q[3];
rz(0.13959985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
