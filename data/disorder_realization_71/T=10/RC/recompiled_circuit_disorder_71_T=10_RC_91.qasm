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
rz(-1.6439438) q[0];
sx q[0];
rz(-0.82984501) q[0];
rz(-2.3614376) q[1];
sx q[1];
rz(-1.0649788) q[1];
sx q[1];
rz(0.87632626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(2.602306) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19199065) q[2];
sx q[2];
rz(-2.1477094) q[2];
sx q[2];
rz(0.38919762) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0383366) q[1];
sx q[1];
rz(-1.5858608) q[1];
sx q[1];
rz(-1.2407833) q[1];
x q[2];
rz(-0.34378864) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(2.6478298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(-0.87444011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53423184) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(0.94632728) q[0];
rz(-1.8057683) q[2];
sx q[2];
rz(-0.58832303) q[2];
sx q[2];
rz(0.36662835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0793822) q[1];
sx q[1];
rz(-0.96696889) q[1];
sx q[1];
rz(1.0152597) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3784834) q[3];
sx q[3];
rz(-0.27413878) q[3];
sx q[3];
rz(0.4916693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-1.0748192) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(2.7456465) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9112644) q[0];
sx q[0];
rz(-1.0697782) q[0];
sx q[0];
rz(-0.45341861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6171574) q[2];
sx q[2];
rz(-0.29435396) q[2];
sx q[2];
rz(-2.4183395) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.695895) q[1];
sx q[1];
rz(-0.88765111) q[1];
sx q[1];
rz(-1.5584857) q[1];
rz(2.8860502) q[3];
sx q[3];
rz(-1.5934172) q[3];
sx q[3];
rz(2.8850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-0.81992942) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(-0.23342361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080290681) q[0];
sx q[0];
rz(-1.6675321) q[0];
sx q[0];
rz(-0.060782766) q[0];
rz(0.64220631) q[2];
sx q[2];
rz(-1.7608479) q[2];
sx q[2];
rz(-1.1914636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0049131752) q[1];
sx q[1];
rz(-0.79489743) q[1];
sx q[1];
rz(2.6585048) q[1];
rz(-2.3573973) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(-0.22339544) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8262186) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(1.3099567) q[0];
rz(2.9391187) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(1.3445878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.099247301) q[1];
sx q[1];
rz(-1.595101) q[1];
sx q[1];
rz(-2.821032) q[1];
rz(-0.77869271) q[3];
sx q[3];
rz(-0.38749309) q[3];
sx q[3];
rz(-3.1117698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(0.09207329) q[2];
rz(-0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(0.87310711) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(-0.32593265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317141) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(-0.563234) q[0];
rz(-pi) q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.9982669) q[2];
sx q[2];
rz(-1.3104591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6335878) q[1];
sx q[1];
rz(-1.8691917) q[1];
sx q[1];
rz(1.7119346) q[1];
rz(-pi) q[2];
rz(1.553922) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-0.11925764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8598547) q[0];
sx q[0];
rz(-2.2681232) q[0];
sx q[0];
rz(-1.2233234) q[0];
x q[1];
rz(0.10995933) q[2];
sx q[2];
rz(-1.409515) q[2];
sx q[2];
rz(1.1683299) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99265656) q[1];
sx q[1];
rz(-1.7033556) q[1];
sx q[1];
rz(1.0692783) q[1];
x q[2];
rz(1.9414385) q[3];
sx q[3];
rz(-2.2439085) q[3];
sx q[3];
rz(-0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(2.3826777) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83157241) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6016156) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(0.45666306) q[0];
rz(3.1259414) q[2];
sx q[2];
rz(-2.1502697) q[2];
sx q[2];
rz(1.5755115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9319716) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(0.45100905) q[1];
rz(2.3732244) q[3];
sx q[3];
rz(-2.1185015) q[3];
sx q[3];
rz(2.2636569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-1.1317066) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866859) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-0.19219877) q[0];
x q[1];
rz(0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(1.6795295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3867934) q[1];
sx q[1];
rz(-0.42897412) q[1];
sx q[1];
rz(-0.11744833) q[1];
rz(-pi) q[2];
rz(-0.92000658) q[3];
sx q[3];
rz(-1.2244867) q[3];
sx q[3];
rz(-1.3045834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.9445673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.409317) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(2.2073295) q[0];
rz(-pi) q[1];
rz(-0.12038259) q[2];
sx q[2];
rz(-1.3345846) q[2];
sx q[2];
rz(-1.2506968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5843643) q[1];
sx q[1];
rz(-0.96712501) q[1];
sx q[1];
rz(-0.61324688) q[1];
rz(-pi) q[2];
rz(-2.5246546) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(-0.85655772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-2.7675046) q[2];
sx q[2];
rz(-1.903152) q[2];
sx q[2];
rz(-2.8852035) q[2];
rz(-0.63511499) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
