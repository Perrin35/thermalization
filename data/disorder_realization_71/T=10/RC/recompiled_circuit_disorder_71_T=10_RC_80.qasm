OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(2.2652664) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1286368) q[0];
sx q[0];
rz(-2.6814333) q[0];
sx q[0];
rz(0.53928661) q[0];
x q[1];
rz(-2.1562188) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(2.0656245) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6688924) q[1];
sx q[1];
rz(-1.2408222) q[1];
sx q[1];
rz(3.1256691) q[1];
rz(0.34378864) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(0.49376282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(-0.87444011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16337285) q[0];
sx q[0];
rz(-1.1436497) q[0];
sx q[0];
rz(-2.8073729) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9875056) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(2.4947583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.16972152) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(2.4596618) q[1];
x q[2];
rz(-3.087895) q[3];
sx q[3];
rz(-1.8397545) q[3];
sx q[3];
rz(-0.69124903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(-0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(1.0748192) q[0];
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
rz(-3.031035) q[0];
sx q[0];
rz(-1.1764515) q[0];
sx q[0];
rz(1.0236543) q[0];
rz(-pi) q[1];
rz(0.01404889) q[2];
sx q[2];
rz(-1.8648246) q[2];
sx q[2];
rz(-0.77169466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6763941) q[1];
sx q[1];
rz(-2.4583543) q[1];
sx q[1];
rz(-3.1264683) q[1];
x q[2];
rz(0.25554244) q[3];
sx q[3];
rz(-1.5481755) q[3];
sx q[3];
rz(2.8850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(0.23920693) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(-0.23342361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596286) q[0];
sx q[0];
rz(-0.11419645) q[0];
sx q[0];
rz(1.0114848) q[0];
rz(-2.4993863) q[2];
sx q[2];
rz(-1.7608479) q[2];
sx q[2];
rz(1.9501291) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1366795) q[1];
sx q[1];
rz(-0.79489743) q[1];
sx q[1];
rz(2.6585048) q[1];
rz(0.78419533) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.6915879) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(0.064237021) q[0];
rz(-2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8781389) q[0];
sx q[0];
rz(-1.3100776) q[0];
sx q[0];
rz(-3.1104452) q[0];
x q[1];
rz(-2.9391187) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(1.7970049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(3.0645963) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77869271) q[3];
sx q[3];
rz(-2.7540996) q[3];
sx q[3];
rz(3.1117698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(0.09207329) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(0.87310711) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(2.81566) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.35615) q[0];
sx q[0];
rz(-2.4768562) q[0];
sx q[0];
rz(-2.507693) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6961156) q[2];
sx q[2];
rz(-1.2947086) q[2];
sx q[2];
rz(0.38976994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(2.8403776) q[1];
x q[2];
rz(-1.553922) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(-0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
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
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-0.72189271) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-3.022335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0810453) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(-0.72780769) q[0];
x q[1];
rz(-1.4085521) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(-2.721399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.99265656) q[1];
sx q[1];
rz(-1.438237) q[1];
sx q[1];
rz(2.0723144) q[1];
rz(0.42640949) q[3];
sx q[3];
rz(-2.3873513) q[3];
sx q[3];
rz(-0.38503669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44935903) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(2.2559821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53997707) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(-2.6849296) q[0];
rz(-3.1259414) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(1.5755115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9319716) q[1];
sx q[1];
rz(-2.0721657) q[1];
sx q[1];
rz(-0.45100905) q[1];
rz(2.4212491) q[3];
sx q[3];
rz(-0.90984905) q[3];
sx q[3];
rz(-1.1870445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(1.2148946) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-2.0057604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275883) q[0];
sx q[0];
rz(-1.7363318) q[0];
sx q[0];
rz(1.0323314) q[0];
rz(-0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(1.4620632) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75479924) q[1];
sx q[1];
rz(-2.7126185) q[1];
sx q[1];
rz(3.0241443) q[1];
x q[2];
rz(-0.42580749) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-0.18994722) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68162936) q[0];
sx q[0];
rz(-0.75461331) q[0];
sx q[0];
rz(0.90453903) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3329266) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(-2.7931917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5843643) q[1];
sx q[1];
rz(-2.1744676) q[1];
sx q[1];
rz(0.61324688) q[1];
rz(0.50935575) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(-1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(-2.1879788) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-0.52195436) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-1.9258826) q[2];
sx q[2];
rz(-1.9234895) q[2];
sx q[2];
rz(-1.1870155) q[2];
rz(-2.342631) q[3];
sx q[3];
rz(-2.3286455) q[3];
sx q[3];
rz(0.13959985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
