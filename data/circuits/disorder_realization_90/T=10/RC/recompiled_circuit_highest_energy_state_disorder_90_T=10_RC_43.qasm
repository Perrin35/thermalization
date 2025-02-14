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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(2.1460549) q[0];
rz(-2.0350463) q[1];
sx q[1];
rz(-2.3732329) q[1];
sx q[1];
rz(0.33222693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7357644) q[0];
sx q[0];
rz(-1.9737184) q[0];
sx q[0];
rz(0.71056239) q[0];
x q[1];
rz(-1.840749) q[2];
sx q[2];
rz(-0.47856646) q[2];
sx q[2];
rz(-1.1978483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.886604) q[1];
sx q[1];
rz(-1.1082778) q[1];
sx q[1];
rz(-0.47420331) q[1];
x q[2];
rz(-2.291392) q[3];
sx q[3];
rz(-2.4293141) q[3];
sx q[3];
rz(-0.18092052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60653162) q[2];
sx q[2];
rz(-1.0142832) q[2];
sx q[2];
rz(0.052058546) q[2];
rz(1.082513) q[3];
sx q[3];
rz(-2.1942997) q[3];
sx q[3];
rz(-1.989971) q[3];
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
rz(-1.624619) q[0];
sx q[0];
rz(-0.032051429) q[0];
sx q[0];
rz(-2.8042941) q[0];
rz(2.9889122) q[1];
sx q[1];
rz(-2.3744507) q[1];
sx q[1];
rz(1.5229567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168287) q[0];
sx q[0];
rz(-1.7375907) q[0];
sx q[0];
rz(2.5186064) q[0];
rz(-pi) q[1];
rz(-0.27045336) q[2];
sx q[2];
rz(-0.65545299) q[2];
sx q[2];
rz(1.7846817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2238521) q[1];
sx q[1];
rz(-1.187993) q[1];
sx q[1];
rz(-0.21767016) q[1];
x q[2];
rz(0.60373276) q[3];
sx q[3];
rz(-1.4430178) q[3];
sx q[3];
rz(1.6953682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1107669) q[2];
sx q[2];
rz(-2.775707) q[2];
sx q[2];
rz(0.3863253) q[2];
rz(1.2836766) q[3];
sx q[3];
rz(-1.9648896) q[3];
sx q[3];
rz(-1.9234689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.06685) q[0];
sx q[0];
rz(-1.0800986) q[0];
sx q[0];
rz(2.1225488) q[0];
rz(-0.77541238) q[1];
sx q[1];
rz(-1.2762028) q[1];
sx q[1];
rz(2.6554241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4751525) q[0];
sx q[0];
rz(-1.5678582) q[0];
sx q[0];
rz(-0.03133987) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48392754) q[2];
sx q[2];
rz(-0.84887767) q[2];
sx q[2];
rz(-0.98361042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3348766) q[1];
sx q[1];
rz(-1.3522321) q[1];
sx q[1];
rz(-0.49592917) q[1];
rz(2.8193238) q[3];
sx q[3];
rz(-1.8506088) q[3];
sx q[3];
rz(-1.1424949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39573085) q[2];
sx q[2];
rz(-1.1660601) q[2];
sx q[2];
rz(-0.80500785) q[2];
rz(-0.65195596) q[3];
sx q[3];
rz(-1.1019573) q[3];
sx q[3];
rz(-0.93418795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66251278) q[0];
sx q[0];
rz(-1.8204239) q[0];
sx q[0];
rz(1.3385734) q[0];
rz(-1.592912) q[1];
sx q[1];
rz(-2.4376696) q[1];
sx q[1];
rz(2.7159363) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073624728) q[0];
sx q[0];
rz(-1.5789766) q[0];
sx q[0];
rz(-0.091853022) q[0];
rz(0.034599853) q[2];
sx q[2];
rz(-1.8784754) q[2];
sx q[2];
rz(-0.75242413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4257926) q[1];
sx q[1];
rz(-1.1669901) q[1];
sx q[1];
rz(2.3753662) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2888186) q[3];
sx q[3];
rz(-0.87347841) q[3];
sx q[3];
rz(3.092749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8694596) q[2];
sx q[2];
rz(-1.3674066) q[2];
sx q[2];
rz(-2.7222705) q[2];
rz(1.9648633) q[3];
sx q[3];
rz(-0.71507016) q[3];
sx q[3];
rz(2.5756605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6856573) q[0];
sx q[0];
rz(-2.3929907) q[0];
sx q[0];
rz(-1.6393882) q[0];
rz(-1.4337076) q[1];
sx q[1];
rz(-1.8610443) q[1];
sx q[1];
rz(0.40275231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0549598) q[0];
sx q[0];
rz(-1.827359) q[0];
sx q[0];
rz(2.199928) q[0];
rz(-pi) q[1];
rz(1.6300768) q[2];
sx q[2];
rz(-1.0067954) q[2];
sx q[2];
rz(0.21280542) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13632475) q[1];
sx q[1];
rz(-1.439716) q[1];
sx q[1];
rz(-1.6682427) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6879657) q[3];
sx q[3];
rz(-0.90725431) q[3];
sx q[3];
rz(-2.508916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6413642) q[2];
sx q[2];
rz(-1.0470752) q[2];
sx q[2];
rz(0.25263986) q[2];
rz(-0.80444515) q[3];
sx q[3];
rz(-2.6201456) q[3];
sx q[3];
rz(-0.49853244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.3280403) q[0];
sx q[0];
rz(-0.10843065) q[0];
sx q[0];
rz(0.88388467) q[0];
rz(-0.49993316) q[1];
sx q[1];
rz(-1.4686613) q[1];
sx q[1];
rz(-2.6271741) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0262951) q[0];
sx q[0];
rz(-2.4567488) q[0];
sx q[0];
rz(2.1158873) q[0];
rz(-pi) q[1];
rz(-0.40174257) q[2];
sx q[2];
rz(-2.5726924) q[2];
sx q[2];
rz(0.57050975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7307694) q[1];
sx q[1];
rz(-1.9716773) q[1];
sx q[1];
rz(-2.904179) q[1];
x q[2];
rz(3.1153999) q[3];
sx q[3];
rz(-2.3520654) q[3];
sx q[3];
rz(-1.7408371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9875235) q[2];
sx q[2];
rz(-1.46393) q[2];
sx q[2];
rz(-0.12518159) q[2];
rz(-0.82031885) q[3];
sx q[3];
rz(-1.9671974) q[3];
sx q[3];
rz(-2.7952747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5527375) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(-3.0297025) q[0];
rz(2.0018068) q[1];
sx q[1];
rz(-2.9264989) q[1];
sx q[1];
rz(1.8591759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49408052) q[0];
sx q[0];
rz(-1.2984797) q[0];
sx q[0];
rz(2.168505) q[0];
x q[1];
rz(-2.5865104) q[2];
sx q[2];
rz(-2.8709445) q[2];
sx q[2];
rz(-0.78449434) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0874028) q[1];
sx q[1];
rz(-0.8381745) q[1];
sx q[1];
rz(2.0778632) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5130643) q[3];
sx q[3];
rz(-2.2542076) q[3];
sx q[3];
rz(-2.6635955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45184267) q[2];
sx q[2];
rz(-2.8359783) q[2];
sx q[2];
rz(1.3172733) q[2];
rz(-0.02056038) q[3];
sx q[3];
rz(-0.33187425) q[3];
sx q[3];
rz(-2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9356215) q[0];
sx q[0];
rz(-2.1501849) q[0];
sx q[0];
rz(2.8427065) q[0];
rz(1.0609974) q[1];
sx q[1];
rz(-1.3875049) q[1];
sx q[1];
rz(-1.2334197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3372634) q[0];
sx q[0];
rz(-1.4739887) q[0];
sx q[0];
rz(1.9053354) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.100343) q[2];
sx q[2];
rz(-1.4564351) q[2];
sx q[2];
rz(-3.0806266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1202482) q[1];
sx q[1];
rz(-1.2750191) q[1];
sx q[1];
rz(2.2463754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66090195) q[3];
sx q[3];
rz(-1.9993625) q[3];
sx q[3];
rz(1.0180343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3491106) q[2];
sx q[2];
rz(-1.0976378) q[2];
sx q[2];
rz(0.20723542) q[2];
rz(-2.4810897) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(-1.2787261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.58885634) q[0];
sx q[0];
rz(-1.4016466) q[0];
sx q[0];
rz(1.2218342) q[0];
rz(-2.6264722) q[1];
sx q[1];
rz(-0.11038596) q[1];
sx q[1];
rz(2.5882904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56263262) q[0];
sx q[0];
rz(-1.4427358) q[0];
sx q[0];
rz(-1.0109896) q[0];
x q[1];
rz(1.7153928) q[2];
sx q[2];
rz(-0.42486496) q[2];
sx q[2];
rz(-1.0185084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8749145) q[1];
sx q[1];
rz(-1.5648702) q[1];
sx q[1];
rz(-1.3937217) q[1];
x q[2];
rz(3.1114486) q[3];
sx q[3];
rz(-0.6693535) q[3];
sx q[3];
rz(-3.1140259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8052266) q[2];
sx q[2];
rz(-1.1104501) q[2];
sx q[2];
rz(0.52555788) q[2];
rz(-1.4793652) q[3];
sx q[3];
rz(-0.95366228) q[3];
sx q[3];
rz(0.76013887) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568327) q[0];
sx q[0];
rz(-0.79462093) q[0];
sx q[0];
rz(2.7123465) q[0];
rz(1.5224573) q[1];
sx q[1];
rz(-1.2702962) q[1];
sx q[1];
rz(-2.2116908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4711205) q[0];
sx q[0];
rz(-1.3567302) q[0];
sx q[0];
rz(-1.9164969) q[0];
x q[1];
rz(1.3512813) q[2];
sx q[2];
rz(-1.4633312) q[2];
sx q[2];
rz(-1.4023413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86204862) q[1];
sx q[1];
rz(-0.32717184) q[1];
sx q[1];
rz(-1.0485093) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0449379) q[3];
sx q[3];
rz(-0.55516637) q[3];
sx q[3];
rz(-0.47278178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4055736) q[2];
sx q[2];
rz(-2.6655727) q[2];
sx q[2];
rz(2.948577) q[2];
rz(-0.080282601) q[3];
sx q[3];
rz(-2.2716227) q[3];
sx q[3];
rz(1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3748462) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(-0.97024067) q[1];
sx q[1];
rz(-1.2087676) q[1];
sx q[1];
rz(-1.2420775) q[1];
rz(2.8497981) q[2];
sx q[2];
rz(-0.84635432) q[2];
sx q[2];
rz(1.9683471) q[2];
rz(0.24701444) q[3];
sx q[3];
rz(-0.09709662) q[3];
sx q[3];
rz(0.63963565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
