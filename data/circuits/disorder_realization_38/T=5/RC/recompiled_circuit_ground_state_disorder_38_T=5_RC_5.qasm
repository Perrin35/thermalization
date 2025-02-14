OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91650668) q[0];
sx q[0];
rz(-2.9986311) q[0];
sx q[0];
rz(1.4885055) q[0];
rz(1.1324963) q[1];
sx q[1];
rz(-2.7096665) q[1];
sx q[1];
rz(2.5052524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778438) q[0];
sx q[0];
rz(-1.1882458) q[0];
sx q[0];
rz(-2.2027459) q[0];
rz(-pi) q[1];
rz(1.9540399) q[2];
sx q[2];
rz(-1.9343209) q[2];
sx q[2];
rz(-2.5853047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.591037) q[1];
sx q[1];
rz(-1.6948957) q[1];
sx q[1];
rz(1.1731338) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8543872) q[3];
sx q[3];
rz(-0.99528377) q[3];
sx q[3];
rz(-2.7137386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0296313) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(0.65043989) q[2];
rz(-1.316831) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0541075) q[0];
sx q[0];
rz(-0.58687812) q[0];
sx q[0];
rz(0.45463872) q[0];
rz(-3.1184323) q[1];
sx q[1];
rz(-1.6975479) q[1];
sx q[1];
rz(0.32618162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9232193) q[0];
sx q[0];
rz(-2.120864) q[0];
sx q[0];
rz(-2.0327507) q[0];
x q[1];
rz(-1.147338) q[2];
sx q[2];
rz(-1.2988161) q[2];
sx q[2];
rz(-0.47951298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4224976) q[1];
sx q[1];
rz(-1.0651673) q[1];
sx q[1];
rz(1.0984332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2315027) q[3];
sx q[3];
rz(-1.4779363) q[3];
sx q[3];
rz(1.6805436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11639103) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(2.1873059) q[2];
rz(-1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8702451) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(1.9648319) q[0];
rz(-0.36508834) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(1.6129859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58816649) q[0];
sx q[0];
rz(-1.5840826) q[0];
sx q[0];
rz(-1.5499328) q[0];
rz(-pi) q[1];
rz(-1.1272548) q[2];
sx q[2];
rz(-1.4993641) q[2];
sx q[2];
rz(0.88692666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6512201) q[1];
sx q[1];
rz(-1.4191322) q[1];
sx q[1];
rz(-3.1107799) q[1];
x q[2];
rz(-3.0729483) q[3];
sx q[3];
rz(-1.513283) q[3];
sx q[3];
rz(2.7781093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8927346) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(-1.9888606) q[2];
rz(-1.24498) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(0.28321701) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98274851) q[0];
sx q[0];
rz(-pi/12) q[0];
sx q[0];
rz(2.4520279) q[0];
rz(-0.025029643) q[1];
sx q[1];
rz(-1.5182779) q[1];
sx q[1];
rz(1.4208581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11224225) q[0];
sx q[0];
rz(-2.1922702) q[0];
sx q[0];
rz(1.6981324) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17949149) q[2];
sx q[2];
rz(-1.1467666) q[2];
sx q[2];
rz(1.9551203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1240439) q[1];
sx q[1];
rz(-2.6376403) q[1];
sx q[1];
rz(-0.057226463) q[1];
rz(-1.8426651) q[3];
sx q[3];
rz(-1.2449578) q[3];
sx q[3];
rz(0.0001903521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.932852) q[2];
sx q[2];
rz(-2.9661861) q[2];
sx q[2];
rz(-1.1669195) q[2];
rz(0.44421998) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(2.8319632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8532448) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(-2.7254768) q[0];
rz(-2.6314645) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(-2.3663734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213971) q[0];
sx q[0];
rz(-0.96548432) q[0];
sx q[0];
rz(1.3277675) q[0];
x q[1];
rz(1.3890319) q[2];
sx q[2];
rz(-1.9700288) q[2];
sx q[2];
rz(0.71208524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0803804) q[1];
sx q[1];
rz(-2.9728372) q[1];
sx q[1];
rz(-0.83909281) q[1];
x q[2];
rz(0.52948496) q[3];
sx q[3];
rz(-2.4765402) q[3];
sx q[3];
rz(2.3163026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9267209) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(2.1057687) q[2];
rz(-0.18946798) q[3];
sx q[3];
rz(-0.47179705) q[3];
sx q[3];
rz(2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8871317) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(-0.3558085) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-2.1864086) q[1];
sx q[1];
rz(-2.0004418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1352859) q[0];
sx q[0];
rz(-1.3606439) q[0];
sx q[0];
rz(-2.2457613) q[0];
rz(-pi) q[1];
rz(-3.1325794) q[2];
sx q[2];
rz(-0.95925943) q[2];
sx q[2];
rz(2.1076941) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.43737632) q[1];
sx q[1];
rz(-0.17568888) q[1];
sx q[1];
rz(-2.499627) q[1];
rz(-pi) q[2];
rz(0.012142556) q[3];
sx q[3];
rz(-2.0046765) q[3];
sx q[3];
rz(-1.7858684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63153875) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(1.4368524) q[2];
rz(2.0884183) q[3];
sx q[3];
rz(-1.6097693) q[3];
sx q[3];
rz(0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.961504) q[0];
sx q[0];
rz(-0.29619521) q[0];
sx q[0];
rz(0.12251138) q[0];
rz(-2.4854614) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(-1.8001385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264474) q[0];
sx q[0];
rz(-1.6065803) q[0];
sx q[0];
rz(-3.0962178) q[0];
rz(-pi) q[1];
rz(0.28320988) q[2];
sx q[2];
rz(-1.6746175) q[2];
sx q[2];
rz(1.4701172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16640284) q[1];
sx q[1];
rz(-0.67968785) q[1];
sx q[1];
rz(1.1389334) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44388598) q[3];
sx q[3];
rz(-1.4680913) q[3];
sx q[3];
rz(-0.041180276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5682257) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(1.4637671) q[2];
rz(0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(1.8370139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952591) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(2.6829868) q[0];
rz(2.1268225) q[1];
sx q[1];
rz(-1.1943123) q[1];
sx q[1];
rz(2.0789304) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39358562) q[0];
sx q[0];
rz(-0.17608041) q[0];
sx q[0];
rz(0.3349456) q[0];
rz(0.52824654) q[2];
sx q[2];
rz(-0.88111231) q[2];
sx q[2];
rz(1.4572342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96982376) q[1];
sx q[1];
rz(-0.14482982) q[1];
sx q[1];
rz(1.1317731) q[1];
rz(-pi) q[2];
rz(-2.0902781) q[3];
sx q[3];
rz(-1.2881345) q[3];
sx q[3];
rz(1.9611028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.7532316) q[2];
sx q[2];
rz(-1.3512705) q[2];
rz(1.9225559) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(2.3082699) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(0.90743995) q[0];
rz(2.9145248) q[1];
sx q[1];
rz(-2.0957969) q[1];
sx q[1];
rz(2.2423832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.718841) q[0];
sx q[0];
rz(-2.1678574) q[0];
sx q[0];
rz(-2.5121252) q[0];
x q[1];
rz(-0.81041281) q[2];
sx q[2];
rz(-3.028026) q[2];
sx q[2];
rz(1.0134987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6688706) q[1];
sx q[1];
rz(-2.5963915) q[1];
sx q[1];
rz(-1.9342058) q[1];
x q[2];
rz(1.6257004) q[3];
sx q[3];
rz(-2.3927352) q[3];
sx q[3];
rz(2.1532173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9565309) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(2.5409307) q[2];
rz(1.0968084) q[3];
sx q[3];
rz(-2.5513702) q[3];
sx q[3];
rz(0.87305951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3908865) q[0];
sx q[0];
rz(-1.0989256) q[0];
sx q[0];
rz(2.1571958) q[0];
rz(0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.9459928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988783) q[0];
sx q[0];
rz(-0.98473583) q[0];
sx q[0];
rz(1.1436966) q[0];
rz(-pi) q[1];
rz(2.6457328) q[2];
sx q[2];
rz(-0.95520077) q[2];
sx q[2];
rz(-0.95814182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25567935) q[1];
sx q[1];
rz(-1.1885841) q[1];
sx q[1];
rz(0.087319386) q[1];
rz(-pi) q[2];
rz(-1.1117842) q[3];
sx q[3];
rz(-2.2892366) q[3];
sx q[3];
rz(1.3415003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0857346) q[2];
sx q[2];
rz(-1.6207638) q[2];
sx q[2];
rz(0.98304191) q[2];
rz(0.25162697) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(-2.1380641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428103) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(1.8819173) q[1];
sx q[1];
rz(-1.7301662) q[1];
sx q[1];
rz(0.87283254) q[1];
rz(-2.2839584) q[2];
sx q[2];
rz(-1.4076172) q[2];
sx q[2];
rz(2.313688) q[2];
rz(2.3771277) q[3];
sx q[3];
rz(-2.6938312) q[3];
sx q[3];
rz(-2.5934283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
