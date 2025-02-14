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
rz(0.1641195) q[0];
sx q[0];
rz(4.1657148) q[0];
sx q[0];
rz(9.9088718) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9021716) q[0];
sx q[0];
rz(-2.4706984) q[0];
sx q[0];
rz(-1.5205199) q[0];
x q[1];
rz(1.820716) q[2];
sx q[2];
rz(-1.5139765) q[2];
sx q[2];
rz(2.394258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.37147) q[1];
sx q[1];
rz(-2.9767163) q[1];
sx q[1];
rz(0.6368963) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1160158) q[3];
sx q[3];
rz(-0.6729799) q[3];
sx q[3];
rz(-2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.310828) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(2.8669299) q[2];
rz(0.87485391) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(0.39875317) q[0];
rz(-0.2440456) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(1.1057378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4717446) q[0];
sx q[0];
rz(-0.5037743) q[0];
sx q[0];
rz(0.48724799) q[0];
rz(-2.2244338) q[2];
sx q[2];
rz(-2.0431402) q[2];
sx q[2];
rz(1.2231959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.058814) q[1];
sx q[1];
rz(-1.5734488) q[1];
sx q[1];
rz(-3.1383731) q[1];
x q[2];
rz(1.8375754) q[3];
sx q[3];
rz(-1.6872129) q[3];
sx q[3];
rz(-0.67177012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70046052) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(2.5221141) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19720896) q[0];
sx q[0];
rz(-1.5329307) q[0];
sx q[0];
rz(0.029649563) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4476851) q[2];
sx q[2];
rz(-1.6173714) q[2];
sx q[2];
rz(1.2706626) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6782376) q[1];
sx q[1];
rz(-1.872552) q[1];
sx q[1];
rz(-2.9185118) q[1];
rz(2.1501174) q[3];
sx q[3];
rz(-2.1001171) q[3];
sx q[3];
rz(0.73562276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8850024) q[2];
sx q[2];
rz(-1.8424748) q[2];
sx q[2];
rz(-3.0510862) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(-2.0335782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3312382) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(2.2099387) q[0];
rz(-2.9725507) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(2.1021252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1202034) q[0];
sx q[0];
rz(-1.8112371) q[0];
sx q[0];
rz(2.405015) q[0];
x q[1];
rz(-2.9157964) q[2];
sx q[2];
rz(-1.8183299) q[2];
sx q[2];
rz(-3.053726) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4770866) q[1];
sx q[1];
rz(-1.5110104) q[1];
sx q[1];
rz(-1.7782446) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5040197) q[3];
sx q[3];
rz(-0.3488144) q[3];
sx q[3];
rz(-1.82774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52183759) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-0.49631611) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(2.1697178) q[0];
rz(-3.019849) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(1.3294539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7065763) q[0];
sx q[0];
rz(-1.6321275) q[0];
sx q[0];
rz(1.2641843) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2586063) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(2.4508053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59898224) q[1];
sx q[1];
rz(-0.81331454) q[1];
sx q[1];
rz(-2.5573362) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0092486898) q[3];
sx q[3];
rz(-1.4917177) q[3];
sx q[3];
rz(-1.6759863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10151265) q[2];
sx q[2];
rz(-1.9209346) q[2];
sx q[2];
rz(2.9902048) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720035) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(1.0714916) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(-0.62613097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75979489) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(-1.2977029) q[0];
rz(-0.59711908) q[2];
sx q[2];
rz(-1.3665939) q[2];
sx q[2];
rz(1.1115505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31589139) q[1];
sx q[1];
rz(-2.4448312) q[1];
sx q[1];
rz(2.6543762) q[1];
rz(0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(0.56624352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(2.1273071) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(-2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(-3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(1.1309518) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3944954) q[0];
sx q[0];
rz(-2.1685026) q[0];
sx q[0];
rz(-1.5776442) q[0];
x q[1];
rz(2.5860687) q[2];
sx q[2];
rz(-1.0706524) q[2];
sx q[2];
rz(1.9490567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9662094) q[1];
sx q[1];
rz(-1.4384603) q[1];
sx q[1];
rz(-0.90465178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6252296) q[3];
sx q[3];
rz(-1.5011906) q[3];
sx q[3];
rz(-2.7854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(-0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(-2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(2.193702) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(1.9416169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24058293) q[0];
sx q[0];
rz(-1.173061) q[0];
sx q[0];
rz(2.564604) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84917111) q[2];
sx q[2];
rz(-2.2361922) q[2];
sx q[2];
rz(-2.9862491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7941798) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.7296687) q[1];
rz(-pi) q[2];
rz(3.0407716) q[3];
sx q[3];
rz(-1.3885048) q[3];
sx q[3];
rz(-1.3829447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.348749) q[2];
rz(1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(-2.0674904) q[0];
rz(2.7117924) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(-0.46357402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(3.0449165) q[0];
rz(-pi) q[1];
rz(-0.65166574) q[2];
sx q[2];
rz(-2.8164688) q[2];
sx q[2];
rz(0.59949694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1164828) q[1];
sx q[1];
rz(-1.1949202) q[1];
sx q[1];
rz(-0.89565887) q[1];
rz(-0.14318569) q[3];
sx q[3];
rz(-2.5783263) q[3];
sx q[3];
rz(1.3384502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27434906) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(2.3331433) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(-3.1032108) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(-0.29676944) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0050186) q[0];
sx q[0];
rz(-1.2941735) q[0];
sx q[0];
rz(1.0199976) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.692524) q[2];
sx q[2];
rz(-1.9757604) q[2];
sx q[2];
rz(-1.8835889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3546875) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(-2.812318) q[1];
rz(-pi) q[2];
rz(2.2567883) q[3];
sx q[3];
rz(-1.2594885) q[3];
sx q[3];
rz(1.6303568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9639637) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(2.5753944) q[2];
rz(1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.1894777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(-2.879907) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-1.5952806) q[2];
sx q[2];
rz(-1.4888121) q[2];
sx q[2];
rz(-0.76406995) q[2];
rz(0.64601267) q[3];
sx q[3];
rz(-2.1727242) q[3];
sx q[3];
rz(0.77632191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
