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
rz(-2.5858606) q[0];
sx q[0];
rz(-1.2795804) q[0];
sx q[0];
rz(0.32787856) q[0];
rz(-2.9887587) q[1];
sx q[1];
rz(-2.6522377) q[1];
sx q[1];
rz(-1.0110923) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3837102) q[0];
sx q[0];
rz(-2.5284736) q[0];
sx q[0];
rz(-0.28045537) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7462872) q[2];
sx q[2];
rz(-2.7031941) q[2];
sx q[2];
rz(-0.47847462) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0786405) q[1];
sx q[1];
rz(-0.83940369) q[1];
sx q[1];
rz(2.6432829) q[1];
rz(-pi) q[2];
rz(1.9888844) q[3];
sx q[3];
rz(-1.7888513) q[3];
sx q[3];
rz(2.5421028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8622387) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(-1.5834825) q[2];
rz(2.8067348) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(0.28086942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25664499) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(0.33194342) q[0];
rz(0.360082) q[1];
sx q[1];
rz(-1.304345) q[1];
sx q[1];
rz(-2.8538381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5482564) q[0];
sx q[0];
rz(-1.3704482) q[0];
sx q[0];
rz(0.25357539) q[0];
rz(-pi) q[1];
rz(-1.8910725) q[2];
sx q[2];
rz(-0.8647635) q[2];
sx q[2];
rz(-1.6635739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2079434) q[1];
sx q[1];
rz(-0.48743409) q[1];
sx q[1];
rz(-2.3323564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7340363) q[3];
sx q[3];
rz(-1.3426011) q[3];
sx q[3];
rz(1.850224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(-1.7506556) q[2];
rz(-0.85919291) q[3];
sx q[3];
rz(-1.8682559) q[3];
sx q[3];
rz(2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.745382) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(-0.23042738) q[0];
rz(0.51741171) q[1];
sx q[1];
rz(-1.0942065) q[1];
sx q[1];
rz(-2.8527625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6073501) q[0];
sx q[0];
rz(-1.5114044) q[0];
sx q[0];
rz(-1.9797822) q[0];
x q[1];
rz(1.6452698) q[2];
sx q[2];
rz(-1.1033325) q[2];
sx q[2];
rz(-0.68231586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7728665) q[1];
sx q[1];
rz(-1.1613559) q[1];
sx q[1];
rz(0.057842908) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90694859) q[3];
sx q[3];
rz(-0.28093279) q[3];
sx q[3];
rz(-2.8317766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0383703) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(3.1008516) q[2];
rz(0.1013969) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(-2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(2.907584) q[0];
rz(0.19800828) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3458684) q[0];
sx q[0];
rz(-1.0875174) q[0];
sx q[0];
rz(0.048075284) q[0];
rz(-pi) q[1];
rz(-2.6279672) q[2];
sx q[2];
rz(-1.6503929) q[2];
sx q[2];
rz(-2.3964376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9069104) q[1];
sx q[1];
rz(-0.89212067) q[1];
sx q[1];
rz(-0.80295103) q[1];
rz(-pi) q[2];
rz(-1.4928905) q[3];
sx q[3];
rz(-2.7512105) q[3];
sx q[3];
rz(1.7254894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733751) q[2];
sx q[2];
rz(-0.28926352) q[2];
sx q[2];
rz(2.3667228) q[2];
rz(1.3261999) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.73959094) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(3.1091029) q[0];
rz(1.2292817) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(-0.083018735) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1047528) q[0];
sx q[0];
rz(-2.6304465) q[0];
sx q[0];
rz(2.6698861) q[0];
rz(-pi) q[1];
rz(1.7402224) q[2];
sx q[2];
rz(-1.8352785) q[2];
sx q[2];
rz(-1.6319909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28369812) q[1];
sx q[1];
rz(-2.2304728) q[1];
sx q[1];
rz(-0.96925756) q[1];
rz(-pi) q[2];
rz(1.5918077) q[3];
sx q[3];
rz(-2.5046231) q[3];
sx q[3];
rz(1.6629499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7897537) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(2.9808673) q[2];
rz(-1.9488581) q[3];
sx q[3];
rz(-2.9294117) q[3];
sx q[3];
rz(2.3737657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6650894) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(2.8506668) q[0];
rz(1.7581958) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(0.49984041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2384278) q[0];
sx q[0];
rz(-3.0808093) q[0];
sx q[0];
rz(-2.3709036) q[0];
rz(-pi) q[1];
rz(2.1792214) q[2];
sx q[2];
rz(-2.314724) q[2];
sx q[2];
rz(2.8679071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0840869) q[1];
sx q[1];
rz(-2.018714) q[1];
sx q[1];
rz(2.491076) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96792241) q[3];
sx q[3];
rz(-2.335066) q[3];
sx q[3];
rz(1.4462808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9342039) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(0.70154166) q[2];
rz(-0.97405854) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(2.2402703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3153673) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(1.9024128) q[1];
sx q[1];
rz(-1.0737123) q[1];
sx q[1];
rz(2.0674131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095431134) q[0];
sx q[0];
rz(-3.1218596) q[0];
sx q[0];
rz(2.31953) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4137109) q[2];
sx q[2];
rz(-0.54833503) q[2];
sx q[2];
rz(0.20073433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3201808) q[1];
sx q[1];
rz(-1.3548676) q[1];
sx q[1];
rz(-1.126299) q[1];
rz(-2.3775565) q[3];
sx q[3];
rz(-1.2437399) q[3];
sx q[3];
rz(0.59248176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(2.5568753) q[2];
rz(-1.0662474) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(-0.40759459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(-2.6608652) q[0];
rz(1.2113384) q[1];
sx q[1];
rz(-1.5296661) q[1];
sx q[1];
rz(0.617625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0636198) q[0];
sx q[0];
rz(-1.5665861) q[0];
sx q[0];
rz(1.2456889) q[0];
x q[1];
rz(-2.4489538) q[2];
sx q[2];
rz(-0.73441457) q[2];
sx q[2];
rz(1.278233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.86194) q[1];
sx q[1];
rz(-1.2766925) q[1];
sx q[1];
rz(-2.9169464) q[1];
rz(-0.17980735) q[3];
sx q[3];
rz(-0.19528611) q[3];
sx q[3];
rz(-2.2168468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9628613) q[2];
sx q[2];
rz(-1.7469254) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(1.0293055) q[3];
sx q[3];
rz(-1.7512713) q[3];
sx q[3];
rz(1.5503957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7903098) q[0];
sx q[0];
rz(-0.93515486) q[0];
sx q[0];
rz(-1.9679605) q[0];
rz(-1.1890746) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(0.48114166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8011337) q[0];
sx q[0];
rz(-1.0164355) q[0];
sx q[0];
rz(1.8340684) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6716154) q[2];
sx q[2];
rz(-1.013213) q[2];
sx q[2];
rz(-2.4876311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6509241) q[1];
sx q[1];
rz(-1.7177525) q[1];
sx q[1];
rz(-0.96687324) q[1];
x q[2];
rz(-2.1012382) q[3];
sx q[3];
rz(-1.181349) q[3];
sx q[3];
rz(0.77270179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22566191) q[2];
sx q[2];
rz(-1.2292726) q[2];
sx q[2];
rz(-1.4813598) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.1888844) q[3];
sx q[3];
rz(-1.9143298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34761053) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(-3.0718497) q[0];
rz(-0.28930411) q[1];
sx q[1];
rz(-1.2846839) q[1];
sx q[1];
rz(-1.0522254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5890516) q[0];
sx q[0];
rz(-0.49363401) q[0];
sx q[0];
rz(0.73666739) q[0];
rz(-pi) q[1];
rz(1.7830816) q[2];
sx q[2];
rz(-0.67321482) q[2];
sx q[2];
rz(-0.5897943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6695646) q[1];
sx q[1];
rz(-1.9306364) q[1];
sx q[1];
rz(1.1675952) q[1];
rz(-0.77731384) q[3];
sx q[3];
rz(-2.5717426) q[3];
sx q[3];
rz(-0.22138813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5436486) q[2];
sx q[2];
rz(-1.1755377) q[2];
sx q[2];
rz(3.0808466) q[2];
rz(-0.28901035) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(-1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48988265) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(-1.0501077) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(-0.15090987) q[2];
sx q[2];
rz(-2.0357292) q[2];
sx q[2];
rz(-0.72307088) q[2];
rz(-0.28239863) q[3];
sx q[3];
rz(-1.1544607) q[3];
sx q[3];
rz(-0.37731597) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
