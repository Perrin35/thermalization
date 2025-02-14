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
rz(-1.7718908) q[0];
sx q[0];
rz(-2.7112024) q[0];
sx q[0];
rz(0.16464591) q[0];
rz(-1.6904263) q[1];
sx q[1];
rz(-0.37835205) q[1];
sx q[1];
rz(-2.4317256) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043308) q[0];
sx q[0];
rz(-1.834615) q[0];
sx q[0];
rz(3.0065048) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9989292) q[2];
sx q[2];
rz(-1.1096508) q[2];
sx q[2];
rz(-2.3335339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30395092) q[1];
sx q[1];
rz(-1.5216989) q[1];
sx q[1];
rz(-1.9242679) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23386896) q[3];
sx q[3];
rz(-1.4483671) q[3];
sx q[3];
rz(-0.32952884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6115438) q[2];
sx q[2];
rz(-2.7908466) q[2];
sx q[2];
rz(-1.9225527) q[2];
rz(-0.48398316) q[3];
sx q[3];
rz(-0.84570208) q[3];
sx q[3];
rz(-1.774196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71500635) q[0];
sx q[0];
rz(-2.8026717) q[0];
sx q[0];
rz(1.4434848) q[0];
rz(1.8679856) q[1];
sx q[1];
rz(-2.1111635) q[1];
sx q[1];
rz(-0.66676203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0979157) q[0];
sx q[0];
rz(-0.34404342) q[0];
sx q[0];
rz(0.52283575) q[0];
x q[1];
rz(1.6826434) q[2];
sx q[2];
rz(-1.4750655) q[2];
sx q[2];
rz(2.7494631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0280062) q[1];
sx q[1];
rz(-1.1328814) q[1];
sx q[1];
rz(-2.6922845) q[1];
rz(2.784235) q[3];
sx q[3];
rz(-2.2893421) q[3];
sx q[3];
rz(-2.0524764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71191177) q[2];
sx q[2];
rz(-1.1744171) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(-2.1176254) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-2.021324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73151082) q[0];
sx q[0];
rz(-0.54786587) q[0];
sx q[0];
rz(1.6743073) q[0];
rz(3.1027377) q[1];
sx q[1];
rz(-1.7170186) q[1];
sx q[1];
rz(-2.255596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666834) q[0];
sx q[0];
rz(-1.0499448) q[0];
sx q[0];
rz(-2.0684) q[0];
rz(-0.90068659) q[2];
sx q[2];
rz(-1.613613) q[2];
sx q[2];
rz(0.17806554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6941144) q[1];
sx q[1];
rz(-1.2130108) q[1];
sx q[1];
rz(-1.5286137) q[1];
x q[2];
rz(1.4913849) q[3];
sx q[3];
rz(-1.6326427) q[3];
sx q[3];
rz(1.1140149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92750612) q[2];
sx q[2];
rz(-1.6897886) q[2];
sx q[2];
rz(2.5778594) q[2];
rz(-0.8756513) q[3];
sx q[3];
rz(-2.0893658) q[3];
sx q[3];
rz(-1.0912857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617274) q[0];
sx q[0];
rz(-0.848333) q[0];
sx q[0];
rz(2.0735829) q[0];
rz(0.39455286) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(1.669917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739597) q[0];
sx q[0];
rz(-1.9955705) q[0];
sx q[0];
rz(1.1921117) q[0];
rz(0.11278048) q[2];
sx q[2];
rz(-1.28714) q[2];
sx q[2];
rz(3.1007161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.010943451) q[1];
sx q[1];
rz(-0.92045451) q[1];
sx q[1];
rz(-2.5992812) q[1];
rz(-pi) q[2];
rz(0.25217523) q[3];
sx q[3];
rz(-2.7311595) q[3];
sx q[3];
rz(-1.3595734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7431405) q[2];
sx q[2];
rz(-1.8299711) q[2];
sx q[2];
rz(-3.0121646) q[2];
rz(-2.5982924) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.3885385) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1035136) q[0];
sx q[0];
rz(-2.1383801) q[0];
sx q[0];
rz(-0.45573086) q[0];
rz(2.575846) q[1];
sx q[1];
rz(-2.5028298) q[1];
sx q[1];
rz(-0.21557132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8889897) q[0];
sx q[0];
rz(-0.71160331) q[0];
sx q[0];
rz(-1.3093185) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59026123) q[2];
sx q[2];
rz(-1.0059085) q[2];
sx q[2];
rz(1.60162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0255007) q[1];
sx q[1];
rz(-1.3057723) q[1];
sx q[1];
rz(3.070036) q[1];
rz(-pi) q[2];
rz(1.2858538) q[3];
sx q[3];
rz(-1.4654034) q[3];
sx q[3];
rz(-1.2440497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4244708) q[2];
sx q[2];
rz(-2.8346546) q[2];
sx q[2];
rz(1.4324987) q[2];
rz(-2.0004382) q[3];
sx q[3];
rz(-1.9211831) q[3];
sx q[3];
rz(-0.46877638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6345217) q[0];
sx q[0];
rz(-2.2880726) q[0];
sx q[0];
rz(-0.61990196) q[0];
rz(1.2076521) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(0.29144731) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675596) q[0];
sx q[0];
rz(-0.58692044) q[0];
sx q[0];
rz(-1.8700325) q[0];
rz(0.37967988) q[2];
sx q[2];
rz(-0.34647339) q[2];
sx q[2];
rz(2.068067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48780729) q[1];
sx q[1];
rz(-1.8099306) q[1];
sx q[1];
rz(-0.59822786) q[1];
rz(-pi) q[2];
rz(-2.484887) q[3];
sx q[3];
rz(-1.5791594) q[3];
sx q[3];
rz(1.4073296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42346272) q[2];
sx q[2];
rz(-2.4726157) q[2];
sx q[2];
rz(-0.97406975) q[2];
rz(1.2232716) q[3];
sx q[3];
rz(-1.7134824) q[3];
sx q[3];
rz(-1.0380925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28737268) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(-0.30782345) q[0];
rz(0.27990118) q[1];
sx q[1];
rz(-0.24903909) q[1];
sx q[1];
rz(1.8797967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9591589) q[0];
sx q[0];
rz(-1.7049346) q[0];
sx q[0];
rz(-3.1073551) q[0];
rz(-1.9597739) q[2];
sx q[2];
rz(-2.2643467) q[2];
sx q[2];
rz(0.77893585) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9692291) q[1];
sx q[1];
rz(-1.5820272) q[1];
sx q[1];
rz(-1.8658691) q[1];
rz(-pi) q[2];
x q[2];
rz(0.070075017) q[3];
sx q[3];
rz(-2.0349906) q[3];
sx q[3];
rz(-1.7572559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5032924) q[2];
sx q[2];
rz(-1.1887447) q[2];
sx q[2];
rz(-0.80257455) q[2];
rz(-1.5447626) q[3];
sx q[3];
rz(-2.7464726) q[3];
sx q[3];
rz(1.6667268) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0751727) q[0];
sx q[0];
rz(-1.5479227) q[0];
sx q[0];
rz(-0.61298925) q[0];
rz(1.0523484) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(-1.8863511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57277623) q[0];
sx q[0];
rz(-0.63268393) q[0];
sx q[0];
rz(0.79214545) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0381211) q[2];
sx q[2];
rz(-1.7591068) q[2];
sx q[2];
rz(1.4092768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8415642) q[1];
sx q[1];
rz(-2.2662913) q[1];
sx q[1];
rz(2.5072875) q[1];
rz(-pi) q[2];
rz(-1.7057214) q[3];
sx q[3];
rz(-1.8308911) q[3];
sx q[3];
rz(1.8985727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5455948) q[2];
sx q[2];
rz(-2.2277446) q[2];
sx q[2];
rz(2.9929898) q[2];
rz(0.654486) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(-1.7053846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504836) q[0];
sx q[0];
rz(-1.5330667) q[0];
sx q[0];
rz(-0.0012375687) q[0];
rz(-1.2507863) q[1];
sx q[1];
rz(-2.2092399) q[1];
sx q[1];
rz(-0.95157448) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5630659) q[0];
sx q[0];
rz(-2.9349083) q[0];
sx q[0];
rz(-2.5220968) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1733093) q[2];
sx q[2];
rz(-0.70455019) q[2];
sx q[2];
rz(0.3144484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9597541) q[1];
sx q[1];
rz(-1.5262456) q[1];
sx q[1];
rz(-0.049330508) q[1];
x q[2];
rz(0.4001049) q[3];
sx q[3];
rz(-1.4830638) q[3];
sx q[3];
rz(-2.5468605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9963659) q[2];
sx q[2];
rz(-1.9731015) q[2];
sx q[2];
rz(-2.1053947) q[2];
rz(-0.869831) q[3];
sx q[3];
rz(-1.6706322) q[3];
sx q[3];
rz(-0.74365348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0319801) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(-1.1085283) q[0];
rz(-0.96500665) q[1];
sx q[1];
rz(-1.7644278) q[1];
sx q[1];
rz(-2.5722497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66273615) q[0];
sx q[0];
rz(-1.7883739) q[0];
sx q[0];
rz(2.6170501) q[0];
rz(-2.2340745) q[2];
sx q[2];
rz(-0.56625) q[2];
sx q[2];
rz(0.50390255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20726897) q[1];
sx q[1];
rz(-2.6187463) q[1];
sx q[1];
rz(-1.5542514) q[1];
x q[2];
rz(1.4042013) q[3];
sx q[3];
rz(-2.3218465) q[3];
sx q[3];
rz(-2.6566237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0137332) q[2];
sx q[2];
rz(-0.49488417) q[2];
sx q[2];
rz(-3.0950756) q[2];
rz(2.8401996) q[3];
sx q[3];
rz(-1.2208341) q[3];
sx q[3];
rz(-2.9197599) q[3];
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
rz(-2.2409869) q[0];
sx q[0];
rz(-1.931668) q[0];
sx q[0];
rz(1.3264309) q[0];
rz(-1.2251414) q[1];
sx q[1];
rz(-2.6381208) q[1];
sx q[1];
rz(-2.0874964) q[1];
rz(1.7635119) q[2];
sx q[2];
rz(-2.1753934) q[2];
sx q[2];
rz(2.8080804) q[2];
rz(-0.88171885) q[3];
sx q[3];
rz(-0.93986012) q[3];
sx q[3];
rz(-2.8840251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
