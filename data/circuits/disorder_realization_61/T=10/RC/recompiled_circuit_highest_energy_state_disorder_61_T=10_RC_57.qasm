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
rz(1.4690118) q[0];
sx q[0];
rz(5.3634085) q[0];
sx q[0];
rz(10.07977) q[0];
rz(1.6545777) q[1];
sx q[1];
rz(4.5643212) q[1];
sx q[1];
rz(8.1017452) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56795492) q[0];
sx q[0];
rz(-1.0229737) q[0];
sx q[0];
rz(1.2090888) q[0];
rz(-2.5147314) q[2];
sx q[2];
rz(-2.755609) q[2];
sx q[2];
rz(0.71039334) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2046584) q[1];
sx q[1];
rz(-1.8465252) q[1];
sx q[1];
rz(-0.10358056) q[1];
x q[2];
rz(0.58889525) q[3];
sx q[3];
rz(-0.82020226) q[3];
sx q[3];
rz(-0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(0.65832552) q[2];
rz(-0.63411921) q[3];
sx q[3];
rz(-1.6987957) q[3];
sx q[3];
rz(1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.051801) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(-3.064503) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.8010151) q[1];
sx q[1];
rz(1.9416521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18835078) q[0];
sx q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(-1.9882697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.029424981) q[2];
sx q[2];
rz(-1.4064854) q[2];
sx q[2];
rz(-1.2677594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89948326) q[1];
sx q[1];
rz(-0.98190183) q[1];
sx q[1];
rz(-0.45581006) q[1];
rz(-pi) q[2];
x q[2];
rz(0.097278519) q[3];
sx q[3];
rz(-2.5803498) q[3];
sx q[3];
rz(-1.8986072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2526907) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-2.4661031) q[2];
rz(-3.1318956) q[3];
sx q[3];
rz(-0.24295013) q[3];
sx q[3];
rz(-3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552396) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(0.010628788) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(2.1307814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1026238) q[0];
sx q[0];
rz(-1.083923) q[0];
sx q[0];
rz(-2.7953933) q[0];
rz(0.028378475) q[2];
sx q[2];
rz(-0.88936964) q[2];
sx q[2];
rz(0.703237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39937035) q[1];
sx q[1];
rz(-2.3934747) q[1];
sx q[1];
rz(-3.0518603) q[1];
rz(-pi) q[2];
rz(-0.63186462) q[3];
sx q[3];
rz(-2.2711828) q[3];
sx q[3];
rz(-0.66955843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(-2.6348616) q[2];
rz(1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(-2.8633964) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5469359) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(1.1154255) q[0];
rz(1.2660654) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(0.87475264) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.311888) q[0];
sx q[0];
rz(-0.26712298) q[0];
sx q[0];
rz(1.9127089) q[0];
x q[1];
rz(-1.5953996) q[2];
sx q[2];
rz(-0.65381351) q[2];
sx q[2];
rz(3.1179414) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.03434788) q[1];
sx q[1];
rz(-2.6470091) q[1];
sx q[1];
rz(-0.58652189) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86438365) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(-0.19696008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(-1.9258707) q[2];
rz(-1.7440965) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(-1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.27595156) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(-1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.7313622) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263162) q[0];
sx q[0];
rz(-1.363376) q[0];
sx q[0];
rz(-2.0735969) q[0];
x q[1];
rz(0.67362154) q[2];
sx q[2];
rz(-1.3811908) q[2];
sx q[2];
rz(3.0233235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23996767) q[1];
sx q[1];
rz(-2.8291467) q[1];
sx q[1];
rz(0.0052253763) q[1];
rz(-1.1377119) q[3];
sx q[3];
rz(-1.8942465) q[3];
sx q[3];
rz(2.5057305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1269647) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(-1.8178168) q[2];
rz(1.3592367) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0303665) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(-0.02221814) q[0];
rz(-2.4413595) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(2.1846695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9807463) q[0];
sx q[0];
rz(-0.87781802) q[0];
sx q[0];
rz(-1.7599544) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21778743) q[2];
sx q[2];
rz(-1.3931257) q[2];
sx q[2];
rz(-1.2051932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0007504) q[1];
sx q[1];
rz(-2.6472808) q[1];
sx q[1];
rz(2.1782081) q[1];
x q[2];
rz(0.92659183) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(-2.7972972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2172829) q[2];
sx q[2];
rz(-1.1223015) q[2];
sx q[2];
rz(1.4420606) q[2];
rz(2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42665136) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(-0.46698025) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(-0.39042815) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4680088) q[0];
sx q[0];
rz(-1.7112186) q[0];
sx q[0];
rz(-1.6037081) q[0];
rz(-pi) q[1];
rz(3.0340334) q[2];
sx q[2];
rz(-1.2778408) q[2];
sx q[2];
rz(-2.7358472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(2.5659849) q[1];
rz(-pi) q[2];
rz(0.14777811) q[3];
sx q[3];
rz(-2.1782603) q[3];
sx q[3];
rz(-2.1097418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(-2.1245655) q[2];
rz(2.2022066) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(-0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(1.2128879) q[0];
rz(-0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903666) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(1.4506819) q[0];
rz(-0.076831623) q[2];
sx q[2];
rz(-1.4203826) q[2];
sx q[2];
rz(2.7422649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.054476995) q[1];
sx q[1];
rz(-1.3605788) q[1];
sx q[1];
rz(1.3940548) q[1];
rz(0.16984197) q[3];
sx q[3];
rz(-1.5780996) q[3];
sx q[3];
rz(0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6508871) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(-0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(-1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438743) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(-3.0978715) q[0];
rz(1.9678736) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(-0.79183212) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4872015) q[0];
sx q[0];
rz(-2.374361) q[0];
sx q[0];
rz(0.61875288) q[0];
rz(-pi) q[1];
rz(-0.20251198) q[2];
sx q[2];
rz(-1.8648476) q[2];
sx q[2];
rz(-2.685315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14230141) q[1];
sx q[1];
rz(-2.5901234) q[1];
sx q[1];
rz(-0.89889185) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6641162) q[3];
sx q[3];
rz(-1.9100185) q[3];
sx q[3];
rz(-0.71996688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23058471) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(-1.9045551) q[2];
rz(0.48219314) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1173387) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(1.3379958) q[0];
rz(-2.2894739) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.710708) q[0];
sx q[0];
rz(-2.0938494) q[0];
sx q[0];
rz(2.5780748) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.030376) q[2];
sx q[2];
rz(-1.2555946) q[2];
sx q[2];
rz(2.1325796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4070421) q[1];
sx q[1];
rz(-1.7104539) q[1];
sx q[1];
rz(-1.9662844) q[1];
rz(-pi) q[2];
rz(-2.1086333) q[3];
sx q[3];
rz(-2.4696484) q[3];
sx q[3];
rz(2.5047461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.368025) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(2.8311938) q[2];
rz(2.5144905) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8695759) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(1.9238453) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(1.942684) q[2];
sx q[2];
rz(-0.52996427) q[2];
sx q[2];
rz(0.90373273) q[2];
rz(2.8195856) q[3];
sx q[3];
rz(-1.6290355) q[3];
sx q[3];
rz(1.6215948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
