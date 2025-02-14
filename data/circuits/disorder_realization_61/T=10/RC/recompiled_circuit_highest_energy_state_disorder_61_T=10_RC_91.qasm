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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(2.4866009) q[0];
rz(1.6545777) q[1];
sx q[1];
rz(-1.7188641) q[1];
sx q[1];
rz(-1.3230327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5736377) q[0];
sx q[0];
rz(-1.0229737) q[0];
sx q[0];
rz(-1.2090888) q[0];
rz(-pi) q[1];
rz(-2.8236515) q[2];
sx q[2];
rz(-1.3481209) q[2];
sx q[2];
rz(0.26938619) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93693426) q[1];
sx q[1];
rz(-1.8465252) q[1];
sx q[1];
rz(0.10358056) q[1];
rz(-0.58889525) q[3];
sx q[3];
rz(-2.3213904) q[3];
sx q[3];
rz(-0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46372867) q[2];
sx q[2];
rz(-0.68631309) q[2];
sx q[2];
rz(2.4832671) q[2];
rz(-0.63411921) q[3];
sx q[3];
rz(-1.6987957) q[3];
sx q[3];
rz(1.3707976) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089791678) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(0.077089699) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(1.9416521) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.752992) q[0];
sx q[0];
rz(-1.7677099) q[0];
sx q[0];
rz(1.1042751) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3951833) q[2];
sx q[2];
rz(-2.9746911) q[2];
sx q[2];
rz(1.4457955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6226095) q[1];
sx q[1];
rz(-2.4138192) q[1];
sx q[1];
rz(-0.98811291) q[1];
rz(-pi) q[2];
rz(-1.6317815) q[3];
sx q[3];
rz(-2.1290696) q[3];
sx q[3];
rz(-1.7838441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8889019) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(-0.67548951) q[2];
rz(-3.1318956) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(-0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49552396) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(2.38696) q[0];
rz(0.010628788) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-1.0108112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6957798) q[0];
sx q[0];
rz(-0.58923972) q[0];
sx q[0];
rz(1.0007829) q[0];
x q[1];
rz(-2.25242) q[2];
sx q[2];
rz(-1.592836) q[2];
sx q[2];
rz(-2.2561548) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39937035) q[1];
sx q[1];
rz(-0.74811799) q[1];
sx q[1];
rz(-0.089732371) q[1];
rz(-pi) q[2];
rz(2.378024) q[3];
sx q[3];
rz(-1.1022304) q[3];
sx q[3];
rz(2.6811622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2406771) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(0.50673103) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.3576018) q[3];
sx q[3];
rz(0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5469359) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(-2.0261672) q[0];
rz(-1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-2.26684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0698358) q[0];
sx q[0];
rz(-1.6594145) q[0];
sx q[0];
rz(1.3184692) q[0];
rz(-1.546193) q[2];
sx q[2];
rz(-2.4877791) q[2];
sx q[2];
rz(3.1179414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.03434788) q[1];
sx q[1];
rz(-2.6470091) q[1];
sx q[1];
rz(-2.5550708) q[1];
rz(-0.86438365) q[3];
sx q[3];
rz(-1.7076075) q[3];
sx q[3];
rz(0.19696008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(-1.215722) q[2];
rz(-1.3974961) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.6708466) q[0];
rz(-1.7631081) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.4102304) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61527644) q[0];
sx q[0];
rz(-1.363376) q[0];
sx q[0];
rz(-1.0679958) q[0];
rz(0.29844923) q[2];
sx q[2];
rz(-2.4458234) q[2];
sx q[2];
rz(-1.6844105) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3258563) q[1];
sx q[1];
rz(-1.5724025) q[1];
sx q[1];
rz(2.8291507) q[1];
x q[2];
rz(-1.1377119) q[3];
sx q[3];
rz(-1.2473462) q[3];
sx q[3];
rz(0.63586217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1269647) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(-1.3592367) q[3];
sx q[3];
rz(-1.8424572) q[3];
sx q[3];
rz(0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(3.1193745) q[0];
rz(0.70023316) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(-0.95692316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6895804) q[0];
sx q[0];
rz(-2.4274024) q[0];
sx q[0];
rz(-0.22269188) q[0];
x q[1];
rz(2.9238052) q[2];
sx q[2];
rz(-1.3931257) q[2];
sx q[2];
rz(-1.9363994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80920389) q[1];
sx q[1];
rz(-1.1706377) q[1];
sx q[1];
rz(-0.2984115) q[1];
rz(-pi) q[2];
rz(2.0104353) q[3];
sx q[3];
rz(-1.2613858) q[3];
sx q[3];
rz(0.65195665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2172829) q[2];
sx q[2];
rz(-1.1223015) q[2];
sx q[2];
rz(1.699532) q[2];
rz(1.1066655) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(1.1520011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(2.6746124) q[0];
rz(1.3106208) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4425499) q[0];
sx q[0];
rz(-0.14420284) q[0];
sx q[0];
rz(-0.22871916) q[0];
x q[1];
rz(3.0340334) q[2];
sx q[2];
rz(-1.8637519) q[2];
sx q[2];
rz(-0.40574542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3876593) q[1];
sx q[1];
rz(-1.1186592) q[1];
sx q[1];
rz(-2.2958295) q[1];
rz(-0.14777811) q[3];
sx q[3];
rz(-2.1782603) q[3];
sx q[3];
rz(2.1097418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9499669) q[2];
sx q[2];
rz(-1.2167296) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(-2.2022066) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(-2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(-1.2128879) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-1.1319356) q[1];
sx q[1];
rz(-2.718198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0807665) q[0];
sx q[0];
rz(-0.87196022) q[0];
sx q[0];
rz(-0.10159512) q[0];
x q[1];
rz(1.7216484) q[2];
sx q[2];
rz(-1.6467588) q[2];
sx q[2];
rz(1.9585889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5535721) q[1];
sx q[1];
rz(-1.3979853) q[1];
sx q[1];
rz(-0.2134448) q[1];
rz(-pi) q[2];
rz(1.5782062) q[3];
sx q[3];
rz(-1.4009589) q[3];
sx q[3];
rz(-1.3337204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49070552) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(-2.1631964) q[2];
rz(2.4466416) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438743) q[0];
sx q[0];
rz(-2.9053423) q[0];
sx q[0];
rz(-0.043721113) q[0];
rz(-1.1737191) q[1];
sx q[1];
rz(-2.3269188) q[1];
sx q[1];
rz(-2.3497605) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2671473) q[0];
sx q[0];
rz(-2.1717779) q[0];
sx q[0];
rz(-2.0807666) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2709684) q[2];
sx q[2];
rz(-1.3770896) q[2];
sx q[2];
rz(1.1739588) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9992912) q[1];
sx q[1];
rz(-2.5901234) q[1];
sx q[1];
rz(0.89889185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1926226) q[3];
sx q[3];
rz(-2.0190051) q[3];
sx q[3];
rz(-2.1202655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(1.2370375) q[2];
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
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.024254) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(1.8035969) q[0];
rz(-2.2894739) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8088732) q[0];
sx q[0];
rz(-0.74902636) q[0];
sx q[0];
rz(2.3179884) q[0];
rz(-pi) q[1];
rz(0.11121662) q[2];
sx q[2];
rz(-1.885998) q[2];
sx q[2];
rz(1.0090131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.627189) q[1];
sx q[1];
rz(-2.7233988) q[1];
sx q[1];
rz(1.2209284) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0329594) q[3];
sx q[3];
rz(-0.67194429) q[3];
sx q[3];
rz(-0.6368466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77356768) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(0.31039882) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(-1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695759) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(-1.2177474) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(-0.2097585) q[2];
sx q[2];
rz(-1.0804313) q[2];
sx q[2];
rz(-2.6624138) q[2];
rz(-0.18219215) q[3];
sx q[3];
rz(-2.8145418) q[3];
sx q[3];
rz(-2.9180632) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
