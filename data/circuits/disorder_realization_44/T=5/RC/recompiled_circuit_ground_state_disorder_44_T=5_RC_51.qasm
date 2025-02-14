OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(-2.1679572) q[1];
sx q[1];
rz(-1.4373625) q[1];
sx q[1];
rz(0.91926423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0231173) q[0];
sx q[0];
rz(-3.1065515) q[0];
sx q[0];
rz(-0.50993292) q[0];
rz(-pi) q[1];
rz(-1.6617352) q[2];
sx q[2];
rz(-1.6906066) q[2];
sx q[2];
rz(0.38209846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8265958) q[1];
sx q[1];
rz(-0.8657786) q[1];
sx q[1];
rz(-0.12048529) q[1];
x q[2];
rz(0.42114847) q[3];
sx q[3];
rz(-2.0055565) q[3];
sx q[3];
rz(-2.8041294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(-0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(0.10281674) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8973812) q[0];
sx q[0];
rz(-2.057071) q[0];
sx q[0];
rz(2.4308423) q[0];
rz(-0.76849014) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(2.1220727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7026414) q[0];
sx q[0];
rz(-1.1239237) q[0];
sx q[0];
rz(1.5269482) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0654991) q[2];
sx q[2];
rz(-2.00092) q[2];
sx q[2];
rz(-0.96904749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.63377127) q[1];
sx q[1];
rz(-1.1999167) q[1];
sx q[1];
rz(0.033286496) q[1];
x q[2];
rz(-3.0583303) q[3];
sx q[3];
rz(-2.1234649) q[3];
sx q[3];
rz(-3.0577615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(-1.2126806) q[3];
sx q[3];
rz(-1.0168394) q[3];
sx q[3];
rz(-2.4299664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(-1.9236176) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(-2.2191494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3271752) q[0];
sx q[0];
rz(-1.6285768) q[0];
sx q[0];
rz(1.6114651) q[0];
rz(0.72539665) q[2];
sx q[2];
rz(-1.9624406) q[2];
sx q[2];
rz(-1.2466516) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6092012) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(-2.3228541) q[1];
rz(-0.32914583) q[3];
sx q[3];
rz(-1.7308047) q[3];
sx q[3];
rz(0.015586675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1232274) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(2.5943622) q[3];
sx q[3];
rz(-0.42357835) q[3];
sx q[3];
rz(-2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0860586) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-2.8841282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2672294) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(0.49212014) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6370378) q[2];
sx q[2];
rz(-0.77755723) q[2];
sx q[2];
rz(-2.24868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3066669) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(-2.5132781) q[1];
rz(-pi) q[2];
rz(2.2063755) q[3];
sx q[3];
rz(-0.87163371) q[3];
sx q[3];
rz(-2.2548667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3782392) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(-2.4371448) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(-2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064297) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.7154891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97010342) q[0];
sx q[0];
rz(-1.3521776) q[0];
sx q[0];
rz(-1.4755558) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9272396) q[2];
sx q[2];
rz(-0.39696908) q[2];
sx q[2];
rz(0.036338417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1864422) q[1];
sx q[1];
rz(-1.3480061) q[1];
sx q[1];
rz(2.8563315) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9687443) q[3];
sx q[3];
rz(-1.8573055) q[3];
sx q[3];
rz(-2.2423173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9466729) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(1.8348414) q[2];
rz(1.9715747) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65866798) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(2.7401127) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(-0.85404095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14388785) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(-1.7884939) q[0];
rz(-pi) q[1];
rz(-1.7207022) q[2];
sx q[2];
rz(-1.7542363) q[2];
sx q[2];
rz(2.4114087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5062949) q[1];
sx q[1];
rz(-0.98165252) q[1];
sx q[1];
rz(1.6904669) q[1];
rz(-pi) q[2];
rz(0.60131945) q[3];
sx q[3];
rz(-2.4375705) q[3];
sx q[3];
rz(2.6773334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(1.1775449) q[2];
rz(2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(-0.83561713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(1.1516736) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(0.91845671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8480523) q[0];
sx q[0];
rz(-1.9013202) q[0];
sx q[0];
rz(-2.8577515) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8807441) q[2];
sx q[2];
rz(-2.2493304) q[2];
sx q[2];
rz(1.4229753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53808182) q[1];
sx q[1];
rz(-0.0099364837) q[1];
sx q[1];
rz(2.4007829) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2270532) q[3];
sx q[3];
rz(-1.731428) q[3];
sx q[3];
rz(2.4317222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.7983961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4485432) q[0];
sx q[0];
rz(-2.9149084) q[0];
sx q[0];
rz(1.6402798) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(2.0223845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6014746) q[0];
sx q[0];
rz(-2.1261423) q[0];
sx q[0];
rz(2.0106273) q[0];
rz(-pi) q[1];
rz(-2.881024) q[2];
sx q[2];
rz(-2.0379279) q[2];
sx q[2];
rz(-1.9595944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70485605) q[1];
sx q[1];
rz(-2.4697127) q[1];
sx q[1];
rz(-2.9465972) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7469278) q[3];
sx q[3];
rz(-1.695444) q[3];
sx q[3];
rz(3.0277071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98562733) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(-0.30820942) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(-1.6109899) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823328) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(2.7985213) q[0];
rz(-1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(-1.7128568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7633782) q[0];
sx q[0];
rz(-1.175048) q[0];
sx q[0];
rz(1.8431435) q[0];
rz(-2.3104383) q[2];
sx q[2];
rz(-0.88696431) q[2];
sx q[2];
rz(-1.2338975) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45537585) q[1];
sx q[1];
rz(-1.1955402) q[1];
sx q[1];
rz(0.065836716) q[1];
rz(-2.1200646) q[3];
sx q[3];
rz(-0.77199751) q[3];
sx q[3];
rz(-2.271351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(1.3915871) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(-1.3579177) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(2.3416669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185266) q[0];
sx q[0];
rz(-0.79341054) q[0];
sx q[0];
rz(1.8961468) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8736035) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-0.15632665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8854196) q[1];
sx q[1];
rz(-0.90682632) q[1];
sx q[1];
rz(-1.3732984) q[1];
rz(-pi) q[2];
rz(-2.0560451) q[3];
sx q[3];
rz(-1.3486514) q[3];
sx q[3];
rz(-2.9911161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7633729) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(-2.8857005) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(-1.488744) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8681317) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(1.801626) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(-2.432178) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(-2.1382016) q[3];
sx q[3];
rz(-2.1324674) q[3];
sx q[3];
rz(2.4218925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
