OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(1.7680661) q[0];
sx q[0];
rz(11.058523) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190174) q[0];
sx q[0];
rz(-0.28261533) q[0];
sx q[0];
rz(2.9119888) q[0];
rz(-2.9688641) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(1.8634862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3589188) q[1];
sx q[1];
rz(-2.3896857) q[1];
sx q[1];
rz(2.9208675) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1148881) q[3];
sx q[3];
rz(-1.6741802) q[3];
sx q[3];
rz(-1.0552693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(-2.1388334) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8589477) q[0];
sx q[0];
rz(-0.98416057) q[0];
sx q[0];
rz(-1.0413917) q[0];
rz(-pi) q[1];
rz(0.1728671) q[2];
sx q[2];
rz(-1.3860518) q[2];
sx q[2];
rz(2.6181521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36675378) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(1.8444091) q[1];
rz(2.5428883) q[3];
sx q[3];
rz(-1.2748534) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(-0.15094748) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944091) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(0.39788525) q[0];
x q[1];
rz(2.8475464) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(-1.2750212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5933697) q[1];
sx q[1];
rz(1.2858461) q[1];
rz(2.482445) q[3];
sx q[3];
rz(-1.4775606) q[3];
sx q[3];
rz(2.6574709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(1.6741333) q[0];
x q[1];
rz(-0.30087154) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(3.0645264) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(-2.1469841) q[1];
x q[2];
rz(-0.21869603) q[3];
sx q[3];
rz(-0.23467206) q[3];
sx q[3];
rz(-2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-2.2263288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(-1.6508474) q[0];
rz(-3.0643164) q[2];
sx q[2];
rz(-0.52646598) q[2];
sx q[2];
rz(0.72409814) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(0.032226493) q[1];
x q[2];
rz(-2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81178938) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(-1.1008218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47607143) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(-1.6518041) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(-pi) q[2];
rz(1.3271689) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(-2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(1.126359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267245) q[0];
sx q[0];
rz(-1.65979) q[0];
sx q[0];
rz(0.29863775) q[0];
rz(0.14234219) q[2];
sx q[2];
rz(-2.0747613) q[2];
sx q[2];
rz(-0.0056643639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(1.7794442) q[1];
x q[2];
rz(1.2231636) q[3];
sx q[3];
rz(-1.6602483) q[3];
sx q[3];
rz(-1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927133) q[0];
sx q[0];
rz(-1.8939051) q[0];
sx q[0];
rz(0.26922853) q[0];
rz(-pi) q[1];
rz(0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(-1.966114) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2331088) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(-0.064518708) q[1];
x q[2];
rz(-1.950374) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(0.4656725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-0.48535767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94333121) q[0];
sx q[0];
rz(-1.4364463) q[0];
sx q[0];
rz(-1.6249379) q[0];
rz(2.543407) q[2];
sx q[2];
rz(-1.7028114) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-2.328518) q[1];
rz(-pi) q[2];
rz(-0.63391179) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(0.5982582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0043250672) q[0];
sx q[0];
rz(-2.8398872) q[0];
sx q[0];
rz(0.76122491) q[0];
x q[1];
rz(0.21226378) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-3.0002909) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(-1.2050864) q[1];
rz(-0.15828295) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-1.4478222) q[2];
sx q[2];
rz(-0.95477827) q[2];
sx q[2];
rz(1.9670602) q[2];
rz(2.7491309) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
