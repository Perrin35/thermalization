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
rz(-1.3735266) q[0];
sx q[0];
rz(-1.6337448) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(3.9197796) q[1];
sx q[1];
rz(9.9240886) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12257523) q[0];
sx q[0];
rz(-0.28261533) q[0];
sx q[0];
rz(0.22960381) q[0];
rz(1.7668031) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(2.1240049) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3589188) q[1];
sx q[1];
rz(-0.75190699) q[1];
sx q[1];
rz(0.2207252) q[1];
rz(-pi) q[2];
rz(1.6742168) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(-2.6288222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-1.0027592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8589477) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(-2.100201) q[0];
rz(-pi) q[1];
rz(0.82680306) q[2];
sx q[2];
rz(-0.25233341) q[2];
sx q[2];
rz(1.8581055) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4304632) q[1];
sx q[1];
rz(-1.4158447) q[1];
sx q[1];
rz(-0.98053812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9243794) q[3];
sx q[3];
rz(-2.1401322) q[3];
sx q[3];
rz(-0.47488892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-0.77899581) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(-1.7494739) q[0];
rz(-pi) q[1];
rz(-2.1951139) q[2];
sx q[2];
rz(-2.6636071) q[2];
sx q[2];
rz(-1.1849272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(1.6509389) q[1];
x q[2];
rz(1.6885353) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(-1.1586787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69458354) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(1.4674594) q[0];
rz(-pi) q[1];
rz(2.8407211) q[2];
sx q[2];
rz(-0.30756018) q[2];
sx q[2];
rz(-0.077066271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6039227) q[1];
sx q[1];
rz(-0.99891716) q[1];
sx q[1];
rz(0.13725431) q[1];
x q[2];
rz(1.5189734) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(0.94829544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(0.91526389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-2.1992654) q[0];
x q[1];
rz(-2.6164242) q[2];
sx q[2];
rz(-1.5319954) q[2];
sx q[2];
rz(-0.9135439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60587864) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(1.72623) q[1];
rz(-pi) q[2];
rz(-0.74113412) q[3];
sx q[3];
rz(-2.1949265) q[3];
sx q[3];
rz(-0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-2.129014) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81178938) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(-2.0407709) q[0];
x q[1];
rz(1.8814538) q[2];
sx q[2];
rz(-2.0271745) q[2];
sx q[2];
rz(3.0820456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.035316) q[1];
sx q[1];
rz(-0.66674399) q[1];
sx q[1];
rz(-1.7229401) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(-2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-1.126359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(-2.8471332) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3190086) q[2];
sx q[2];
rz(-0.52201027) q[2];
sx q[2];
rz(0.28282794) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2396444) q[1];
sx q[1];
rz(-2.6522954) q[1];
sx q[1];
rz(-0.40892618) q[1];
rz(-1.3133615) q[3];
sx q[3];
rz(-0.35850393) q[3];
sx q[3];
rz(0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(2.9879925) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.37384) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54887933) q[0];
sx q[0];
rz(-1.2476876) q[0];
sx q[0];
rz(-2.8723641) q[0];
rz(-pi) q[1];
rz(0.44185408) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(1.2043124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.9883518) q[1];
x q[2];
rz(-1.350358) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(2.2409852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982614) q[0];
sx q[0];
rz(-1.4364463) q[0];
sx q[0];
rz(1.5166548) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7301241) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(1.2194022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8977412) q[1];
sx q[1];
rz(-2.284639) q[1];
sx q[1];
rz(-0.41390151) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58166196) q[3];
sx q[3];
rz(-1.9546486) q[3];
sx q[3];
rz(1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-2.5433345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.7822595) q[0];
rz(2.5355849) q[2];
sx q[2];
rz(-2.8851644) q[2];
sx q[2];
rz(2.3026349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2559402) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(3.0401858) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-0.17157208) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-1.0929575) q[3];
sx q[3];
rz(-0.7328877) q[3];
sx q[3];
rz(2.4536798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];