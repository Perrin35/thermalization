OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(3.9197796) q[1];
sx q[1];
rz(9.9240886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190174) q[0];
sx q[0];
rz(-2.8589773) q[0];
sx q[0];
rz(-2.9119888) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2917065) q[2];
sx q[2];
rz(-1.4404785) q[2];
sx q[2];
rz(-0.40638129) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0807304) q[1];
sx q[1];
rz(-0.84134358) q[1];
sx q[1];
rz(-1.7727477) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8226837) q[3];
sx q[3];
rz(-3.0348274) q[3];
sx q[3];
rz(-1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(-1.0144368) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6111345) q[0];
sx q[0];
rz(-0.76871745) q[0];
sx q[0];
rz(-2.4918633) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3147896) q[2];
sx q[2];
rz(-2.8892592) q[2];
sx q[2];
rz(-1.2834872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1042852) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(-0.18584713) q[1];
x q[2];
rz(2.6456804) q[3];
sx q[3];
rz(-0.65973385) q[3];
sx q[3];
rz(2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(0.88358203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040366216) q[0];
sx q[0];
rz(-1.4063615) q[0];
sx q[0];
rz(0.40584392) q[0];
rz(0.94647879) q[2];
sx q[2];
rz(-0.4779856) q[2];
sx q[2];
rz(-1.9566655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5933697) q[1];
sx q[1];
rz(1.2858461) q[1];
x q[2];
rz(-0.15150841) q[3];
sx q[3];
rz(-2.476859) q[3];
sx q[3];
rz(-0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-0.42923129) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4470091) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(1.6741333) q[0];
x q[1];
rz(1.6646531) q[2];
sx q[2];
rz(-1.2774733) q[2];
sx q[2];
rz(-2.9039127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1830814) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(-0.99460852) q[1];
x q[2];
rz(1.6226193) q[3];
sx q[3];
rz(-1.341815) q[3];
sx q[3];
rz(-2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1137431) q[0];
sx q[0];
rz(-2.5115974) q[0];
sx q[0];
rz(-1.6508474) q[0];
x q[1];
rz(-3.0643164) q[2];
sx q[2];
rz(-0.52646598) q[2];
sx q[2];
rz(0.72409814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(0.032226493) q[1];
x q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(-1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(2.6780224) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148465) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(0.39658305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81178938) q[0];
sx q[0];
rz(-1.3178696) q[0];
sx q[0];
rz(-2.0407709) q[0];
rz(-pi) q[1];
rz(2.5846892) q[2];
sx q[2];
rz(-0.54585005) q[2];
sx q[2];
rz(-0.57005537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(0.1187101) q[1];
rz(-pi) q[2];
rz(1.8144238) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(-0.43587303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(1.5863824) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-1.126359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(2.8471332) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.4462496) q[2];
sx q[2];
rz(-1.6342271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4451051) q[1];
sx q[1];
rz(-1.7587874) q[1];
sx q[1];
rz(0.45447116) q[1];
x q[2];
rz(-1.8282312) q[3];
sx q[3];
rz(-2.7830887) q[3];
sx q[3];
rz(-3.0646216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(3.0287108) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-3.1088366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(2.2420922) q[0];
rz(-pi) q[1];
rz(-0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(-1.1754787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0667463) q[1];
sx q[1];
rz(-2.7201338) q[1];
sx q[1];
rz(1.4261817) q[1];
rz(-3.0524391) q[3];
sx q[3];
rz(-1.9490064) q[3];
sx q[3];
rz(-1.13812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(2.5194871) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(0.39045236) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63472414) q[0];
sx q[0];
rz(-1.5171432) q[0];
sx q[0];
rz(3.0070478) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4114686) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(-1.2194022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5345926) q[1];
sx q[1];
rz(-1.2619164) q[1];
sx q[1];
rz(-0.81307462) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0210578) q[3];
sx q[3];
rz(-1.0362719) q[3];
sx q[3];
rz(0.083732097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(0.1594485) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-2.5433345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(3.0002909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4926589) q[1];
sx q[1];
rz(-1.6655386) q[1];
sx q[1];
rz(-1.2050864) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(-1.602802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-1.1317071) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(0.17157208) q[2];
sx q[2];
rz(-2.5149859) q[2];
sx q[2];
rz(1.7563216) q[2];
rz(-0.89647722) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];