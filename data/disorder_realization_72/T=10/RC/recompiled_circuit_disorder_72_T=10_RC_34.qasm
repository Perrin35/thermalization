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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36131418) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(-1.5048024) q[0];
rz(-pi) q[1];
rz(-1.3747896) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(2.1240049) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62568134) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(-2.4019269) q[1];
rz(-pi) q[2];
rz(1.4673759) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9177861) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-1.0043253) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-2.1388334) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53045814) q[0];
sx q[0];
rz(-0.76871745) q[0];
sx q[0];
rz(-2.4918633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7582714) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(-1.0152917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7111295) q[1];
sx q[1];
rz(-1.4158447) q[1];
sx q[1];
rz(-2.1610545) q[1];
x q[2];
rz(-2.6456804) q[3];
sx q[3];
rz(-0.65973385) q[3];
sx q[3];
rz(-2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484444) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(1.3921188) q[0];
x q[1];
rz(0.29404624) q[2];
sx q[2];
rz(-1.1883192) q[2];
sx q[2];
rz(1.8665714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.98595847) q[1];
sx q[1];
rz(-1.8556719) q[1];
sx q[1];
rz(-0.023521544) q[1];
rz(-0.15150841) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(-2.1745149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(-2.5298932) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.600425) q[0];
sx q[0];
rz(-2.4006002) q[0];
sx q[0];
rz(-3.0279972) q[0];
x q[1];
rz(0.29454622) q[2];
sx q[2];
rz(-1.4809594) q[2];
sx q[2];
rz(1.3603269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1830814) q[1];
sx q[1];
rz(-1.6861048) q[1];
sx q[1];
rz(2.1469841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6226193) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(-0.94829544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(-0.11166212) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(-1.6508474) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6164242) q[2];
sx q[2];
rz(-1.6095973) q[2];
sx q[2];
rz(2.2280488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(0.032226493) q[1];
x q[2];
rz(-0.74113412) q[3];
sx q[3];
rz(-2.1949265) q[3];
sx q[3];
rz(2.4853064) q[3];
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
rz(0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148465) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-0.39658305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81178938) q[0];
sx q[0];
rz(-1.3178696) q[0];
sx q[0];
rz(-1.1008218) q[0];
rz(2.5846892) q[2];
sx q[2];
rz(-2.5957426) q[2];
sx q[2];
rz(0.57005537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(-3.0228826) q[1];
x q[2];
rz(0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(2.356785) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2130177) q[0];
sx q[0];
rz(-1.8682161) q[0];
sx q[0];
rz(-1.6638882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8225841) q[2];
sx q[2];
rz(-0.52201027) q[2];
sx q[2];
rz(-2.8587647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095110006) q[3];
sx q[3];
rz(-1.9169807) q[3];
sx q[3];
rz(-2.9444875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(2.9879925) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-3.1088366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93452867) q[0];
sx q[0];
rz(-1.8257739) q[0];
sx q[0];
rz(1.9051139) q[0];
rz(-pi) q[1];
rz(0.44185408) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(1.2043124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9084839) q[1];
sx q[1];
rz(-1.9875803) q[1];
sx q[1];
rz(-3.0770739) q[1];
rz(-pi) q[2];
x q[2];
rz(1.950374) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(-0.48535767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63472414) q[0];
sx q[0];
rz(-1.6244495) q[0];
sx q[0];
rz(3.0070478) q[0];
rz(-1.4114686) q[2];
sx q[2];
rz(-0.97852856) q[2];
sx q[2];
rz(1.9221905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34708729) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(-1.1361213) q[1];
rz(2.5599307) q[3];
sx q[3];
rz(-1.9546486) q[3];
sx q[3];
rz(1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(1.4032646) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.52206802) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-2.5433345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132968) q[0];
sx q[0];
rz(-1.7772355) q[0];
sx q[0];
rz(2.9199827) q[0];
rz(-2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-3.0002909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.977539) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(-1.3110728) q[1];
rz(-pi) q[2];
rz(2.9833097) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(-0.76009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-2.0098856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-1.4478222) q[2];
sx q[2];
rz(-0.95477827) q[2];
sx q[2];
rz(1.9670602) q[2];
rz(2.2451154) q[3];
sx q[3];
rz(-1.258068) q[3];
sx q[3];
rz(0.51545959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
