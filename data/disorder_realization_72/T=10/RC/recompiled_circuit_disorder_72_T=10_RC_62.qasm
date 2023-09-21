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
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(1.5048024) q[0];
x q[1];
rz(-0.17272858) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(-1.8634862) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7826739) q[1];
sx q[1];
rz(-0.75190699) q[1];
sx q[1];
rz(0.2207252) q[1];
rz(-pi) q[2];
rz(1.4673759) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402055) q[0];
sx q[0];
rz(-2.0048855) q[0];
sx q[0];
rz(-2.4853404) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3147896) q[2];
sx q[2];
rz(-0.25233341) q[2];
sx q[2];
rz(1.8581055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4304632) q[1];
sx q[1];
rz(-1.7257479) q[1];
sx q[1];
rz(0.98053812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9243794) q[3];
sx q[3];
rz(-1.0014605) q[3];
sx q[3];
rz(-2.6667037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-0.99439159) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484444) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(-2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(-1.7494739) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8475464) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(1.8665714) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90244251) q[1];
sx q[1];
rz(-2.855774) q[1];
sx q[1];
rz(-1.4906537) q[1];
rz(-pi) q[2];
rz(2.9900842) q[3];
sx q[3];
rz(-2.476859) q[3];
sx q[3];
rz(-0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(0.99003506) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(-2.591419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94565369) q[0];
sx q[0];
rz(-1.6473856) q[0];
sx q[0];
rz(2.4038195) q[0];
rz(2.8470464) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(-1.7812658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53766996) q[1];
sx q[1];
rz(-0.99891716) q[1];
sx q[1];
rz(-0.13725431) q[1];
rz(-2.9228966) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(0.72363879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-2.8660529) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-0.91526389) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(-1.6508474) q[0];
rz(-pi) q[1];
rz(-2.6164242) q[2];
sx q[2];
rz(-1.5319954) q[2];
sx q[2];
rz(2.2280488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.535714) q[1];
sx q[1];
rz(-1.6026346) q[1];
sx q[1];
rz(1.4153626) q[1];
rz(-pi) q[2];
rz(2.3239273) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(-2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.0625431) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81178938) q[0];
sx q[0];
rz(-1.3178696) q[0];
sx q[0];
rz(-2.0407709) q[0];
rz(-pi) q[1];
rz(0.47607143) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(1.4897886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10627667) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.7229401) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7820884) q[3];
sx q[3];
rz(-2.5264969) q[3];
sx q[3];
rz(2.2744327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(-1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(2.0152337) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(-0.29445946) q[0];
rz(-2.0790714) q[2];
sx q[2];
rz(-1.4462496) q[2];
sx q[2];
rz(1.5073656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78325242) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(1.7794442) q[1];
rz(-pi) q[2];
rz(3.0464826) q[3];
sx q[3];
rz(-1.2246119) q[3];
sx q[3];
rz(2.9444875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.37384) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(0.032756068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.207064) q[0];
sx q[0];
rz(-1.3158187) q[0];
sx q[0];
rz(1.9051139) q[0];
rz(-2.2135127) q[2];
sx q[2];
rz(-2.4734481) q[2];
sx q[2];
rz(-1.1754787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2331088) q[1];
sx q[1];
rz(-1.1540124) q[1];
sx q[1];
rz(3.0770739) q[1];
rz(0.0891536) q[3];
sx q[3];
rz(-1.1925863) q[3];
sx q[3];
rz(1.13812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
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
rz(-0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.5266248) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(0.48535767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94333121) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(1.5166548) q[0];
x q[1];
rz(2.910026) q[2];
sx q[2];
rz(-0.61083691) q[2];
sx q[2];
rz(1.6419186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-0.81307462) q[1];
rz(-pi) q[2];
rz(-0.63391179) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-0.1594485) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.7074701) q[0];
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
rz(-0.78865096) q[0];
sx q[0];
rz(-1.3539679) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(2.9293289) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-0.14130172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64893374) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.2050864) q[1];
rz(-pi) q[2];
rz(-1.4062905) q[3];
sx q[3];
rz(-1.4146155) q[3];
sx q[3];
rz(-0.7848878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(0.65199488) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(0.17157208) q[2];
sx q[2];
rz(-2.5149859) q[2];
sx q[2];
rz(1.7563216) q[2];
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
