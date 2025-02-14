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
rz(2.0107438) q[0];
sx q[0];
rz(-1.7717489) q[0];
sx q[0];
rz(-2.9983591) q[0];
rz(3.1205966) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(2.6666759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.582481) q[0];
sx q[0];
rz(-1.6266291) q[0];
sx q[0];
rz(-3.0654967) q[0];
rz(-pi) q[1];
rz(-0.61491809) q[2];
sx q[2];
rz(-1.5496829) q[2];
sx q[2];
rz(-2.9647765) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.43496599) q[1];
sx q[1];
rz(-0.7966704) q[1];
sx q[1];
rz(0.15310751) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1056312) q[3];
sx q[3];
rz(-1.1316) q[3];
sx q[3];
rz(0.10094563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1249866) q[2];
sx q[2];
rz(-1.6281444) q[2];
sx q[2];
rz(2.5561257) q[2];
rz(0.79750195) q[3];
sx q[3];
rz(-1.5471285) q[3];
sx q[3];
rz(-2.7378979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0320691) q[0];
sx q[0];
rz(-1.7970947) q[0];
sx q[0];
rz(-2.4865785) q[0];
rz(-0.0034927448) q[1];
sx q[1];
rz(-2.4066636) q[1];
sx q[1];
rz(-0.014911501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61136041) q[0];
sx q[0];
rz(-1.6515628) q[0];
sx q[0];
rz(-1.5967036) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0833561) q[2];
sx q[2];
rz(-0.76448802) q[2];
sx q[2];
rz(-0.22138324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88483459) q[1];
sx q[1];
rz(-1.2169098) q[1];
sx q[1];
rz(-2.3164151) q[1];
rz(-pi) q[2];
rz(3.0767617) q[3];
sx q[3];
rz(-0.81063834) q[3];
sx q[3];
rz(-3.1136581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1917176) q[2];
sx q[2];
rz(-1.8052552) q[2];
sx q[2];
rz(0.28676644) q[2];
rz(-3.1112572) q[3];
sx q[3];
rz(-1.560863) q[3];
sx q[3];
rz(-1.102977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.2716118) q[0];
sx q[0];
rz(-2.5219707) q[0];
sx q[0];
rz(-1.9770589) q[0];
rz(-2.8097235) q[1];
sx q[1];
rz(-1.7121366) q[1];
sx q[1];
rz(1.9934995) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6004922) q[0];
sx q[0];
rz(-2.9229188) q[0];
sx q[0];
rz(-1.650901) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9752047) q[2];
sx q[2];
rz(-1.0999318) q[2];
sx q[2];
rz(1.3883049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0169605) q[1];
sx q[1];
rz(-1.3784882) q[1];
sx q[1];
rz(-1.8796199) q[1];
rz(-0.049552187) q[3];
sx q[3];
rz(-1.3930071) q[3];
sx q[3];
rz(-0.73977284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1295192) q[2];
sx q[2];
rz(-1.4612863) q[2];
sx q[2];
rz(0.47895437) q[2];
rz(2.3524513) q[3];
sx q[3];
rz(-2.456587) q[3];
sx q[3];
rz(-1.7710549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0399465) q[0];
sx q[0];
rz(-0.78831035) q[0];
sx q[0];
rz(-0.76860651) q[0];
rz(0.45306122) q[1];
sx q[1];
rz(-2.1078347) q[1];
sx q[1];
rz(0.57910848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6957618) q[0];
sx q[0];
rz(-0.33798744) q[0];
sx q[0];
rz(-0.72622444) q[0];
rz(0.3884811) q[2];
sx q[2];
rz(-1.0916096) q[2];
sx q[2];
rz(0.84333429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6542289) q[1];
sx q[1];
rz(-1.0768862) q[1];
sx q[1];
rz(-3.0904675) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7594537) q[3];
sx q[3];
rz(-1.9629729) q[3];
sx q[3];
rz(0.79526633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1733178) q[2];
sx q[2];
rz(-1.3710794) q[2];
sx q[2];
rz(1.0260065) q[2];
rz(2.6722028) q[3];
sx q[3];
rz(-1.0204851) q[3];
sx q[3];
rz(-2.6494086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-0.39407179) q[0];
sx q[0];
rz(-1.515027) q[0];
sx q[0];
rz(2.1409905) q[0];
rz(-1.3032777) q[1];
sx q[1];
rz(-2.3293827) q[1];
sx q[1];
rz(-0.8808879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6495889) q[0];
sx q[0];
rz(-1.0179217) q[0];
sx q[0];
rz(-2.1472424) q[0];
rz(-pi) q[1];
rz(-0.32994907) q[2];
sx q[2];
rz(-0.68701982) q[2];
sx q[2];
rz(-2.6181412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92426244) q[1];
sx q[1];
rz(-1.2479022) q[1];
sx q[1];
rz(2.4004637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39805746) q[3];
sx q[3];
rz(-2.0144793) q[3];
sx q[3];
rz(-1.8026082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96054968) q[2];
sx q[2];
rz(-2.0281823) q[2];
sx q[2];
rz(1.95365) q[2];
rz(-0.12399593) q[3];
sx q[3];
rz(-1.8143727) q[3];
sx q[3];
rz(-2.7441062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7749629) q[0];
sx q[0];
rz(-2.9149945) q[0];
sx q[0];
rz(-0.15283787) q[0];
rz(2.9951908) q[1];
sx q[1];
rz(-1.791879) q[1];
sx q[1];
rz(-0.70405594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2938151) q[0];
sx q[0];
rz(-3.1070947) q[0];
sx q[0];
rz(0.52963799) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33284171) q[2];
sx q[2];
rz(-2.2031076) q[2];
sx q[2];
rz(1.2472356) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0378677) q[1];
sx q[1];
rz(-1.6462391) q[1];
sx q[1];
rz(1.451205) q[1];
rz(-pi) q[2];
rz(-2.0364159) q[3];
sx q[3];
rz(-2.4676968) q[3];
sx q[3];
rz(0.020581882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.081508962) q[2];
sx q[2];
rz(-2.1174105) q[2];
sx q[2];
rz(0.3788968) q[2];
rz(2.1156408) q[3];
sx q[3];
rz(-2.43695) q[3];
sx q[3];
rz(2.0254693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7260471) q[0];
sx q[0];
rz(-1.8141831) q[0];
sx q[0];
rz(-2.5782247) q[0];
rz(2.7401961) q[1];
sx q[1];
rz(-1.0029663) q[1];
sx q[1];
rz(0.67515236) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440985) q[0];
sx q[0];
rz(-1.5821274) q[0];
sx q[0];
rz(-1.5392153) q[0];
rz(-0.60826081) q[2];
sx q[2];
rz(-1.741) q[2];
sx q[2];
rz(0.93675429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71175486) q[1];
sx q[1];
rz(-0.40374175) q[1];
sx q[1];
rz(-0.92947361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3067857) q[3];
sx q[3];
rz(-1.2865598) q[3];
sx q[3];
rz(-0.72407297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99951619) q[2];
sx q[2];
rz(-0.88366214) q[2];
sx q[2];
rz(2.7189972) q[2];
rz(-0.26250592) q[3];
sx q[3];
rz(-1.0650977) q[3];
sx q[3];
rz(2.0601823) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.066684) q[0];
sx q[0];
rz(-0.96195641) q[0];
sx q[0];
rz(0.99895507) q[0];
rz(-1.2133489) q[1];
sx q[1];
rz(-1.9478925) q[1];
sx q[1];
rz(0.34128571) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0851819) q[0];
sx q[0];
rz(-0.23427948) q[0];
sx q[0];
rz(-2.5368669) q[0];
x q[1];
rz(0.26087324) q[2];
sx q[2];
rz(-0.91756454) q[2];
sx q[2];
rz(3.0034371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5708551) q[1];
sx q[1];
rz(-2.8767817) q[1];
sx q[1];
rz(-2.6874506) q[1];
rz(1.8496277) q[3];
sx q[3];
rz(-1.5145455) q[3];
sx q[3];
rz(1.0644399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21218941) q[2];
sx q[2];
rz(-1.5985649) q[2];
sx q[2];
rz(2.8821778) q[2];
rz(2.8560396) q[3];
sx q[3];
rz(-2.4694337) q[3];
sx q[3];
rz(1.0478421) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015273) q[0];
sx q[0];
rz(-2.6492388) q[0];
sx q[0];
rz(-2.0557851) q[0];
rz(1.0049413) q[1];
sx q[1];
rz(-1.137038) q[1];
sx q[1];
rz(0.55652943) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5421931) q[0];
sx q[0];
rz(-1.43723) q[0];
sx q[0];
rz(-0.041102701) q[0];
rz(-2.25726) q[2];
sx q[2];
rz(-2.2811675) q[2];
sx q[2];
rz(-2.7820003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2697325) q[1];
sx q[1];
rz(-0.29127866) q[1];
sx q[1];
rz(-0.32637432) q[1];
rz(-pi) q[2];
rz(0.50332467) q[3];
sx q[3];
rz(-1.1098813) q[3];
sx q[3];
rz(3.1129376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0050547) q[2];
sx q[2];
rz(-1.7359066) q[2];
sx q[2];
rz(1.8704937) q[2];
rz(-2.2942885) q[3];
sx q[3];
rz(-1.3264791) q[3];
sx q[3];
rz(0.29720753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.51410455) q[0];
sx q[0];
rz(-1.6922373) q[0];
sx q[0];
rz(2.3741212) q[0];
rz(-0.96495676) q[1];
sx q[1];
rz(-1.283353) q[1];
sx q[1];
rz(-2.8585785) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8022158) q[0];
sx q[0];
rz(-0.87516038) q[0];
sx q[0];
rz(-2.8921769) q[0];
rz(-0.04472132) q[2];
sx q[2];
rz(-2.7095717) q[2];
sx q[2];
rz(0.17145874) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9486299) q[1];
sx q[1];
rz(-1.9300736) q[1];
sx q[1];
rz(2.6032843) q[1];
rz(-1.0177294) q[3];
sx q[3];
rz(-0.8700287) q[3];
sx q[3];
rz(1.7297945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9877801) q[2];
sx q[2];
rz(-0.38975468) q[2];
sx q[2];
rz(-1.4022931) q[2];
rz(0.66579372) q[3];
sx q[3];
rz(-1.2716764) q[3];
sx q[3];
rz(-1.3027035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6233728) q[0];
sx q[0];
rz(-1.2517396) q[0];
sx q[0];
rz(-1.804833) q[0];
rz(-1.549859) q[1];
sx q[1];
rz(-1.5668329) q[1];
sx q[1];
rz(-0.11650539) q[1];
rz(2.5575175) q[2];
sx q[2];
rz(-1.1933977) q[2];
sx q[2];
rz(-2.615821) q[2];
rz(-1.3605098) q[3];
sx q[3];
rz(-2.3514699) q[3];
sx q[3];
rz(-2.0024553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
