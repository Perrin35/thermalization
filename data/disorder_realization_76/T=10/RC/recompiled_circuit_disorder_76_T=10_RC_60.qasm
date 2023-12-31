OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(5.4260317) q[0];
sx q[0];
rz(9.5572588) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(0.69256988) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(1.7696487) q[0];
rz(-pi) q[1];
rz(-0.64451005) q[2];
sx q[2];
rz(-1.1877726) q[2];
sx q[2];
rz(-0.884998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2107271) q[1];
sx q[1];
rz(-0.14666808) q[1];
sx q[1];
rz(1.055483) q[1];
rz(-pi) q[2];
rz(-1.0971783) q[3];
sx q[3];
rz(-1.7930822) q[3];
sx q[3];
rz(1.3003295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(1.6050603) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55727977) q[0];
sx q[0];
rz(-1.6788388) q[0];
sx q[0];
rz(0.31038196) q[0];
rz(-pi) q[1];
rz(2.5035628) q[2];
sx q[2];
rz(-1.3943765) q[2];
sx q[2];
rz(1.0002913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7932574) q[1];
sx q[1];
rz(-2.1264592) q[1];
sx q[1];
rz(-0.72850119) q[1];
x q[2];
rz(0.35238102) q[3];
sx q[3];
rz(-2.0003194) q[3];
sx q[3];
rz(-2.9126715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4124174) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(-2.2632329) q[2];
rz(-0.39204028) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(2.5382606) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6648401) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(0.77366775) q[0];
rz(0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(3.1087648) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87356991) q[0];
sx q[0];
rz(-0.39186726) q[0];
sx q[0];
rz(0.64143945) q[0];
x q[1];
rz(2.122934) q[2];
sx q[2];
rz(-1.3214006) q[2];
sx q[2];
rz(0.92555911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5579538) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(-2.3418531) q[1];
rz(-1.0333943) q[3];
sx q[3];
rz(-2.5279547) q[3];
sx q[3];
rz(0.19526853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(2.8642505) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-0.69916454) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(-0.016419134) q[0];
rz(-0.1044354) q[2];
sx q[2];
rz(-2.0586788) q[2];
sx q[2];
rz(1.3654396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82476393) q[1];
sx q[1];
rz(-2.0889805) q[1];
sx q[1];
rz(1.1512685) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7951489) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(0.88691521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9757441) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(1.7011401) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(0.17000155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6572896) q[0];
sx q[0];
rz(-2.347725) q[0];
sx q[0];
rz(-0.33052175) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0001569) q[2];
sx q[2];
rz(-0.51157198) q[2];
sx q[2];
rz(0.19829743) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.070378455) q[1];
sx q[1];
rz(-2.4134679) q[1];
sx q[1];
rz(2.1945164) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65808987) q[3];
sx q[3];
rz(-0.98539017) q[3];
sx q[3];
rz(-0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99469441) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(-1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(0.68516723) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165669) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(2.2216703) q[0];
rz(-pi) q[1];
rz(2.3643656) q[2];
sx q[2];
rz(-1.6738552) q[2];
sx q[2];
rz(-1.4003786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7428776) q[1];
sx q[1];
rz(-1.0814953) q[1];
sx q[1];
rz(0.87280416) q[1];
x q[2];
rz(0.018304304) q[3];
sx q[3];
rz(-2.1566475) q[3];
sx q[3];
rz(0.5154807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3926065) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.1475295) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80619752) q[0];
sx q[0];
rz(-1.9837539) q[0];
sx q[0];
rz(-1.4892682) q[0];
rz(-2.7717934) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(0.72593867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7123375) q[1];
sx q[1];
rz(-1.0865679) q[1];
sx q[1];
rz(1.6171364) q[1];
rz(-pi) q[2];
x q[2];
rz(1.958193) q[3];
sx q[3];
rz(-1.1068739) q[3];
sx q[3];
rz(-1.0138318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-2.2275887) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99326837) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.3759026) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(-0.39852279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14197037) q[0];
sx q[0];
rz(-2.6451783) q[0];
sx q[0];
rz(1.2312141) q[0];
rz(-2.987791) q[2];
sx q[2];
rz(-0.6066583) q[2];
sx q[2];
rz(0.80468824) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7290013) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(0.97357915) q[1];
rz(-pi) q[2];
rz(1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(-1.8438967) q[2];
rz(2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(-0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(-1.7522316) q[0];
rz(-1.5147491) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(1.0983889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669945) q[0];
sx q[0];
rz(-2.5962788) q[0];
sx q[0];
rz(-1.1101515) q[0];
x q[1];
rz(-0.82010014) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(-1.0320013) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2213649) q[1];
sx q[1];
rz(-2.5623294) q[1];
sx q[1];
rz(3.1406162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.8212574) q[3];
sx q[3];
rz(-2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(-1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(2.675132) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254612) q[0];
sx q[0];
rz(-2.7636508) q[0];
sx q[0];
rz(-1.7037237) q[0];
x q[1];
rz(1.2538818) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(-0.72310477) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.65044636) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(-2.6994929) q[1];
rz(1.915515) q[3];
sx q[3];
rz(-0.88092062) q[3];
sx q[3];
rz(-2.7310839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8978867) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(-2.0824599) q[2];
rz(-0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8582936) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(2.4181096) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(-1.4678636) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
