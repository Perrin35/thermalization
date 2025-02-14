OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(-2.13621) q[0];
sx q[0];
rz(-0.96860743) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(-0.029835116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56108367) q[0];
sx q[0];
rz(-2.5445815) q[0];
sx q[0];
rz(-0.42177864) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.070646324) q[2];
sx q[2];
rz(-2.321817) q[2];
sx q[2];
rz(-1.1653314) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0152143) q[1];
sx q[1];
rz(-2.3087037) q[1];
sx q[1];
rz(-0.53453858) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4814867) q[3];
sx q[3];
rz(-1.4353283) q[3];
sx q[3];
rz(0.40218807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6385103) q[2];
sx q[2];
rz(-0.57380399) q[2];
sx q[2];
rz(-0.32162515) q[2];
rz(-0.62421787) q[3];
sx q[3];
rz(-1.7523242) q[3];
sx q[3];
rz(-1.483884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6931848) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-0.93628991) q[0];
rz(1.2721277) q[1];
sx q[1];
rz(-1.5969758) q[1];
sx q[1];
rz(-1.3880656) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084685582) q[0];
sx q[0];
rz(-1.3515679) q[0];
sx q[0];
rz(3.0830383) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3642367) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(-2.0344025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2876828) q[1];
sx q[1];
rz(-2.1016995) q[1];
sx q[1];
rz(2.4413051) q[1];
rz(-0.8901435) q[3];
sx q[3];
rz(-2.806744) q[3];
sx q[3];
rz(2.072842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0442514) q[2];
sx q[2];
rz(-0.96434957) q[2];
sx q[2];
rz(-0.67330366) q[2];
rz(1.3168969) q[3];
sx q[3];
rz(-1.8484867) q[3];
sx q[3];
rz(0.84002686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8116542) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(0.04976186) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.3866837) q[1];
sx q[1];
rz(-1.6516364) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508036) q[0];
sx q[0];
rz(-0.62701462) q[0];
sx q[0];
rz(-2.9360082) q[0];
rz(-0.65138731) q[2];
sx q[2];
rz(-2.4804401) q[2];
sx q[2];
rz(-2.9686454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4088285) q[1];
sx q[1];
rz(-1.8547426) q[1];
sx q[1];
rz(0.27393053) q[1];
rz(1.1127031) q[3];
sx q[3];
rz(-1.8761523) q[3];
sx q[3];
rz(-0.28766838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9516307) q[2];
sx q[2];
rz(-0.54533521) q[2];
sx q[2];
rz(3.087431) q[2];
rz(-1.9870029) q[3];
sx q[3];
rz(-1.8535987) q[3];
sx q[3];
rz(2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40256777) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(-2.2836852) q[0];
rz(-0.67445451) q[1];
sx q[1];
rz(-0.94991389) q[1];
sx q[1];
rz(1.3309006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20786634) q[0];
sx q[0];
rz(-0.27369341) q[0];
sx q[0];
rz(0.52723186) q[0];
rz(-pi) q[1];
rz(1.0687549) q[2];
sx q[2];
rz(-2.7981361) q[2];
sx q[2];
rz(2.1844027) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96283276) q[1];
sx q[1];
rz(-0.58738999) q[1];
sx q[1];
rz(1.2580832) q[1];
x q[2];
rz(0.0033802948) q[3];
sx q[3];
rz(-1.9900609) q[3];
sx q[3];
rz(1.5730891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3746419) q[2];
sx q[2];
rz(-2.9810814) q[2];
sx q[2];
rz(3.0373419) q[2];
rz(0.57779622) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(2.2195063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.6660026) q[0];
sx q[0];
rz(2.1971028) q[0];
rz(-1.7780444) q[1];
sx q[1];
rz(-1.8133769) q[1];
sx q[1];
rz(1.69467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76857012) q[0];
sx q[0];
rz(-1.3813586) q[0];
sx q[0];
rz(-0.85590881) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4291006) q[2];
sx q[2];
rz(-0.25280127) q[2];
sx q[2];
rz(1.6293874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8268298) q[1];
sx q[1];
rz(-2.3896425) q[1];
sx q[1];
rz(0.58677267) q[1];
rz(-pi) q[2];
rz(-0.99926104) q[3];
sx q[3];
rz(-2.0700698) q[3];
sx q[3];
rz(1.5972135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1030582) q[2];
sx q[2];
rz(-0.47163042) q[2];
sx q[2];
rz(0.72059694) q[2];
rz(3.0125812) q[3];
sx q[3];
rz(-1.7804264) q[3];
sx q[3];
rz(-1.2989929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1593889) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(0.61656117) q[0];
rz(2.189134) q[1];
sx q[1];
rz(-1.9472232) q[1];
sx q[1];
rz(-1.8985101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12981249) q[0];
sx q[0];
rz(-1.3115378) q[0];
sx q[0];
rz(-3.0002322) q[0];
rz(-pi) q[1];
rz(-0.86402969) q[2];
sx q[2];
rz(-1.8908769) q[2];
sx q[2];
rz(-0.79369407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3042708) q[1];
sx q[1];
rz(-2.4189286) q[1];
sx q[1];
rz(0.091462084) q[1];
rz(-pi) q[2];
x q[2];
rz(0.074148103) q[3];
sx q[3];
rz(-0.95063248) q[3];
sx q[3];
rz(2.1886203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9055966) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(3.0671885) q[2];
rz(-2.060804) q[3];
sx q[3];
rz(-1.6483043) q[3];
sx q[3];
rz(1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3022795) q[0];
sx q[0];
rz(-1.1458719) q[0];
sx q[0];
rz(-3.0795414) q[0];
rz(-1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(-3.0516023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71848559) q[0];
sx q[0];
rz(-0.88750091) q[0];
sx q[0];
rz(2.8124269) q[0];
x q[1];
rz(-1.3700468) q[2];
sx q[2];
rz(-2.4436969) q[2];
sx q[2];
rz(3.0332886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1931928) q[1];
sx q[1];
rz(-0.71454778) q[1];
sx q[1];
rz(2.1195861) q[1];
rz(-pi) q[2];
rz(-0.84109938) q[3];
sx q[3];
rz(-0.78783082) q[3];
sx q[3];
rz(0.5238409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0061079582) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(-2.4583859) q[2];
rz(-1.3225383) q[3];
sx q[3];
rz(-2.1893978) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.6739864) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(-2.6737387) q[0];
rz(0.6894919) q[1];
sx q[1];
rz(-0.71527022) q[1];
sx q[1];
rz(-0.47324866) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5734542) q[0];
sx q[0];
rz(-1.7466508) q[0];
sx q[0];
rz(3.0699845) q[0];
rz(-pi) q[1];
rz(0.44395776) q[2];
sx q[2];
rz(-1.4251764) q[2];
sx q[2];
rz(0.99010003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93750226) q[1];
sx q[1];
rz(-0.023462208) q[1];
sx q[1];
rz(-0.84659441) q[1];
rz(-2.8611818) q[3];
sx q[3];
rz(-1.3254998) q[3];
sx q[3];
rz(0.51285686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0143934) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(0.42352208) q[2];
rz(-1.315377) q[3];
sx q[3];
rz(-2.6500621) q[3];
sx q[3];
rz(-1.8486842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2143283) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(-1.3182053) q[1];
sx q[1];
rz(-0.91000906) q[1];
sx q[1];
rz(3.1277025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5197138) q[0];
sx q[0];
rz(-1.5962999) q[0];
sx q[0];
rz(-1.9344781) q[0];
x q[1];
rz(-0.40851302) q[2];
sx q[2];
rz(-2.7613104) q[2];
sx q[2];
rz(2.399596) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74456066) q[1];
sx q[1];
rz(-0.63493997) q[1];
sx q[1];
rz(-1.145361) q[1];
x q[2];
rz(0.12356476) q[3];
sx q[3];
rz(-2.2964604) q[3];
sx q[3];
rz(0.35620505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5615329) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(2.0884936) q[2];
rz(-2.2140391) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(2.659306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574361) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(0.60047737) q[0];
rz(-0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(-1.3788266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10256448) q[0];
sx q[0];
rz(-1.2192654) q[0];
sx q[0];
rz(2.8536056) q[0];
rz(-2.4587833) q[2];
sx q[2];
rz(-0.46418846) q[2];
sx q[2];
rz(-0.95303553) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8178508) q[1];
sx q[1];
rz(-1.3881589) q[1];
sx q[1];
rz(2.707213) q[1];
rz(0.61609488) q[3];
sx q[3];
rz(-1.9368251) q[3];
sx q[3];
rz(0.6835365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.025734162) q[2];
sx q[2];
rz(-0.33418843) q[2];
sx q[2];
rz(0.96246976) q[2];
rz(0.95364237) q[3];
sx q[3];
rz(-2.0680659) q[3];
sx q[3];
rz(-1.3941221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3772603) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(1.3537021) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(0.66966343) q[2];
sx q[2];
rz(-1.0367994) q[2];
sx q[2];
rz(-1.8112469) q[2];
rz(-1.7103473) q[3];
sx q[3];
rz(-2.0736251) q[3];
sx q[3];
rz(0.74794378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
