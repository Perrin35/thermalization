OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27788568) q[0];
sx q[0];
rz(3.570896) q[0];
sx q[0];
rz(12.255393) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(-2.0452979) q[1];
sx q[1];
rz(-0.57408339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1659083) q[0];
sx q[0];
rz(-2.568779) q[0];
sx q[0];
rz(0.84337916) q[0];
rz(-pi) q[1];
rz(-0.034717807) q[2];
sx q[2];
rz(-2.1534734) q[2];
sx q[2];
rz(-2.5947844) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.61389) q[1];
sx q[1];
rz(-0.37606323) q[1];
sx q[1];
rz(1.8164053) q[1];
rz(-pi) q[2];
rz(-1.4101683) q[3];
sx q[3];
rz(-2.7607252) q[3];
sx q[3];
rz(-2.4695652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8671888) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(1.3192419) q[2];
rz(3.0684209) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(0.81899548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0935593) q[0];
sx q[0];
rz(-1.1807384) q[0];
sx q[0];
rz(2.4617885) q[0];
rz(0.80348429) q[1];
sx q[1];
rz(-1.3619224) q[1];
sx q[1];
rz(0.63978535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6006192) q[0];
sx q[0];
rz(-2.2007211) q[0];
sx q[0];
rz(-2.7472904) q[0];
rz(-pi) q[1];
rz(-3.0306973) q[2];
sx q[2];
rz(-2.408658) q[2];
sx q[2];
rz(-1.1711677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3181139) q[1];
sx q[1];
rz(-2.1770855) q[1];
sx q[1];
rz(0.73474523) q[1];
x q[2];
rz(-1.7946258) q[3];
sx q[3];
rz(-1.8719684) q[3];
sx q[3];
rz(-2.6909242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5378319) q[2];
sx q[2];
rz(-2.64309) q[2];
sx q[2];
rz(0.70408386) q[2];
rz(3.0294561) q[3];
sx q[3];
rz(-1.6928558) q[3];
sx q[3];
rz(1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0391487) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(2.802134) q[0];
rz(-2.0231694) q[1];
sx q[1];
rz(-0.92670852) q[1];
sx q[1];
rz(-0.035645398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729926) q[0];
sx q[0];
rz(-0.22449271) q[0];
sx q[0];
rz(-0.39546449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8730803) q[2];
sx q[2];
rz(-1.2460104) q[2];
sx q[2];
rz(-0.12476441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0782203) q[1];
sx q[1];
rz(-2.2148892) q[1];
sx q[1];
rz(-1.8404585) q[1];
rz(-pi) q[2];
rz(0.44575739) q[3];
sx q[3];
rz(-0.66346079) q[3];
sx q[3];
rz(-1.9395246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3493335) q[2];
sx q[2];
rz(-0.68816853) q[2];
sx q[2];
rz(-2.2815857) q[2];
rz(2.3490014) q[3];
sx q[3];
rz(-0.20321295) q[3];
sx q[3];
rz(1.0205166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-1.0122796) q[0];
sx q[0];
rz(-1.7565933) q[0];
sx q[0];
rz(2.5087575) q[0];
rz(1.1526147) q[1];
sx q[1];
rz(-0.8747789) q[1];
sx q[1];
rz(-1.6713743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36258478) q[0];
sx q[0];
rz(-0.92155313) q[0];
sx q[0];
rz(3.1253689) q[0];
rz(-pi) q[1];
rz(2.987902) q[2];
sx q[2];
rz(-0.77506232) q[2];
sx q[2];
rz(2.3714921) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1296334) q[1];
sx q[1];
rz(-2.2599038) q[1];
sx q[1];
rz(2.8065475) q[1];
rz(-pi) q[2];
rz(-0.65514647) q[3];
sx q[3];
rz(-1.0628502) q[3];
sx q[3];
rz(2.9367713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4319438) q[2];
sx q[2];
rz(-2.180763) q[2];
sx q[2];
rz(-2.7107837) q[2];
rz(2.4591947) q[3];
sx q[3];
rz(-1.6513446) q[3];
sx q[3];
rz(0.94990134) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7948941) q[0];
sx q[0];
rz(-1.2656724) q[0];
sx q[0];
rz(1.9516113) q[0];
rz(0.31040141) q[1];
sx q[1];
rz(-1.0810532) q[1];
sx q[1];
rz(-2.6365872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520262) q[0];
sx q[0];
rz(-2.1329693) q[0];
sx q[0];
rz(-0.52345353) q[0];
rz(-pi) q[1];
rz(3.0104899) q[2];
sx q[2];
rz(-1.982568) q[2];
sx q[2];
rz(0.51136651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4601535) q[1];
sx q[1];
rz(-2.6362754) q[1];
sx q[1];
rz(2.9605335) q[1];
rz(-pi) q[2];
rz(-2.4159031) q[3];
sx q[3];
rz(-1.5066875) q[3];
sx q[3];
rz(0.025246092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7348822) q[2];
sx q[2];
rz(-2.4411185) q[2];
sx q[2];
rz(-2.238838) q[2];
rz(-1.0550176) q[3];
sx q[3];
rz(-2.2746634) q[3];
sx q[3];
rz(-2.6796851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91586739) q[0];
sx q[0];
rz(-2.4287455) q[0];
sx q[0];
rz(2.0676887) q[0];
rz(-1.7621) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(0.097537907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65501761) q[0];
sx q[0];
rz(-1.7518105) q[0];
sx q[0];
rz(-0.54365309) q[0];
rz(-pi) q[1];
rz(-0.16867192) q[2];
sx q[2];
rz(-1.9430117) q[2];
sx q[2];
rz(-1.8255359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40332023) q[1];
sx q[1];
rz(-3.0636859) q[1];
sx q[1];
rz(-2.2408443) q[1];
rz(-pi) q[2];
rz(-1.0076526) q[3];
sx q[3];
rz(-2.030367) q[3];
sx q[3];
rz(0.28467049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11233106) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(-0.41109273) q[2];
rz(1.4136275) q[3];
sx q[3];
rz(-2.3634383) q[3];
sx q[3];
rz(-1.5467862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77808648) q[0];
sx q[0];
rz(-0.030487617) q[0];
sx q[0];
rz(1.7331069) q[0];
rz(-1.015649) q[1];
sx q[1];
rz(-1.3280832) q[1];
sx q[1];
rz(-1.4782864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7789981) q[0];
sx q[0];
rz(-1.2752879) q[0];
sx q[0];
rz(-0.63531718) q[0];
rz(2.4934216) q[2];
sx q[2];
rz(-1.0888466) q[2];
sx q[2];
rz(-0.75343695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95132213) q[1];
sx q[1];
rz(-0.45004216) q[1];
sx q[1];
rz(1.7461036) q[1];
rz(-0.13161259) q[3];
sx q[3];
rz(-1.0743047) q[3];
sx q[3];
rz(-0.59959664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0279072) q[2];
sx q[2];
rz(-0.62782851) q[2];
sx q[2];
rz(-0.17024635) q[2];
rz(1.8611192) q[3];
sx q[3];
rz(-1.4743285) q[3];
sx q[3];
rz(-1.9835651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.668642) q[0];
sx q[0];
rz(-0.96422115) q[0];
sx q[0];
rz(2.6131795) q[0];
rz(-0.93327418) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(2.5859213) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1707204) q[0];
sx q[0];
rz(-1.375388) q[0];
sx q[0];
rz(-2.1989345) q[0];
rz(-pi) q[1];
rz(-3.0425983) q[2];
sx q[2];
rz(-2.3470391) q[2];
sx q[2];
rz(-1.0459448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0289149) q[1];
sx q[1];
rz(-2.2820447) q[1];
sx q[1];
rz(1.5284786) q[1];
rz(2.4497568) q[3];
sx q[3];
rz(-1.8134549) q[3];
sx q[3];
rz(1.7586046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8133424) q[2];
sx q[2];
rz(-0.80208653) q[2];
sx q[2];
rz(-2.8343406) q[2];
rz(-0.79636374) q[3];
sx q[3];
rz(-0.73085228) q[3];
sx q[3];
rz(0.097631924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873213) q[0];
sx q[0];
rz(-1.048943) q[0];
sx q[0];
rz(1.3294504) q[0];
rz(0.21996552) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(-1.3185917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47596915) q[0];
sx q[0];
rz(-2.7950056) q[0];
sx q[0];
rz(-2.2524538) q[0];
rz(-1.6873785) q[2];
sx q[2];
rz(-2.0676842) q[2];
sx q[2];
rz(-2.1379444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4107549) q[1];
sx q[1];
rz(-1.2304092) q[1];
sx q[1];
rz(1.0761258) q[1];
rz(2.044146) q[3];
sx q[3];
rz(-0.90765491) q[3];
sx q[3];
rz(-1.1569661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69507504) q[2];
sx q[2];
rz(-1.7138285) q[2];
sx q[2];
rz(1.8565149) q[2];
rz(2.255693) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(1.6133962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000243) q[0];
sx q[0];
rz(-2.4055241) q[0];
sx q[0];
rz(-0.31684434) q[0];
rz(0.44081229) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(0.37030927) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.646831) q[0];
sx q[0];
rz(-1.5398105) q[0];
sx q[0];
rz(1.9520743) q[0];
rz(-pi) q[1];
rz(0.84443386) q[2];
sx q[2];
rz(-1.053868) q[2];
sx q[2];
rz(-2.5518565) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.036344254) q[1];
sx q[1];
rz(-0.40648983) q[1];
sx q[1];
rz(2.2986733) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0256151) q[3];
sx q[3];
rz(-0.11546455) q[3];
sx q[3];
rz(-2.786557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73654282) q[2];
sx q[2];
rz(-2.6506347) q[2];
sx q[2];
rz(-1.5092108) q[2];
rz(-2.3368808) q[3];
sx q[3];
rz(-1.6377621) q[3];
sx q[3];
rz(0.10778431) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0263154) q[0];
sx q[0];
rz(-1.4125217) q[0];
sx q[0];
rz(-1.9052292) q[0];
rz(0.23030494) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(0.7028107) q[2];
sx q[2];
rz(-0.79844777) q[2];
sx q[2];
rz(-2.993882) q[2];
rz(-2.1896918) q[3];
sx q[3];
rz(-1.3018082) q[3];
sx q[3];
rz(-0.34886532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
