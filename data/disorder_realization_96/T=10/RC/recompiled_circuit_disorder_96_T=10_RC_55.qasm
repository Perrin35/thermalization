OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49657208) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(-2.7215331) q[0];
rz(3.009216) q[2];
sx q[2];
rz(-1.0356324) q[2];
sx q[2];
rz(-1.7420499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.28509724) q[1];
sx q[1];
rz(-2.3772117) q[1];
sx q[1];
rz(0.38715036) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31906268) q[3];
sx q[3];
rz(-0.82108077) q[3];
sx q[3];
rz(1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(3.0875207) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0756404) q[0];
sx q[0];
rz(-1.4464805) q[0];
sx q[0];
rz(0.9888222) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29157721) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(-2.4183194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9745001) q[1];
sx q[1];
rz(-1.8247316) q[1];
sx q[1];
rz(2.7212935) q[1];
x q[2];
rz(1.4278533) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-3.1157852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255071) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(0.020629701) q[0];
rz(-pi) q[1];
rz(1.8823207) q[2];
sx q[2];
rz(-1.4853962) q[2];
sx q[2];
rz(2.4740919) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3734696) q[1];
sx q[1];
rz(-0.79263955) q[1];
sx q[1];
rz(-2.5440689) q[1];
x q[2];
rz(0.13173007) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6484084) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(-1.0113082) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6279814) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(2.016071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0828447) q[1];
sx q[1];
rz(-1.9472329) q[1];
sx q[1];
rz(0.67614268) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67646497) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(1.7959309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2252786) q[0];
sx q[0];
rz(-1.7597223) q[0];
sx q[0];
rz(-2.180045) q[0];
x q[1];
rz(2.8565065) q[2];
sx q[2];
rz(-1.0694155) q[2];
sx q[2];
rz(-1.8269055) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4700714) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(-0.77312153) q[1];
x q[2];
rz(2.352036) q[3];
sx q[3];
rz(-1.4649179) q[3];
sx q[3];
rz(1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(-2.7894003) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6569865) q[0];
sx q[0];
rz(-1.779042) q[0];
sx q[0];
rz(-0.44021846) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52144737) q[2];
sx q[2];
rz(-1.2505184) q[2];
sx q[2];
rz(-2.5495868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2211654) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(1.2616874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3771463) q[3];
sx q[3];
rz(-0.87246694) q[3];
sx q[3];
rz(-2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(-1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-2.2999433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8182897) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(-1.3099758) q[0];
rz(0.26288962) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(0.4292683) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6517666) q[1];
sx q[1];
rz(-3.0511599) q[1];
sx q[1];
rz(-1.0033146) q[1];
rz(-pi) q[2];
rz(-1.8089201) q[3];
sx q[3];
rz(-0.75546414) q[3];
sx q[3];
rz(-1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-3.0916396) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-2.9522827) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4293489) q[0];
sx q[0];
rz(-1.1375543) q[0];
sx q[0];
rz(-2.594125) q[0];
rz(-pi) q[1];
rz(-1.2049963) q[2];
sx q[2];
rz(-1.4636283) q[2];
sx q[2];
rz(2.0049713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8466134) q[1];
sx q[1];
rz(-2.0366922) q[1];
sx q[1];
rz(-2.2094775) q[1];
rz(-pi) q[2];
rz(1.9440218) q[3];
sx q[3];
rz(-1.8456568) q[3];
sx q[3];
rz(-1.9762135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-0.94690698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98012053) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(1.253771) q[0];
rz(-pi) q[1];
rz(-1.5309179) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(-3.1181042) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7864101) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(1.4817609) q[1];
x q[2];
rz(-1.1831207) q[3];
sx q[3];
rz(-2.6689853) q[3];
sx q[3];
rz(-1.2848867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(-2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-0.64430976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(0.17382646) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23004736) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(-0.41781296) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5035142) q[1];
sx q[1];
rz(-2.6942309) q[1];
sx q[1];
rz(2.4380986) q[1];
rz(-pi) q[2];
rz(0.44390042) q[3];
sx q[3];
rz(-2.0022941) q[3];
sx q[3];
rz(-2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(-2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29466378) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-1.3953801) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(-0.39137822) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];