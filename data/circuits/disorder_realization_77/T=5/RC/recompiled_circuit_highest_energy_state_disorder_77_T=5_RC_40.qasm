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
rz(-1.0722941) q[0];
sx q[0];
rz(-0.18360734) q[0];
sx q[0];
rz(2.973383) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(1.1163196) q[1];
sx q[1];
rz(7.4012227) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.428225) q[0];
sx q[0];
rz(-1.2452188) q[0];
sx q[0];
rz(-1.6859562) q[0];
rz(-pi) q[1];
rz(-0.69324915) q[2];
sx q[2];
rz(-1.8680671) q[2];
sx q[2];
rz(1.3849486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.11169) q[1];
sx q[1];
rz(-2.3808834) q[1];
sx q[1];
rz(-3.1055443) q[1];
rz(-2.0327031) q[3];
sx q[3];
rz(-2.8402764) q[3];
sx q[3];
rz(1.0404943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2844124) q[2];
sx q[2];
rz(-2.1511011) q[2];
sx q[2];
rz(1.277479) q[2];
rz(-2.5340581) q[3];
sx q[3];
rz(-1.7888864) q[3];
sx q[3];
rz(-1.7551306) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4472537) q[0];
sx q[0];
rz(-2.4632833) q[0];
sx q[0];
rz(0.69764486) q[0];
rz(0.33277008) q[1];
sx q[1];
rz(-2.0336626) q[1];
sx q[1];
rz(2.8093801) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282917) q[0];
sx q[0];
rz(-1.6915404) q[0];
sx q[0];
rz(3.0101454) q[0];
rz(-1.9273913) q[2];
sx q[2];
rz(-2.5348672) q[2];
sx q[2];
rz(-1.0068829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1134791) q[1];
sx q[1];
rz(-1.81908) q[1];
sx q[1];
rz(-1.1400272) q[1];
rz(-pi) q[2];
rz(-1.9929797) q[3];
sx q[3];
rz(-2.2530972) q[3];
sx q[3];
rz(2.9158031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.701) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(-2.8458703) q[2];
rz(0.010585636) q[3];
sx q[3];
rz(-2.1518095) q[3];
sx q[3];
rz(0.73339614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877141) q[0];
sx q[0];
rz(-2.6831477) q[0];
sx q[0];
rz(0.53471765) q[0];
rz(1.1010822) q[1];
sx q[1];
rz(-1.705876) q[1];
sx q[1];
rz(-0.91317493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3252439) q[0];
sx q[0];
rz(-1.2486098) q[0];
sx q[0];
rz(-1.0464099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1276563) q[2];
sx q[2];
rz(-1.3512313) q[2];
sx q[2];
rz(2.1923421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.53131858) q[1];
sx q[1];
rz(-0.90435678) q[1];
sx q[1];
rz(0.6998723) q[1];
x q[2];
rz(1.6524773) q[3];
sx q[3];
rz(-0.73058999) q[3];
sx q[3];
rz(-3.0944892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59213006) q[2];
sx q[2];
rz(-0.89207804) q[2];
sx q[2];
rz(-2.7590052) q[2];
rz(1.4866746) q[3];
sx q[3];
rz(-0.29522172) q[3];
sx q[3];
rz(-2.2216643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53052467) q[0];
sx q[0];
rz(-1.684573) q[0];
sx q[0];
rz(-1.0462421) q[0];
rz(-2.5776082) q[1];
sx q[1];
rz(-1.0020741) q[1];
sx q[1];
rz(-1.8246338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038669) q[0];
sx q[0];
rz(-0.90529862) q[0];
sx q[0];
rz(0.72030495) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91805075) q[2];
sx q[2];
rz(-1.8260644) q[2];
sx q[2];
rz(1.7328615) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35933381) q[1];
sx q[1];
rz(-2.8247807) q[1];
sx q[1];
rz(0.13765361) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.149156) q[3];
sx q[3];
rz(-0.8332592) q[3];
sx q[3];
rz(-1.6374388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6947393) q[2];
sx q[2];
rz(-1.6971735) q[2];
sx q[2];
rz(2.1430338) q[2];
rz(-0.57991943) q[3];
sx q[3];
rz(-1.4812352) q[3];
sx q[3];
rz(-1.3045605) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9272083) q[0];
sx q[0];
rz(-1.2283607) q[0];
sx q[0];
rz(2.4217915) q[0];
rz(0.7431227) q[1];
sx q[1];
rz(-1.2179951) q[1];
sx q[1];
rz(-3.0845089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8419619) q[0];
sx q[0];
rz(-1.7052545) q[0];
sx q[0];
rz(2.7706783) q[0];
x q[1];
rz(2.7774335) q[2];
sx q[2];
rz(-0.79150598) q[2];
sx q[2];
rz(3.0730545) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.067043153) q[1];
sx q[1];
rz(-0.32575575) q[1];
sx q[1];
rz(0.52149741) q[1];
rz(-pi) q[2];
rz(-2.1716398) q[3];
sx q[3];
rz(-0.68884736) q[3];
sx q[3];
rz(-0.057202489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4382978) q[2];
sx q[2];
rz(-1.1763828) q[2];
sx q[2];
rz(-1.8717742) q[2];
rz(-2.5979009) q[3];
sx q[3];
rz(-0.68259493) q[3];
sx q[3];
rz(-1.3683176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.066034) q[0];
sx q[0];
rz(-0.37519255) q[0];
sx q[0];
rz(-0.22848465) q[0];
rz(-2.1560419) q[1];
sx q[1];
rz(-1.5551714) q[1];
sx q[1];
rz(1.9353345) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.120283) q[0];
sx q[0];
rz(-2.7643959) q[0];
sx q[0];
rz(-0.46964236) q[0];
rz(-1.2689077) q[2];
sx q[2];
rz(-1.3573555) q[2];
sx q[2];
rz(-0.53415307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4691733) q[1];
sx q[1];
rz(-0.9658199) q[1];
sx q[1];
rz(-2.9939326) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4643741) q[3];
sx q[3];
rz(-1.0297349) q[3];
sx q[3];
rz(0.18782561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1365405) q[2];
sx q[2];
rz(-1.8661934) q[2];
sx q[2];
rz(-1.2460111) q[2];
rz(-0.15240845) q[3];
sx q[3];
rz(-2.992575) q[3];
sx q[3];
rz(-2.3685031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8361255) q[0];
sx q[0];
rz(-1.2141328) q[0];
sx q[0];
rz(-2.8866696) q[0];
rz(2.1615324) q[1];
sx q[1];
rz(-1.3420811) q[1];
sx q[1];
rz(2.9068388) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4027883) q[0];
sx q[0];
rz(-2.1128034) q[0];
sx q[0];
rz(-2.8554585) q[0];
x q[1];
rz(2.9696941) q[2];
sx q[2];
rz(-1.8020639) q[2];
sx q[2];
rz(-2.8794546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67803176) q[1];
sx q[1];
rz(-1.1593736) q[1];
sx q[1];
rz(-1.2829343) q[1];
rz(-pi) q[2];
rz(-1.717155) q[3];
sx q[3];
rz(-0.68181935) q[3];
sx q[3];
rz(0.54659971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0945101) q[2];
sx q[2];
rz(-2.7213056) q[2];
sx q[2];
rz(-1.1959929) q[2];
rz(0.99810537) q[3];
sx q[3];
rz(-1.866303) q[3];
sx q[3];
rz(0.78817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.81414476) q[0];
sx q[0];
rz(-0.87004167) q[0];
sx q[0];
rz(2.5092464) q[0];
rz(1.5517722) q[1];
sx q[1];
rz(-1.7951782) q[1];
sx q[1];
rz(-1.8428892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0576026) q[0];
sx q[0];
rz(-0.51483908) q[0];
sx q[0];
rz(2.0979416) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.024237337) q[2];
sx q[2];
rz(-2.4757984) q[2];
sx q[2];
rz(2.0074304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82567333) q[1];
sx q[1];
rz(-2.6330067) q[1];
sx q[1];
rz(0.085476106) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15141204) q[3];
sx q[3];
rz(-1.5431343) q[3];
sx q[3];
rz(0.20229483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87045264) q[2];
sx q[2];
rz(-2.2282034) q[2];
sx q[2];
rz(-0.75922981) q[2];
rz(-0.40902725) q[3];
sx q[3];
rz(-1.7274011) q[3];
sx q[3];
rz(0.11212382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8522393) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(0.36561832) q[0];
rz(-0.002937142) q[1];
sx q[1];
rz(-1.5298839) q[1];
sx q[1];
rz(1.0681577) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6474175) q[0];
sx q[0];
rz(-1.023698) q[0];
sx q[0];
rz(-1.2450594) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67068213) q[2];
sx q[2];
rz(-0.95206958) q[2];
sx q[2];
rz(-2.3681792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6456071) q[1];
sx q[1];
rz(-1.9580972) q[1];
sx q[1];
rz(3.1003968) q[1];
x q[2];
rz(3.0084988) q[3];
sx q[3];
rz(-0.20340098) q[3];
sx q[3];
rz(-1.4730367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3785765) q[2];
sx q[2];
rz(-0.91332674) q[2];
sx q[2];
rz(-0.98706377) q[2];
rz(-0.8118363) q[3];
sx q[3];
rz(-0.87849179) q[3];
sx q[3];
rz(-1.2442376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064747485) q[0];
sx q[0];
rz(-2.5223314) q[0];
sx q[0];
rz(1.6774119) q[0];
rz(0.96425104) q[1];
sx q[1];
rz(-1.827927) q[1];
sx q[1];
rz(1.7589232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2428873) q[0];
sx q[0];
rz(-0.90685493) q[0];
sx q[0];
rz(-2.6580145) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5141508) q[2];
sx q[2];
rz(-0.44763264) q[2];
sx q[2];
rz(1.0120376) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23338887) q[1];
sx q[1];
rz(-1.0393794) q[1];
sx q[1];
rz(-2.6446093) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2887705) q[3];
sx q[3];
rz(-2.3896165) q[3];
sx q[3];
rz(0.2096006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1796403) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(1.4947653) q[2];
rz(2.0252114) q[3];
sx q[3];
rz(-2.3514082) q[3];
sx q[3];
rz(-2.5377929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2434926) q[0];
sx q[0];
rz(-1.6742764) q[0];
sx q[0];
rz(1.8754638) q[0];
rz(2.948214) q[1];
sx q[1];
rz(-1.0693751) q[1];
sx q[1];
rz(1.3042915) q[1];
rz(0.73430772) q[2];
sx q[2];
rz(-1.3398017) q[2];
sx q[2];
rz(-0.2504713) q[2];
rz(-1.1361271) q[3];
sx q[3];
rz(-2.2129314) q[3];
sx q[3];
rz(1.5066838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
