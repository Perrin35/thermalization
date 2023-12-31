OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(-2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5380733) q[0];
sx q[0];
rz(-0.78342122) q[0];
sx q[0];
rz(2.6253683) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4835303) q[2];
sx q[2];
rz(-2.1076638) q[2];
sx q[2];
rz(-1.3047578) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4989657) q[1];
sx q[1];
rz(-1.9470012) q[1];
sx q[1];
rz(2.7514003) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3124247) q[3];
sx q[3];
rz(-1.4989304) q[3];
sx q[3];
rz(1.3952599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26596507) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(-1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(2.6665376) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(1.0377201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1514725) q[0];
sx q[0];
rz(-1.3511786) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-1.0449045) q[2];
sx q[2];
rz(-1.0962588) q[2];
sx q[2];
rz(-1.431312) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2983919) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(-3.1133679) q[1];
x q[2];
rz(-0.59229895) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4804046) q[2];
sx q[2];
rz(-1.8030689) q[2];
sx q[2];
rz(0.084687106) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-0.27733222) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-0.95570046) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(2.5684165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5983551) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(1.4198562) q[0];
x q[1];
rz(1.8884044) q[2];
sx q[2];
rz(-2.3177958) q[2];
sx q[2];
rz(-0.93712805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55858559) q[1];
sx q[1];
rz(-1.188394) q[1];
sx q[1];
rz(-1.8797727) q[1];
rz(2.0315941) q[3];
sx q[3];
rz(-1.2013544) q[3];
sx q[3];
rz(0.24584578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1406143) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8594584) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-2.5909246) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(2.7788924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3804647) q[0];
sx q[0];
rz(-1.8109545) q[0];
sx q[0];
rz(2.2855177) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5393125) q[2];
sx q[2];
rz(-1.5523947) q[2];
sx q[2];
rz(-0.035426332) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6994233) q[1];
sx q[1];
rz(-1.6740834) q[1];
sx q[1];
rz(0.63421722) q[1];
rz(0.94621559) q[3];
sx q[3];
rz(-1.6295625) q[3];
sx q[3];
rz(-2.0277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-0.0030227946) q[2];
rz(0.65888843) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(0.80250424) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(-1.4594706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0844903) q[0];
sx q[0];
rz(-1.717289) q[0];
sx q[0];
rz(1.811972) q[0];
rz(2.2541788) q[2];
sx q[2];
rz(-1.1537342) q[2];
sx q[2];
rz(-2.8487157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7534605) q[1];
sx q[1];
rz(-1.3452483) q[1];
sx q[1];
rz(1.3688341) q[1];
rz(2.1305389) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(-1.2572446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83546272) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275948) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(1.9110514) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-3.0029283) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9753871) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-3.1211981) q[0];
rz(-pi) q[1];
rz(1.4939098) q[2];
sx q[2];
rz(-1.7711519) q[2];
sx q[2];
rz(0.69918699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87673346) q[1];
sx q[1];
rz(-2.2899592) q[1];
sx q[1];
rz(2.369146) q[1];
rz(-pi) q[2];
rz(0.45913978) q[3];
sx q[3];
rz(-2.3448179) q[3];
sx q[3];
rz(-1.9285942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3665294) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5360864) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-0.53652525) q[0];
rz(2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(2.5172863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7621967) q[0];
sx q[0];
rz(-1.7805829) q[0];
sx q[0];
rz(2.4136153) q[0];
x q[1];
rz(-2.1937218) q[2];
sx q[2];
rz(-1.9588753) q[2];
sx q[2];
rz(0.58448234) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.045947) q[1];
sx q[1];
rz(-2.1150757) q[1];
sx q[1];
rz(-0.73927684) q[1];
rz(-pi) q[2];
rz(2.1950486) q[3];
sx q[3];
rz(-2.4157899) q[3];
sx q[3];
rz(0.040369999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(0.13993046) q[0];
rz(1.6775999) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-3.0126742) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3710204) q[0];
sx q[0];
rz(-0.41045529) q[0];
sx q[0];
rz(-1.8380941) q[0];
rz(2.2419937) q[2];
sx q[2];
rz(-1.1102144) q[2];
sx q[2];
rz(-2.7923498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.376437) q[1];
sx q[1];
rz(-0.51424485) q[1];
sx q[1];
rz(-2.5625485) q[1];
x q[2];
rz(-0.017459004) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(-0.12714126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(2.5174985) q[2];
rz(0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.8918442) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6791145) q[0];
sx q[0];
rz(-2.9318641) q[0];
sx q[0];
rz(1.0366584) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5596041) q[2];
sx q[2];
rz(-0.66170035) q[2];
sx q[2];
rz(2.2616507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13547922) q[1];
sx q[1];
rz(-1.0891501) q[1];
sx q[1];
rz(2.3856132) q[1];
rz(-pi) q[2];
rz(1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-0.11432153) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(-0.21324883) q[0];
rz(0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0471668) q[0];
sx q[0];
rz(-1.3494274) q[0];
sx q[0];
rz(0.84035994) q[0];
rz(0.28995138) q[2];
sx q[2];
rz(-2.6486514) q[2];
sx q[2];
rz(1.314756) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(1.6092369) q[1];
rz(1.7095079) q[3];
sx q[3];
rz(-1.7943534) q[3];
sx q[3];
rz(0.15299882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0032349) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(-0.004301087) q[2];
rz(-2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(-0.44395631) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(1.2658723) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(1.0989582) q[3];
sx q[3];
rz(-1.2755339) q[3];
sx q[3];
rz(-0.76225029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
