OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(-0.34531265) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(3.4808573) q[1];
sx q[1];
rz(9.1453287) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61235185) q[0];
sx q[0];
rz(-1.2623598) q[0];
sx q[0];
rz(0.2150857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9136393) q[2];
sx q[2];
rz(-2.3228085) q[2];
sx q[2];
rz(2.0057099) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2742259) q[1];
sx q[1];
rz(-0.19051954) q[1];
sx q[1];
rz(-1.6412559) q[1];
rz(-pi) q[2];
rz(1.6104326) q[3];
sx q[3];
rz(-2.2124024) q[3];
sx q[3];
rz(-1.8455781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(-0.41369307) q[0];
rz(1.3445688) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-0.63562524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7660852) q[0];
sx q[0];
rz(-1.5264395) q[0];
sx q[0];
rz(-1.4384598) q[0];
x q[1];
rz(-0.046819709) q[2];
sx q[2];
rz(-1.702582) q[2];
sx q[2];
rz(1.5891967) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0932514) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(-1.3202207) q[1];
x q[2];
rz(0.91931822) q[3];
sx q[3];
rz(-2.161536) q[3];
sx q[3];
rz(-0.35349333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.7522493) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(-2.6129369) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5813331) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(-0.43446508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70179825) q[0];
sx q[0];
rz(-0.5895624) q[0];
sx q[0];
rz(-0.88296417) q[0];
rz(0.065504727) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(2.5158109) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3420521) q[1];
sx q[1];
rz(-0.5040579) q[1];
sx q[1];
rz(0.44517681) q[1];
rz(-0.64819077) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(-1.0685514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(3.0857962) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(0.14053024) q[0];
rz(-2.9699516) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(0.26352873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22334252) q[0];
sx q[0];
rz(-1.6726057) q[0];
sx q[0];
rz(2.5725767) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0657004) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(-0.87841735) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9159689) q[1];
sx q[1];
rz(-1.2265424) q[1];
sx q[1];
rz(-2.4179789) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9043546) q[3];
sx q[3];
rz(-1.3106489) q[3];
sx q[3];
rz(-3.1317657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.2496703) q[2];
rz(1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12954457) q[0];
sx q[0];
rz(-0.49783266) q[0];
sx q[0];
rz(1.7853907) q[0];
rz(-2.8721696) q[2];
sx q[2];
rz(-1.7413483) q[2];
sx q[2];
rz(1.5829057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8935834) q[1];
sx q[1];
rz(-1.6343405) q[1];
sx q[1];
rz(-2.7081972) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6167496) q[3];
sx q[3];
rz(-2.3224324) q[3];
sx q[3];
rz(-0.51845779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(-1.1462071) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(-2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-0.78654003) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(-0.4424817) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1174283) q[0];
sx q[0];
rz(-2.6167343) q[0];
sx q[0];
rz(2.3470122) q[0];
rz(1.8779249) q[2];
sx q[2];
rz(-1.0182292) q[2];
sx q[2];
rz(-0.27586684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1861371) q[1];
sx q[1];
rz(-1.8147787) q[1];
sx q[1];
rz(-1.9013491) q[1];
rz(0.093591452) q[3];
sx q[3];
rz(-1.8851265) q[3];
sx q[3];
rz(-2.4623354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68142146) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(0.9712514) q[2];
rz(0.69532895) q[3];
sx q[3];
rz(-0.23614241) q[3];
sx q[3];
rz(-2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5044395) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(-1.0094281) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(1.1632464) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.163994) q[0];
sx q[0];
rz(-1.7118651) q[0];
sx q[0];
rz(1.903957) q[0];
rz(-1.8613569) q[2];
sx q[2];
rz(-0.72939789) q[2];
sx q[2];
rz(1.9733719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35112652) q[1];
sx q[1];
rz(-2.8664221) q[1];
sx q[1];
rz(0.64063425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2054687) q[3];
sx q[3];
rz(-2.4118773) q[3];
sx q[3];
rz(-0.038289379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(-3.0187507) q[0];
rz(-3.0154862) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.2164446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2467263) q[0];
sx q[0];
rz(-2.3658532) q[0];
sx q[0];
rz(0.6070444) q[0];
rz(-pi) q[1];
rz(-0.33151303) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(-2.8216459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.053902102) q[1];
sx q[1];
rz(-1.9551827) q[1];
sx q[1];
rz(1.3575421) q[1];
rz(-pi) q[2];
rz(0.27463953) q[3];
sx q[3];
rz(-0.60987771) q[3];
sx q[3];
rz(-0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2609666) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(0.20877008) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(1.1245022) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1632657) q[0];
sx q[0];
rz(-1.7074589) q[0];
sx q[0];
rz(-0.13701963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3051946) q[2];
sx q[2];
rz(-2.2679066) q[2];
sx q[2];
rz(-0.12649378) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75301718) q[1];
sx q[1];
rz(-1.57197) q[1];
sx q[1];
rz(-1.7041824) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87499683) q[3];
sx q[3];
rz(-1.9584624) q[3];
sx q[3];
rz(-2.8029203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889325) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(-1.0349405) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4231739) q[0];
sx q[0];
rz(-0.75773865) q[0];
sx q[0];
rz(-0.54990479) q[0];
rz(-pi) q[1];
rz(-1.5871928) q[2];
sx q[2];
rz(-1.2842442) q[2];
sx q[2];
rz(2.3404442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9518937) q[1];
sx q[1];
rz(-0.67283291) q[1];
sx q[1];
rz(-2.2646907) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3224162) q[3];
sx q[3];
rz(-2.2349149) q[3];
sx q[3];
rz(-2.3770203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-2.7330772) q[2];
rz(0.92489964) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395441) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.3760024) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(0.31870141) q[2];
sx q[2];
rz(-0.8209043) q[2];
sx q[2];
rz(2.2742297) q[2];
rz(-0.61990191) q[3];
sx q[3];
rz(-1.3719659) q[3];
sx q[3];
rz(-0.93056783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];