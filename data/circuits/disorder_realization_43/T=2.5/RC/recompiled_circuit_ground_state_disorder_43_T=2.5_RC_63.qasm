OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(-2.9147122) q[0];
sx q[0];
rz(1.3751295) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(-2.8606666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048427933) q[0];
sx q[0];
rz(-0.68138382) q[0];
sx q[0];
rz(-1.3171893) q[0];
x q[1];
rz(1.9380577) q[2];
sx q[2];
rz(-2.8312771) q[2];
sx q[2];
rz(-0.65904891) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0620131) q[1];
sx q[1];
rz(-1.307319) q[1];
sx q[1];
rz(-2.0017712) q[1];
x q[2];
rz(-2.6219764) q[3];
sx q[3];
rz(-2.5027788) q[3];
sx q[3];
rz(-2.6878302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1137587) q[2];
sx q[2];
rz(-1.7366624) q[2];
sx q[2];
rz(0.11715451) q[2];
rz(0.33200085) q[3];
sx q[3];
rz(-2.3947075) q[3];
sx q[3];
rz(-2.3565256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(-2.7681328) q[0];
rz(-2.7442878) q[1];
sx q[1];
rz(-1.9970857) q[1];
sx q[1];
rz(-2.4041046) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0720754) q[0];
sx q[0];
rz(-0.89651075) q[0];
sx q[0];
rz(-0.86541318) q[0];
rz(-pi) q[1];
rz(-2.3829997) q[2];
sx q[2];
rz(-0.70655381) q[2];
sx q[2];
rz(-1.8667081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0651689) q[1];
sx q[1];
rz(-2.5701414) q[1];
sx q[1];
rz(2.854611) q[1];
x q[2];
rz(2.3839398) q[3];
sx q[3];
rz(-1.831358) q[3];
sx q[3];
rz(-2.9716932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.069933683) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(0.37754479) q[2];
rz(-2.8318882) q[3];
sx q[3];
rz(-1.0891424) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51652235) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(-1.2155493) q[0];
rz(-2.1274321) q[1];
sx q[1];
rz(-1.7182257) q[1];
sx q[1];
rz(-2.1713712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.559186) q[0];
sx q[0];
rz(-1.0716039) q[0];
sx q[0];
rz(3.0954719) q[0];
rz(0.54527905) q[2];
sx q[2];
rz(-0.98731326) q[2];
sx q[2];
rz(2.3346221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.045779) q[1];
sx q[1];
rz(-1.9111553) q[1];
sx q[1];
rz(0.18758054) q[1];
rz(-pi) q[2];
rz(-1.4554475) q[3];
sx q[3];
rz(-1.9206502) q[3];
sx q[3];
rz(-3.030341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0479451) q[2];
sx q[2];
rz(-0.84257546) q[2];
sx q[2];
rz(0.45316163) q[2];
rz(1.9122745) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969084) q[0];
sx q[0];
rz(-2.4561645) q[0];
sx q[0];
rz(-2.7759283) q[0];
rz(2.2293495) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(-0.40547392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5580229) q[0];
sx q[0];
rz(-1.6896392) q[0];
sx q[0];
rz(-1.7773377) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7833961) q[2];
sx q[2];
rz(-2.3153164) q[2];
sx q[2];
rz(-0.43897334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1846611) q[1];
sx q[1];
rz(-1.4029364) q[1];
sx q[1];
rz(0.54001684) q[1];
rz(-pi) q[2];
rz(-2.7248091) q[3];
sx q[3];
rz(-2.458771) q[3];
sx q[3];
rz(-2.5522751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0696062) q[2];
sx q[2];
rz(-1.4402086) q[2];
sx q[2];
rz(-1.5427422) q[2];
rz(-2.767848) q[3];
sx q[3];
rz(-1.475324) q[3];
sx q[3];
rz(1.6811949) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986346) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(2.4454818) q[0];
rz(1.2593345) q[1];
sx q[1];
rz(-1.101661) q[1];
sx q[1];
rz(2.7764244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8445963) q[0];
sx q[0];
rz(-1.1410603) q[0];
sx q[0];
rz(0.55083042) q[0];
x q[1];
rz(-2.8234473) q[2];
sx q[2];
rz(-1.9650302) q[2];
sx q[2];
rz(2.1925558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7988749) q[1];
sx q[1];
rz(-1.5694071) q[1];
sx q[1];
rz(-0.54552127) q[1];
rz(-pi) q[2];
rz(2.5138084) q[3];
sx q[3];
rz(-1.2300228) q[3];
sx q[3];
rz(-2.9140811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68088561) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(2.9648901) q[2];
rz(-1.3867311) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(-2.7669014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.52727592) q[0];
sx q[0];
rz(-0.85955954) q[0];
sx q[0];
rz(-1.543462) q[0];
rz(1.8371001) q[1];
sx q[1];
rz(-2.0871128) q[1];
sx q[1];
rz(2.3354882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62063187) q[0];
sx q[0];
rz(-1.3343108) q[0];
sx q[0];
rz(0.76115258) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97053846) q[2];
sx q[2];
rz(-2.0036864) q[2];
sx q[2];
rz(2.8899756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7215868) q[1];
sx q[1];
rz(-1.7750106) q[1];
sx q[1];
rz(-3.0908683) q[1];
x q[2];
rz(-1.0433572) q[3];
sx q[3];
rz(-2.0813848) q[3];
sx q[3];
rz(0.012618806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.97258893) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(0.58498597) q[3];
sx q[3];
rz(-1.4002742) q[3];
sx q[3];
rz(1.5144279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1503898) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(-2.4124131) q[0];
rz(1.179262) q[1];
sx q[1];
rz(-2.3919892) q[1];
sx q[1];
rz(0.29456219) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23239947) q[0];
sx q[0];
rz(-0.56944427) q[0];
sx q[0];
rz(0.48869407) q[0];
x q[1];
rz(-0.11004098) q[2];
sx q[2];
rz(-1.1400643) q[2];
sx q[2];
rz(-0.34039341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.197261) q[1];
sx q[1];
rz(-2.3880868) q[1];
sx q[1];
rz(-1.7427518) q[1];
rz(-pi) q[2];
rz(-2.0414646) q[3];
sx q[3];
rz(-2.4146955) q[3];
sx q[3];
rz(-1.5125546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4906759) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(0.73224625) q[2];
rz(-2.2722774) q[3];
sx q[3];
rz(-1.1704051) q[3];
sx q[3];
rz(1.4066345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.231584) q[0];
sx q[0];
rz(-1.5554447) q[0];
sx q[0];
rz(2.3439132) q[0];
rz(2.8186467) q[1];
sx q[1];
rz(-2.1452417) q[1];
sx q[1];
rz(-1.4083883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080264576) q[0];
sx q[0];
rz(-1.3664075) q[0];
sx q[0];
rz(-2.5811282) q[0];
x q[1];
rz(1.1118719) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(1.2101419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8127784) q[1];
sx q[1];
rz(-1.3359038) q[1];
sx q[1];
rz(2.3702904) q[1];
x q[2];
rz(-1.3878294) q[3];
sx q[3];
rz(-1.8865693) q[3];
sx q[3];
rz(1.0435728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.306376) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(2.9442673) q[2];
rz(2.3399682) q[3];
sx q[3];
rz(-2.1620731) q[3];
sx q[3];
rz(2.7200123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590969) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(0.91484797) q[0];
rz(0.21028701) q[1];
sx q[1];
rz(-2.5631914) q[1];
sx q[1];
rz(-2.4519144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.698477) q[0];
sx q[0];
rz(-1.7962828) q[0];
sx q[0];
rz(1.4754292) q[0];
x q[1];
rz(-1.6946227) q[2];
sx q[2];
rz(-3.0231907) q[2];
sx q[2];
rz(-1.8121383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9079118) q[1];
sx q[1];
rz(-2.4123976) q[1];
sx q[1];
rz(-0.75320525) q[1];
rz(1.0178863) q[3];
sx q[3];
rz(-2.5431271) q[3];
sx q[3];
rz(-2.726647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0476394) q[2];
sx q[2];
rz(-0.84838212) q[2];
sx q[2];
rz(2.8116255) q[2];
rz(-1.0432358) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(-2.6250093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107776) q[0];
sx q[0];
rz(-1.8106221) q[0];
sx q[0];
rz(2.0844841) q[0];
rz(0.2541751) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(-0.70456299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1833558) q[0];
sx q[0];
rz(-0.57794774) q[0];
sx q[0];
rz(-2.1653752) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0800663) q[2];
sx q[2];
rz(-1.7545934) q[2];
sx q[2];
rz(-1.4169324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5479991) q[1];
sx q[1];
rz(-1.750046) q[1];
sx q[1];
rz(-1.1186734) q[1];
rz(-pi) q[2];
rz(0.20829717) q[3];
sx q[3];
rz(-1.7476255) q[3];
sx q[3];
rz(0.014315072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5083984) q[2];
sx q[2];
rz(-1.1722379) q[2];
sx q[2];
rz(-2.634826) q[2];
rz(-0.1861598) q[3];
sx q[3];
rz(-2.9045744) q[3];
sx q[3];
rz(-0.53562927) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7645466) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(-2.2346732) q[1];
sx q[1];
rz(-2.4005371) q[1];
sx q[1];
rz(-1.3225318) q[1];
rz(-2.604031) q[2];
sx q[2];
rz(-2.1765709) q[2];
sx q[2];
rz(-1.4449262) q[2];
rz(-1.3033397) q[3];
sx q[3];
rz(-1.8831913) q[3];
sx q[3];
rz(2.5358653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
