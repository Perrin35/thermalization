OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7795774) q[0];
sx q[0];
rz(-0.22688046) q[0];
sx q[0];
rz(1.7664631) q[0];
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
rz(1.8178806) q[0];
sx q[0];
rz(-1.4120988) q[0];
sx q[0];
rz(-2.236332) q[0];
rz(-pi) q[1];
rz(1.9380577) q[2];
sx q[2];
rz(-2.8312771) q[2];
sx q[2];
rz(-0.65904891) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5136265) q[1];
sx q[1];
rz(-1.1556323) q[1];
sx q[1];
rz(-0.28860753) q[1];
rz(2.6219764) q[3];
sx q[3];
rz(-2.5027788) q[3];
sx q[3];
rz(2.6878302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0278339) q[2];
sx q[2];
rz(-1.4049302) q[2];
sx q[2];
rz(-0.11715451) q[2];
rz(-0.33200085) q[3];
sx q[3];
rz(-2.3947075) q[3];
sx q[3];
rz(-0.78506708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431045) q[0];
sx q[0];
rz(-0.76258689) q[0];
sx q[0];
rz(-0.37345988) q[0];
rz(-2.7442878) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(-0.73748803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0720754) q[0];
sx q[0];
rz(-2.2450819) q[0];
sx q[0];
rz(-2.2761795) q[0];
x q[1];
rz(-1.039871) q[2];
sx q[2];
rz(-2.061452) q[2];
sx q[2];
rz(0.38035989) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4026248) q[1];
sx q[1];
rz(-2.1161656) q[1];
sx q[1];
rz(-1.3907554) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37009671) q[3];
sx q[3];
rz(-0.79269275) q[3];
sx q[3];
rz(-2.0065789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.069933683) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(2.7640479) q[2];
rz(2.8318882) q[3];
sx q[3];
rz(-2.0524502) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6250703) q[0];
sx q[0];
rz(-2.3023038) q[0];
sx q[0];
rz(-1.9260433) q[0];
rz(-1.0141605) q[1];
sx q[1];
rz(-1.7182257) q[1];
sx q[1];
rz(-0.97022143) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5824066) q[0];
sx q[0];
rz(-1.0716039) q[0];
sx q[0];
rz(-0.046120709) q[0];
x q[1];
rz(-0.91328158) q[2];
sx q[2];
rz(-2.0184419) q[2];
sx q[2];
rz(-2.7003082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7298895) q[1];
sx q[1];
rz(-1.7474993) q[1];
sx q[1];
rz(-1.9167625) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7895864) q[3];
sx q[3];
rz(-1.6791293) q[3];
sx q[3];
rz(1.7217404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0479451) q[2];
sx q[2];
rz(-0.84257546) q[2];
sx q[2];
rz(2.688431) q[2];
rz(1.9122745) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1969084) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(2.7759283) q[0];
rz(2.2293495) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(2.7361187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.643754) q[0];
sx q[0];
rz(-0.23786834) q[0];
sx q[0];
rz(-2.0980623) q[0];
rz(-pi) q[1];
rz(2.3857152) q[2];
sx q[2];
rz(-1.7265994) q[2];
sx q[2];
rz(2.1550117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7995389) q[1];
sx q[1];
rz(-0.56302445) q[1];
sx q[1];
rz(0.31837007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2527189) q[3];
sx q[3];
rz(-2.1858099) q[3];
sx q[3];
rz(-1.1078887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0696062) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(1.5988505) q[2];
rz(-0.37374464) q[3];
sx q[3];
rz(-1.6662686) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986346) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(-2.4454818) q[0];
rz(1.2593345) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(-2.7764244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1183305) q[0];
sx q[0];
rz(-1.0749165) q[0];
sx q[0];
rz(1.0773247) q[0];
x q[1];
rz(0.31814534) q[2];
sx q[2];
rz(-1.1765624) q[2];
sx q[2];
rz(0.94903681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2303673) q[1];
sx q[1];
rz(-0.54552286) q[1];
sx q[1];
rz(-0.0026774252) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9837512) q[3];
sx q[3];
rz(-2.1574508) q[3];
sx q[3];
rz(1.1053305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68088561) q[2];
sx q[2];
rz(-1.1628217) q[2];
sx q[2];
rz(-0.17670259) q[2];
rz(-1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(-0.37469125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52727592) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(-1.5981307) q[0];
rz(-1.8371001) q[1];
sx q[1];
rz(-2.0871128) q[1];
sx q[1];
rz(-2.3354882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62063187) q[0];
sx q[0];
rz(-1.8072819) q[0];
sx q[0];
rz(2.3804401) q[0];
rz(-pi) q[1];
rz(2.2564933) q[2];
sx q[2];
rz(-2.4174567) q[2];
sx q[2];
rz(-2.372218) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1404952) q[1];
sx q[1];
rz(-1.6204658) q[1];
sx q[1];
rz(-1.7752663) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0433572) q[3];
sx q[3];
rz(-2.0813848) q[3];
sx q[3];
rz(3.1289738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97258893) q[2];
sx q[2];
rz(-2.6469595) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(-0.58498597) q[3];
sx q[3];
rz(-1.7413185) q[3];
sx q[3];
rz(1.5144279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99120283) q[0];
sx q[0];
rz(-1.3911893) q[0];
sx q[0];
rz(2.4124131) q[0];
rz(-1.9623307) q[1];
sx q[1];
rz(-0.74960342) q[1];
sx q[1];
rz(2.8470305) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23239947) q[0];
sx q[0];
rz(-0.56944427) q[0];
sx q[0];
rz(-2.6528986) q[0];
rz(-0.11004098) q[2];
sx q[2];
rz(-1.1400643) q[2];
sx q[2];
rz(-0.34039341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.197261) q[1];
sx q[1];
rz(-0.75350584) q[1];
sx q[1];
rz(1.7427518) q[1];
rz(-pi) q[2];
rz(0.90055777) q[3];
sx q[3];
rz(-1.2646741) q[3];
sx q[3];
rz(0.30508074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4906759) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(2.4093464) q[2];
rz(2.2722774) q[3];
sx q[3];
rz(-1.9711875) q[3];
sx q[3];
rz(1.4066345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9100087) q[0];
sx q[0];
rz(-1.5554447) q[0];
sx q[0];
rz(2.3439132) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-2.1452417) q[1];
sx q[1];
rz(-1.4083883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080264576) q[0];
sx q[0];
rz(-1.3664075) q[0];
sx q[0];
rz(-2.5811282) q[0];
rz(-pi) q[1];
rz(-0.23741053) q[2];
sx q[2];
rz(-1.123073) q[2];
sx q[2];
rz(2.6765649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6646386) q[1];
sx q[1];
rz(-0.79915291) q[1];
sx q[1];
rz(-2.8108912) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8208169) q[3];
sx q[3];
rz(-1.7446221) q[3];
sx q[3];
rz(2.6717693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8352167) q[2];
sx q[2];
rz(-2.0312467) q[2];
sx q[2];
rz(-2.9442673) q[2];
rz(-0.80162445) q[3];
sx q[3];
rz(-0.97951952) q[3];
sx q[3];
rz(0.42158034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9824958) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(0.91484797) q[0];
rz(0.21028701) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-0.68967825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.698477) q[0];
sx q[0];
rz(-1.3453099) q[0];
sx q[0];
rz(1.4754292) q[0];
rz(1.6946227) q[2];
sx q[2];
rz(-3.0231907) q[2];
sx q[2];
rz(1.8121383) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2336809) q[1];
sx q[1];
rz(-0.72919508) q[1];
sx q[1];
rz(-0.75320525) q[1];
rz(-1.0178863) q[3];
sx q[3];
rz(-2.5431271) q[3];
sx q[3];
rz(-0.41494568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0476394) q[2];
sx q[2];
rz(-0.84838212) q[2];
sx q[2];
rz(-0.32996714) q[2];
rz(-2.0983569) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(-0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030815) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(2.0844841) q[0];
rz(2.8874176) q[1];
sx q[1];
rz(-1.9994241) q[1];
sx q[1];
rz(-0.70456299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097226133) q[0];
sx q[0];
rz(-1.8818047) q[0];
sx q[0];
rz(-1.0754271) q[0];
x q[1];
rz(-1.754934) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(-2.9764701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59359351) q[1];
sx q[1];
rz(-1.750046) q[1];
sx q[1];
rz(2.0229193) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7514486) q[3];
sx q[3];
rz(-1.3657939) q[3];
sx q[3];
rz(1.6222749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5083984) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(0.50676662) q[2];
rz(-0.1861598) q[3];
sx q[3];
rz(-0.23701826) q[3];
sx q[3];
rz(-2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7645466) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(-0.90691943) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(-2.2494153) q[2];
sx q[2];
rz(-2.0052675) q[2];
sx q[2];
rz(0.45309767) q[2];
rz(2.4558057) q[3];
sx q[3];
rz(-2.7332173) q[3];
sx q[3];
rz(0.12242534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
