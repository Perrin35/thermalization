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
rz(0.22688046) q[0];
sx q[0];
rz(14.332834) q[0];
rz(-0.094376266) q[1];
sx q[1];
rz(4.0552858) q[1];
sx q[1];
rz(6.0022592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710815) q[0];
sx q[0];
rz(-0.91508741) q[0];
sx q[0];
rz(-0.20072584) q[0];
x q[1];
rz(1.2035349) q[2];
sx q[2];
rz(-2.8312771) q[2];
sx q[2];
rz(0.65904891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0795796) q[1];
sx q[1];
rz(-1.307319) q[1];
sx q[1];
rz(1.1398214) q[1];
rz(1.2174843) q[3];
sx q[3];
rz(-2.114776) q[3];
sx q[3];
rz(-2.0687188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0278339) q[2];
sx q[2];
rz(-1.4049302) q[2];
sx q[2];
rz(3.0244381) q[2];
rz(-2.8095918) q[3];
sx q[3];
rz(-2.3947075) q[3];
sx q[3];
rz(-2.3565256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431045) q[0];
sx q[0];
rz(-0.76258689) q[0];
sx q[0];
rz(2.7681328) q[0];
rz(-0.39730486) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(-2.4041046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0076548) q[0];
sx q[0];
rz(-2.2077401) q[0];
sx q[0];
rz(-0.68151125) q[0];
x q[1];
rz(-0.55464427) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(2.2210768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8792845) q[1];
sx q[1];
rz(-1.4170987) q[1];
sx q[1];
rz(-0.55262312) q[1];
rz(2.3839398) q[3];
sx q[3];
rz(-1.3102346) q[3];
sx q[3];
rz(-0.16989947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(-2.7640479) q[2];
rz(2.8318882) q[3];
sx q[3];
rz(-2.0524502) q[3];
sx q[3];
rz(1.496605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6250703) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(-1.9260433) q[0];
rz(-2.1274321) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(2.1713712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6785203) q[0];
sx q[0];
rz(-0.50113916) q[0];
sx q[0];
rz(1.6551514) q[0];
rz(-2.2367291) q[2];
sx q[2];
rz(-0.7762802) q[2];
sx q[2];
rz(1.5011476) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6128224) q[1];
sx q[1];
rz(-0.38685683) q[1];
sx q[1];
rz(-2.0554789) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6861451) q[3];
sx q[3];
rz(-1.9206502) q[3];
sx q[3];
rz(3.030341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0936475) q[2];
sx q[2];
rz(-0.84257546) q[2];
sx q[2];
rz(-2.688431) q[2];
rz(1.2293182) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(-0.17190988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.8625926) q[1];
sx q[1];
rz(0.40547392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835698) q[0];
sx q[0];
rz(-1.6896392) q[0];
sx q[0];
rz(-1.364255) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75587749) q[2];
sx q[2];
rz(-1.4149932) q[2];
sx q[2];
rz(-2.1550117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1846611) q[1];
sx q[1];
rz(-1.7386562) q[1];
sx q[1];
rz(-2.6015758) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41678352) q[3];
sx q[3];
rz(-0.68282167) q[3];
sx q[3];
rz(-0.58931755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0696062) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(1.5427422) q[2];
rz(2.767848) q[3];
sx q[3];
rz(-1.6662686) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.54295802) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(0.69611088) q[0];
rz(1.2593345) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(0.36516821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8445963) q[0];
sx q[0];
rz(-1.1410603) q[0];
sx q[0];
rz(0.55083042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92607195) q[2];
sx q[2];
rz(-0.50130166) q[2];
sx q[2];
rz(0.24025606) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9112253) q[1];
sx q[1];
rz(-0.54552286) q[1];
sx q[1];
rz(3.1389152) q[1];
rz(-2.5138084) q[3];
sx q[3];
rz(-1.2300228) q[3];
sx q[3];
rz(2.9140811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68088561) q[2];
sx q[2];
rz(-1.1628217) q[2];
sx q[2];
rz(-2.9648901) q[2];
rz(-1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(2.7669014) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143167) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(-1.543462) q[0];
rz(-1.8371001) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(2.3354882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70901727) q[0];
sx q[0];
rz(-0.78993778) q[0];
sx q[0];
rz(0.33613251) q[0];
rz(-pi) q[1];
rz(2.2564933) q[2];
sx q[2];
rz(-2.4174567) q[2];
sx q[2];
rz(0.76937461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4762905) q[1];
sx q[1];
rz(-2.9312583) q[1];
sx q[1];
rz(-1.3306985) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0433572) q[3];
sx q[3];
rz(-1.0602078) q[3];
sx q[3];
rz(0.012618806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1690037) q[2];
sx q[2];
rz(-0.4946332) q[2];
sx q[2];
rz(0.063610323) q[2];
rz(0.58498597) q[3];
sx q[3];
rz(-1.4002742) q[3];
sx q[3];
rz(1.5144279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99120283) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(-2.4124131) q[0];
rz(-1.179262) q[1];
sx q[1];
rz(-0.74960342) q[1];
sx q[1];
rz(-2.8470305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23239947) q[0];
sx q[0];
rz(-0.56944427) q[0];
sx q[0];
rz(0.48869407) q[0];
rz(3.0315517) q[2];
sx q[2];
rz(-2.0015284) q[2];
sx q[2];
rz(0.34039341) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7105508) q[1];
sx q[1];
rz(-2.3105544) q[1];
sx q[1];
rz(-2.9824216) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38334966) q[3];
sx q[3];
rz(-2.2046996) q[3];
sx q[3];
rz(-1.0312361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6509167) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(-0.73224625) q[2];
rz(-0.8693153) q[3];
sx q[3];
rz(-1.9711875) q[3];
sx q[3];
rz(-1.7349582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3638316) q[0];
sx q[0];
rz(-2.1182502) q[0];
sx q[0];
rz(-1.3307894) q[0];
rz(-pi) q[1];
rz(2.0261835) q[2];
sx q[2];
rz(-0.50297996) q[2];
sx q[2];
rz(-3.0969381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8127784) q[1];
sx q[1];
rz(-1.3359038) q[1];
sx q[1];
rz(0.77130227) q[1];
rz(-pi) q[2];
rz(-1.3878294) q[3];
sx q[3];
rz(-1.2550233) q[3];
sx q[3];
rz(2.0980199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8352167) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(2.9442673) q[2];
rz(-2.3399682) q[3];
sx q[3];
rz(-0.97951952) q[3];
sx q[3];
rz(-0.42158034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.9313056) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-2.4519144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0352958) q[0];
sx q[0];
rz(-1.4778504) q[0];
sx q[0];
rz(-0.22648099) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6883001) q[2];
sx q[2];
rz(-1.5562061) q[2];
sx q[2];
rz(2.7772825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0090078) q[1];
sx q[1];
rz(-1.0632391) q[1];
sx q[1];
rz(1.0222438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7977169) q[3];
sx q[3];
rz(-2.0708041) q[3];
sx q[3];
rz(-0.22658928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0476394) q[2];
sx q[2];
rz(-2.2932105) q[2];
sx q[2];
rz(-0.32996714) q[2];
rz(2.0983569) q[3];
sx q[3];
rz(-1.4179351) q[3];
sx q[3];
rz(-0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030815) q[0];
sx q[0];
rz(-1.3309706) q[0];
sx q[0];
rz(-1.0571085) q[0];
rz(0.2541751) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(-0.70456299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0443665) q[0];
sx q[0];
rz(-1.8818047) q[0];
sx q[0];
rz(2.0661656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8902337) q[2];
sx q[2];
rz(-2.9478812) q[2];
sx q[2];
rz(1.0918231) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8125638) q[1];
sx q[1];
rz(-0.48406752) q[1];
sx q[1];
rz(-1.1776393) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71263616) q[3];
sx q[3];
rz(-2.8691926) q[3];
sx q[3];
rz(-2.2505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6331943) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(2.634826) q[2];
rz(2.9554328) q[3];
sx q[3];
rz(-0.23701826) q[3];
sx q[3];
rz(-2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3770461) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(-0.90691943) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(0.93420784) q[2];
sx q[2];
rz(-2.3547966) q[2];
sx q[2];
rz(2.5049868) q[2];
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
