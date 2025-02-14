OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(4.0816981) q[0];
sx q[0];
rz(8.8844086) q[0];
rz(0.49343935) q[1];
sx q[1];
rz(-2.4210338) q[1];
sx q[1];
rz(2.5277353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9073528) q[0];
sx q[0];
rz(-2.1932903) q[0];
sx q[0];
rz(-1.7982152) q[0];
x q[1];
rz(-1.5510606) q[2];
sx q[2];
rz(-0.29183772) q[2];
sx q[2];
rz(-0.29668754) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93537027) q[1];
sx q[1];
rz(-2.3435278) q[1];
sx q[1];
rz(-0.3978637) q[1];
x q[2];
rz(-0.24718376) q[3];
sx q[3];
rz(-1.5498115) q[3];
sx q[3];
rz(-2.0387797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8453688) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(-2.5073012) q[2];
rz(0.93572179) q[3];
sx q[3];
rz(-0.24565419) q[3];
sx q[3];
rz(-0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3212386) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(-0.4441922) q[0];
rz(2.3815637) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(2.1601423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8151756) q[0];
sx q[0];
rz(-1.278864) q[0];
sx q[0];
rz(-0.53262226) q[0];
rz(-pi) q[1];
rz(1.3992127) q[2];
sx q[2];
rz(-1.7497471) q[2];
sx q[2];
rz(1.5092261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1552298) q[1];
sx q[1];
rz(-1.6587509) q[1];
sx q[1];
rz(-0.4390688) q[1];
rz(-pi) q[2];
rz(-0.53175521) q[3];
sx q[3];
rz(-0.49177836) q[3];
sx q[3];
rz(-0.39612202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11357073) q[2];
sx q[2];
rz(-0.61442033) q[2];
sx q[2];
rz(-1.2051955) q[2];
rz(-0.53660721) q[3];
sx q[3];
rz(-1.9034932) q[3];
sx q[3];
rz(1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2236915) q[0];
sx q[0];
rz(-1.2811998) q[0];
sx q[0];
rz(-2.3028288) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(-0.020523358) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588226) q[0];
sx q[0];
rz(-0.8626482) q[0];
sx q[0];
rz(2.3282611) q[0];
x q[1];
rz(1.1032365) q[2];
sx q[2];
rz(-0.9976495) q[2];
sx q[2];
rz(1.3064885) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1908326) q[1];
sx q[1];
rz(-1.5632707) q[1];
sx q[1];
rz(-2.7286367) q[1];
x q[2];
rz(2.0666615) q[3];
sx q[3];
rz(-2.819546) q[3];
sx q[3];
rz(-2.5888458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3402349) q[2];
sx q[2];
rz(-1.5284208) q[2];
sx q[2];
rz(-0.94432962) q[2];
rz(-1.0229735) q[3];
sx q[3];
rz(-1.9187656) q[3];
sx q[3];
rz(1.1378707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2706547) q[0];
sx q[0];
rz(-1.174467) q[0];
sx q[0];
rz(1.9842072) q[0];
rz(1.4986787) q[1];
sx q[1];
rz(-1.4721556) q[1];
sx q[1];
rz(1.9532983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0152215) q[0];
sx q[0];
rz(-2.043521) q[0];
sx q[0];
rz(-2.0822099) q[0];
x q[1];
rz(-1.525773) q[2];
sx q[2];
rz(-1.0226215) q[2];
sx q[2];
rz(-0.5612824) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77515652) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-2.3597673) q[1];
x q[2];
rz(-1.2318939) q[3];
sx q[3];
rz(-1.7247685) q[3];
sx q[3];
rz(0.13974259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1184065) q[2];
sx q[2];
rz(-1.964317) q[2];
sx q[2];
rz(-0.6905306) q[2];
rz(2.0265419) q[3];
sx q[3];
rz(-2.3670022) q[3];
sx q[3];
rz(1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(-2.114356) q[0];
rz(1.3014303) q[1];
sx q[1];
rz(-2.4899028) q[1];
sx q[1];
rz(2.8177736) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80211783) q[0];
sx q[0];
rz(-0.81011745) q[0];
sx q[0];
rz(-1.0056061) q[0];
rz(-2.1236046) q[2];
sx q[2];
rz(-1.5828642) q[2];
sx q[2];
rz(2.7928305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71621543) q[1];
sx q[1];
rz(-1.0310804) q[1];
sx q[1];
rz(-1.9464689) q[1];
x q[2];
rz(1.2965167) q[3];
sx q[3];
rz(-2.6135332) q[3];
sx q[3];
rz(-1.5590645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41574898) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(-2.6591163) q[2];
rz(1.9735362) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(0.67679685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93550682) q[0];
sx q[0];
rz(-2.2914903) q[0];
sx q[0];
rz(-2.2555943) q[0];
rz(-0.63367263) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(-0.69127965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1544428) q[0];
sx q[0];
rz(-1.6073174) q[0];
sx q[0];
rz(-2.4958688) q[0];
rz(2.4497767) q[2];
sx q[2];
rz(-1.6199281) q[2];
sx q[2];
rz(-0.092158801) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9554841) q[1];
sx q[1];
rz(-0.82064522) q[1];
sx q[1];
rz(3.0226743) q[1];
rz(-0.78929928) q[3];
sx q[3];
rz(-1.3914445) q[3];
sx q[3];
rz(-1.8478249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0014235) q[2];
sx q[2];
rz(-0.836335) q[2];
sx q[2];
rz(-0.26958618) q[2];
rz(-0.94830281) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(-1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7798994) q[0];
sx q[0];
rz(-2.7044856) q[0];
sx q[0];
rz(-2.1345188) q[0];
rz(-0.40564793) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(-0.55353037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53109518) q[0];
sx q[0];
rz(-2.880504) q[0];
sx q[0];
rz(1.3469264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.361176) q[2];
sx q[2];
rz(-0.71418205) q[2];
sx q[2];
rz(1.2250021) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0900768) q[1];
sx q[1];
rz(-2.0231658) q[1];
sx q[1];
rz(-0.69454792) q[1];
x q[2];
rz(-1.5958435) q[3];
sx q[3];
rz(-1.950693) q[3];
sx q[3];
rz(0.27835007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3125399) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(-1.9276169) q[2];
rz(-2.35516) q[3];
sx q[3];
rz(-1.395547) q[3];
sx q[3];
rz(3.1184375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9455652) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(-1.8096402) q[0];
rz(2.5281483) q[1];
sx q[1];
rz(-1.0204126) q[1];
sx q[1];
rz(0.44949284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7004823) q[0];
sx q[0];
rz(-1.4876517) q[0];
sx q[0];
rz(2.0594199) q[0];
rz(-1.3276991) q[2];
sx q[2];
rz(-0.34533325) q[2];
sx q[2];
rz(-1.9497046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7498988) q[1];
sx q[1];
rz(-0.79872455) q[1];
sx q[1];
rz(-2.9500089) q[1];
rz(0.2576377) q[3];
sx q[3];
rz(-1.6272568) q[3];
sx q[3];
rz(0.27654058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.31676644) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(-0.60834926) q[2];
rz(1.9781205) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(2.5782862) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1083531) q[0];
sx q[0];
rz(-2.2305363) q[0];
sx q[0];
rz(-2.5392927) q[0];
rz(0.96744084) q[1];
sx q[1];
rz(-2.2106705) q[1];
sx q[1];
rz(-2.0379351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3924343) q[0];
sx q[0];
rz(-1.9604248) q[0];
sx q[0];
rz(-0.97831877) q[0];
rz(2.5050312) q[2];
sx q[2];
rz(-1.5825854) q[2];
sx q[2];
rz(1.4135041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97776088) q[1];
sx q[1];
rz(-1.1196616) q[1];
sx q[1];
rz(-3.1402863) q[1];
x q[2];
rz(-1.4673442) q[3];
sx q[3];
rz(-2.561063) q[3];
sx q[3];
rz(0.79642297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.935282) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-2.804011) q[2];
rz(-2.857699) q[3];
sx q[3];
rz(-0.68113911) q[3];
sx q[3];
rz(-2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5934481) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(3.0737851) q[0];
rz(-1.1495122) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(-0.4745208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15856811) q[0];
sx q[0];
rz(-1.3320384) q[0];
sx q[0];
rz(0.013834133) q[0];
rz(-pi) q[1];
rz(-1.7065918) q[2];
sx q[2];
rz(-1.8095922) q[2];
sx q[2];
rz(0.39993024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.671512) q[1];
sx q[1];
rz(-1.5639515) q[1];
sx q[1];
rz(-1.6575302) q[1];
x q[2];
rz(-2.2978333) q[3];
sx q[3];
rz(-0.90137989) q[3];
sx q[3];
rz(0.17217522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6600251) q[2];
sx q[2];
rz(-0.93067545) q[2];
sx q[2];
rz(1.0732667) q[2];
rz(-0.16452161) q[3];
sx q[3];
rz(-1.3273032) q[3];
sx q[3];
rz(0.32065121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58707033) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(1.5578237) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(-2.6804994) q[2];
sx q[2];
rz(-1.8965707) q[2];
sx q[2];
rz(-2.0988219) q[2];
rz(1.6526374) q[3];
sx q[3];
rz(-1.428953) q[3];
sx q[3];
rz(-0.088464213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
