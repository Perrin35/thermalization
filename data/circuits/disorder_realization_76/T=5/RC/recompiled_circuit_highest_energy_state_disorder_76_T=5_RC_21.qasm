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
rz(0.70145488) q[0];
sx q[0];
rz(-2.6534046) q[0];
sx q[0];
rz(-0.39946431) q[0];
rz(2.0692628) q[1];
sx q[1];
rz(-0.89858276) q[1];
sx q[1];
rz(-2.8314765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7880023) q[0];
sx q[0];
rz(-1.912477) q[0];
sx q[0];
rz(2.7284639) q[0];
rz(2.4422855) q[2];
sx q[2];
rz(-1.0808873) q[2];
sx q[2];
rz(1.8037947) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9609144) q[1];
sx q[1];
rz(-2.6506458) q[1];
sx q[1];
rz(0.21638201) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0324494) q[3];
sx q[3];
rz(-1.9535084) q[3];
sx q[3];
rz(0.34527895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.255379) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(-1.0211241) q[2];
rz(-1.6518263) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(0.79045734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0806231) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(-0.27394295) q[0];
rz(0.25320369) q[1];
sx q[1];
rz(-2.4174523) q[1];
sx q[1];
rz(-0.23987548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0027837) q[0];
sx q[0];
rz(-2.0820419) q[0];
sx q[0];
rz(-0.89936043) q[0];
x q[1];
rz(0.48440964) q[2];
sx q[2];
rz(-0.61867117) q[2];
sx q[2];
rz(-1.8687488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4644313) q[1];
sx q[1];
rz(-1.7993235) q[1];
sx q[1];
rz(0.34739574) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20407014) q[3];
sx q[3];
rz(-0.7925472) q[3];
sx q[3];
rz(-2.8784405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46951744) q[2];
sx q[2];
rz(-1.8263475) q[2];
sx q[2];
rz(-1.391927) q[2];
rz(-0.64432708) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.367584) q[0];
sx q[0];
rz(-2.3180361) q[0];
sx q[0];
rz(2.9665663) q[0];
rz(1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-0.79644901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957297) q[0];
sx q[0];
rz(-1.8289164) q[0];
sx q[0];
rz(1.2392736) q[0];
rz(1.8548707) q[2];
sx q[2];
rz(-0.98978087) q[2];
sx q[2];
rz(0.43339455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1323582) q[1];
sx q[1];
rz(-1.9446484) q[1];
sx q[1];
rz(0.49970766) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6692596) q[3];
sx q[3];
rz(-1.0659127) q[3];
sx q[3];
rz(1.7652924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.986787) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(-3.0955443) q[2];
rz(2.0264528) q[3];
sx q[3];
rz(-1.0067078) q[3];
sx q[3];
rz(2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030877) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(0.16511551) q[0];
rz(-2.0021527) q[1];
sx q[1];
rz(-2.20729) q[1];
sx q[1];
rz(-1.6868235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44581977) q[0];
sx q[0];
rz(-0.26560387) q[0];
sx q[0];
rz(-0.37047343) q[0];
x q[1];
rz(0.43920774) q[2];
sx q[2];
rz(-1.7602709) q[2];
sx q[2];
rz(-2.2486692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6858721) q[1];
sx q[1];
rz(-0.8511976) q[1];
sx q[1];
rz(0.058579727) q[1];
x q[2];
rz(2.8533972) q[3];
sx q[3];
rz(-1.981498) q[3];
sx q[3];
rz(-2.2626143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4399061) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-0.17882027) q[2];
rz(1.6245101) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521859) q[0];
sx q[0];
rz(-1.4302197) q[0];
sx q[0];
rz(-1.7328316) q[0];
rz(2.5694555) q[1];
sx q[1];
rz(-1.9913543) q[1];
sx q[1];
rz(2.9170654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8978724) q[0];
sx q[0];
rz(-2.2099054) q[0];
sx q[0];
rz(2.0452818) q[0];
rz(-1.9064181) q[2];
sx q[2];
rz(-2.3773509) q[2];
sx q[2];
rz(-2.3847876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71345879) q[1];
sx q[1];
rz(-1.1468219) q[1];
sx q[1];
rz(-0.78935539) q[1];
rz(-1.8711617) q[3];
sx q[3];
rz(-1.1704418) q[3];
sx q[3];
rz(2.3107236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81395927) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(2.853788) q[2];
rz(-0.52590251) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(0.56720081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.7228058) q[0];
sx q[0];
rz(-0.64544183) q[0];
sx q[0];
rz(-0.01471113) q[0];
rz(-2.9171433) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(0.86123484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447198) q[0];
sx q[0];
rz(-1.6137436) q[0];
sx q[0];
rz(1.5728024) q[0];
x q[1];
rz(-1.6887356) q[2];
sx q[2];
rz(-0.33706306) q[2];
sx q[2];
rz(-1.3729505) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76892454) q[1];
sx q[1];
rz(-0.75809352) q[1];
sx q[1];
rz(-0.43129403) q[1];
rz(2.2179763) q[3];
sx q[3];
rz(-1.7024534) q[3];
sx q[3];
rz(0.090196808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8031249) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(3.0306446) q[2];
rz(1.126193) q[3];
sx q[3];
rz(-1.5092827) q[3];
sx q[3];
rz(2.4748928) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768537) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.8524843) q[0];
rz(1.3339169) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(-1.2510373) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90207267) q[0];
sx q[0];
rz(-0.7131359) q[0];
sx q[0];
rz(-1.1807185) q[0];
rz(0.47908135) q[2];
sx q[2];
rz(-0.63036171) q[2];
sx q[2];
rz(0.97409526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2062155) q[1];
sx q[1];
rz(-0.70997974) q[1];
sx q[1];
rz(1.1919271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2716633) q[3];
sx q[3];
rz(-1.342648) q[3];
sx q[3];
rz(0.46534416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0879849) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(1.4245707) q[2];
rz(2.9234431) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(-2.654259) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3394534) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(0.46808991) q[0];
rz(0.83174902) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(-0.98240486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8881852) q[0];
sx q[0];
rz(-1.5063852) q[0];
sx q[0];
rz(3.0712434) q[0];
rz(2.1929412) q[2];
sx q[2];
rz(-1.4013883) q[2];
sx q[2];
rz(-0.51727766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0259798) q[1];
sx q[1];
rz(-1.8048688) q[1];
sx q[1];
rz(-1.8744517) q[1];
rz(1.8194852) q[3];
sx q[3];
rz(-0.2137972) q[3];
sx q[3];
rz(0.71960654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9437647) q[2];
sx q[2];
rz(-2.6524537) q[2];
sx q[2];
rz(-0.93734199) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.849544) q[3];
sx q[3];
rz(-0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.26405239) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(-0.75793761) q[0];
rz(-1.0559878) q[1];
sx q[1];
rz(-0.8989369) q[1];
sx q[1];
rz(-0.24076375) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42007459) q[0];
sx q[0];
rz(-1.2596247) q[0];
sx q[0];
rz(-2.5779186) q[0];
rz(1.1409615) q[2];
sx q[2];
rz(-2.3464964) q[2];
sx q[2];
rz(-1.5305331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71320885) q[1];
sx q[1];
rz(-1.2782103) q[1];
sx q[1];
rz(2.2702515) q[1];
rz(-pi) q[2];
rz(-2.3189697) q[3];
sx q[3];
rz(-0.47815505) q[3];
sx q[3];
rz(-1.9229696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9869438) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(-1.1692952) q[2];
rz(0.5430921) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(0.57524601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6146522) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(0.27969435) q[0];
rz(-1.1138629) q[1];
sx q[1];
rz(-1.056347) q[1];
sx q[1];
rz(2.9934771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046127) q[0];
sx q[0];
rz(-1.3805423) q[0];
sx q[0];
rz(-2.4664761) q[0];
rz(0.19377562) q[2];
sx q[2];
rz(-0.69952088) q[2];
sx q[2];
rz(1.7278838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.95947) q[1];
sx q[1];
rz(-2.39678) q[1];
sx q[1];
rz(-0.710979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68365109) q[3];
sx q[3];
rz(-1.3908852) q[3];
sx q[3];
rz(-0.87547567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1088045) q[2];
sx q[2];
rz(-0.84454909) q[2];
sx q[2];
rz(-2.8631701) q[2];
rz(-1.6935211) q[3];
sx q[3];
rz(-1.4689987) q[3];
sx q[3];
rz(1.1650803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096810452) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(2.0432368) q[1];
sx q[1];
rz(-1.8041246) q[1];
sx q[1];
rz(0.95380797) q[1];
rz(-1.7395217) q[2];
sx q[2];
rz(-1.0015971) q[2];
sx q[2];
rz(0.88605273) q[2];
rz(0.050508113) q[3];
sx q[3];
rz(-1.4359513) q[3];
sx q[3];
rz(1.1953228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
