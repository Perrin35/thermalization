OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(0.90149108) q[0];
rz(1.3540406) q[1];
sx q[1];
rz(-1.5259589) q[1];
sx q[1];
rz(1.807133) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62518277) q[0];
sx q[0];
rz(-1.3526655) q[0];
sx q[0];
rz(-1.5152009) q[0];
x q[1];
rz(0.15057474) q[2];
sx q[2];
rz(-1.5158487) q[2];
sx q[2];
rz(-0.51412941) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1062231) q[1];
sx q[1];
rz(-1.5179885) q[1];
sx q[1];
rz(1.04671) q[1];
rz(2.9525312) q[3];
sx q[3];
rz(-2.8950518) q[3];
sx q[3];
rz(-2.2061359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2287075) q[2];
sx q[2];
rz(-0.097147377) q[2];
sx q[2];
rz(2.482282) q[2];
rz(2.3653476) q[3];
sx q[3];
rz(-3.1228437) q[3];
sx q[3];
rz(-0.68500486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.9126251) q[0];
sx q[0];
rz(-1.2094867) q[0];
sx q[0];
rz(1.9483161) q[0];
rz(0.044366447) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(-0.22656974) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18767522) q[0];
sx q[0];
rz(-1.5724475) q[0];
sx q[0];
rz(1.5673914) q[0];
rz(-pi) q[1];
rz(0.042828538) q[2];
sx q[2];
rz(-0.37938243) q[2];
sx q[2];
rz(1.6216175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6084131) q[1];
sx q[1];
rz(-0.4371818) q[1];
sx q[1];
rz(2.9730148) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3188521) q[3];
sx q[3];
rz(-0.74352194) q[3];
sx q[3];
rz(-1.3321259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6988397) q[2];
sx q[2];
rz(-1.5696462) q[2];
sx q[2];
rz(-1.6119831) q[2];
rz(2.2085341) q[3];
sx q[3];
rz(-1.6589087) q[3];
sx q[3];
rz(-0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642692) q[0];
sx q[0];
rz(-3.1260335) q[0];
sx q[0];
rz(2.9413057) q[0];
rz(3.1411723) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(3.1281085) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6967675) q[0];
sx q[0];
rz(-0.42967859) q[0];
sx q[0];
rz(-1.4151156) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57245596) q[2];
sx q[2];
rz(-0.079266313) q[2];
sx q[2];
rz(-2.5545189) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36253769) q[1];
sx q[1];
rz(-1.7064306) q[1];
sx q[1];
rz(-0.070222994) q[1];
rz(-1.5192658) q[3];
sx q[3];
rz(-1.7331707) q[3];
sx q[3];
rz(0.80013093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9424332) q[2];
sx q[2];
rz(-1.556267) q[2];
sx q[2];
rz(1.5419434) q[2];
rz(2.13983) q[3];
sx q[3];
rz(-2.9250513) q[3];
sx q[3];
rz(-0.17222968) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7672985) q[0];
sx q[0];
rz(-2.9758487) q[0];
sx q[0];
rz(2.7941008) q[0];
rz(-2.5047498) q[1];
sx q[1];
rz(-0.0056191365) q[1];
sx q[1];
rz(-1.1905131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325) q[0];
sx q[0];
rz(-0.14483368) q[0];
sx q[0];
rz(-1.1454789) q[0];
rz(-pi) q[1];
rz(1.4585481) q[2];
sx q[2];
rz(-0.069321037) q[2];
sx q[2];
rz(1.4510509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0034297) q[1];
sx q[1];
rz(-0.82050475) q[1];
sx q[1];
rz(1.3123039) q[1];
rz(-pi) q[2];
rz(2.6116294) q[3];
sx q[3];
rz(-2.7254482) q[3];
sx q[3];
rz(-0.79485369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5570598) q[2];
sx q[2];
rz(-0.0396885) q[2];
sx q[2];
rz(1.7812799) q[2];
rz(1.6654061) q[3];
sx q[3];
rz(-1.5740266) q[3];
sx q[3];
rz(0.40569693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.0535468) q[0];
sx q[0];
rz(-0.73943728) q[0];
sx q[0];
rz(-3.1355701) q[0];
rz(1.7209523) q[1];
sx q[1];
rz(-3.0907478) q[1];
sx q[1];
rz(0.086070148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94347826) q[0];
sx q[0];
rz(-3.0331342) q[0];
sx q[0];
rz(-0.15546851) q[0];
rz(-pi) q[1];
rz(-0.055069607) q[2];
sx q[2];
rz(-1.5447642) q[2];
sx q[2];
rz(0.30773417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9169315) q[1];
sx q[1];
rz(-1.6078533) q[1];
sx q[1];
rz(0.091450973) q[1];
rz(0.045343355) q[3];
sx q[3];
rz(-1.6312977) q[3];
sx q[3];
rz(-0.44024548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-0.72282183) q[2];
sx q[2];
rz(-1.7576199) q[2];
rz(0.39517394) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(-1.9743732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60642099) q[0];
sx q[0];
rz(-2.9223154) q[0];
sx q[0];
rz(2.1411335) q[0];
rz(-2.4011627) q[1];
sx q[1];
rz(-0.43559566) q[1];
sx q[1];
rz(-2.7620517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1943037) q[0];
sx q[0];
rz(-1.6662046) q[0];
sx q[0];
rz(-1.5964776) q[0];
rz(-pi) q[1];
rz(1.2984811) q[2];
sx q[2];
rz(-2.0305995) q[2];
sx q[2];
rz(2.1260171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4011456) q[1];
sx q[1];
rz(-1.5441431) q[1];
sx q[1];
rz(-1.0668287) q[1];
x q[2];
rz(-0.11594144) q[3];
sx q[3];
rz(-2.017147) q[3];
sx q[3];
rz(1.6801113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.015811054) q[2];
sx q[2];
rz(-0.28818473) q[2];
sx q[2];
rz(3.06456) q[2];
rz(-0.020126255) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(-2.3441815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035148419) q[0];
sx q[0];
rz(-3.1290717) q[0];
sx q[0];
rz(1.4323819) q[0];
rz(-0.2969946) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(3.0800379) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765574) q[0];
sx q[0];
rz(-1.9752968) q[0];
sx q[0];
rz(-0.96488953) q[0];
rz(-pi) q[1];
rz(1.5764357) q[2];
sx q[2];
rz(-1.3562849) q[2];
sx q[2];
rz(2.0129865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0377215) q[1];
sx q[1];
rz(-1.5024517) q[1];
sx q[1];
rz(0.29199021) q[1];
x q[2];
rz(1.6084303) q[3];
sx q[3];
rz(-2.0052344) q[3];
sx q[3];
rz(-1.5805336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84294549) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(-1.7227777) q[2];
rz(-1.336054) q[3];
sx q[3];
rz(-3.1137443) q[3];
sx q[3];
rz(1.7347887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0272738) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(-1.5007716) q[0];
rz(-1.1227135) q[1];
sx q[1];
rz(-2.8148459) q[1];
sx q[1];
rz(1.2935125) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0121078) q[0];
sx q[0];
rz(-0.92381682) q[0];
sx q[0];
rz(0.39374521) q[0];
x q[1];
rz(-2.4475936) q[2];
sx q[2];
rz(-1.8376164) q[2];
sx q[2];
rz(2.9129183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97630097) q[1];
sx q[1];
rz(-1.8433851) q[1];
sx q[1];
rz(1.3957681) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82610749) q[3];
sx q[3];
rz(-0.87658054) q[3];
sx q[3];
rz(-2.7339767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5571755) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(-1.7822251) q[2];
rz(2.6946097) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(2.9744448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8166703) q[0];
sx q[0];
rz(-1.336038) q[0];
sx q[0];
rz(-2.0409806) q[0];
rz(-1.0758411) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(-0.67695391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8146801) q[0];
sx q[0];
rz(-1.956358) q[0];
sx q[0];
rz(1.7837693) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4482037) q[2];
sx q[2];
rz(-1.4559146) q[2];
sx q[2];
rz(-2.5455395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0648411) q[1];
sx q[1];
rz(-3.1411618) q[1];
sx q[1];
rz(1.3382124) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6864945) q[3];
sx q[3];
rz(-1.7069723) q[3];
sx q[3];
rz(1.0973339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1415652) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(0.72762093) q[2];
rz(-1.8992807) q[3];
sx q[3];
rz(-3.1051903) q[3];
sx q[3];
rz(1.9696994) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2321371) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(-2.452028) q[0];
rz(-1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(0.15588674) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0077654) q[0];
sx q[0];
rz(-1.7031914) q[0];
sx q[0];
rz(-0.65806972) q[0];
rz(2.8901398) q[2];
sx q[2];
rz(-1.4551468) q[2];
sx q[2];
rz(2.419099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.879133) q[1];
sx q[1];
rz(-0.004460881) q[1];
sx q[1];
rz(-0.35426472) q[1];
rz(-pi) q[2];
rz(1.4149551) q[3];
sx q[3];
rz(-1.6644125) q[3];
sx q[3];
rz(-0.35767143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9547687) q[2];
sx q[2];
rz(-3.0195152) q[2];
sx q[2];
rz(2.1464777) q[2];
rz(2.9156445) q[3];
sx q[3];
rz(-3.093231) q[3];
sx q[3];
rz(2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36289287) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(-1.4317935) q[1];
sx q[1];
rz(-1.8451537) q[1];
sx q[1];
rz(0.61737212) q[1];
rz(0.68516784) q[2];
sx q[2];
rz(-1.6772267) q[2];
sx q[2];
rz(-0.75135372) q[2];
rz(1.2479242) q[3];
sx q[3];
rz(-0.40798305) q[3];
sx q[3];
rz(1.8508607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
