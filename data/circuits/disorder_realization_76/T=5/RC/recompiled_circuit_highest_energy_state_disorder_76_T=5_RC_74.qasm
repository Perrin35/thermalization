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
rz(5.2108555) q[1];
sx q[1];
rz(4.0401754) q[1];
sx q[1];
rz(5.9730692) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3535904) q[0];
sx q[0];
rz(-1.912477) q[0];
sx q[0];
rz(0.41312878) q[0];
rz(0.96220509) q[2];
sx q[2];
rz(-0.96675379) q[2];
sx q[2];
rz(-0.60985987) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.063608544) q[1];
sx q[1];
rz(-1.0922752) q[1];
sx q[1];
rz(-1.4565181) q[1];
rz(-pi) q[2];
rz(3.0324494) q[3];
sx q[3];
rz(-1.1880842) q[3];
sx q[3];
rz(-2.7963137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8862137) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(2.1204685) q[2];
rz(1.4897664) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(0.79045734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0806231) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(0.27394295) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-2.4174523) q[1];
sx q[1];
rz(-2.9017172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0027837) q[0];
sx q[0];
rz(-1.0595508) q[0];
sx q[0];
rz(0.89936043) q[0];
rz(2.657183) q[2];
sx q[2];
rz(-2.5229215) q[2];
sx q[2];
rz(-1.8687488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1170711) q[1];
sx q[1];
rz(-1.2327984) q[1];
sx q[1];
rz(1.8132957) q[1];
rz(0.78206324) q[3];
sx q[3];
rz(-1.4259699) q[3];
sx q[3];
rz(1.1633671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46951744) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(-1.7496656) q[2];
rz(-2.4972656) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(0.91359514) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.367584) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(-0.17502633) q[0];
rz(1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-0.79644901) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0458629) q[0];
sx q[0];
rz(-1.3126763) q[0];
sx q[0];
rz(-1.2392736) q[0];
x q[1];
rz(-2.7381684) q[2];
sx q[2];
rz(-0.63948389) q[2];
sx q[2];
rz(0.055482023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15115498) q[1];
sx q[1];
rz(-0.61444047) q[1];
sx q[1];
rz(-0.68604705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9655592) q[3];
sx q[3];
rz(-0.51357944) q[3];
sx q[3];
rz(1.5638417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1548057) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(-0.046048306) q[2];
rz(2.0264528) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(1.0435638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107156) q[0];
sx q[0];
rz(-0.59030384) q[0];
sx q[0];
rz(0.16511551) q[0];
rz(-1.1394399) q[1];
sx q[1];
rz(-2.20729) q[1];
sx q[1];
rz(-1.4547691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0784521) q[0];
sx q[0];
rz(-1.3236031) q[0];
sx q[0];
rz(1.4726223) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7179081) q[2];
sx q[2];
rz(-0.47587816) q[2];
sx q[2];
rz(-2.0824473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6858721) q[1];
sx q[1];
rz(-0.8511976) q[1];
sx q[1];
rz(-3.0830129) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28819542) q[3];
sx q[3];
rz(-1.1600947) q[3];
sx q[3];
rz(0.87897838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4399061) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(0.17882027) q[2];
rz(-1.5170826) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(-1.4823683) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4894067) q[0];
sx q[0];
rz(-1.4302197) q[0];
sx q[0];
rz(1.7328316) q[0];
rz(-0.57213712) q[1];
sx q[1];
rz(-1.1502384) q[1];
sx q[1];
rz(-2.9170654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8978724) q[0];
sx q[0];
rz(-0.93168726) q[0];
sx q[0];
rz(-2.0452818) q[0];
rz(-2.8357887) q[2];
sx q[2];
rz(-0.85875466) q[2];
sx q[2];
rz(-0.30669566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.677331) q[1];
sx q[1];
rz(-2.2744997) q[1];
sx q[1];
rz(-2.1407077) q[1];
rz(-1.8711617) q[3];
sx q[3];
rz(-1.9711509) q[3];
sx q[3];
rz(-2.3107236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3276334) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(2.853788) q[2];
rz(0.52590251) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(-0.56720081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41878685) q[0];
sx q[0];
rz(-0.64544183) q[0];
sx q[0];
rz(-0.01471113) q[0];
rz(-0.22444935) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(2.2803578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4968729) q[0];
sx q[0];
rz(-1.6137436) q[0];
sx q[0];
rz(-1.5728024) q[0];
x q[1];
rz(-1.6887356) q[2];
sx q[2];
rz(-2.8045296) q[2];
sx q[2];
rz(1.3729505) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76892454) q[1];
sx q[1];
rz(-0.75809352) q[1];
sx q[1];
rz(-0.43129403) q[1];
rz(1.786993) q[3];
sx q[3];
rz(-0.65854581) q[3];
sx q[3];
rz(1.6525846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8031249) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(3.0306446) q[2];
rz(-1.126193) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(-0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647389) q[0];
sx q[0];
rz(-0.56346098) q[0];
sx q[0];
rz(-1.8524843) q[0];
rz(-1.3339169) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(-1.8905554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.23952) q[0];
sx q[0];
rz(-2.4284568) q[0];
sx q[0];
rz(-1.9608742) q[0];
rz(2.5669615) q[2];
sx q[2];
rz(-1.845965) q[2];
sx q[2];
rz(-0.99400101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.928628) q[1];
sx q[1];
rz(-1.3273094) q[1];
sx q[1];
rz(-0.89694556) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9034702) q[3];
sx q[3];
rz(-0.37411896) q[3];
sx q[3];
rz(-1.4033409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.053607792) q[2];
sx q[2];
rz(-1.3259652) q[2];
sx q[2];
rz(-1.717022) q[2];
rz(-2.9234431) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(-0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80213928) q[0];
sx q[0];
rz(-0.55291432) q[0];
sx q[0];
rz(-2.6735027) q[0];
rz(0.83174902) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(2.1591878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7187854) q[0];
sx q[0];
rz(-3.0462461) q[0];
sx q[0];
rz(0.74246343) q[0];
rz(-0.94865145) q[2];
sx q[2];
rz(-1.7402044) q[2];
sx q[2];
rz(0.51727766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.959572) q[1];
sx q[1];
rz(-0.38117711) q[1];
sx q[1];
rz(0.89761727) q[1];
x q[2];
rz(0.053388552) q[3];
sx q[3];
rz(-1.777919) q[3];
sx q[3];
rz(-2.6762251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19782797) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8775403) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(-0.75793761) q[0];
rz(2.0856048) q[1];
sx q[1];
rz(-0.8989369) q[1];
sx q[1];
rz(2.9008289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42007459) q[0];
sx q[0];
rz(-1.2596247) q[0];
sx q[0];
rz(0.56367409) q[0];
rz(1.1409615) q[2];
sx q[2];
rz(-0.79509622) q[2];
sx q[2];
rz(1.5305331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5220753) q[1];
sx q[1];
rz(-0.9065827) q[1];
sx q[1];
rz(0.37503506) q[1];
rz(0.33896995) q[3];
sx q[3];
rz(-1.2267988) q[3];
sx q[3];
rz(-1.1153592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15464887) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(1.1692952) q[2];
rz(0.5430921) q[3];
sx q[3];
rz(-1.8648059) q[3];
sx q[3];
rz(2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6146522) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(-0.27969435) q[0];
rz(-1.1138629) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(-2.9934771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046127) q[0];
sx q[0];
rz(-1.7610504) q[0];
sx q[0];
rz(2.4664761) q[0];
rz(2.4513638) q[2];
sx q[2];
rz(-1.6950995) q[2];
sx q[2];
rz(-3.1335434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.95947) q[1];
sx q[1];
rz(-2.39678) q[1];
sx q[1];
rz(-2.4306137) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.801228) q[3];
sx q[3];
rz(-2.2413669) q[3];
sx q[3];
rz(2.5910395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0327882) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(-1.4480715) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(1.1650803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096810452) q[0];
sx q[0];
rz(-0.7881931) q[0];
sx q[0];
rz(-1.0544554) q[0];
rz(1.0983559) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(-0.57571147) q[2];
sx q[2];
rz(-1.7127219) q[2];
sx q[2];
rz(-0.59319152) q[2];
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
