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
rz(3.6297807) q[0];
sx q[0];
rz(9.0253137) q[0];
rz(2.0692628) q[1];
sx q[1];
rz(-0.89858276) q[1];
sx q[1];
rz(-2.8314765) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7785491) q[0];
sx q[0];
rz(-1.1828711) q[0];
sx q[0];
rz(-1.9411731) q[0];
rz(-pi) q[1];
rz(0.96220509) q[2];
sx q[2];
rz(-0.96675379) q[2];
sx q[2];
rz(-0.60985987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5599903) q[1];
sx q[1];
rz(-1.4694012) q[1];
sx q[1];
rz(-2.6603917) q[1];
rz(-pi) q[2];
rz(1.9555803) q[3];
sx q[3];
rz(-1.6720155) q[3];
sx q[3];
rz(1.1846194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8862137) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(1.0211241) q[2];
rz(1.6518263) q[3];
sx q[3];
rz(-1.3083873) q[3];
sx q[3];
rz(-2.3511353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0609695) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(2.8676497) q[0];
rz(2.888389) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(2.9017172) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0027837) q[0];
sx q[0];
rz(-1.0595508) q[0];
sx q[0];
rz(0.89936043) q[0];
x q[1];
rz(0.56218688) q[2];
sx q[2];
rz(-1.8442684) q[2];
sx q[2];
rz(-2.4386465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.024521527) q[1];
sx q[1];
rz(-1.2327984) q[1];
sx q[1];
rz(-1.8132957) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3595294) q[3];
sx q[3];
rz(-1.7156228) q[3];
sx q[3];
rz(-1.9782255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6720752) q[2];
sx q[2];
rz(-1.8263475) q[2];
sx q[2];
rz(-1.391927) q[2];
rz(-2.4972656) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(0.91359514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(2.9665663) q[0];
rz(1.7063829) q[1];
sx q[1];
rz(-1.274546) q[1];
sx q[1];
rz(-0.79644901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957297) q[0];
sx q[0];
rz(-1.3126763) q[0];
sx q[0];
rz(1.902319) q[0];
rz(0.4034243) q[2];
sx q[2];
rz(-0.63948389) q[2];
sx q[2];
rz(-3.0861106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5063616) q[1];
sx q[1];
rz(-2.0331675) q[1];
sx q[1];
rz(-1.1504786) q[1];
x q[2];
rz(-2.6346509) q[3];
sx q[3];
rz(-1.6569417) q[3];
sx q[3];
rz(-2.994842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1548057) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(3.0955443) q[2];
rz(-1.1151399) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(-2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1107156) q[0];
sx q[0];
rz(-0.59030384) q[0];
sx q[0];
rz(2.9764771) q[0];
rz(1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(1.6868235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4835614) q[0];
sx q[0];
rz(-1.665977) q[0];
sx q[0];
rz(-2.8932518) q[0];
rz(-pi) q[1];
rz(1.3619992) q[2];
sx q[2];
rz(-2.0016124) q[2];
sx q[2];
rz(-0.58961678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1537097) q[1];
sx q[1];
rz(-1.6148414) q[1];
sx q[1];
rz(2.2912461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9971354) q[3];
sx q[3];
rz(-1.3071663) q[3];
sx q[3];
rz(-2.3319649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4399061) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-2.9627724) q[2];
rz(-1.6245101) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(1.4823683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521859) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(-1.408761) q[0];
rz(0.57213712) q[1];
sx q[1];
rz(-1.1502384) q[1];
sx q[1];
rz(-0.22452721) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1869769) q[0];
sx q[0];
rz(-0.77565926) q[0];
sx q[0];
rz(-0.55121724) q[0];
rz(-2.3064086) q[2];
sx q[2];
rz(-1.3408644) q[2];
sx q[2];
rz(1.060677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.677331) q[1];
sx q[1];
rz(-0.867093) q[1];
sx q[1];
rz(1.000885) q[1];
rz(-2.5314674) q[3];
sx q[3];
rz(-0.49558276) q[3];
sx q[3];
rz(-0.15935824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3276334) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(0.28780469) q[2];
rz(2.6156901) q[3];
sx q[3];
rz(-1.9764427) q[3];
sx q[3];
rz(2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41878685) q[0];
sx q[0];
rz(-0.64544183) q[0];
sx q[0];
rz(-0.01471113) q[0];
rz(-2.9171433) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(-2.2803578) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5435641) q[0];
sx q[0];
rz(-3.0985986) q[0];
sx q[0];
rz(-0.046648101) q[0];
rz(-pi) q[1];
rz(1.4528571) q[2];
sx q[2];
rz(-2.8045296) q[2];
sx q[2];
rz(-1.7686421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3726681) q[1];
sx q[1];
rz(-2.3834991) q[1];
sx q[1];
rz(-2.7102986) q[1];
rz(0.16448824) q[3];
sx q[3];
rz(-0.93014088) q[3];
sx q[3];
rz(1.759884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8031249) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(-3.0306446) q[2];
rz(1.126193) q[3];
sx q[3];
rz(-1.5092827) q[3];
sx q[3];
rz(2.4748928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2647389) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.8524843) q[0];
rz(-1.8076757) q[1];
sx q[1];
rz(-1.3197118) q[1];
sx q[1];
rz(1.2510373) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97020075) q[0];
sx q[0];
rz(-1.8222061) q[0];
sx q[0];
rz(-0.89604495) q[0];
x q[1];
rz(-0.47908135) q[2];
sx q[2];
rz(-0.63036171) q[2];
sx q[2];
rz(-0.97409526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45197645) q[1];
sx q[1];
rz(-2.2213081) q[1];
sx q[1];
rz(-2.8337949) q[1];
rz(-pi) q[2];
rz(-2.9032307) q[3];
sx q[3];
rz(-1.8619478) q[3];
sx q[3];
rz(-2.1057768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0879849) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(1.717022) q[2];
rz(0.21814957) q[3];
sx q[3];
rz(-1.679922) q[3];
sx q[3];
rz(0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394534) q[0];
sx q[0];
rz(-0.55291432) q[0];
sx q[0];
rz(0.46808991) q[0];
rz(0.83174902) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(2.1591878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7187854) q[0];
sx q[0];
rz(-0.095346538) q[0];
sx q[0];
rz(2.3991292) q[0];
x q[1];
rz(-1.2853107) q[2];
sx q[2];
rz(-2.4997556) q[2];
sx q[2];
rz(-1.2844207) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6693263) q[1];
sx q[1];
rz(-1.865918) q[1];
sx q[1];
rz(-2.8967316) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7782061) q[3];
sx q[3];
rz(-1.6230427) q[3];
sx q[3];
rz(2.047153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9437647) q[2];
sx q[2];
rz(-0.48913893) q[2];
sx q[2];
rz(0.93734199) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(-2.9924801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8775403) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(2.383655) q[0];
rz(2.0856048) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(0.24076375) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4419542) q[0];
sx q[0];
rz(-0.6356217) q[0];
sx q[0];
rz(-0.54185303) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1409615) q[2];
sx q[2];
rz(-0.79509622) q[2];
sx q[2];
rz(1.5305331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61951738) q[1];
sx q[1];
rz(-0.9065827) q[1];
sx q[1];
rz(0.37503506) q[1];
x q[2];
rz(2.8026227) q[3];
sx q[3];
rz(-1.2267988) q[3];
sx q[3];
rz(-2.0262335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15464887) q[2];
sx q[2];
rz(-1.1661531) q[2];
sx q[2];
rz(1.9722975) q[2];
rz(2.5985006) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(2.5663466) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6146522) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(-0.27969435) q[0];
rz(-2.0277297) q[1];
sx q[1];
rz(-1.056347) q[1];
sx q[1];
rz(-2.9934771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0954656) q[0];
sx q[0];
rz(-1.3805423) q[0];
sx q[0];
rz(0.67511651) q[0];
x q[1];
rz(-2.4513638) q[2];
sx q[2];
rz(-1.6950995) q[2];
sx q[2];
rz(-0.0080492654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.95947) q[1];
sx q[1];
rz(-0.74481267) q[1];
sx q[1];
rz(-0.710979) q[1];
rz(2.8612265) q[3];
sx q[3];
rz(-0.70322817) q[3];
sx q[3];
rz(2.230068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1088045) q[2];
sx q[2];
rz(-0.84454909) q[2];
sx q[2];
rz(-0.27842251) q[2];
rz(1.6935211) q[3];
sx q[3];
rz(-1.4689987) q[3];
sx q[3];
rz(1.9765123) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0447822) q[0];
sx q[0];
rz(-0.7881931) q[0];
sx q[0];
rz(-1.0544554) q[0];
rz(1.0983559) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(-2.8849309) q[2];
sx q[2];
rz(-0.59102249) q[2];
sx q[2];
rz(-1.9494117) q[2];
rz(1.7058115) q[3];
sx q[3];
rz(-1.6208456) q[3];
sx q[3];
rz(2.7593233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
