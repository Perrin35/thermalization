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
rz(-1.7880023) q[0];
sx q[0];
rz(-1.2291156) q[0];
sx q[0];
rz(0.41312878) q[0];
rz(-pi) q[1];
rz(-2.1793876) q[2];
sx q[2];
rz(-2.1748389) q[2];
sx q[2];
rz(0.60985987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5599903) q[1];
sx q[1];
rz(-1.6721915) q[1];
sx q[1];
rz(0.48120099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1860124) q[3];
sx q[3];
rz(-1.6720155) q[3];
sx q[3];
rz(-1.1846194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.255379) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(-1.0211241) q[2];
rz(-1.6518263) q[3];
sx q[3];
rz(-1.3083873) q[3];
sx q[3];
rz(2.3511353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0609695) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(0.27394295) q[0];
rz(2.888389) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(2.9017172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1254933) q[0];
sx q[0];
rz(-2.3224127) q[0];
sx q[0];
rz(-2.3045834) q[0];
rz(-0.48440964) q[2];
sx q[2];
rz(-0.61867117) q[2];
sx q[2];
rz(-1.2728438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1170711) q[1];
sx q[1];
rz(-1.2327984) q[1];
sx q[1];
rz(1.8132957) q[1];
rz(-pi) q[2];
rz(0.20407014) q[3];
sx q[3];
rz(-0.7925472) q[3];
sx q[3];
rz(-0.26315215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6720752) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(1.7496656) q[2];
rz(0.64432708) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(-2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740087) q[0];
sx q[0];
rz(-2.3180361) q[0];
sx q[0];
rz(-0.17502633) q[0];
rz(1.7063829) q[1];
sx q[1];
rz(-1.274546) q[1];
sx q[1];
rz(2.3451436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0957297) q[0];
sx q[0];
rz(-1.3126763) q[0];
sx q[0];
rz(-1.2392736) q[0];
x q[1];
rz(1.286722) q[2];
sx q[2];
rz(-0.98978087) q[2];
sx q[2];
rz(-0.43339455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9904377) q[1];
sx q[1];
rz(-0.61444047) q[1];
sx q[1];
rz(2.4555456) q[1];
rz(-pi) q[2];
rz(1.6692596) q[3];
sx q[3];
rz(-2.0756799) q[3];
sx q[3];
rz(-1.7652924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1548057) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(-0.046048306) q[2];
rz(1.1151399) q[3];
sx q[3];
rz(-1.0067078) q[3];
sx q[3];
rz(-2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107156) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(-0.16511551) q[0];
rz(1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(-1.4547691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063140537) q[0];
sx q[0];
rz(-1.3236031) q[0];
sx q[0];
rz(-1.4726223) q[0];
rz(-pi) q[1];
rz(-2.7179081) q[2];
sx q[2];
rz(-0.47587816) q[2];
sx q[2];
rz(-2.0824473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.36697179) q[1];
sx q[1];
rz(-0.72155385) q[1];
sx q[1];
rz(-1.5040892) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99253486) q[3];
sx q[3];
rz(-0.49697486) q[3];
sx q[3];
rz(0.24028905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4399061) q[2];
sx q[2];
rz(-0.88752282) q[2];
sx q[2];
rz(-0.17882027) q[2];
rz(1.6245101) q[3];
sx q[3];
rz(-2.8330018) q[3];
sx q[3];
rz(-1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521859) q[0];
sx q[0];
rz(-1.4302197) q[0];
sx q[0];
rz(1.408761) q[0];
rz(2.5694555) q[1];
sx q[1];
rz(-1.1502384) q[1];
sx q[1];
rz(0.22452721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029789797) q[0];
sx q[0];
rz(-1.9462612) q[0];
sx q[0];
rz(0.69598868) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3064086) q[2];
sx q[2];
rz(-1.8007282) q[2];
sx q[2];
rz(2.0809157) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8967486) q[1];
sx q[1];
rz(-0.87376423) q[1];
sx q[1];
rz(2.5752707) q[1];
rz(-pi) q[2];
rz(-2.5314674) q[3];
sx q[3];
rz(-2.6460099) q[3];
sx q[3];
rz(0.15935824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.81395927) q[2];
sx q[2];
rz(-1.9093134) q[2];
sx q[2];
rz(-0.28780469) q[2];
rz(0.52590251) q[3];
sx q[3];
rz(-1.9764427) q[3];
sx q[3];
rz(0.56720081) q[3];
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
rz(-pi/2) q[0];
x q[1];
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
rz(0.22444935) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(-2.2803578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0677554) q[0];
sx q[0];
rz(-1.5687921) q[0];
sx q[0];
rz(-3.0986453) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4528571) q[2];
sx q[2];
rz(-2.8045296) q[2];
sx q[2];
rz(1.3729505) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47937313) q[1];
sx q[1];
rz(-1.279261) q[1];
sx q[1];
rz(0.71034224) q[1];
rz(-pi) q[2];
rz(-1.3545996) q[3];
sx q[3];
rz(-0.65854581) q[3];
sx q[3];
rz(1.6525846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33846778) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(-0.11094805) q[2];
rz(2.0153996) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(-0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.2647389) q[0];
sx q[0];
rz(-0.56346098) q[0];
sx q[0];
rz(-1.8524843) q[0];
rz(-1.3339169) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(1.2510373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1713919) q[0];
sx q[0];
rz(-1.8222061) q[0];
sx q[0];
rz(-2.2455477) q[0];
rz(-pi) q[1];
rz(2.6625113) q[2];
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
rz(-pi) q[0];
x q[0];
rz(0.928628) q[1];
sx q[1];
rz(-1.8142833) q[1];
sx q[1];
rz(2.2446471) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.053607792) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(1.4245707) q[2];
rz(-2.9234431) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(-0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.3394534) q[0];
sx q[0];
rz(-0.55291432) q[0];
sx q[0];
rz(-2.6735027) q[0];
rz(-2.3098436) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(-0.98240486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2534075) q[0];
sx q[0];
rz(-1.6352075) q[0];
sx q[0];
rz(-3.0712434) q[0];
rz(-2.1929412) q[2];
sx q[2];
rz(-1.4013883) q[2];
sx q[2];
rz(-2.624315) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0259798) q[1];
sx q[1];
rz(-1.3367238) q[1];
sx q[1];
rz(1.267141) q[1];
rz(-pi) q[2];
rz(1.7782061) q[3];
sx q[3];
rz(-1.5185499) q[3];
sx q[3];
rz(-1.0944397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19782797) q[2];
sx q[2];
rz(-0.48913893) q[2];
sx q[2];
rz(2.2042507) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(-2.9924801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.0559878) q[1];
sx q[1];
rz(-0.8989369) q[1];
sx q[1];
rz(2.9008289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4419542) q[0];
sx q[0];
rz(-2.505971) q[0];
sx q[0];
rz(2.5997396) q[0];
rz(2.7398213) q[2];
sx q[2];
rz(-0.86454287) q[2];
sx q[2];
rz(-2.1101949) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61951738) q[1];
sx q[1];
rz(-2.23501) q[1];
sx q[1];
rz(-2.7665576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2077755) q[3];
sx q[3];
rz(-1.8891834) q[3];
sx q[3];
rz(-0.33708474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9869438) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(-1.1692952) q[2];
rz(-0.5430921) q[3];
sx q[3];
rz(-1.8648059) q[3];
sx q[3];
rz(-2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52694046) q[0];
sx q[0];
rz(-1.4450547) q[0];
sx q[0];
rz(0.27969435) q[0];
rz(-1.1138629) q[1];
sx q[1];
rz(-1.056347) q[1];
sx q[1];
rz(2.9934771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24337846) q[0];
sx q[0];
rz(-2.4442456) q[0];
sx q[0];
rz(2.8426857) q[0];
x q[1];
rz(2.4513638) q[2];
sx q[2];
rz(-1.6950995) q[2];
sx q[2];
rz(-3.1335434) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.95947) q[1];
sx q[1];
rz(-2.39678) q[1];
sx q[1];
rz(0.710979) q[1];
x q[2];
rz(-2.4579416) q[3];
sx q[3];
rz(-1.3908852) q[3];
sx q[3];
rz(-0.87547567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1088045) q[2];
sx q[2];
rz(-0.84454909) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(1.4480715) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.0432368) q[1];
sx q[1];
rz(-1.8041246) q[1];
sx q[1];
rz(0.95380797) q[1];
rz(0.57571147) q[2];
sx q[2];
rz(-1.4288708) q[2];
sx q[2];
rz(2.5484011) q[2];
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
