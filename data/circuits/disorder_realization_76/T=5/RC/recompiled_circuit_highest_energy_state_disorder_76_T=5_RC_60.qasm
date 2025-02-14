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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3535904) q[0];
sx q[0];
rz(-1.912477) q[0];
sx q[0];
rz(-2.7284639) q[0];
rz(2.1793876) q[2];
sx q[2];
rz(-2.1748389) q[2];
sx q[2];
rz(-0.60985987) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18067828) q[1];
sx q[1];
rz(-2.6506458) q[1];
sx q[1];
rz(-2.9252106) q[1];
rz(-3.0324494) q[3];
sx q[3];
rz(-1.1880842) q[3];
sx q[3];
rz(-0.34527895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8862137) q[2];
sx q[2];
rz(-1.9510521) q[2];
sx q[2];
rz(2.1204685) q[2];
rz(1.6518263) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(2.3511353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0609695) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(2.8676497) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-2.4174523) q[1];
sx q[1];
rz(0.23987548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016099361) q[0];
sx q[0];
rz(-0.81917995) q[0];
sx q[0];
rz(-0.83700928) q[0];
x q[1];
rz(0.56218688) q[2];
sx q[2];
rz(-1.2973243) q[2];
sx q[2];
rz(-0.70294619) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6771614) q[1];
sx q[1];
rz(-1.7993235) q[1];
sx q[1];
rz(2.7941969) q[1];
x q[2];
rz(-2.3595294) q[3];
sx q[3];
rz(-1.7156228) q[3];
sx q[3];
rz(-1.1633671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46951744) q[2];
sx q[2];
rz(-1.8263475) q[2];
sx q[2];
rz(1.7496656) q[2];
rz(2.4972656) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(0.17502633) q[0];
rz(1.7063829) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-2.3451436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0284517) q[0];
sx q[0];
rz(-2.7243834) q[0];
sx q[0];
rz(0.88930486) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5416614) q[2];
sx q[2];
rz(-1.8072727) q[2];
sx q[2];
rz(-1.8452871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9904377) q[1];
sx q[1];
rz(-0.61444047) q[1];
sx q[1];
rz(0.68604705) q[1];
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
rz(-0.046048306) q[2];
rz(-1.1151399) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(1.0435638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030877) q[0];
sx q[0];
rz(-0.59030384) q[0];
sx q[0];
rz(-2.9764771) q[0];
rz(1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(1.6868235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6957729) q[0];
sx q[0];
rz(-0.26560387) q[0];
sx q[0];
rz(0.37047343) q[0];
rz(-pi) q[1];
rz(0.42368453) q[2];
sx q[2];
rz(-0.47587816) q[2];
sx q[2];
rz(1.0591454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6858721) q[1];
sx q[1];
rz(-0.8511976) q[1];
sx q[1];
rz(3.0830129) q[1];
rz(2.8533972) q[3];
sx q[3];
rz(-1.1600947) q[3];
sx q[3];
rz(-0.87897838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7016865) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-2.9627724) q[2];
rz(1.5170826) q[3];
sx q[3];
rz(-2.8330018) q[3];
sx q[3];
rz(1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521859) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(-1.7328316) q[0];
rz(2.5694555) q[1];
sx q[1];
rz(-1.9913543) q[1];
sx q[1];
rz(2.9170654) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8978724) q[0];
sx q[0];
rz(-0.93168726) q[0];
sx q[0];
rz(1.0963109) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3064086) q[2];
sx q[2];
rz(-1.8007282) q[2];
sx q[2];
rz(2.0809157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.677331) q[1];
sx q[1];
rz(-2.2744997) q[1];
sx q[1];
rz(-1.000885) q[1];
rz(-2.5314674) q[3];
sx q[3];
rz(-0.49558276) q[3];
sx q[3];
rz(-0.15935824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3276334) q[2];
sx q[2];
rz(-1.9093134) q[2];
sx q[2];
rz(-2.853788) q[2];
rz(-0.52590251) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(-2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7228058) q[0];
sx q[0];
rz(-2.4961508) q[0];
sx q[0];
rz(0.01471113) q[0];
rz(-0.22444935) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(-0.86123484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447198) q[0];
sx q[0];
rz(-1.6137436) q[0];
sx q[0];
rz(-1.5687902) q[0];
x q[1];
rz(-1.235903) q[2];
sx q[2];
rz(-1.5318724) q[2];
sx q[2];
rz(-0.30920497) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3337832) q[1];
sx q[1];
rz(-0.89620608) q[1];
sx q[1];
rz(1.9476931) q[1];
rz(-pi) q[2];
rz(-1.786993) q[3];
sx q[3];
rz(-0.65854581) q[3];
sx q[3];
rz(1.4890081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33846778) q[2];
sx q[2];
rz(-0.71102342) q[2];
sx q[2];
rz(-0.11094805) q[2];
rz(-2.0153996) q[3];
sx q[3];
rz(-1.5092827) q[3];
sx q[3];
rz(2.4748928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8768537) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.2891084) q[0];
rz(1.8076757) q[1];
sx q[1];
rz(-1.3197118) q[1];
sx q[1];
rz(-1.2510373) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1713919) q[0];
sx q[0];
rz(-1.3193865) q[0];
sx q[0];
rz(-2.2455477) q[0];
rz(-pi) q[1];
rz(1.8952605) q[2];
sx q[2];
rz(-1.0203385) q[2];
sx q[2];
rz(2.7389527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93537718) q[1];
sx q[1];
rz(-0.70997974) q[1];
sx q[1];
rz(1.9496656) q[1];
x q[2];
rz(0.9034702) q[3];
sx q[3];
rz(-0.37411896) q[3];
sx q[3];
rz(-1.7382517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053607792) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(-1.4245707) q[2];
rz(2.9234431) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80213928) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(0.46808991) q[0];
rz(-2.3098436) q[1];
sx q[1];
rz(-2.2305198) q[1];
sx q[1];
rz(0.98240486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2534075) q[0];
sx q[0];
rz(-1.5063852) q[0];
sx q[0];
rz(-3.0712434) q[0];
rz(0.20745785) q[2];
sx q[2];
rz(-2.1827123) q[2];
sx q[2];
rz(0.93319172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.959572) q[1];
sx q[1];
rz(-2.7604155) q[1];
sx q[1];
rz(-0.89761727) q[1];
rz(-pi) q[2];
rz(-3.0882041) q[3];
sx q[3];
rz(-1.777919) q[3];
sx q[3];
rz(0.46536756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19782797) q[2];
sx q[2];
rz(-2.6524537) q[2];
sx q[2];
rz(2.2042507) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.26405239) q[0];
sx q[0];
rz(-0.5101246) q[0];
sx q[0];
rz(0.75793761) q[0];
rz(2.0856048) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(-2.9008289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997076) q[0];
sx q[0];
rz(-2.1043964) q[0];
sx q[0];
rz(-1.207229) q[0];
rz(2.3182436) q[2];
sx q[2];
rz(-1.8728796) q[2];
sx q[2];
rz(-2.8712517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1878464) q[1];
sx q[1];
rz(-2.3930379) q[1];
sx q[1];
rz(1.1331609) q[1];
rz(-2.8026227) q[3];
sx q[3];
rz(-1.9147938) q[3];
sx q[3];
rz(-2.0262335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15464887) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(-1.1692952) q[2];
rz(-0.5430921) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(2.8618983) q[0];
rz(-2.0277297) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(2.9934771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8982142) q[0];
sx q[0];
rz(-0.69734708) q[0];
sx q[0];
rz(2.8426857) q[0];
x q[1];
rz(1.4101548) q[2];
sx q[2];
rz(-2.2546708) q[2];
sx q[2];
rz(-1.4768254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.095299) q[1];
sx q[1];
rz(-2.110184) q[1];
sx q[1];
rz(2.1124243) q[1];
rz(-pi) q[2];
rz(-2.4579416) q[3];
sx q[3];
rz(-1.7507075) q[3];
sx q[3];
rz(0.87547567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1088045) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(1.6935211) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(-1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.096810452) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(-1.0983559) q[1];
sx q[1];
rz(-1.8041246) q[1];
sx q[1];
rz(0.95380797) q[1];
rz(-0.25666176) q[2];
sx q[2];
rz(-2.5505702) q[2];
sx q[2];
rz(1.1921809) q[2];
rz(-1.7058115) q[3];
sx q[3];
rz(-1.5207471) q[3];
sx q[3];
rz(-0.38226939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
