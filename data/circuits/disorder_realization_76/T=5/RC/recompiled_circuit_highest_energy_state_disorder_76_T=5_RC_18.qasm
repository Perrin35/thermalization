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
rz(-1.0723298) q[1];
sx q[1];
rz(-2.2430099) q[1];
sx q[1];
rz(2.8314765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7785491) q[0];
sx q[0];
rz(-1.1828711) q[0];
sx q[0];
rz(1.9411731) q[0];
rz(-pi) q[1];
rz(-0.69930716) q[2];
sx q[2];
rz(-2.0607053) q[2];
sx q[2];
rz(-1.8037947) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0779841) q[1];
sx q[1];
rz(-2.0493174) q[1];
sx q[1];
rz(-1.6850745) q[1];
x q[2];
rz(-1.8350527) q[3];
sx q[3];
rz(-0.39723165) q[3];
sx q[3];
rz(-0.63069944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8862137) q[2];
sx q[2];
rz(-1.9510521) q[2];
sx q[2];
rz(1.0211241) q[2];
rz(1.6518263) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(-0.79045734) q[3];
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
rz(-pi/2) q[3];
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
rz(1.0609695) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(-2.8676497) q[0];
rz(0.25320369) q[1];
sx q[1];
rz(-2.4174523) q[1];
sx q[1];
rz(2.9017172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93879733) q[0];
sx q[0];
rz(-0.99730856) q[0];
sx q[0];
rz(2.5198562) q[0];
x q[1];
rz(-1.8909177) q[2];
sx q[2];
rz(-2.1097398) q[2];
sx q[2];
rz(-0.69931617) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66544724) q[1];
sx q[1];
rz(-2.7283333) q[1];
sx q[1];
rz(-0.59935092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3595294) q[3];
sx q[3];
rz(-1.4259699) q[3];
sx q[3];
rz(-1.1633671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46951744) q[2];
sx q[2];
rz(-1.8263475) q[2];
sx q[2];
rz(1.391927) q[2];
rz(-0.64432708) q[3];
sx q[3];
rz(-1.143012) q[3];
sx q[3];
rz(0.91359514) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(-2.9665663) q[0];
rz(-1.4352098) q[1];
sx q[1];
rz(-1.274546) q[1];
sx q[1];
rz(-0.79644901) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0284517) q[0];
sx q[0];
rz(-2.7243834) q[0];
sx q[0];
rz(2.2522878) q[0];
rz(1.286722) q[2];
sx q[2];
rz(-0.98978087) q[2];
sx q[2];
rz(2.7081981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0092344) q[1];
sx q[1];
rz(-1.1969442) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1548057) q[2];
sx q[2];
rz(-3.0974168) q[2];
sx q[2];
rz(0.046048306) q[2];
rz(2.0264528) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(-2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.030877) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(0.16511551) q[0];
rz(-1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(1.4547691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063140537) q[0];
sx q[0];
rz(-1.8179896) q[0];
sx q[0];
rz(1.4726223) q[0];
rz(-0.42368453) q[2];
sx q[2];
rz(-2.6657145) q[2];
sx q[2];
rz(1.0591454) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7746209) q[1];
sx q[1];
rz(-2.4200388) q[1];
sx q[1];
rz(-1.5040892) q[1];
rz(-pi) q[2];
rz(-2.1490578) q[3];
sx q[3];
rz(-2.6446178) q[3];
sx q[3];
rz(2.9013036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7016865) q[2];
sx q[2];
rz(-0.88752282) q[2];
sx q[2];
rz(2.9627724) q[2];
rz(-1.6245101) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(1.4823683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4894067) q[0];
sx q[0];
rz(-1.4302197) q[0];
sx q[0];
rz(1.408761) q[0];
rz(0.57213712) q[1];
sx q[1];
rz(-1.9913543) q[1];
sx q[1];
rz(-2.9170654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8978724) q[0];
sx q[0];
rz(-2.2099054) q[0];
sx q[0];
rz(1.0963109) q[0];
x q[1];
rz(2.8357887) q[2];
sx q[2];
rz(-2.282838) q[2];
sx q[2];
rz(2.834897) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4281339) q[1];
sx q[1];
rz(-1.1468219) q[1];
sx q[1];
rz(-2.3522373) q[1];
rz(-pi) q[2];
rz(1.8711617) q[3];
sx q[3];
rz(-1.1704418) q[3];
sx q[3];
rz(0.8308691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81395927) q[2];
sx q[2];
rz(-1.9093134) q[2];
sx q[2];
rz(0.28780469) q[2];
rz(0.52590251) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(-0.56720081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41878685) q[0];
sx q[0];
rz(-0.64544183) q[0];
sx q[0];
rz(-0.01471113) q[0];
rz(2.9171433) q[1];
sx q[1];
rz(-1.4597471) q[1];
sx q[1];
rz(0.86123484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4968729) q[0];
sx q[0];
rz(-1.5278491) q[0];
sx q[0];
rz(1.5687902) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1003816) q[2];
sx q[2];
rz(-1.2361666) q[2];
sx q[2];
rz(-1.24805) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3337832) q[1];
sx q[1];
rz(-2.2453866) q[1];
sx q[1];
rz(-1.9476931) q[1];
rz(-pi) q[2];
rz(-0.16448824) q[3];
sx q[3];
rz(-0.93014088) q[3];
sx q[3];
rz(1.3817087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33846778) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(0.11094805) q[2];
rz(2.0153996) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(2.4748928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647389) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.2891084) q[0];
rz(1.3339169) q[1];
sx q[1];
rz(-1.3197118) q[1];
sx q[1];
rz(-1.8905554) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7374501) q[0];
sx q[0];
rz(-0.92094983) q[0];
sx q[0];
rz(2.823816) q[0];
x q[1];
rz(1.8952605) q[2];
sx q[2];
rz(-2.1212541) q[2];
sx q[2];
rz(0.40263995) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2062155) q[1];
sx q[1];
rz(-2.4316129) q[1];
sx q[1];
rz(1.9496656) q[1];
rz(-pi) q[2];
rz(0.23836191) q[3];
sx q[3];
rz(-1.8619478) q[3];
sx q[3];
rz(1.0358159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.053607792) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(-1.717022) q[2];
rz(0.21814957) q[3];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394534) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(-0.46808991) q[0];
rz(2.3098436) q[1];
sx q[1];
rz(-2.2305198) q[1];
sx q[1];
rz(2.1591878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32192445) q[0];
sx q[0];
rz(-1.6409994) q[0];
sx q[0];
rz(1.6353668) q[0];
rz(0.94865145) q[2];
sx q[2];
rz(-1.4013883) q[2];
sx q[2];
rz(-2.624315) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0259798) q[1];
sx q[1];
rz(-1.8048688) q[1];
sx q[1];
rz(1.8744517) q[1];
x q[2];
rz(-1.8194852) q[3];
sx q[3];
rz(-0.2137972) q[3];
sx q[3];
rz(-0.71960654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19782797) q[2];
sx q[2];
rz(-0.48913893) q[2];
sx q[2];
rz(-0.93734199) q[2];
rz(1.7249974) q[3];
sx q[3];
rz(-1.849544) q[3];
sx q[3];
rz(-2.9924801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(2.9008289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69963843) q[0];
sx q[0];
rz(-2.505971) q[0];
sx q[0];
rz(2.5997396) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0006311) q[2];
sx q[2];
rz(-0.79509622) q[2];
sx q[2];
rz(-1.5305331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4283838) q[1];
sx q[1];
rz(-1.2782103) q[1];
sx q[1];
rz(-2.2702515) q[1];
rz(-0.82262294) q[3];
sx q[3];
rz(-2.6634376) q[3];
sx q[3];
rz(-1.9229696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9869438) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(-1.1692952) q[2];
rz(-2.5985006) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(-2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(2.0277297) q[1];
sx q[1];
rz(-1.056347) q[1];
sx q[1];
rz(-0.14811555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8982142) q[0];
sx q[0];
rz(-2.4442456) q[0];
sx q[0];
rz(-0.298907) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4101548) q[2];
sx q[2];
rz(-0.88692188) q[2];
sx q[2];
rz(-1.6647673) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9657987) q[1];
sx q[1];
rz(-2.0290004) q[1];
sx q[1];
rz(-2.5318164) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.801228) q[3];
sx q[3];
rz(-0.90022579) q[3];
sx q[3];
rz(-2.5910395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0327882) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(-0.27842251) q[2];
rz(-1.6935211) q[3];
sx q[3];
rz(-1.4689987) q[3];
sx q[3];
rz(-1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447822) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(-1.0983559) q[1];
sx q[1];
rz(-1.8041246) q[1];
sx q[1];
rz(0.95380797) q[1];
rz(-2.5658812) q[2];
sx q[2];
rz(-1.4288708) q[2];
sx q[2];
rz(2.5484011) q[2];
rz(1.9270509) q[3];
sx q[3];
rz(-2.9976525) q[3];
sx q[3];
rz(0.83565686) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
