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
rz(-2.4401378) q[0];
sx q[0];
rz(-0.48818809) q[0];
sx q[0];
rz(0.39946431) q[0];
rz(2.0692628) q[1];
sx q[1];
rz(-0.89858276) q[1];
sx q[1];
rz(-2.8314765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7061569) q[0];
sx q[0];
rz(-2.6118267) q[0];
sx q[0];
rz(0.72491531) q[0];
x q[1];
rz(0.69930716) q[2];
sx q[2];
rz(-2.0607053) q[2];
sx q[2];
rz(1.8037947) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9609144) q[1];
sx q[1];
rz(-0.49094683) q[1];
sx q[1];
rz(0.21638201) q[1];
x q[2];
rz(-1.8350527) q[3];
sx q[3];
rz(-0.39723165) q[3];
sx q[3];
rz(-0.63069944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.255379) q[2];
sx q[2];
rz(-1.9510521) q[2];
sx q[2];
rz(-1.0211241) q[2];
rz(-1.4897664) q[3];
sx q[3];
rz(-1.3083873) q[3];
sx q[3];
rz(-2.3511353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0609695) q[0];
sx q[0];
rz(-1.0401833) q[0];
sx q[0];
rz(2.8676497) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(-0.23987548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2027953) q[0];
sx q[0];
rz(-0.99730856) q[0];
sx q[0];
rz(0.62173642) q[0];
rz(-pi) q[1];
rz(-2.657183) q[2];
sx q[2];
rz(-2.5229215) q[2];
sx q[2];
rz(1.8687488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66544724) q[1];
sx q[1];
rz(-2.7283333) q[1];
sx q[1];
rz(-2.5422417) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20407014) q[3];
sx q[3];
rz(-2.3490454) q[3];
sx q[3];
rz(-2.8784405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46951744) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(-1.391927) q[2];
rz(-0.64432708) q[3];
sx q[3];
rz(-1.143012) q[3];
sx q[3];
rz(-2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7740087) q[0];
sx q[0];
rz(-2.3180361) q[0];
sx q[0];
rz(-2.9665663) q[0];
rz(-1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-2.3451436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38742313) q[0];
sx q[0];
rz(-1.8909374) q[0];
sx q[0];
rz(-2.8693136) q[0];
rz(-pi) q[1];
rz(1.8548707) q[2];
sx q[2];
rz(-2.1518118) q[2];
sx q[2];
rz(-0.43339455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63523102) q[1];
sx q[1];
rz(-2.0331675) q[1];
sx q[1];
rz(1.1504786) q[1];
x q[2];
rz(-0.50694179) q[3];
sx q[3];
rz(-1.6569417) q[3];
sx q[3];
rz(2.994842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.986787) q[2];
sx q[2];
rz(-3.0974168) q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1107156) q[0];
sx q[0];
rz(-0.59030384) q[0];
sx q[0];
rz(2.9764771) q[0];
rz(-1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(1.4547691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4835614) q[0];
sx q[0];
rz(-1.4756157) q[0];
sx q[0];
rz(0.24834085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3619992) q[2];
sx q[2];
rz(-2.0016124) q[2];
sx q[2];
rz(0.58961678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45572051) q[1];
sx q[1];
rz(-0.8511976) q[1];
sx q[1];
rz(3.0830129) q[1];
x q[2];
rz(0.28819542) q[3];
sx q[3];
rz(-1.1600947) q[3];
sx q[3];
rz(0.87897838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4399061) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-2.9627724) q[2];
rz(-1.5170826) q[3];
sx q[3];
rz(-0.30859083) q[3];
sx q[3];
rz(-1.4823683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.1502384) q[1];
sx q[1];
rz(0.22452721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1118029) q[0];
sx q[0];
rz(-1.9462612) q[0];
sx q[0];
rz(-2.445604) q[0];
rz(-0.30580394) q[2];
sx q[2];
rz(-0.85875466) q[2];
sx q[2];
rz(0.30669566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8967486) q[1];
sx q[1];
rz(-2.2678284) q[1];
sx q[1];
rz(2.5752707) q[1];
rz(-pi) q[2];
rz(-1.270431) q[3];
sx q[3];
rz(-1.1704418) q[3];
sx q[3];
rz(0.8308691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3276334) q[2];
sx q[2];
rz(-1.9093134) q[2];
sx q[2];
rz(-0.28780469) q[2];
rz(-2.6156901) q[3];
sx q[3];
rz(-1.9764427) q[3];
sx q[3];
rz(-2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.7228058) q[0];
sx q[0];
rz(-2.4961508) q[0];
sx q[0];
rz(0.01471113) q[0];
rz(2.9171433) q[1];
sx q[1];
rz(-1.4597471) q[1];
sx q[1];
rz(0.86123484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5435641) q[0];
sx q[0];
rz(-0.042994067) q[0];
sx q[0];
rz(-0.046648101) q[0];
rz(-pi) q[1];
rz(-1.4528571) q[2];
sx q[2];
rz(-2.8045296) q[2];
sx q[2];
rz(1.7686421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3726681) q[1];
sx q[1];
rz(-0.75809352) q[1];
sx q[1];
rz(-2.7102986) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3545996) q[3];
sx q[3];
rz(-2.4830468) q[3];
sx q[3];
rz(-1.4890081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8031249) q[2];
sx q[2];
rz(-0.71102342) q[2];
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
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768537) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.8524843) q[0];
rz(1.8076757) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(-1.8905554) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40414251) q[0];
sx q[0];
rz(-0.92094983) q[0];
sx q[0];
rz(-2.823816) q[0];
rz(-pi) q[1];
rz(-2.5669615) q[2];
sx q[2];
rz(-1.2956276) q[2];
sx q[2];
rz(-0.99400101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45197645) q[1];
sx q[1];
rz(-0.9202846) q[1];
sx q[1];
rz(2.8337949) q[1];
rz(-pi) q[2];
rz(-1.2716633) q[3];
sx q[3];
rz(-1.7989446) q[3];
sx q[3];
rz(-2.6762485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053607792) q[2];
sx q[2];
rz(-1.3259652) q[2];
sx q[2];
rz(-1.717022) q[2];
rz(2.9234431) q[3];
sx q[3];
rz(-1.679922) q[3];
sx q[3];
rz(-0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3394534) q[0];
sx q[0];
rz(-0.55291432) q[0];
sx q[0];
rz(-0.46808991) q[0];
rz(-0.83174902) q[1];
sx q[1];
rz(-0.91107285) q[1];
sx q[1];
rz(-2.1591878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8196682) q[0];
sx q[0];
rz(-1.5005932) q[0];
sx q[0];
rz(-1.5062259) q[0];
rz(2.1929412) q[2];
sx q[2];
rz(-1.4013883) q[2];
sx q[2];
rz(2.624315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.959572) q[1];
sx q[1];
rz(-2.7604155) q[1];
sx q[1];
rz(-0.89761727) q[1];
rz(-pi) q[2];
rz(-1.8194852) q[3];
sx q[3];
rz(-0.2137972) q[3];
sx q[3];
rz(2.4219861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19782797) q[2];
sx q[2];
rz(-2.6524537) q[2];
sx q[2];
rz(0.93734199) q[2];
rz(1.7249974) q[3];
sx q[3];
rz(-1.849544) q[3];
sx q[3];
rz(0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.8775403) q[0];
sx q[0];
rz(-0.5101246) q[0];
sx q[0];
rz(0.75793761) q[0];
rz(-1.0559878) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(-2.9008289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42007459) q[0];
sx q[0];
rz(-1.2596247) q[0];
sx q[0];
rz(-2.5779186) q[0];
x q[1];
rz(-2.3182436) q[2];
sx q[2];
rz(-1.2687131) q[2];
sx q[2];
rz(0.27034098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1878464) q[1];
sx q[1];
rz(-0.74855474) q[1];
sx q[1];
rz(-2.0084318) q[1];
rz(-0.82262294) q[3];
sx q[3];
rz(-2.6634376) q[3];
sx q[3];
rz(1.218623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15464887) q[2];
sx q[2];
rz(-1.1661531) q[2];
sx q[2];
rz(1.9722975) q[2];
rz(-2.5985006) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(0.57524601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.52694046) q[0];
sx q[0];
rz(-1.4450547) q[0];
sx q[0];
rz(0.27969435) q[0];
rz(1.1138629) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(2.9934771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0954656) q[0];
sx q[0];
rz(-1.3805423) q[0];
sx q[0];
rz(-2.4664761) q[0];
rz(-2.947817) q[2];
sx q[2];
rz(-0.69952088) q[2];
sx q[2];
rz(1.7278838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9657987) q[1];
sx q[1];
rz(-1.1125922) q[1];
sx q[1];
rz(0.60977625) q[1];
rz(1.3403647) q[3];
sx q[3];
rz(-2.2413669) q[3];
sx q[3];
rz(2.5910395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1088045) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(2.8631701) q[2];
rz(1.6935211) q[3];
sx q[3];
rz(-1.4689987) q[3];
sx q[3];
rz(1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096810452) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(-2.0432368) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(-1.4020709) q[2];
sx q[2];
rz(-2.1399956) q[2];
sx q[2];
rz(-2.2555399) q[2];
rz(1.2145417) q[3];
sx q[3];
rz(-0.14394017) q[3];
sx q[3];
rz(-2.3059358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
