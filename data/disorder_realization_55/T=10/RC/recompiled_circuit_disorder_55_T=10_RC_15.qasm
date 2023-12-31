OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(0.2642785) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(2.4063453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3448479) q[0];
sx q[0];
rz(-1.42294) q[0];
sx q[0];
rz(-0.66230886) q[0];
x q[1];
rz(0.60230435) q[2];
sx q[2];
rz(-2.3659083) q[2];
sx q[2];
rz(1.3210981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.211261) q[1];
sx q[1];
rz(-0.32713612) q[1];
sx q[1];
rz(2.0698018) q[1];
rz(-pi) q[2];
rz(2.5611476) q[3];
sx q[3];
rz(-1.5942897) q[3];
sx q[3];
rz(1.1593727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(1.1038587) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(-0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274149) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(0.46288681) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-0.26611051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-0.87478144) q[0];
sx q[0];
rz(-0.50177411) q[0];
rz(-0.99533178) q[2];
sx q[2];
rz(-1.4333712) q[2];
sx q[2];
rz(2.1910138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38763603) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(2.3073879) q[1];
x q[2];
rz(-0.18272419) q[3];
sx q[3];
rz(-0.97191873) q[3];
sx q[3];
rz(0.11174186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(0.51149386) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-2.8539343) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(-1.7354895) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12362326) q[0];
sx q[0];
rz(-2.6051913) q[0];
sx q[0];
rz(0.53039741) q[0];
rz(-1.0727097) q[2];
sx q[2];
rz(-0.32215298) q[2];
sx q[2];
rz(-2.1307532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82145065) q[1];
sx q[1];
rz(-0.67678932) q[1];
sx q[1];
rz(-2.0435145) q[1];
x q[2];
rz(0.8637572) q[3];
sx q[3];
rz(-1.6861526) q[3];
sx q[3];
rz(2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(-1.5220801) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(2.6456397) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494444) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(0.55975634) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20630079) q[0];
sx q[0];
rz(-1.9959404) q[0];
sx q[0];
rz(1.4817119) q[0];
x q[1];
rz(-0.83300029) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(-2.8766362) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2019129) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(1.4918785) q[1];
x q[2];
rz(-2.9066554) q[3];
sx q[3];
rz(-0.80312356) q[3];
sx q[3];
rz(-2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-2.7752303) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(2.7979134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77728358) q[0];
sx q[0];
rz(-2.1099578) q[0];
sx q[0];
rz(2.3794412) q[0];
x q[1];
rz(-2.2291549) q[2];
sx q[2];
rz(-2.3235333) q[2];
sx q[2];
rz(-0.91615265) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1829454) q[1];
sx q[1];
rz(-1.6445541) q[1];
sx q[1];
rz(2.069509) q[1];
x q[2];
rz(-1.5991391) q[3];
sx q[3];
rz(-0.64405555) q[3];
sx q[3];
rz(1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3917824) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(-0.053744944) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4869726) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(-0.25973928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1007337) q[0];
sx q[0];
rz(-1.20964) q[0];
sx q[0];
rz(2.8054537) q[0];
rz(-pi) q[1];
rz(0.01979205) q[2];
sx q[2];
rz(-1.2345825) q[2];
sx q[2];
rz(-2.3078231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4252792) q[1];
sx q[1];
rz(-1.3448161) q[1];
sx q[1];
rz(1.2334361) q[1];
rz(-pi) q[2];
rz(1.1214244) q[3];
sx q[3];
rz(-1.5546038) q[3];
sx q[3];
rz(1.3390954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-0.10061131) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(-1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(-2.8314262) q[0];
rz(0.50225964) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(-0.60595864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46442407) q[0];
sx q[0];
rz(-0.66093984) q[0];
sx q[0];
rz(0.37207694) q[0];
rz(-pi) q[1];
rz(-0.59553869) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(-0.66875848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.052918606) q[1];
sx q[1];
rz(-0.62401906) q[1];
sx q[1];
rz(-2.1799929) q[1];
rz(-pi) q[2];
rz(-2.1171655) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(1.446561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(0.85285464) q[2];
rz(1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(-2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.461331) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(3.0775552) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-2.1770283) q[0];
sx q[0];
rz(-0.22867815) q[0];
rz(-pi) q[1];
rz(-2.0766826) q[2];
sx q[2];
rz(-1.3008899) q[2];
sx q[2];
rz(3.0907061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36754164) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(-0.35518412) q[1];
rz(-pi) q[2];
rz(-3.1021032) q[3];
sx q[3];
rz(-1.2590623) q[3];
sx q[3];
rz(1.5453651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5861661) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-2.9973307) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(0.56888467) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(-1.127839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147946) q[0];
sx q[0];
rz(-0.1495805) q[0];
sx q[0];
rz(3.0339255) q[0];
rz(-pi) q[1];
rz(-2.1986507) q[2];
sx q[2];
rz(-2.0161511) q[2];
sx q[2];
rz(-3.1415591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7833264) q[1];
sx q[1];
rz(-1.8124541) q[1];
sx q[1];
rz(3.1296455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0117709) q[3];
sx q[3];
rz(-1.2563007) q[3];
sx q[3];
rz(2.3712456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(0.17364994) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-0.5724268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6515598) q[0];
sx q[0];
rz(-2.3736554) q[0];
sx q[0];
rz(2.4134273) q[0];
rz(-2.0296713) q[2];
sx q[2];
rz(-1.8177114) q[2];
sx q[2];
rz(1.0690451) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13092566) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(-2.8972577) q[1];
x q[2];
rz(1.278986) q[3];
sx q[3];
rz(-2.2581873) q[3];
sx q[3];
rz(0.22112267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(1.6981381) q[2];
sx q[2];
rz(-0.6586532) q[2];
sx q[2];
rz(-2.3011343) q[2];
rz(2.0588856) q[3];
sx q[3];
rz(-1.516021) q[3];
sx q[3];
rz(1.0513023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
