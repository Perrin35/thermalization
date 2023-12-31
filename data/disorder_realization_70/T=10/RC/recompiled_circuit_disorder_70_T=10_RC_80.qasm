OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6887309) q[0];
sx q[0];
rz(-2.9714669) q[0];
sx q[0];
rz(-2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85672985) q[0];
sx q[0];
rz(-2.3861109) q[0];
sx q[0];
rz(-1.612624) q[0];
rz(-pi) q[1];
rz(1.6526821) q[2];
sx q[2];
rz(-0.6559283) q[2];
sx q[2];
rz(1.3002849) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4418728) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(-0.91200478) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(-2.9633629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(2.5105387) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9843922) q[0];
sx q[0];
rz(-1.7565787) q[0];
sx q[0];
rz(1.1093344) q[0];
x q[1];
rz(1.2220076) q[2];
sx q[2];
rz(-1.6752599) q[2];
sx q[2];
rz(-1.8609357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088614956) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(2.8739724) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6392194) q[3];
sx q[3];
rz(-1.4851735) q[3];
sx q[3];
rz(-1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9235886) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(2.5533) q[2];
rz(0.44899392) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56851971) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(0.75769889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.889567) q[0];
sx q[0];
rz(-0.55253367) q[0];
rz(-pi) q[1];
rz(-2.9341142) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(-0.66397882) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.12049) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
rz(-pi) q[2];
rz(-2.5030701) q[3];
sx q[3];
rz(-1.4811852) q[3];
sx q[3];
rz(1.5479969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(2.3051252) q[2];
rz(-1.6992735) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(0.077443667) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(-1.3285332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.3550718) q[0];
x q[1];
rz(-2.3782016) q[2];
sx q[2];
rz(-1.1995458) q[2];
sx q[2];
rz(3.0766746) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9424142) q[1];
sx q[1];
rz(-1.5702827) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(-0.10947157) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040314019) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(-1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-0.50450528) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0230334) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(-2.8708007) q[0];
rz(1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(-1.2010241) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5704755) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(1.8829324) q[1];
rz(2.4818146) q[3];
sx q[3];
rz(-1.0228844) q[3];
sx q[3];
rz(1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.4939235) q[2];
rz(-1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(-1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81297368) q[0];
sx q[0];
rz(-0.44437528) q[0];
sx q[0];
rz(-2.5286753) q[0];
x q[1];
rz(1.2994453) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(-0.041610418) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69953883) q[1];
sx q[1];
rz(-2.8901254) q[1];
sx q[1];
rz(0.26775189) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8998428) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(-1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.3791929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1349072) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(1.8968614) q[0];
x q[1];
rz(1.9786644) q[2];
sx q[2];
rz(-1.9991572) q[2];
sx q[2];
rz(-2.5025764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52723253) q[1];
sx q[1];
rz(-0.34253201) q[1];
sx q[1];
rz(0.50369461) q[1];
rz(-pi) q[2];
rz(-2.561065) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(2.011727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.8483298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6587257) q[0];
sx q[0];
rz(-1.3728766) q[0];
sx q[0];
rz(2.0107962) q[0];
x q[1];
rz(0.60111945) q[2];
sx q[2];
rz(-1.4919315) q[2];
sx q[2];
rz(-2.2476946) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9570436) q[1];
sx q[1];
rz(-0.42889412) q[1];
sx q[1];
rz(-1.9098319) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.771365) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.2506437) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8563961) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(-1.517209) q[0];
rz(-pi) q[1];
rz(0.48607562) q[2];
sx q[2];
rz(-1.7239778) q[2];
sx q[2];
rz(-1.4391522) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(-2.8833564) q[1];
x q[2];
rz(2.3171114) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(2.0064158) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(0.36718711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627535) q[0];
sx q[0];
rz(-0.05719962) q[0];
sx q[0];
rz(0.57533933) q[0];
rz(-pi) q[1];
rz(2.514421) q[2];
sx q[2];
rz(-1.7998724) q[2];
sx q[2];
rz(-1.0692182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2314184) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(2.477596) q[1];
rz(-pi) q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(-1.184954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8463678) q[0];
sx q[0];
rz(-1.5867148) q[0];
sx q[0];
rz(-1.5750194) q[0];
rz(-0.29905839) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(-1.4916186) q[2];
sx q[2];
rz(-2.3186602) q[2];
sx q[2];
rz(3.0849948) q[2];
rz(-0.84898938) q[3];
sx q[3];
rz(-1.6860387) q[3];
sx q[3];
rz(0.54308346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
