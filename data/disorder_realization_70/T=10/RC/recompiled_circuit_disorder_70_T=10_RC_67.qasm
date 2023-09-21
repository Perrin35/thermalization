OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2848628) q[0];
sx q[0];
rz(-2.3861109) q[0];
sx q[0];
rz(1.5289686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2251031) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-2.8061342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69971985) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(-0.91200478) q[1];
rz(1.5942469) q[3];
sx q[3];
rz(-1.1707077) q[3];
sx q[3];
rz(-1.7398906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.263164) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(-2.9204869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1572004) q[0];
sx q[0];
rz(-1.7565787) q[0];
sx q[0];
rz(-2.0322582) q[0];
x q[1];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(-0.32804104) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0529777) q[1];
sx q[1];
rz(-2.6664475) q[1];
sx q[1];
rz(0.26762025) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(2.6925987) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-2.0670481) q[1];
sx q[1];
rz(0.75769889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66540668) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(0.55253367) q[0];
rz(-0.65501113) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(-2.8734145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.12049) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-0.91336577) q[1];
x q[2];
rz(1.4593616) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(-1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-0.20382717) q[3];
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
rz(pi/2) q[2];
sx q[2];
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
rz(-0.42235342) q[1];
sx q[1];
rz(-1.8130594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0275832) q[0];
sx q[0];
rz(-2.8767577) q[0];
sx q[0];
rz(-0.94075216) q[0];
x q[1];
rz(0.76339108) q[2];
sx q[2];
rz(-1.9420468) q[2];
sx q[2];
rz(0.064918092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1991785) q[1];
sx q[1];
rz(-1.5702827) q[1];
sx q[1];
rz(-1.3039939) q[1];
rz(1.7926746) q[3];
sx q[3];
rz(-1.6776049) q[3];
sx q[3];
rz(-0.91339236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(2.034534) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040314019) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(-0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-0.50450528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550068) q[0];
sx q[0];
rz(-1.803777) q[0];
sx q[0];
rz(1.0250807) q[0];
rz(-pi) q[1];
rz(1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(-1.2010241) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5704755) q[1];
sx q[1];
rz(-1.8091396) q[1];
sx q[1];
rz(-1.2586602) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2284331) q[3];
sx q[3];
rz(-2.1214161) q[3];
sx q[3];
rz(0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.7593613) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(2.0507623) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-0.18879034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15121962) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(1.3034526) q[0];
x q[1];
rz(3.0374755) q[2];
sx q[2];
rz(-1.8407485) q[2];
sx q[2];
rz(-1.5013258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61154762) q[1];
sx q[1];
rz(-1.6366742) q[1];
sx q[1];
rz(-2.8987315) q[1];
rz(-pi) q[2];
rz(-1.3879697) q[3];
sx q[3];
rz(-1.332924) q[3];
sx q[3];
rz(-0.49926234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-2.3110456) q[0];
sx q[0];
rz(2.8175957) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24810476) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(-1.8155314) q[0];
rz(-pi) q[1];
rz(2.4263779) q[2];
sx q[2];
rz(-0.5826125) q[2];
sx q[2];
rz(-1.4441393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5768347) q[1];
sx q[1];
rz(-1.733629) q[1];
sx q[1];
rz(0.30270438) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2077683) q[3];
sx q[3];
rz(-1.0855506) q[3];
sx q[3];
rz(-2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.8483298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9613567) q[0];
sx q[0];
rz(-2.0016252) q[0];
sx q[0];
rz(0.21813099) q[0];
rz(3.0027577) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(2.5790737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5544719) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(-0.15092571) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44869081) q[3];
sx q[3];
rz(-2.7033227) q[3];
sx q[3];
rz(1.1718307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3747037) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(2.5602706) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.2506437) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8563961) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(1.517209) q[0];
x q[1];
rz(1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(-2.929504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7848283) q[1];
sx q[1];
rz(-0.67912662) q[1];
sx q[1];
rz(-1.2374415) q[1];
x q[2];
rz(-2.3171114) q[3];
sx q[3];
rz(-1.5319676) q[3];
sx q[3];
rz(-2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-1.0725718) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(-1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(-2.7744055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8242278) q[0];
sx q[0];
rz(-1.5396848) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7634002) q[2];
sx q[2];
rz(-2.4792255) q[2];
sx q[2];
rz(2.3364002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91017427) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(-0.66399666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(-1.1993053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24511589) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(1.4916186) q[2];
sx q[2];
rz(-0.82293246) q[2];
sx q[2];
rz(-0.056597829) q[2];
rz(2.9885837) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
