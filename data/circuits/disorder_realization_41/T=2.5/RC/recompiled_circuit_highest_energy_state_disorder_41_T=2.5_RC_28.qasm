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
rz(-0.37201878) q[0];
sx q[0];
rz(-2.7899185) q[0];
sx q[0];
rz(-0.055211842) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(-0.13394314) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3833544) q[0];
sx q[0];
rz(-1.3516434) q[0];
sx q[0];
rz(1.4683649) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8419398) q[2];
sx q[2];
rz(-2.5404394) q[2];
sx q[2];
rz(-0.047182949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3250934) q[1];
sx q[1];
rz(-1.2383922) q[1];
sx q[1];
rz(1.2219882) q[1];
rz(3.1024801) q[3];
sx q[3];
rz(-2.0958423) q[3];
sx q[3];
rz(-1.0421747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.967531) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(-2.3270712) q[2];
rz(2.9461765) q[3];
sx q[3];
rz(-2.312909) q[3];
sx q[3];
rz(2.101208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8622417) q[0];
sx q[0];
rz(-0.26222721) q[0];
sx q[0];
rz(2.1694515) q[0];
rz(-0.60130087) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(-1.7960637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1336254) q[0];
sx q[0];
rz(-1.04038) q[0];
sx q[0];
rz(-2.6625457) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53671543) q[2];
sx q[2];
rz(-1.3379352) q[2];
sx q[2];
rz(0.95271207) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3516772) q[1];
sx q[1];
rz(-0.43444217) q[1];
sx q[1];
rz(-1.2277568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18109326) q[3];
sx q[3];
rz(-1.5646184) q[3];
sx q[3];
rz(0.97038499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1874275) q[2];
sx q[2];
rz(-1.2792055) q[2];
sx q[2];
rz(-2.1698451) q[2];
rz(-1.0176954) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(-3.0903604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.63075066) q[0];
sx q[0];
rz(-2.172281) q[0];
sx q[0];
rz(0.72072679) q[0];
rz(-3.0156056) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(-2.7350977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46716786) q[0];
sx q[0];
rz(-0.95646383) q[0];
sx q[0];
rz(-2.7697639) q[0];
x q[1];
rz(-2.3977621) q[2];
sx q[2];
rz(-1.1421428) q[2];
sx q[2];
rz(-0.33996087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6449922) q[1];
sx q[1];
rz(-0.82372249) q[1];
sx q[1];
rz(-0.83622593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2062827) q[3];
sx q[3];
rz(-1.3008833) q[3];
sx q[3];
rz(-2.4382537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5410109) q[2];
sx q[2];
rz(-1.7966248) q[2];
sx q[2];
rz(-2.7715032) q[2];
rz(-0.078941405) q[3];
sx q[3];
rz(-0.97349662) q[3];
sx q[3];
rz(0.48648849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9321891) q[0];
sx q[0];
rz(-1.3511193) q[0];
sx q[0];
rz(-2.6655647) q[0];
rz(-2.027482) q[1];
sx q[1];
rz(-2.1962491) q[1];
sx q[1];
rz(1.5155189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1080804) q[0];
sx q[0];
rz(-2.7191504) q[0];
sx q[0];
rz(-2.6080934) q[0];
rz(-pi) q[1];
rz(2.5427688) q[2];
sx q[2];
rz(-2.235417) q[2];
sx q[2];
rz(-1.40203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0994708) q[1];
sx q[1];
rz(-1.9212706) q[1];
sx q[1];
rz(-1.6609236) q[1];
rz(-2.2882624) q[3];
sx q[3];
rz(-1.322041) q[3];
sx q[3];
rz(-2.0847818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6473306) q[2];
sx q[2];
rz(-1.6202972) q[2];
sx q[2];
rz(-2.1935513) q[2];
rz(2.7028132) q[3];
sx q[3];
rz(-0.56883562) q[3];
sx q[3];
rz(2.2426864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(-2.7137252) q[0];
rz(-0.95871344) q[1];
sx q[1];
rz(-1.9045279) q[1];
sx q[1];
rz(-1.0520891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31715441) q[0];
sx q[0];
rz(-2.053932) q[0];
sx q[0];
rz(0.89886959) q[0];
rz(-pi) q[1];
rz(-2.1053931) q[2];
sx q[2];
rz(-1.2367168) q[2];
sx q[2];
rz(0.47185791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8278049) q[1];
sx q[1];
rz(-2.1422285) q[1];
sx q[1];
rz(3.0302583) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73486272) q[3];
sx q[3];
rz(-2.4331193) q[3];
sx q[3];
rz(2.3264309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.170257) q[2];
sx q[2];
rz(-1.3054138) q[2];
sx q[2];
rz(0.35856426) q[2];
rz(-1.9606918) q[3];
sx q[3];
rz(-0.87555331) q[3];
sx q[3];
rz(-2.6023279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81826687) q[0];
sx q[0];
rz(-1.8805255) q[0];
sx q[0];
rz(0.64875025) q[0];
rz(0.58748856) q[1];
sx q[1];
rz(-1.3222539) q[1];
sx q[1];
rz(2.9041451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64024583) q[0];
sx q[0];
rz(-3.0398439) q[0];
sx q[0];
rz(2.8215088) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.468979) q[2];
sx q[2];
rz(-1.8378496) q[2];
sx q[2];
rz(1.36324) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9623649) q[1];
sx q[1];
rz(-0.74848524) q[1];
sx q[1];
rz(-0.039123936) q[1];
rz(-pi) q[2];
rz(2.1319207) q[3];
sx q[3];
rz(-2.3497407) q[3];
sx q[3];
rz(0.71969024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6610403) q[2];
sx q[2];
rz(-2.48017) q[2];
sx q[2];
rz(-2.1858369) q[2];
rz(-1.3794948) q[3];
sx q[3];
rz(-2.5199514) q[3];
sx q[3];
rz(-0.1563589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697923) q[0];
sx q[0];
rz(-0.0026230165) q[0];
sx q[0];
rz(-0.045259137) q[0];
rz(-2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(-2.9387567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58205523) q[0];
sx q[0];
rz(-1.4180611) q[0];
sx q[0];
rz(2.5688897) q[0];
rz(-pi) q[1];
rz(1.3812106) q[2];
sx q[2];
rz(-0.84824296) q[2];
sx q[2];
rz(-1.0740395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5109628) q[1];
sx q[1];
rz(-1.4893526) q[1];
sx q[1];
rz(0.54334531) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84879227) q[3];
sx q[3];
rz(-1.5986852) q[3];
sx q[3];
rz(1.3031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63991919) q[2];
sx q[2];
rz(-1.4036274) q[2];
sx q[2];
rz(-0.61140927) q[2];
rz(2.0407138) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(-2.8618405) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6708577) q[0];
sx q[0];
rz(-1.9453456) q[0];
sx q[0];
rz(1.083495) q[0];
rz(2.1776543) q[1];
sx q[1];
rz(-1.0716535) q[1];
sx q[1];
rz(0.49577698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475233) q[0];
sx q[0];
rz(-0.076193245) q[0];
sx q[0];
rz(1.2890362) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18270418) q[2];
sx q[2];
rz(-1.2983381) q[2];
sx q[2];
rz(0.24531103) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3246877) q[1];
sx q[1];
rz(-0.73313112) q[1];
sx q[1];
rz(-2.2156961) q[1];
rz(2.1585309) q[3];
sx q[3];
rz(-1.5957629) q[3];
sx q[3];
rz(2.4189396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44020161) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(-0.68897796) q[2];
rz(1.5135328) q[3];
sx q[3];
rz(-2.3115034) q[3];
sx q[3];
rz(1.0015063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643519) q[0];
sx q[0];
rz(-1.5761201) q[0];
sx q[0];
rz(-1.0871357) q[0];
rz(-0.76308909) q[1];
sx q[1];
rz(-1.3975846) q[1];
sx q[1];
rz(0.22689247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339461) q[0];
sx q[0];
rz(-1.3469704) q[0];
sx q[0];
rz(-0.86854078) q[0];
x q[1];
rz(-0.82590286) q[2];
sx q[2];
rz(-2.7159198) q[2];
sx q[2];
rz(0.22944006) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5530738) q[1];
sx q[1];
rz(-1.5344193) q[1];
sx q[1];
rz(0.22471551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1710728) q[3];
sx q[3];
rz(-1.9991989) q[3];
sx q[3];
rz(0.71479931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6202966) q[2];
sx q[2];
rz(-2.2038348) q[2];
sx q[2];
rz(-2.1916981) q[2];
rz(1.0319483) q[3];
sx q[3];
rz(-2.2622006) q[3];
sx q[3];
rz(-2.7040645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49143219) q[0];
sx q[0];
rz(-3.037945) q[0];
sx q[0];
rz(-1.9895122) q[0];
rz(3.0488455) q[1];
sx q[1];
rz(-0.68789613) q[1];
sx q[1];
rz(0.50382096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9418056) q[0];
sx q[0];
rz(-1.5212987) q[0];
sx q[0];
rz(-1.0440396) q[0];
rz(-2.8724818) q[2];
sx q[2];
rz(-1.7935441) q[2];
sx q[2];
rz(2.3151223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42119869) q[1];
sx q[1];
rz(-2.6662146) q[1];
sx q[1];
rz(-0.092751547) q[1];
x q[2];
rz(-1.1931879) q[3];
sx q[3];
rz(-1.2758299) q[3];
sx q[3];
rz(1.1136536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4656226) q[2];
sx q[2];
rz(-1.114782) q[2];
sx q[2];
rz(-2.5176804) q[2];
rz(-0.55656773) q[3];
sx q[3];
rz(-2.4534295) q[3];
sx q[3];
rz(2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9482166) q[0];
sx q[0];
rz(-0.93127903) q[0];
sx q[0];
rz(-2.2122526) q[0];
rz(2.9931862) q[1];
sx q[1];
rz(-1.8971309) q[1];
sx q[1];
rz(2.94577) q[1];
rz(0.23870809) q[2];
sx q[2];
rz(-2.6663824) q[2];
sx q[2];
rz(2.0486365) q[2];
rz(1.2986167) q[3];
sx q[3];
rz(-1.6485572) q[3];
sx q[3];
rz(0.89383517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
