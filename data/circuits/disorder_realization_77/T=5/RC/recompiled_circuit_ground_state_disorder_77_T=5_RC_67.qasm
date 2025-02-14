OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(2.0030608) q[0];
sx q[0];
rz(8.3264405) q[0];
rz(-1.9810454) q[1];
sx q[1];
rz(-2.0055973) q[1];
sx q[1];
rz(-2.5392037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576958) q[0];
sx q[0];
rz(-2.0838782) q[0];
sx q[0];
rz(-2.0747572) q[0];
rz(-pi) q[1];
rz(-2.4279934) q[2];
sx q[2];
rz(-1.6125792) q[2];
sx q[2];
rz(-3.0609727) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5307926) q[1];
sx q[1];
rz(-1.6850753) q[1];
sx q[1];
rz(2.0487993) q[1];
x q[2];
rz(-1.4887047) q[3];
sx q[3];
rz(-0.42969504) q[3];
sx q[3];
rz(1.2712417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7263553) q[2];
sx q[2];
rz(-2.0550315) q[2];
sx q[2];
rz(-1.6437257) q[2];
rz(3.030153) q[3];
sx q[3];
rz(-1.8569511) q[3];
sx q[3];
rz(2.1845412) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1794671) q[0];
sx q[0];
rz(-2.2158556) q[0];
sx q[0];
rz(-2.977648) q[0];
rz(2.1666849) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.499739) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3079118) q[0];
sx q[0];
rz(-1.320086) q[0];
sx q[0];
rz(-2.3019019) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4663839) q[2];
sx q[2];
rz(-0.9281425) q[2];
sx q[2];
rz(2.2810138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30399366) q[1];
sx q[1];
rz(-1.6483867) q[1];
sx q[1];
rz(-2.4291888) q[1];
rz(-pi) q[2];
rz(-1.896656) q[3];
sx q[3];
rz(-1.3658524) q[3];
sx q[3];
rz(2.4166493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6770596) q[2];
sx q[2];
rz(-2.634282) q[2];
sx q[2];
rz(0.23059174) q[2];
rz(0.683189) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(-0.75146875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36000073) q[0];
sx q[0];
rz(-1.3104023) q[0];
sx q[0];
rz(0.42136425) q[0];
rz(1.5158117) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(-2.1741507) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7712647) q[0];
sx q[0];
rz(-1.2441085) q[0];
sx q[0];
rz(-2.8730434) q[0];
rz(-pi) q[1];
rz(1.2363628) q[2];
sx q[2];
rz(-0.49116557) q[2];
sx q[2];
rz(0.14310357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6306259) q[1];
sx q[1];
rz(-0.31019638) q[1];
sx q[1];
rz(1.0794182) q[1];
x q[2];
rz(2.580242) q[3];
sx q[3];
rz(-1.8551984) q[3];
sx q[3];
rz(-1.4903414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4030054) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(2.9834413) q[2];
rz(2.3824597) q[3];
sx q[3];
rz(-1.4594376) q[3];
sx q[3];
rz(2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247093) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(1.9547113) q[0];
rz(3.0184556) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(1.311696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13611159) q[0];
sx q[0];
rz(-1.4887047) q[0];
sx q[0];
rz(1.7186706) q[0];
rz(2.4995828) q[2];
sx q[2];
rz(-3.0138956) q[2];
sx q[2];
rz(-2.0754432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1748743) q[1];
sx q[1];
rz(-2.5481564) q[1];
sx q[1];
rz(0.75812854) q[1];
rz(-pi) q[2];
rz(-1.7100641) q[3];
sx q[3];
rz(-3.1374806) q[3];
sx q[3];
rz(2.508956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19795236) q[2];
sx q[2];
rz(-1.3200878) q[2];
sx q[2];
rz(0.8521592) q[2];
rz(0.85396829) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(-1.0796237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78344408) q[0];
sx q[0];
rz(-1.7572948) q[0];
sx q[0];
rz(-3.070991) q[0];
rz(1.0380896) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(2.3322754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22915086) q[0];
sx q[0];
rz(-2.1437867) q[0];
sx q[0];
rz(-1.4848764) q[0];
x q[1];
rz(2.2373166) q[2];
sx q[2];
rz(-0.94221409) q[2];
sx q[2];
rz(-1.6055799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0496724) q[1];
sx q[1];
rz(-0.61245239) q[1];
sx q[1];
rz(2.6785707) q[1];
x q[2];
rz(-1.7820939) q[3];
sx q[3];
rz(-0.70970067) q[3];
sx q[3];
rz(-3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(1.9522379) q[2];
rz(-0.14063028) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(-2.752059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73914948) q[0];
sx q[0];
rz(-2.2279255) q[0];
sx q[0];
rz(1.6074578) q[0];
rz(2.2588579) q[1];
sx q[1];
rz(-2.2969756) q[1];
sx q[1];
rz(-0.61558634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12703757) q[0];
sx q[0];
rz(-1.1786228) q[0];
sx q[0];
rz(-1.6704301) q[0];
rz(0.83258038) q[2];
sx q[2];
rz(-2.2924149) q[2];
sx q[2];
rz(-0.20038062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57762417) q[1];
sx q[1];
rz(-1.3805439) q[1];
sx q[1];
rz(-3.0012194) q[1];
rz(2.1911931) q[3];
sx q[3];
rz(-1.5614034) q[3];
sx q[3];
rz(-0.91419807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1659871) q[2];
sx q[2];
rz(-2.1773982) q[2];
sx q[2];
rz(-1.1791505) q[2];
rz(-2.9278582) q[3];
sx q[3];
rz(-1.390018) q[3];
sx q[3];
rz(0.24553044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24446503) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(0.044064673) q[0];
rz(0.38912082) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(-0.74180952) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76352966) q[0];
sx q[0];
rz(-2.6776095) q[0];
sx q[0];
rz(2.4876809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35105437) q[2];
sx q[2];
rz(-0.85509713) q[2];
sx q[2];
rz(2.9911016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1434113) q[1];
sx q[1];
rz(-2.8074968) q[1];
sx q[1];
rz(2.8715474) q[1];
rz(-0.70176418) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(-3.0879741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47702181) q[2];
sx q[2];
rz(-1.3571813) q[2];
sx q[2];
rz(0.50194293) q[2];
rz(-1.7160412) q[3];
sx q[3];
rz(-2.612096) q[3];
sx q[3];
rz(0.79290032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75941706) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(0.39696804) q[0];
rz(-1.5396897) q[1];
sx q[1];
rz(-0.27286467) q[1];
sx q[1];
rz(-2.4698965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8803589) q[0];
sx q[0];
rz(-1.4443928) q[0];
sx q[0];
rz(2.0194299) q[0];
rz(-pi) q[1];
rz(-0.88550466) q[2];
sx q[2];
rz(-0.54965082) q[2];
sx q[2];
rz(3.0424712) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.063350141) q[1];
sx q[1];
rz(-2.1890609) q[1];
sx q[1];
rz(1.9552014) q[1];
rz(-2.3109452) q[3];
sx q[3];
rz(-1.1201356) q[3];
sx q[3];
rz(0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.73528618) q[2];
sx q[2];
rz(-2.6267509) q[2];
sx q[2];
rz(-1.3675281) q[2];
rz(-2.1788518) q[3];
sx q[3];
rz(-1.5906426) q[3];
sx q[3];
rz(2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3108567) q[0];
sx q[0];
rz(-1.5565358) q[0];
sx q[0];
rz(0.32064015) q[0];
rz(2.2036208) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(1.7373614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2187248) q[0];
sx q[0];
rz(-1.467112) q[0];
sx q[0];
rz(-1.6746878) q[0];
x q[1];
rz(-3.0530372) q[2];
sx q[2];
rz(-0.27805304) q[2];
sx q[2];
rz(1.1403569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93130805) q[1];
sx q[1];
rz(-3.0979994) q[1];
sx q[1];
rz(-2.0491854) q[1];
rz(-0.25038428) q[3];
sx q[3];
rz(-2.3634849) q[3];
sx q[3];
rz(0.66533662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9414901) q[2];
sx q[2];
rz(-2.3788171) q[2];
sx q[2];
rz(0.74083677) q[2];
rz(-2.0790993) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(-1.9443996) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9366539) q[0];
sx q[0];
rz(-0.42977253) q[0];
sx q[0];
rz(2.3858261) q[0];
rz(-1.3142122) q[1];
sx q[1];
rz(-1.6270437) q[1];
sx q[1];
rz(-2.4155713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7814815) q[0];
sx q[0];
rz(-1.6051634) q[0];
sx q[0];
rz(-0.78256677) q[0];
rz(-pi) q[1];
rz(1.4439031) q[2];
sx q[2];
rz(-0.68918258) q[2];
sx q[2];
rz(2.5886332) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0418813) q[1];
sx q[1];
rz(-1.9181632) q[1];
sx q[1];
rz(2.9798286) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9922759) q[3];
sx q[3];
rz(-2.0155613) q[3];
sx q[3];
rz(2.5721511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3296335) q[2];
sx q[2];
rz(-0.52874955) q[2];
sx q[2];
rz(0.94919666) q[2];
rz(0.1720998) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(-0.43898359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028845499) q[0];
sx q[0];
rz(-1.5486568) q[0];
sx q[0];
rz(-1.5966709) q[0];
rz(-1.159809) q[1];
sx q[1];
rz(-1.3216959) q[1];
sx q[1];
rz(-1.6285223) q[1];
rz(0.25111689) q[2];
sx q[2];
rz(-2.579183) q[2];
sx q[2];
rz(-0.05280799) q[2];
rz(2.0598026) q[3];
sx q[3];
rz(-1.5086969) q[3];
sx q[3];
rz(2.8668565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
