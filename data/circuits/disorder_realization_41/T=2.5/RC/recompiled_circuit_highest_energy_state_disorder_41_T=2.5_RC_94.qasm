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
rz(3.4932669) q[0];
sx q[0];
rz(9.3695661) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(3.0076495) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3165749) q[0];
sx q[0];
rz(-2.900032) q[0];
sx q[0];
rz(2.7112083) q[0];
rz(-pi) q[1];
rz(0.29965286) q[2];
sx q[2];
rz(-0.60115325) q[2];
sx q[2];
rz(-3.0944097) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81649924) q[1];
sx q[1];
rz(-1.9032005) q[1];
sx q[1];
rz(1.9196045) q[1];
x q[2];
rz(1.6381959) q[3];
sx q[3];
rz(-2.615228) q[3];
sx q[3];
rz(-1.1200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17406164) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(0.81452149) q[2];
rz(-0.19541611) q[3];
sx q[3];
rz(-0.82868367) q[3];
sx q[3];
rz(-2.101208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.279351) q[0];
sx q[0];
rz(-0.26222721) q[0];
sx q[0];
rz(2.1694515) q[0];
rz(2.5402918) q[1];
sx q[1];
rz(-0.96527946) q[1];
sx q[1];
rz(1.7960637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0079673) q[0];
sx q[0];
rz(-1.04038) q[0];
sx q[0];
rz(-0.47904695) q[0];
rz(-pi) q[1];
rz(1.8400606) q[2];
sx q[2];
rz(-1.0500638) q[2];
sx q[2];
rz(0.7545287) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.094191782) q[1];
sx q[1];
rz(-1.4287474) q[1];
sx q[1];
rz(-1.9827651) q[1];
rz(0.18109326) q[3];
sx q[3];
rz(-1.5769742) q[3];
sx q[3];
rz(2.1712077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95416516) q[2];
sx q[2];
rz(-1.2792055) q[2];
sx q[2];
rz(2.1698451) q[2];
rz(1.0176954) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(3.0903604) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63075066) q[0];
sx q[0];
rz(-2.172281) q[0];
sx q[0];
rz(2.4208659) q[0];
rz(-3.0156056) q[1];
sx q[1];
rz(-0.69460136) q[1];
sx q[1];
rz(2.7350977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.061977) q[0];
sx q[0];
rz(-2.4361389) q[0];
sx q[0];
rz(2.0464226) q[0];
rz(-pi) q[1];
rz(2.1265246) q[2];
sx q[2];
rz(-2.2342618) q[2];
sx q[2];
rz(-2.2760454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5169591) q[1];
sx q[1];
rz(-2.0849094) q[1];
sx q[1];
rz(0.89526432) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2062827) q[3];
sx q[3];
rz(-1.3008833) q[3];
sx q[3];
rz(0.70333896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6005818) q[2];
sx q[2];
rz(-1.3449679) q[2];
sx q[2];
rz(-0.37008944) q[2];
rz(-0.078941405) q[3];
sx q[3];
rz(-2.168096) q[3];
sx q[3];
rz(-0.48648849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2094035) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(-2.6655647) q[0];
rz(-1.1141106) q[1];
sx q[1];
rz(-2.1962491) q[1];
sx q[1];
rz(-1.5155189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0335122) q[0];
sx q[0];
rz(-0.42244222) q[0];
sx q[0];
rz(0.53349924) q[0];
x q[1];
rz(2.1944207) q[2];
sx q[2];
rz(-2.2786254) q[2];
sx q[2];
rz(0.90367452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5596931) q[1];
sx q[1];
rz(-1.6554313) q[1];
sx q[1];
rz(-2.7898048) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8533303) q[3];
sx q[3];
rz(-1.8195517) q[3];
sx q[3];
rz(-2.0847818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49426207) q[2];
sx q[2];
rz(-1.6202972) q[2];
sx q[2];
rz(2.1935513) q[2];
rz(2.7028132) q[3];
sx q[3];
rz(-2.572757) q[3];
sx q[3];
rz(0.89890629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(0.4278675) q[0];
rz(0.95871344) q[1];
sx q[1];
rz(-1.2370647) q[1];
sx q[1];
rz(2.0895035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89969544) q[0];
sx q[0];
rz(-2.1546531) q[0];
sx q[0];
rz(-2.5510699) q[0];
x q[1];
rz(0.97276997) q[2];
sx q[2];
rz(-0.62168078) q[2];
sx q[2];
rz(-1.6046815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8278049) q[1];
sx q[1];
rz(-0.99936411) q[1];
sx q[1];
rz(3.0302583) q[1];
x q[2];
rz(1.0493199) q[3];
sx q[3];
rz(-1.0670033) q[3];
sx q[3];
rz(-1.6871444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.170257) q[2];
sx q[2];
rz(-1.8361788) q[2];
sx q[2];
rz(2.7830284) q[2];
rz(-1.9606918) q[3];
sx q[3];
rz(-0.87555331) q[3];
sx q[3];
rz(-2.6023279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233258) q[0];
sx q[0];
rz(-1.2610672) q[0];
sx q[0];
rz(0.64875025) q[0];
rz(-2.5541041) q[1];
sx q[1];
rz(-1.8193388) q[1];
sx q[1];
rz(0.23744753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61201216) q[0];
sx q[0];
rz(-1.5388312) q[0];
sx q[0];
rz(3.0449788) q[0];
x q[1];
rz(0.35570972) q[2];
sx q[2];
rz(-2.8562219) q[2];
sx q[2];
rz(1.4089546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7213504) q[1];
sx q[1];
rz(-1.5974177) q[1];
sx q[1];
rz(0.74810352) q[1];
rz(-pi) q[2];
rz(2.2797188) q[3];
sx q[3];
rz(-1.1824058) q[3];
sx q[3];
rz(-2.7062922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6610403) q[2];
sx q[2];
rz(-2.48017) q[2];
sx q[2];
rz(-2.1858369) q[2];
rz(1.3794948) q[3];
sx q[3];
rz(-0.62164128) q[3];
sx q[3];
rz(2.9852338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5718004) q[0];
sx q[0];
rz(-0.0026230165) q[0];
sx q[0];
rz(-3.0963335) q[0];
rz(-2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(-2.9387567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2506382) q[0];
sx q[0];
rz(-2.136011) q[0];
sx q[0];
rz(-1.3896451) q[0];
rz(-pi) q[1];
rz(0.21059307) q[2];
sx q[2];
rz(-2.3989429) q[2];
sx q[2];
rz(-2.3499678) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5109628) q[1];
sx q[1];
rz(-1.65224) q[1];
sx q[1];
rz(-2.5982473) q[1];
x q[2];
rz(-1.6129812) q[3];
sx q[3];
rz(-0.72244553) q[3];
sx q[3];
rz(-2.9056321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5016735) q[2];
sx q[2];
rz(-1.4036274) q[2];
sx q[2];
rz(2.5301834) q[2];
rz(-2.0407138) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(2.8618405) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9993047) q[0];
sx q[0];
rz(-1.54963) q[0];
sx q[0];
rz(-1.643996) q[0];
x q[1];
rz(2.1473608) q[2];
sx q[2];
rz(-0.3267757) q[2];
sx q[2];
rz(0.35636679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.878022) q[1];
sx q[1];
rz(-1.1568034) q[1];
sx q[1];
rz(2.1946226) q[1];
rz(-pi) q[2];
rz(3.1115948) q[3];
sx q[3];
rz(-0.98326937) q[3];
sx q[3];
rz(-0.8647747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44020161) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(-0.68897796) q[2];
rz(-1.5135328) q[3];
sx q[3];
rz(-2.3115034) q[3];
sx q[3];
rz(2.1400863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643519) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(-2.054457) q[0];
rz(-0.76308909) q[1];
sx q[1];
rz(-1.7440081) q[1];
sx q[1];
rz(-0.22689247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339461) q[0];
sx q[0];
rz(-1.7946222) q[0];
sx q[0];
rz(-2.2730519) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8434038) q[2];
sx q[2];
rz(-1.2623566) q[2];
sx q[2];
rz(2.1205663) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5530738) q[1];
sx q[1];
rz(-1.5344193) q[1];
sx q[1];
rz(0.22471551) q[1];
x q[2];
rz(-2.6361106) q[3];
sx q[3];
rz(-2.1104321) q[3];
sx q[3];
rz(2.0085378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6202966) q[2];
sx q[2];
rz(-2.2038348) q[2];
sx q[2];
rz(2.1916981) q[2];
rz(-2.1096443) q[3];
sx q[3];
rz(-0.87939206) q[3];
sx q[3];
rz(-0.43752813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501605) q[0];
sx q[0];
rz(-3.037945) q[0];
sx q[0];
rz(1.9895122) q[0];
rz(-3.0488455) q[1];
sx q[1];
rz(-2.4536965) q[1];
sx q[1];
rz(0.50382096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71386218) q[0];
sx q[0];
rz(-0.52885884) q[0];
sx q[0];
rz(1.4725757) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26911084) q[2];
sx q[2];
rz(-1.3480486) q[2];
sx q[2];
rz(2.3151223) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31695932) q[1];
sx q[1];
rz(-1.0976296) q[1];
sx q[1];
rz(-1.5231569) q[1];
rz(-2.2600365) q[3];
sx q[3];
rz(-0.47475749) q[3];
sx q[3];
rz(2.9662544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6759701) q[2];
sx q[2];
rz(-2.0268107) q[2];
sx q[2];
rz(2.5176804) q[2];
rz(-2.5850249) q[3];
sx q[3];
rz(-0.68816319) q[3];
sx q[3];
rz(2.9437959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1933761) q[0];
sx q[0];
rz(-2.2103136) q[0];
sx q[0];
rz(0.92934004) q[0];
rz(0.14840645) q[1];
sx q[1];
rz(-1.2444617) q[1];
sx q[1];
rz(-0.1958227) q[1];
rz(1.6918626) q[2];
sx q[2];
rz(-1.1101223) q[2];
sx q[2];
rz(2.3157673) q[2];
rz(-0.080720223) q[3];
sx q[3];
rz(-1.29946) q[3];
sx q[3];
rz(2.4429532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
