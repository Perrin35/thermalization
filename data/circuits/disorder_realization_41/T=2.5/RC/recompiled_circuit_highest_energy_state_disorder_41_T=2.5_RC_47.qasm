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
rz(2.7695739) q[0];
sx q[0];
rz(-0.3516742) q[0];
sx q[0];
rz(-3.0863808) q[0];
rz(1.4456324) q[1];
sx q[1];
rz(-0.90336019) q[1];
sx q[1];
rz(-0.13394314) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3066912) q[0];
sx q[0];
rz(-1.4708232) q[0];
sx q[0];
rz(2.9213219) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8419398) q[2];
sx q[2];
rz(-0.60115325) q[2];
sx q[2];
rz(3.0944097) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81649924) q[1];
sx q[1];
rz(-1.2383922) q[1];
sx q[1];
rz(1.2219882) q[1];
rz(-pi) q[2];
rz(-1.5033967) q[3];
sx q[3];
rz(-2.615228) q[3];
sx q[3];
rz(2.0215066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.967531) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(-2.3270712) q[2];
rz(2.9461765) q[3];
sx q[3];
rz(-2.312909) q[3];
sx q[3];
rz(-1.0403847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.279351) q[0];
sx q[0];
rz(-0.26222721) q[0];
sx q[0];
rz(0.97214118) q[0];
rz(0.60130087) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(1.7960637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1336254) q[0];
sx q[0];
rz(-1.04038) q[0];
sx q[0];
rz(-2.6625457) q[0];
x q[1];
rz(-0.43429476) q[2];
sx q[2];
rz(-0.58048297) q[2];
sx q[2];
rz(-2.8934997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0474009) q[1];
sx q[1];
rz(-1.4287474) q[1];
sx q[1];
rz(-1.1588276) q[1];
x q[2];
rz(1.5645157) q[3];
sx q[3];
rz(-1.7518861) q[3];
sx q[3];
rz(-2.5400502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95416516) q[2];
sx q[2];
rz(-1.8623872) q[2];
sx q[2];
rz(2.1698451) q[2];
rz(1.0176954) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(-0.051232256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.510842) q[0];
sx q[0];
rz(-2.172281) q[0];
sx q[0];
rz(-0.72072679) q[0];
rz(3.0156056) q[1];
sx q[1];
rz(-0.69460136) q[1];
sx q[1];
rz(0.40649498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2590721) q[0];
sx q[0];
rz(-1.2693624) q[0];
sx q[0];
rz(-0.92275019) q[0];
rz(0.7438306) q[2];
sx q[2];
rz(-1.1421428) q[2];
sx q[2];
rz(2.8016318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49660044) q[1];
sx q[1];
rz(-2.3178702) q[1];
sx q[1];
rz(2.3053667) q[1];
x q[2];
rz(-0.28789374) q[3];
sx q[3];
rz(-1.9215309) q[3];
sx q[3];
rz(2.3755297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5410109) q[2];
sx q[2];
rz(-1.3449679) q[2];
sx q[2];
rz(-0.37008944) q[2];
rz(3.0626512) q[3];
sx q[3];
rz(-0.97349662) q[3];
sx q[3];
rz(-2.6551042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2094035) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(2.6655647) q[0];
rz(-2.027482) q[1];
sx q[1];
rz(-0.94534355) q[1];
sx q[1];
rz(-1.5155189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5334271) q[0];
sx q[0];
rz(-1.2100056) q[0];
sx q[0];
rz(1.7955304) q[0];
rz(2.3298181) q[2];
sx q[2];
rz(-1.1110795) q[2];
sx q[2];
rz(0.22961337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9263401) q[1];
sx q[1];
rz(-2.7801792) q[1];
sx q[1];
rz(0.24140668) q[1];
rz(2.8164316) q[3];
sx q[3];
rz(-0.87988461) q[3];
sx q[3];
rz(-2.8392217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6473306) q[2];
sx q[2];
rz(-1.6202972) q[2];
sx q[2];
rz(-2.1935513) q[2];
rz(0.4387795) q[3];
sx q[3];
rz(-0.56883562) q[3];
sx q[3];
rz(-2.2426864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(0.4278675) q[0];
rz(2.1828792) q[1];
sx q[1];
rz(-1.2370647) q[1];
sx q[1];
rz(-2.0895035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89969544) q[0];
sx q[0];
rz(-0.98693958) q[0];
sx q[0];
rz(-2.5510699) q[0];
rz(-pi) q[1];
rz(-0.97276997) q[2];
sx q[2];
rz(-2.5199119) q[2];
sx q[2];
rz(-1.6046815) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8241938) q[1];
sx q[1];
rz(-1.4772067) q[1];
sx q[1];
rz(-0.99653523) q[1];
x q[2];
rz(-2.4067299) q[3];
sx q[3];
rz(-2.4331193) q[3];
sx q[3];
rz(0.81516176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9713356) q[2];
sx q[2];
rz(-1.3054138) q[2];
sx q[2];
rz(2.7830284) q[2];
rz(-1.1809008) q[3];
sx q[3];
rz(-0.87555331) q[3];
sx q[3];
rz(-0.53926474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81826687) q[0];
sx q[0];
rz(-1.8805255) q[0];
sx q[0];
rz(2.4928424) q[0];
rz(-2.5541041) q[1];
sx q[1];
rz(-1.3222539) q[1];
sx q[1];
rz(2.9041451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96188155) q[0];
sx q[0];
rz(-1.6673606) q[0];
sx q[0];
rz(-1.6029111) q[0];
rz(-0.26837792) q[2];
sx q[2];
rz(-1.4726) q[2];
sx q[2];
rz(-2.9609916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9623649) q[1];
sx q[1];
rz(-2.3931074) q[1];
sx q[1];
rz(0.039123936) q[1];
rz(0.86187382) q[3];
sx q[3];
rz(-1.1824058) q[3];
sx q[3];
rz(-0.43530048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6610403) q[2];
sx q[2];
rz(-0.66142267) q[2];
sx q[2];
rz(0.95575571) q[2];
rz(1.7620979) q[3];
sx q[3];
rz(-0.62164128) q[3];
sx q[3];
rz(-2.9852338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5718004) q[0];
sx q[0];
rz(-3.1389696) q[0];
sx q[0];
rz(0.045259137) q[0];
rz(2.3472002) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(2.9387567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9211333) q[0];
sx q[0];
rz(-2.5510802) q[0];
sx q[0];
rz(-2.8648225) q[0];
rz(-pi) q[1];
rz(-2.9309996) q[2];
sx q[2];
rz(-0.74264975) q[2];
sx q[2];
rz(-0.79162486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0675373) q[1];
sx q[1];
rz(-0.54880868) q[1];
sx q[1];
rz(2.9850053) q[1];
rz(-1.5286115) q[3];
sx q[3];
rz(-0.72244553) q[3];
sx q[3];
rz(2.9056321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5016735) q[2];
sx q[2];
rz(-1.7379652) q[2];
sx q[2];
rz(-0.61140927) q[2];
rz(2.0407138) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(0.27975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6708577) q[0];
sx q[0];
rz(-1.1962471) q[0];
sx q[0];
rz(-1.083495) q[0];
rz(-2.1776543) q[1];
sx q[1];
rz(-2.0699392) q[1];
sx q[1];
rz(0.49577698) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9993047) q[0];
sx q[0];
rz(-1.5919627) q[0];
sx q[0];
rz(1.4975966) q[0];
rz(2.1473608) q[2];
sx q[2];
rz(-2.814817) q[2];
sx q[2];
rz(2.7852259) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3246877) q[1];
sx q[1];
rz(-2.4084615) q[1];
sx q[1];
rz(-2.2156961) q[1];
x q[2];
rz(-3.1115948) q[3];
sx q[3];
rz(-0.98326937) q[3];
sx q[3];
rz(0.8647747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44020161) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(0.68897796) q[2];
rz(1.5135328) q[3];
sx q[3];
rz(-0.83008927) q[3];
sx q[3];
rz(-1.0015063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5772407) q[0];
sx q[0];
rz(-1.5761201) q[0];
sx q[0];
rz(-2.054457) q[0];
rz(0.76308909) q[1];
sx q[1];
rz(-1.7440081) q[1];
sx q[1];
rz(-2.9147002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.034982) q[0];
sx q[0];
rz(-0.7312432) q[0];
sx q[0];
rz(1.2319706) q[0];
x q[1];
rz(0.29818887) q[2];
sx q[2];
rz(-1.2623566) q[2];
sx q[2];
rz(-2.1205663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1676281) q[1];
sx q[1];
rz(-1.7953606) q[1];
sx q[1];
rz(-1.533482) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97051981) q[3];
sx q[3];
rz(-1.9991989) q[3];
sx q[3];
rz(-2.4267933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6202966) q[2];
sx q[2];
rz(-0.93775788) q[2];
sx q[2];
rz(2.1916981) q[2];
rz(1.0319483) q[3];
sx q[3];
rz(-2.2622006) q[3];
sx q[3];
rz(-2.7040645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501605) q[0];
sx q[0];
rz(-0.10364769) q[0];
sx q[0];
rz(1.1520804) q[0];
rz(-0.092747124) q[1];
sx q[1];
rz(-2.4536965) q[1];
sx q[1];
rz(2.6377717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9418056) q[0];
sx q[0];
rz(-1.620294) q[0];
sx q[0];
rz(2.097553) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4359792) q[2];
sx q[2];
rz(-0.34760568) q[2];
sx q[2];
rz(0.068989601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42119869) q[1];
sx q[1];
rz(-2.6662146) q[1];
sx q[1];
rz(3.0488411) q[1];
rz(-pi) q[2];
rz(-1.9484048) q[3];
sx q[3];
rz(-1.8657627) q[3];
sx q[3];
rz(-2.0279391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6759701) q[2];
sx q[2];
rz(-1.114782) q[2];
sx q[2];
rz(-0.62391227) q[2];
rz(2.5850249) q[3];
sx q[3];
rz(-0.68816319) q[3];
sx q[3];
rz(0.19779675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1933761) q[0];
sx q[0];
rz(-2.2103136) q[0];
sx q[0];
rz(0.92934004) q[0];
rz(-0.14840645) q[1];
sx q[1];
rz(-1.8971309) q[1];
sx q[1];
rz(2.94577) q[1];
rz(-0.23870809) q[2];
sx q[2];
rz(-0.47521026) q[2];
sx q[2];
rz(-1.0929562) q[2];
rz(-1.2986167) q[3];
sx q[3];
rz(-1.4930354) q[3];
sx q[3];
rz(-2.2477575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
