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
rz(2.0692985) q[0];
sx q[0];
rz(-2.9579853) q[0];
sx q[0];
rz(-2.973383) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(-2.0252731) q[1];
sx q[1];
rz(-1.1180374) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71336761) q[0];
sx q[0];
rz(-1.8963739) q[0];
sx q[0];
rz(1.6859562) q[0];
rz(-pi) q[1];
rz(0.44702304) q[2];
sx q[2];
rz(-0.7444844) q[2];
sx q[2];
rz(0.15310619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5745821) q[1];
sx q[1];
rz(-1.5956465) q[1];
sx q[1];
rz(-0.76038469) q[1];
rz(-1.1088896) q[3];
sx q[3];
rz(-2.8402764) q[3];
sx q[3];
rz(-1.0404943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2844124) q[2];
sx q[2];
rz(-2.1511011) q[2];
sx q[2];
rz(-1.277479) q[2];
rz(-0.60753456) q[3];
sx q[3];
rz(-1.7888864) q[3];
sx q[3];
rz(1.7551306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4472537) q[0];
sx q[0];
rz(-0.67830938) q[0];
sx q[0];
rz(-2.4439478) q[0];
rz(-2.8088226) q[1];
sx q[1];
rz(-2.0336626) q[1];
sx q[1];
rz(-0.3322126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32658122) q[0];
sx q[0];
rz(-1.4403116) q[0];
sx q[0];
rz(1.4490118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1474481) q[2];
sx q[2];
rz(-1.3704164) q[2];
sx q[2];
rz(2.8746614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0337341) q[1];
sx q[1];
rz(-0.49328557) q[1];
sx q[1];
rz(-2.1164338) q[1];
x q[2];
rz(2.6744889) q[3];
sx q[3];
rz(-0.78416608) q[3];
sx q[3];
rz(-0.84475603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4405926) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(-2.8458703) q[2];
rz(-3.131007) q[3];
sx q[3];
rz(-0.9897832) q[3];
sx q[3];
rz(2.4081965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6538786) q[0];
sx q[0];
rz(-2.6831477) q[0];
sx q[0];
rz(-0.53471765) q[0];
rz(1.1010822) q[1];
sx q[1];
rz(-1.705876) q[1];
sx q[1];
rz(-0.91317493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3951626) q[0];
sx q[0];
rz(-2.5341153) q[0];
sx q[0];
rz(-0.98275252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0139364) q[2];
sx q[2];
rz(-1.3512313) q[2];
sx q[2];
rz(0.94925052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5194381) q[1];
sx q[1];
rz(-1.039912) q[1];
sx q[1];
rz(-0.77150788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.072973386) q[3];
sx q[3];
rz(-0.84318959) q[3];
sx q[3];
rz(-0.062372717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59213006) q[2];
sx q[2];
rz(-2.2495146) q[2];
sx q[2];
rz(-2.7590052) q[2];
rz(-1.654918) q[3];
sx q[3];
rz(-2.8463709) q[3];
sx q[3];
rz(2.2216643) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611068) q[0];
sx q[0];
rz(-1.684573) q[0];
sx q[0];
rz(-2.0953505) q[0];
rz(0.56398448) q[1];
sx q[1];
rz(-2.1395186) q[1];
sx q[1];
rz(-1.3169588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038669) q[0];
sx q[0];
rz(-0.90529862) q[0];
sx q[0];
rz(-2.4212877) q[0];
x q[1];
rz(-2.2235419) q[2];
sx q[2];
rz(-1.8260644) q[2];
sx q[2];
rz(-1.7328615) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0610132) q[1];
sx q[1];
rz(-1.6135585) q[1];
sx q[1];
rz(0.31400915) q[1];
x q[2];
rz(-2.3583007) q[3];
sx q[3];
rz(-1.2630594) q[3];
sx q[3];
rz(-2.7819992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4468533) q[2];
sx q[2];
rz(-1.4444192) q[2];
sx q[2];
rz(-2.1430338) q[2];
rz(-2.5616732) q[3];
sx q[3];
rz(-1.6603575) q[3];
sx q[3];
rz(1.8370321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21438433) q[0];
sx q[0];
rz(-1.2283607) q[0];
sx q[0];
rz(0.71980113) q[0];
rz(-2.39847) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(-0.057083759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.538495) q[0];
sx q[0];
rz(-2.7481226) q[0];
sx q[0];
rz(-0.35719494) q[0];
rz(-pi) q[1];
rz(1.916831) q[2];
sx q[2];
rz(-0.84362307) q[2];
sx q[2];
rz(-2.7131701) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.067043153) q[1];
sx q[1];
rz(-0.32575575) q[1];
sx q[1];
rz(-2.6200952) q[1];
rz(-0.96995285) q[3];
sx q[3];
rz(-2.4527453) q[3];
sx q[3];
rz(3.0843902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4382978) q[2];
sx q[2];
rz(-1.1763828) q[2];
sx q[2];
rz(1.2698184) q[2];
rz(-2.5979009) q[3];
sx q[3];
rz(-0.68259493) q[3];
sx q[3];
rz(-1.3683176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.066034) q[0];
sx q[0];
rz(-2.7664001) q[0];
sx q[0];
rz(-2.913108) q[0];
rz(2.1560419) q[1];
sx q[1];
rz(-1.5551714) q[1];
sx q[1];
rz(1.2062581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0096480308) q[0];
sx q[0];
rz(-1.7382657) q[0];
sx q[0];
rz(2.8020049) q[0];
rz(-pi) q[1];
rz(1.2689077) q[2];
sx q[2];
rz(-1.3573555) q[2];
sx q[2];
rz(0.53415307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6724193) q[1];
sx q[1];
rz(-2.1757728) q[1];
sx q[1];
rz(-0.14766001) q[1];
x q[2];
rz(2.5980232) q[3];
sx q[3];
rz(-1.6619716) q[3];
sx q[3];
rz(1.4379355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0050521684) q[2];
sx q[2];
rz(-1.2753992) q[2];
sx q[2];
rz(1.2460111) q[2];
rz(0.15240845) q[3];
sx q[3];
rz(-0.14901769) q[3];
sx q[3];
rz(-2.3685031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8361255) q[0];
sx q[0];
rz(-1.2141328) q[0];
sx q[0];
rz(2.8866696) q[0];
rz(0.98006025) q[1];
sx q[1];
rz(-1.7995116) q[1];
sx q[1];
rz(2.9068388) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31862751) q[0];
sx q[0];
rz(-1.8150094) q[0];
sx q[0];
rz(2.1313214) q[0];
rz(-pi) q[1];
rz(-2.1990216) q[2];
sx q[2];
rz(-2.8543575) q[2];
sx q[2];
rz(-2.7554905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67803176) q[1];
sx q[1];
rz(-1.982219) q[1];
sx q[1];
rz(-1.8586584) q[1];
rz(2.2473621) q[3];
sx q[3];
rz(-1.4787592) q[3];
sx q[3];
rz(0.91023723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0945101) q[2];
sx q[2];
rz(-2.7213056) q[2];
sx q[2];
rz(-1.1959929) q[2];
rz(2.1434873) q[3];
sx q[3];
rz(-1.866303) q[3];
sx q[3];
rz(-0.78817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.81414476) q[0];
sx q[0];
rz(-0.87004167) q[0];
sx q[0];
rz(0.63234627) q[0];
rz(1.5517722) q[1];
sx q[1];
rz(-1.7951782) q[1];
sx q[1];
rz(-1.8428892) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.08399) q[0];
sx q[0];
rz(-0.51483908) q[0];
sx q[0];
rz(-2.0979416) q[0];
rz(-pi) q[1];
rz(-0.66565158) q[2];
sx q[2];
rz(-1.5857664) q[2];
sx q[2];
rz(-2.7240208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92346159) q[1];
sx q[1];
rz(-2.0773481) q[1];
sx q[1];
rz(1.6183557) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9901806) q[3];
sx q[3];
rz(-1.5431343) q[3];
sx q[3];
rz(-2.9392978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.87045264) q[2];
sx q[2];
rz(-0.91338921) q[2];
sx q[2];
rz(0.75922981) q[2];
rz(-0.40902725) q[3];
sx q[3];
rz(-1.7274011) q[3];
sx q[3];
rz(0.11212382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2893534) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(0.36561832) q[0];
rz(3.1386555) q[1];
sx q[1];
rz(-1.5298839) q[1];
sx q[1];
rz(-2.073435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0442783) q[0];
sx q[0];
rz(-1.8476163) q[0];
sx q[0];
rz(-2.5702049) q[0];
rz(-0.85313646) q[2];
sx q[2];
rz(-0.87867773) q[2];
sx q[2];
rz(-1.4286526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4959855) q[1];
sx q[1];
rz(-1.9580972) q[1];
sx q[1];
rz(0.041195804) q[1];
rz(0.2016507) q[3];
sx q[3];
rz(-1.5976054) q[3];
sx q[3];
rz(-2.9134515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7630161) q[2];
sx q[2];
rz(-2.2282659) q[2];
sx q[2];
rz(-0.98706377) q[2];
rz(-0.8118363) q[3];
sx q[3];
rz(-2.2631009) q[3];
sx q[3];
rz(-1.8973551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768452) q[0];
sx q[0];
rz(-0.61926121) q[0];
sx q[0];
rz(-1.4641807) q[0];
rz(-0.96425104) q[1];
sx q[1];
rz(-1.3136656) q[1];
sx q[1];
rz(-1.3826694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.156517) q[0];
sx q[0];
rz(-1.1958953) q[0];
sx q[0];
rz(0.84701726) q[0];
x q[1];
rz(-1.8455454) q[2];
sx q[2];
rz(-1.2128069) q[2];
sx q[2];
rz(1.4520979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5359491) q[1];
sx q[1];
rz(-1.9943731) q[1];
sx q[1];
rz(0.98137318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85282214) q[3];
sx q[3];
rz(-2.3896165) q[3];
sx q[3];
rz(-0.2096006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1796403) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(1.4947653) q[2];
rz(1.1163813) q[3];
sx q[3];
rz(-0.79018441) q[3];
sx q[3];
rz(-2.5377929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8981001) q[0];
sx q[0];
rz(-1.4673163) q[0];
sx q[0];
rz(-1.2661288) q[0];
rz(-2.948214) q[1];
sx q[1];
rz(-2.0722176) q[1];
sx q[1];
rz(-1.8373012) q[1];
rz(-1.2639574) q[2];
sx q[2];
rz(-2.2813792) q[2];
sx q[2];
rz(1.1165237) q[2];
rz(2.0054655) q[3];
sx q[3];
rz(-2.2129314) q[3];
sx q[3];
rz(1.5066838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
