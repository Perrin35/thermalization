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
rz(-2.3772216) q[0];
sx q[0];
rz(-1.8005014) q[0];
sx q[0];
rz(-2.2361225) q[0];
rz(4.6484923) q[1];
sx q[1];
rz(5.3223106) q[1];
sx q[1];
rz(7.9237908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5643217) q[0];
sx q[0];
rz(-2.4080896) q[0];
sx q[0];
rz(1.8689687) q[0];
rz(2.9653984) q[2];
sx q[2];
rz(-2.004188) q[2];
sx q[2];
rz(-2.7769763) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.562362) q[1];
sx q[1];
rz(-3.0931614) q[1];
sx q[1];
rz(-2.6391451) q[1];
x q[2];
rz(-2.1330058) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(0.62801559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2291439) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(2.8196715) q[2];
rz(-2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(1.9979075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(2.6296997) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(-2.0194676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7417541) q[0];
sx q[0];
rz(-1.6723833) q[0];
sx q[0];
rz(2.8288657) q[0];
rz(1.4830515) q[2];
sx q[2];
rz(-0.43102396) q[2];
sx q[2];
rz(-1.6602248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1074688) q[1];
sx q[1];
rz(-1.9975047) q[1];
sx q[1];
rz(0.36860768) q[1];
rz(-pi) q[2];
rz(0.53776017) q[3];
sx q[3];
rz(-2.3590292) q[3];
sx q[3];
rz(1.7478706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7684218) q[2];
sx q[2];
rz(-1.827652) q[2];
sx q[2];
rz(2.6785417) q[2];
rz(2.566346) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(-3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333882) q[0];
sx q[0];
rz(-2.2938804) q[0];
sx q[0];
rz(0.43011618) q[0];
rz(2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(-2.5045085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67733902) q[0];
sx q[0];
rz(-1.5298944) q[0];
sx q[0];
rz(-0.014057191) q[0];
rz(-pi) q[1];
rz(-0.78539679) q[2];
sx q[2];
rz(-1.5522508) q[2];
sx q[2];
rz(-2.3634499) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7635258) q[1];
sx q[1];
rz(-0.77516205) q[1];
sx q[1];
rz(0.35825348) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5175321) q[3];
sx q[3];
rz(-2.5284323) q[3];
sx q[3];
rz(2.000119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35885262) q[2];
sx q[2];
rz(-2.095486) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(1.8105761) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(-0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.26938874) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(1.4440906) q[0];
rz(0.048642453) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(0.28894249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5015857) q[0];
sx q[0];
rz(-1.3292392) q[0];
sx q[0];
rz(-1.9511186) q[0];
x q[1];
rz(1.4202576) q[2];
sx q[2];
rz(-1.5988013) q[2];
sx q[2];
rz(-1.3893407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1565981) q[1];
sx q[1];
rz(-0.20594507) q[1];
sx q[1];
rz(2.1316281) q[1];
x q[2];
rz(-3.0722202) q[3];
sx q[3];
rz(-2.4751304) q[3];
sx q[3];
rz(1.762076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(-2.5034261) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3087092) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(-2.047245) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(-1.099115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35661429) q[0];
sx q[0];
rz(-1.4334502) q[0];
sx q[0];
rz(2.7524968) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34277041) q[2];
sx q[2];
rz(-1.2452599) q[2];
sx q[2];
rz(2.8936762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.079541072) q[1];
sx q[1];
rz(-1.4126443) q[1];
sx q[1];
rz(1.016721) q[1];
rz(1.7775675) q[3];
sx q[3];
rz(-2.4768157) q[3];
sx q[3];
rz(0.95486508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2431474) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-2.1649427) q[2];
rz(-0.53168932) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(2.4272052) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.8699159) q[0];
rz(-0.42287982) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(-1.5333102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9031808) q[0];
sx q[0];
rz(-0.067771284) q[0];
sx q[0];
rz(-1.3890024) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1321636) q[2];
sx q[2];
rz(-1.4727739) q[2];
sx q[2];
rz(0.67154166) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.803992) q[1];
sx q[1];
rz(-1.873766) q[1];
sx q[1];
rz(1.5125809) q[1];
rz(-pi) q[2];
rz(2.9098784) q[3];
sx q[3];
rz(-2.3149256) q[3];
sx q[3];
rz(0.8207013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61234683) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(-0.38522729) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(-2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92152921) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(-2.7440942) q[0];
rz(-2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(0.0040815512) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3101871) q[0];
sx q[0];
rz(-2.3334529) q[0];
sx q[0];
rz(-0.17669295) q[0];
rz(-pi) q[1];
rz(-1.0881617) q[2];
sx q[2];
rz(-2.3950999) q[2];
sx q[2];
rz(-1.6255524) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4392802) q[1];
sx q[1];
rz(-1.1902307) q[1];
sx q[1];
rz(0.61813942) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8936196) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(0.21305591) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-3.1000948) q[3];
sx q[3];
rz(-1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(1.9785471) q[0];
rz(2.842438) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(-2.1655653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928182) q[0];
sx q[0];
rz(-2.8716757) q[0];
sx q[0];
rz(-0.42645539) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9891122) q[2];
sx q[2];
rz(-1.3950384) q[2];
sx q[2];
rz(2.3703142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8638226) q[1];
sx q[1];
rz(-2.7105717) q[1];
sx q[1];
rz(-2.1716539) q[1];
rz(-1.3221413) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(1.6909042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8396847) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(-1.3236375) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(1.6424204) q[0];
rz(-2.6014853) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(1.3444208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7758862) q[0];
sx q[0];
rz(-1.8593863) q[0];
sx q[0];
rz(-1.9802753) q[0];
rz(-0.70602472) q[2];
sx q[2];
rz(-2.386552) q[2];
sx q[2];
rz(1.4817099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.498913) q[1];
sx q[1];
rz(-1.0084001) q[1];
sx q[1];
rz(-2.845473) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62219859) q[3];
sx q[3];
rz(-2.1902764) q[3];
sx q[3];
rz(3.0752237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6092047) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(1.4538291) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(-0.55117575) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(-0.59250441) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3640219) q[0];
sx q[0];
rz(-0.3705655) q[0];
sx q[0];
rz(0.70623605) q[0];
rz(-0.76982381) q[2];
sx q[2];
rz(-0.69902674) q[2];
sx q[2];
rz(1.5848643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4514102) q[1];
sx q[1];
rz(-1.8333149) q[1];
sx q[1];
rz(-2.7464944) q[1];
rz(-1.5207401) q[3];
sx q[3];
rz(-1.0721237) q[3];
sx q[3];
rz(3.0704481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42980117) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(-0.05833021) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(-0.0037732865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534828) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.7755605) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(0.57353061) q[2];
sx q[2];
rz(-1.7885359) q[2];
sx q[2];
rz(0.89319695) q[2];
rz(-2.1128863) q[3];
sx q[3];
rz(-1.1114612) q[3];
sx q[3];
rz(-2.3041861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
