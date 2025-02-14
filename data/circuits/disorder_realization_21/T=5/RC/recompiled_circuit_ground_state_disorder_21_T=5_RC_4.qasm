OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.8747099) q[0];
sx q[0];
rz(4.2031718) q[0];
sx q[0];
rz(9.940552) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(-0.20144784) q[1];
sx q[1];
rz(2.1405061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2767746) q[0];
sx q[0];
rz(-2.1208115) q[0];
sx q[0];
rz(0.986542) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1273448) q[2];
sx q[2];
rz(-2.7269335) q[2];
sx q[2];
rz(-2.8480107) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42131073) q[1];
sx q[1];
rz(-1.1291468) q[1];
sx q[1];
rz(2.4288105) q[1];
rz(-pi) q[2];
rz(0.14747133) q[3];
sx q[3];
rz(-1.4502954) q[3];
sx q[3];
rz(1.7013753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4150245) q[2];
sx q[2];
rz(-0.51076204) q[2];
sx q[2];
rz(-1.6815574) q[2];
rz(-0.36871746) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(1.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6129795) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(3.0533277) q[0];
rz(2.2093692) q[1];
sx q[1];
rz(-2.014522) q[1];
sx q[1];
rz(2.6384735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4012484) q[0];
sx q[0];
rz(-0.49763864) q[0];
sx q[0];
rz(-0.21792441) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2406949) q[2];
sx q[2];
rz(-1.42282) q[2];
sx q[2];
rz(-1.0690546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14993653) q[1];
sx q[1];
rz(-2.7767065) q[1];
sx q[1];
rz(0.30662341) q[1];
x q[2];
rz(-2.6121268) q[3];
sx q[3];
rz(-0.90116167) q[3];
sx q[3];
rz(2.1807699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3162389) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(0.45315722) q[2];
rz(0.00099269021) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(-0.7411952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2082625) q[0];
sx q[0];
rz(-0.18284155) q[0];
sx q[0];
rz(2.6728447) q[0];
rz(0.58724976) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(-2.8900878) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55401459) q[0];
sx q[0];
rz(-1.9791326) q[0];
sx q[0];
rz(-2.6760882) q[0];
rz(-pi) q[1];
rz(-1.3590888) q[2];
sx q[2];
rz(-1.9518753) q[2];
sx q[2];
rz(0.43109387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1914563) q[1];
sx q[1];
rz(-0.13683441) q[1];
sx q[1];
rz(2.0124643) q[1];
rz(-2.4819314) q[3];
sx q[3];
rz(-0.53852496) q[3];
sx q[3];
rz(-0.3769359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.544439) q[2];
sx q[2];
rz(-2.4033098) q[2];
sx q[2];
rz(0.047133751) q[2];
rz(2.2733287) q[3];
sx q[3];
rz(-1.2406113) q[3];
sx q[3];
rz(-2.6422083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.353001) q[0];
sx q[0];
rz(-2.7356739) q[0];
sx q[0];
rz(1.3141919) q[0];
rz(0.81002533) q[1];
sx q[1];
rz(-1.516073) q[1];
sx q[1];
rz(-0.31585082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0994713) q[0];
sx q[0];
rz(-2.1232623) q[0];
sx q[0];
rz(-0.61037678) q[0];
rz(-0.31762697) q[2];
sx q[2];
rz(-1.7343354) q[2];
sx q[2];
rz(-2.7224772) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12055732) q[1];
sx q[1];
rz(-1.2791954) q[1];
sx q[1];
rz(-0.75052829) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60499936) q[3];
sx q[3];
rz(-1.446786) q[3];
sx q[3];
rz(2.0600788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6541859) q[2];
sx q[2];
rz(-2.5155289) q[2];
sx q[2];
rz(2.249304) q[2];
rz(-1.4687294) q[3];
sx q[3];
rz(-1.059831) q[3];
sx q[3];
rz(-1.0631801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6650218) q[0];
sx q[0];
rz(-2.5044818) q[0];
sx q[0];
rz(-0.88660216) q[0];
rz(2.2992112) q[1];
sx q[1];
rz(-1.8319172) q[1];
sx q[1];
rz(2.0319891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32217596) q[0];
sx q[0];
rz(-2.3482394) q[0];
sx q[0];
rz(-0.94950907) q[0];
rz(2.9363052) q[2];
sx q[2];
rz(-0.73371202) q[2];
sx q[2];
rz(-3.044341) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8362404) q[1];
sx q[1];
rz(-1.0100366) q[1];
sx q[1];
rz(-0.95587041) q[1];
x q[2];
rz(-0.97115406) q[3];
sx q[3];
rz(-0.97929882) q[3];
sx q[3];
rz(-2.9324071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2716918) q[2];
sx q[2];
rz(-2.7503408) q[2];
sx q[2];
rz(-0.58763495) q[2];
rz(2.9247126) q[3];
sx q[3];
rz(-1.8606595) q[3];
sx q[3];
rz(1.3542401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0103535) q[0];
sx q[0];
rz(-2.6709747) q[0];
sx q[0];
rz(1.7577897) q[0];
rz(-2.9673987) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(-1.9536288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9458141) q[0];
sx q[0];
rz(-1.4885689) q[0];
sx q[0];
rz(1.29127) q[0];
rz(-pi) q[1];
rz(-0.50601658) q[2];
sx q[2];
rz(-1.4219473) q[2];
sx q[2];
rz(-2.2300492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7158051) q[1];
sx q[1];
rz(-1.5400229) q[1];
sx q[1];
rz(-1.6795067) q[1];
rz(-pi) q[2];
rz(0.31962304) q[3];
sx q[3];
rz(-0.79053662) q[3];
sx q[3];
rz(-2.9755693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4478628) q[2];
sx q[2];
rz(-1.2749981) q[2];
sx q[2];
rz(1.0542144) q[2];
rz(2.5888455) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(-2.0235846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066910557) q[0];
sx q[0];
rz(-0.56147611) q[0];
sx q[0];
rz(0.018360227) q[0];
rz(1.8439937) q[1];
sx q[1];
rz(-1.8074139) q[1];
sx q[1];
rz(1.1150572) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149688) q[0];
sx q[0];
rz(-1.4415662) q[0];
sx q[0];
rz(-1.5530759) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8283496) q[2];
sx q[2];
rz(-2.0038824) q[2];
sx q[2];
rz(1.05561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1464941) q[1];
sx q[1];
rz(-2.0420688) q[1];
sx q[1];
rz(-2.5033094) q[1];
rz(2.2998718) q[3];
sx q[3];
rz(-0.63243619) q[3];
sx q[3];
rz(3.1068592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4516419) q[2];
sx q[2];
rz(-2.3929598) q[2];
sx q[2];
rz(2.3484223) q[2];
rz(-2.4801109) q[3];
sx q[3];
rz(-1.3756779) q[3];
sx q[3];
rz(2.49559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.20872214) q[0];
sx q[0];
rz(-2.0168309) q[0];
sx q[0];
rz(0.17054184) q[0];
rz(2.1216682) q[1];
sx q[1];
rz(-0.41671697) q[1];
sx q[1];
rz(2.4822617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19457283) q[0];
sx q[0];
rz(-1.6259369) q[0];
sx q[0];
rz(1.5558262) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3691828) q[2];
sx q[2];
rz(-0.67467022) q[2];
sx q[2];
rz(1.9325369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7119321) q[1];
sx q[1];
rz(-1.871456) q[1];
sx q[1];
rz(2.4983404) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.002568) q[3];
sx q[3];
rz(-2.2428277) q[3];
sx q[3];
rz(0.51154414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3063804) q[2];
sx q[2];
rz(-0.39432085) q[2];
sx q[2];
rz(-2.1739056) q[2];
rz(0.59099284) q[3];
sx q[3];
rz(-1.7715745) q[3];
sx q[3];
rz(1.437457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0377334) q[0];
sx q[0];
rz(-0.93872207) q[0];
sx q[0];
rz(2.9528604) q[0];
rz(-2.0813148) q[1];
sx q[1];
rz(-2.4791368) q[1];
sx q[1];
rz(2.2417384) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2162595) q[0];
sx q[0];
rz(-2.1104314) q[0];
sx q[0];
rz(-0.58109053) q[0];
rz(-pi) q[1];
rz(3.138928) q[2];
sx q[2];
rz(-0.84328534) q[2];
sx q[2];
rz(-1.8553986) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.65727216) q[1];
sx q[1];
rz(-1.5725699) q[1];
sx q[1];
rz(0.65048762) q[1];
x q[2];
rz(-1.3542309) q[3];
sx q[3];
rz(-0.73288554) q[3];
sx q[3];
rz(0.93001825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.73682827) q[2];
sx q[2];
rz(-1.2238203) q[2];
sx q[2];
rz(-1.9723816) q[2];
rz(0.87559187) q[3];
sx q[3];
rz(-0.67684042) q[3];
sx q[3];
rz(0.08918795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.5065696) q[0];
sx q[0];
rz(-1.4590141) q[0];
sx q[0];
rz(2.7154679) q[0];
rz(-1.9020853) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(0.32136163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581767) q[0];
sx q[0];
rz(-2.1431648) q[0];
sx q[0];
rz(3.0757927) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05409492) q[2];
sx q[2];
rz(-1.128643) q[2];
sx q[2];
rz(1.7345605) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4453955) q[1];
sx q[1];
rz(-1.9085669) q[1];
sx q[1];
rz(-0.15134751) q[1];
rz(0.71236082) q[3];
sx q[3];
rz(-1.9689631) q[3];
sx q[3];
rz(-1.1172197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1539803) q[2];
sx q[2];
rz(-0.93903792) q[2];
sx q[2];
rz(-0.067848094) q[2];
rz(-2.1765354) q[3];
sx q[3];
rz(-2.8128251) q[3];
sx q[3];
rz(1.4062101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52473749) q[0];
sx q[0];
rz(-2.6829834) q[0];
sx q[0];
rz(0.99880698) q[0];
rz(0.85159341) q[1];
sx q[1];
rz(-1.0706182) q[1];
sx q[1];
rz(0.65382438) q[1];
rz(1.6338255) q[2];
sx q[2];
rz(-1.0152643) q[2];
sx q[2];
rz(2.0282657) q[2];
rz(0.028164955) q[3];
sx q[3];
rz(-0.71816254) q[3];
sx q[3];
rz(-0.95550334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
