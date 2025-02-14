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
rz(-0.93936062) q[0];
sx q[0];
rz(4.2606104) q[0];
sx q[0];
rz(9.4584447) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(-1.6330947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38070883) q[0];
sx q[0];
rz(-0.11822001) q[0];
sx q[0];
rz(1.8272754) q[0];
rz(1.8649419) q[2];
sx q[2];
rz(-1.4987117) q[2];
sx q[2];
rz(-0.26157899) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.6623792) q[1];
sx q[1];
rz(-2.2594927) q[1];
sx q[1];
rz(-1.7583048) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63920984) q[3];
sx q[3];
rz(-0.91851252) q[3];
sx q[3];
rz(-2.2584884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-0.89189947) q[2];
sx q[2];
rz(-0.83667052) q[2];
rz(-0.80638805) q[3];
sx q[3];
rz(-2.216335) q[3];
sx q[3];
rz(-2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760913) q[0];
sx q[0];
rz(-2.0747023) q[0];
sx q[0];
rz(-1.9245032) q[0];
rz(-2.941046) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-2.5977871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4539219) q[0];
sx q[0];
rz(-0.39484307) q[0];
sx q[0];
rz(2.5107491) q[0];
rz(-pi) q[1];
rz(0.5949623) q[2];
sx q[2];
rz(-1.5595072) q[2];
sx q[2];
rz(-1.4608698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3873991) q[1];
sx q[1];
rz(-2.6541944) q[1];
sx q[1];
rz(2.2226366) q[1];
rz(-pi) q[2];
rz(-0.035057706) q[3];
sx q[3];
rz(-0.89588273) q[3];
sx q[3];
rz(1.9652776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66889888) q[2];
sx q[2];
rz(-0.74287477) q[2];
sx q[2];
rz(1.5354068) q[2];
rz(-1.499048) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(1.2523874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1010308) q[0];
sx q[0];
rz(-0.56483785) q[0];
sx q[0];
rz(3.132013) q[0];
rz(-2.1122232) q[1];
sx q[1];
rz(-1.5681489) q[1];
sx q[1];
rz(-1.8052489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9294465) q[0];
sx q[0];
rz(-0.98185724) q[0];
sx q[0];
rz(-2.1959744) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0687683) q[2];
sx q[2];
rz(-1.8548551) q[2];
sx q[2];
rz(1.1570003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.134446) q[1];
sx q[1];
rz(-1.916051) q[1];
sx q[1];
rz(2.9551639) q[1];
rz(-pi) q[2];
rz(-0.16943259) q[3];
sx q[3];
rz(-1.1798254) q[3];
sx q[3];
rz(0.54216551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0336527) q[2];
sx q[2];
rz(-2.630271) q[2];
sx q[2];
rz(-1.8281724) q[2];
rz(-3.1016453) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(2.0136755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0082598) q[0];
sx q[0];
rz(-1.2902211) q[0];
sx q[0];
rz(-0.3966575) q[0];
rz(-0.89504129) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-2.3779714) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067487262) q[0];
sx q[0];
rz(-0.25831902) q[0];
sx q[0];
rz(-1.2506358) q[0];
x q[1];
rz(-0.12053804) q[2];
sx q[2];
rz(-1.5539031) q[2];
sx q[2];
rz(2.8617045) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29429193) q[1];
sx q[1];
rz(-2.9025893) q[1];
sx q[1];
rz(-1.6836749) q[1];
x q[2];
rz(2.850789) q[3];
sx q[3];
rz(-1.5352846) q[3];
sx q[3];
rz(-0.5841271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7274196) q[2];
sx q[2];
rz(-1.6215934) q[2];
sx q[2];
rz(1.1302036) q[2];
rz(-1.1388418) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(-1.9738919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8324757) q[0];
sx q[0];
rz(-0.55823767) q[0];
sx q[0];
rz(-2.6600237) q[0];
rz(0.32749495) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(0.96424261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010134546) q[0];
sx q[0];
rz(-1.1842964) q[0];
sx q[0];
rz(0.23614998) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9498001) q[2];
sx q[2];
rz(-0.99770498) q[2];
sx q[2];
rz(1.7157784) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9157) q[1];
sx q[1];
rz(-1.6927162) q[1];
sx q[1];
rz(1.8173161) q[1];
rz(-pi) q[2];
rz(-0.39059095) q[3];
sx q[3];
rz(-1.3814622) q[3];
sx q[3];
rz(-2.6238497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1714736) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(1.6469275) q[2];
rz(-1.0129048) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2038912) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(-1.0585693) q[0];
rz(-2.2189498) q[1];
sx q[1];
rz(-1.0170931) q[1];
sx q[1];
rz(3.0386818) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.218821) q[0];
sx q[0];
rz(-2.0334683) q[0];
sx q[0];
rz(-0.86313049) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58892624) q[2];
sx q[2];
rz(-1.3385596) q[2];
sx q[2];
rz(-0.015015451) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52572374) q[1];
sx q[1];
rz(-1.7547039) q[1];
sx q[1];
rz(2.1083819) q[1];
rz(-pi) q[2];
rz(-1.5779805) q[3];
sx q[3];
rz(-2.4251591) q[3];
sx q[3];
rz(0.35182692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4395478) q[2];
sx q[2];
rz(-0.44946686) q[2];
sx q[2];
rz(-0.31583819) q[2];
rz(0.99509197) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(2.4163213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.134326) q[0];
sx q[0];
rz(-2.8747929) q[0];
sx q[0];
rz(2.5079978) q[0];
rz(-2.7569356) q[1];
sx q[1];
rz(-1.9622012) q[1];
sx q[1];
rz(2.7124009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0338308) q[0];
sx q[0];
rz(-1.1035447) q[0];
sx q[0];
rz(0.28470566) q[0];
rz(-pi) q[1];
rz(1.6056988) q[2];
sx q[2];
rz(-0.88551025) q[2];
sx q[2];
rz(-1.3250145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8365578) q[1];
sx q[1];
rz(-1.9227131) q[1];
sx q[1];
rz(0.61493997) q[1];
x q[2];
rz(1.0121392) q[3];
sx q[3];
rz(-1.1691165) q[3];
sx q[3];
rz(-0.60810773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.089243285) q[2];
sx q[2];
rz(-0.77755916) q[2];
sx q[2];
rz(-0.26697978) q[2];
rz(3.0698245) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(-1.7223541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(-1.4005533) q[0];
rz(1.1446704) q[1];
sx q[1];
rz(-1.0619699) q[1];
sx q[1];
rz(2.5544419) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42204866) q[0];
sx q[0];
rz(-1.9912155) q[0];
sx q[0];
rz(-2.4499513) q[0];
x q[1];
rz(-1.4540423) q[2];
sx q[2];
rz(-1.9864621) q[2];
sx q[2];
rz(-2.1374201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59881951) q[1];
sx q[1];
rz(-2.2428998) q[1];
sx q[1];
rz(-1.2818976) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8244157) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(2.6871439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8553541) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(-2.722495) q[2];
rz(-1.6176443) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.5986901) q[0];
sx q[0];
rz(-1.9200696) q[0];
sx q[0];
rz(2.8322117) q[0];
rz(1.7912553) q[1];
sx q[1];
rz(-2.5534936) q[1];
sx q[1];
rz(2.5877171) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0623035) q[0];
sx q[0];
rz(-1.252018) q[0];
sx q[0];
rz(-2.7026524) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1872841) q[2];
sx q[2];
rz(-1.7085791) q[2];
sx q[2];
rz(-0.99064186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9123037) q[1];
sx q[1];
rz(-2.4472651) q[1];
sx q[1];
rz(-2.0635384) q[1];
rz(2.7782544) q[3];
sx q[3];
rz(-2.4060892) q[3];
sx q[3];
rz(-0.21061646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2498593) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(-1.7943133) q[2];
rz(2.3740718) q[3];
sx q[3];
rz(-1.1859505) q[3];
sx q[3];
rz(2.2931113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7193741) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(1.2147709) q[0];
rz(-2.1396554) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(-0.92438662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0262605) q[0];
sx q[0];
rz(-0.9965653) q[0];
sx q[0];
rz(-2.6747245) q[0];
rz(3.1141823) q[2];
sx q[2];
rz(-2.1152585) q[2];
sx q[2];
rz(-1.4278864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0707911) q[1];
sx q[1];
rz(-1.7940578) q[1];
sx q[1];
rz(1.5469458) q[1];
rz(0.59238418) q[3];
sx q[3];
rz(-0.96193704) q[3];
sx q[3];
rz(0.99623535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34710106) q[2];
sx q[2];
rz(-0.93239409) q[2];
sx q[2];
rz(-0.36716983) q[2];
rz(-1.9764887) q[3];
sx q[3];
rz(-2.1737183) q[3];
sx q[3];
rz(-2.2754106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1871056) q[0];
sx q[0];
rz(-2.4024873) q[0];
sx q[0];
rz(-3.0229229) q[0];
rz(-2.8754996) q[1];
sx q[1];
rz(-0.91231822) q[1];
sx q[1];
rz(-1.4516861) q[1];
rz(-2.8546147) q[2];
sx q[2];
rz(-1.8331883) q[2];
sx q[2];
rz(2.2856648) q[2];
rz(-0.96919555) q[3];
sx q[3];
rz(-1.7432913) q[3];
sx q[3];
rz(-0.023719214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
