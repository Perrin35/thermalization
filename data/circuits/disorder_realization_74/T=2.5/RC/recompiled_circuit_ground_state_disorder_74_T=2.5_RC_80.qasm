OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7862406) q[0];
sx q[0];
rz(-2.8347637) q[0];
sx q[0];
rz(2.8428349) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(2.5112285) q[1];
sx q[1];
rz(11.093333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10160343) q[0];
sx q[0];
rz(-0.77058661) q[0];
sx q[0];
rz(1.6165074) q[0];
x q[1];
rz(-0.017919964) q[2];
sx q[2];
rz(-1.2700873) q[2];
sx q[2];
rz(-1.199388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2843282) q[1];
sx q[1];
rz(-0.92580399) q[1];
sx q[1];
rz(-0.93289341) q[1];
rz(-2.8657718) q[3];
sx q[3];
rz(-1.6383871) q[3];
sx q[3];
rz(-1.7662314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1851958) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(-1.3107276) q[2];
rz(-3.0500566) q[3];
sx q[3];
rz(-0.59713489) q[3];
sx q[3];
rz(-0.53102791) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(2.6614905) q[0];
rz(-1.9942888) q[1];
sx q[1];
rz(-2.5105748) q[1];
sx q[1];
rz(-2.4400585) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7362721) q[0];
sx q[0];
rz(-2.2715786) q[0];
sx q[0];
rz(2.2474656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.105179) q[2];
sx q[2];
rz(-1.948151) q[2];
sx q[2];
rz(1.5355664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.50276668) q[1];
sx q[1];
rz(-0.78710362) q[1];
sx q[1];
rz(2.2218641) q[1];
x q[2];
rz(-2.5040656) q[3];
sx q[3];
rz(-1.3204323) q[3];
sx q[3];
rz(1.375578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7445765) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(-1.5632632) q[2];
rz(2.8619316) q[3];
sx q[3];
rz(-0.71664387) q[3];
sx q[3];
rz(0.40951148) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(-1.4233587) q[0];
rz(3.0176945) q[1];
sx q[1];
rz(-2.2975477) q[1];
sx q[1];
rz(-1.5843676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4586055) q[0];
sx q[0];
rz(-1.8667954) q[0];
sx q[0];
rz(0.83403765) q[0];
rz(2.2940647) q[2];
sx q[2];
rz(-1.8516527) q[2];
sx q[2];
rz(0.47570634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.023097087) q[1];
sx q[1];
rz(-1.6498936) q[1];
sx q[1];
rz(-1.817607) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4067114) q[3];
sx q[3];
rz(-2.2874444) q[3];
sx q[3];
rz(1.411996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0763756) q[2];
sx q[2];
rz(-0.6535483) q[2];
sx q[2];
rz(-0.93488133) q[2];
rz(0.30385083) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76871753) q[0];
sx q[0];
rz(-0.77943742) q[0];
sx q[0];
rz(-2.8745162) q[0];
rz(3.0067387) q[1];
sx q[1];
rz(-0.57661533) q[1];
sx q[1];
rz(-2.3042302) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7202111) q[0];
sx q[0];
rz(-1.0307587) q[0];
sx q[0];
rz(0.311894) q[0];
rz(-pi) q[1];
rz(-0.83149399) q[2];
sx q[2];
rz(-2.0250426) q[2];
sx q[2];
rz(1.1513833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0935614) q[1];
sx q[1];
rz(-1.8037027) q[1];
sx q[1];
rz(-1.770411) q[1];
x q[2];
rz(3.0947826) q[3];
sx q[3];
rz(-1.6322281) q[3];
sx q[3];
rz(-0.84247473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59552938) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(-0.17892933) q[2];
rz(-1.1692125) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(-2.9235212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1201852) q[0];
sx q[0];
rz(-2.6299801) q[0];
sx q[0];
rz(2.2688493) q[0];
rz(-1.9841638) q[1];
sx q[1];
rz(-2.2667784) q[1];
sx q[1];
rz(-3.1246368) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8166148) q[0];
sx q[0];
rz(-0.8362174) q[0];
sx q[0];
rz(-1.5872383) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5418607) q[2];
sx q[2];
rz(-1.5625192) q[2];
sx q[2];
rz(-3.099583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3185721) q[1];
sx q[1];
rz(-1.8649532) q[1];
sx q[1];
rz(-2.3579954) q[1];
x q[2];
rz(2.8722829) q[3];
sx q[3];
rz(-1.8519173) q[3];
sx q[3];
rz(-0.65015974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0852647) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(-0.31923527) q[2];
rz(0.6790092) q[3];
sx q[3];
rz(-1.7947936) q[3];
sx q[3];
rz(0.80176789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7769258) q[0];
sx q[0];
rz(-2.8209782) q[0];
sx q[0];
rz(2.4989682) q[0];
rz(1.369426) q[1];
sx q[1];
rz(-2.5182928) q[1];
sx q[1];
rz(-2.1646037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95648051) q[0];
sx q[0];
rz(-1.6688884) q[0];
sx q[0];
rz(3.0177659) q[0];
rz(-pi) q[1];
rz(-0.80987038) q[2];
sx q[2];
rz(-2.4634482) q[2];
sx q[2];
rz(-0.43895753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20583892) q[1];
sx q[1];
rz(-0.69200695) q[1];
sx q[1];
rz(1.5490973) q[1];
rz(-0.053835458) q[3];
sx q[3];
rz(-2.1559058) q[3];
sx q[3];
rz(-1.9991635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6151108) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(1.9132445) q[2];
rz(-0.86722106) q[3];
sx q[3];
rz(-0.4129748) q[3];
sx q[3];
rz(-2.373608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7277471) q[0];
sx q[0];
rz(-0.21738805) q[0];
sx q[0];
rz(2.9687498) q[0];
rz(0.8051644) q[1];
sx q[1];
rz(-0.54904896) q[1];
sx q[1];
rz(-0.13596143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14912389) q[0];
sx q[0];
rz(-2.4276016) q[0];
sx q[0];
rz(1.7538422) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0375541) q[2];
sx q[2];
rz(-0.82471961) q[2];
sx q[2];
rz(-2.4989043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6579305) q[1];
sx q[1];
rz(-1.6854981) q[1];
sx q[1];
rz(1.5516267) q[1];
x q[2];
rz(2.7649859) q[3];
sx q[3];
rz(-2.7889502) q[3];
sx q[3];
rz(0.4761179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2148296) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(0.1845486) q[2];
rz(-2.9522225) q[3];
sx q[3];
rz(-0.53818494) q[3];
sx q[3];
rz(0.93809938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1592584) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(-2.2151997) q[0];
rz(0.039208086) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(2.1047986) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77039546) q[0];
sx q[0];
rz(-1.3950893) q[0];
sx q[0];
rz(1.8699339) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43208684) q[2];
sx q[2];
rz(-1.4054023) q[2];
sx q[2];
rz(1.337932) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8421887) q[1];
sx q[1];
rz(-1.328598) q[1];
sx q[1];
rz(2.7224225) q[1];
x q[2];
rz(-2.751334) q[3];
sx q[3];
rz(-1.7446127) q[3];
sx q[3];
rz(2.9040608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89256531) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(2.3426775) q[2];
rz(2.6246081) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(2.5911205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35140458) q[0];
sx q[0];
rz(-3.1365972) q[0];
sx q[0];
rz(1.5193526) q[0];
rz(-1.5478569) q[1];
sx q[1];
rz(-2.3212815) q[1];
sx q[1];
rz(0.69040745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0102745) q[0];
sx q[0];
rz(-1.2206435) q[0];
sx q[0];
rz(-1.8084333) q[0];
x q[1];
rz(1.4682653) q[2];
sx q[2];
rz(-2.3674115) q[2];
sx q[2];
rz(-0.098086327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71296299) q[1];
sx q[1];
rz(-2.3979514) q[1];
sx q[1];
rz(-2.3548467) q[1];
x q[2];
rz(-1.0318831) q[3];
sx q[3];
rz(-1.4461942) q[3];
sx q[3];
rz(0.74927688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(0.969886) q[2];
rz(0.25740933) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(1.359587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53512204) q[0];
sx q[0];
rz(-2.4115925) q[0];
sx q[0];
rz(3.0560793) q[0];
rz(2.8657148) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(1.1891018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333503) q[0];
sx q[0];
rz(-1.2062643) q[0];
sx q[0];
rz(-0.43594283) q[0];
rz(1.5797473) q[2];
sx q[2];
rz(-1.8807966) q[2];
sx q[2];
rz(1.5218889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.738742) q[1];
sx q[1];
rz(-1.1711118) q[1];
sx q[1];
rz(-3.012639) q[1];
rz(-2.889685) q[3];
sx q[3];
rz(-0.35874507) q[3];
sx q[3];
rz(-2.9210966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-0.87916547) q[2];
sx q[2];
rz(-0.48412588) q[2];
rz(3.0020946) q[3];
sx q[3];
rz(-2.3660584) q[3];
sx q[3];
rz(-0.16105306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6763247) q[0];
sx q[0];
rz(-1.9857255) q[0];
sx q[0];
rz(2.3406512) q[0];
rz(-1.6839266) q[1];
sx q[1];
rz(-1.6438345) q[1];
sx q[1];
rz(1.8297292) q[1];
rz(1.7745849) q[2];
sx q[2];
rz(-0.93175722) q[2];
sx q[2];
rz(2.5358806) q[2];
rz(-1.5994208) q[3];
sx q[3];
rz(-0.81394361) q[3];
sx q[3];
rz(1.5725556) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
