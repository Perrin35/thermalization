OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(-0.40894142) q[0];
sx q[0];
rz(-0.9932819) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(-1.7239404) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87365681) q[0];
sx q[0];
rz(-1.9279772) q[0];
sx q[0];
rz(2.3620679) q[0];
rz(-pi) q[1];
rz(1.0101914) q[2];
sx q[2];
rz(-2.1935138) q[2];
sx q[2];
rz(1.6499008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4864145) q[1];
sx q[1];
rz(-1.298212) q[1];
sx q[1];
rz(0.64719871) q[1];
x q[2];
rz(1.3325813) q[3];
sx q[3];
rz(-3.0034667) q[3];
sx q[3];
rz(0.81616831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.950497) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(-1.0547868) q[2];
rz(-2.6705006) q[3];
sx q[3];
rz(-1.4548917) q[3];
sx q[3];
rz(-2.2733222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489814) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(0.23496041) q[0];
rz(-1.8670392) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(-2.4539006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4859568) q[0];
sx q[0];
rz(-1.5036426) q[0];
sx q[0];
rz(-0.32819076) q[0];
rz(2.2830487) q[2];
sx q[2];
rz(-1.1919824) q[2];
sx q[2];
rz(1.0114947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5968273) q[1];
sx q[1];
rz(-3.1278962) q[1];
sx q[1];
rz(-1.9910452) q[1];
x q[2];
rz(-0.95716547) q[3];
sx q[3];
rz(-1.6464697) q[3];
sx q[3];
rz(-1.5694973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71622744) q[2];
sx q[2];
rz(-1.3024412) q[2];
sx q[2];
rz(-0.20935527) q[2];
rz(-2.1685205) q[3];
sx q[3];
rz(-2.2416185) q[3];
sx q[3];
rz(-2.5488241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3000325) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(-5/(11*pi)) q[0];
rz(1.2380098) q[1];
sx q[1];
rz(-0.7716476) q[1];
sx q[1];
rz(0.091781052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973213) q[0];
sx q[0];
rz(-1.034212) q[0];
sx q[0];
rz(-2.7377605) q[0];
x q[1];
rz(3.0726542) q[2];
sx q[2];
rz(-1.4747915) q[2];
sx q[2];
rz(1.4323186) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57942332) q[1];
sx q[1];
rz(-1.3797288) q[1];
sx q[1];
rz(-1.3538989) q[1];
x q[2];
rz(0.067236891) q[3];
sx q[3];
rz(-1.1394258) q[3];
sx q[3];
rz(0.19156415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7564275) q[2];
sx q[2];
rz(-1.6331208) q[2];
sx q[2];
rz(0.37919322) q[2];
rz(2.3181629) q[3];
sx q[3];
rz(-2.7354666) q[3];
sx q[3];
rz(0.2894952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8150811) q[0];
sx q[0];
rz(-2.264475) q[0];
sx q[0];
rz(-2.1429578) q[0];
rz(2.6594992) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(-1.8236209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7319152) q[0];
sx q[0];
rz(-1.240726) q[0];
sx q[0];
rz(0.56185742) q[0];
x q[1];
rz(1.8179719) q[2];
sx q[2];
rz(-0.99946076) q[2];
sx q[2];
rz(-0.94117576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1554679) q[1];
sx q[1];
rz(-1.7805702) q[1];
sx q[1];
rz(-0.54297691) q[1];
x q[2];
rz(-2.34818) q[3];
sx q[3];
rz(-1.1856836) q[3];
sx q[3];
rz(2.957893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92129293) q[2];
sx q[2];
rz(-0.90196323) q[2];
sx q[2];
rz(0.48545066) q[2];
rz(0.38393936) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(-0.78401047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081472814) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(3.1299348) q[0];
rz(-0.13721379) q[1];
sx q[1];
rz(-1.5533181) q[1];
sx q[1];
rz(1.3020017) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0903895) q[0];
sx q[0];
rz(-2.9280991) q[0];
sx q[0];
rz(2.6043686) q[0];
rz(-pi) q[1];
rz(0.88413357) q[2];
sx q[2];
rz(-0.72956402) q[2];
sx q[2];
rz(-1.3183678) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1630262) q[1];
sx q[1];
rz(-2.8103376) q[1];
sx q[1];
rz(1.2322197) q[1];
x q[2];
rz(1.2106895) q[3];
sx q[3];
rz(-1.2146287) q[3];
sx q[3];
rz(0.83009133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7234708) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(-0.29329014) q[2];
rz(0.063118525) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(-2.2466834) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1191331) q[0];
sx q[0];
rz(-2.8526511) q[0];
sx q[0];
rz(0.038507842) q[0];
rz(2.0571845) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(-2.1748621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.181886) q[0];
sx q[0];
rz(-0.27508914) q[0];
sx q[0];
rz(-0.77232124) q[0];
x q[1];
rz(1.6183774) q[2];
sx q[2];
rz(-2.8741925) q[2];
sx q[2];
rz(-0.19879716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3075065) q[1];
sx q[1];
rz(-0.53552578) q[1];
sx q[1];
rz(-0.86043013) q[1];
x q[2];
rz(2.8916675) q[3];
sx q[3];
rz(-2.2062613) q[3];
sx q[3];
rz(-2.9912586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4806369) q[2];
sx q[2];
rz(-2.9065242) q[2];
sx q[2];
rz(0.98141518) q[2];
rz(1.6437982) q[3];
sx q[3];
rz(-1.2621597) q[3];
sx q[3];
rz(1.5952544) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129231) q[0];
sx q[0];
rz(-0.68083119) q[0];
sx q[0];
rz(0.31461) q[0];
rz(0.18868748) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(-2.6729118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.311694) q[0];
sx q[0];
rz(-1.440594) q[0];
sx q[0];
rz(0.24451406) q[0];
rz(-pi) q[1];
rz(1.5596703) q[2];
sx q[2];
rz(-0.16777786) q[2];
sx q[2];
rz(2.0027225) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51408053) q[1];
sx q[1];
rz(-2.2805164) q[1];
sx q[1];
rz(-0.40436097) q[1];
x q[2];
rz(-0.58705892) q[3];
sx q[3];
rz(-2.5958142) q[3];
sx q[3];
rz(-0.097028883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20588747) q[2];
sx q[2];
rz(-0.84172717) q[2];
sx q[2];
rz(0.97071281) q[2];
rz(0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(-0.71603388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.69990528) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(-1.4458789) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(-3.0605002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38302416) q[0];
sx q[0];
rz(-1.8580576) q[0];
sx q[0];
rz(-2.4455435) q[0];
x q[1];
rz(2.4862635) q[2];
sx q[2];
rz(-1.809279) q[2];
sx q[2];
rz(-2.2039889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1782185) q[1];
sx q[1];
rz(-1.6859173) q[1];
sx q[1];
rz(1.8646152) q[1];
x q[2];
rz(-2.567566) q[3];
sx q[3];
rz(-1.9732631) q[3];
sx q[3];
rz(2.5695679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9424092) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(2.8988163) q[2];
rz(1.5593922) q[3];
sx q[3];
rz(-0.42583164) q[3];
sx q[3];
rz(-1.2529713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2389857) q[0];
sx q[0];
rz(-2.1098397) q[0];
sx q[0];
rz(1.2549988) q[0];
rz(1.7651419) q[1];
sx q[1];
rz(-2.1170728) q[1];
sx q[1];
rz(0.66744101) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4723929) q[0];
sx q[0];
rz(-1.0697287) q[0];
sx q[0];
rz(2.7612812) q[0];
x q[1];
rz(1.5282643) q[2];
sx q[2];
rz(-2.0669524) q[2];
sx q[2];
rz(0.68788487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1245728) q[1];
sx q[1];
rz(-1.2806528) q[1];
sx q[1];
rz(0.97495074) q[1];
rz(-pi) q[2];
rz(2.4406901) q[3];
sx q[3];
rz(-0.48279219) q[3];
sx q[3];
rz(2.8822209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5218375) q[2];
sx q[2];
rz(-2.4027368) q[2];
sx q[2];
rz(-0.45836207) q[2];
rz(1.6058263) q[3];
sx q[3];
rz(-2.063844) q[3];
sx q[3];
rz(0.88466907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89852029) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(0.36548734) q[0];
rz(0.99994031) q[1];
sx q[1];
rz(-1.4391856) q[1];
sx q[1];
rz(1.684729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77068771) q[0];
sx q[0];
rz(-1.7720776) q[0];
sx q[0];
rz(-1.007249) q[0];
x q[1];
rz(0.11207726) q[2];
sx q[2];
rz(-1.3044895) q[2];
sx q[2];
rz(-0.2923686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.36023203) q[1];
sx q[1];
rz(-1.7443027) q[1];
sx q[1];
rz(-2.3886015) q[1];
rz(2.6014464) q[3];
sx q[3];
rz(-0.77278256) q[3];
sx q[3];
rz(0.55071044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(2.136039) q[2];
rz(0.65070659) q[3];
sx q[3];
rz(-1.4395827) q[3];
sx q[3];
rz(0.85056359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093813048) q[0];
sx q[0];
rz(-2.35738) q[0];
sx q[0];
rz(2.1398075) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(0.11713709) q[2];
sx q[2];
rz(-1.8979372) q[2];
sx q[2];
rz(-2.7330782) q[2];
rz(-0.10436124) q[3];
sx q[3];
rz(-1.4781393) q[3];
sx q[3];
rz(0.27235336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
