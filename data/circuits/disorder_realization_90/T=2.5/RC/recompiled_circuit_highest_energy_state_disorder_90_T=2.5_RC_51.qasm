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
rz(-2.3174602) q[0];
sx q[0];
rz(-1.1622518) q[0];
sx q[0];
rz(-0.21221575) q[0];
rz(0.68590391) q[1];
sx q[1];
rz(3.8984939) q[1];
sx q[1];
rz(9.6808351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4418418) q[0];
sx q[0];
rz(-0.76053333) q[0];
sx q[0];
rz(2.2470914) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4328622) q[2];
sx q[2];
rz(-0.56521691) q[2];
sx q[2];
rz(-2.3656453) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1009212) q[1];
sx q[1];
rz(-1.7771026) q[1];
sx q[1];
rz(-1.7357769) q[1];
x q[2];
rz(2.0689046) q[3];
sx q[3];
rz(-2.5513322) q[3];
sx q[3];
rz(-0.90567231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.2519007) q[2];
sx q[2];
rz(-0.32521853) q[2];
sx q[2];
rz(-2.0528117) q[2];
rz(2.661656) q[3];
sx q[3];
rz(-1.9710541) q[3];
sx q[3];
rz(-1.1133194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8460409) q[0];
sx q[0];
rz(-0.22838455) q[0];
sx q[0];
rz(-0.26902714) q[0];
rz(-0.61664063) q[1];
sx q[1];
rz(-0.65526217) q[1];
sx q[1];
rz(-2.2356967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3377853) q[0];
sx q[0];
rz(-1.4965987) q[0];
sx q[0];
rz(-2.3548749) q[0];
rz(-2.0935161) q[2];
sx q[2];
rz(-2.0459896) q[2];
sx q[2];
rz(1.1876729) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0720267) q[1];
sx q[1];
rz(-0.2827189) q[1];
sx q[1];
rz(2.7064884) q[1];
rz(-pi) q[2];
rz(-0.57635082) q[3];
sx q[3];
rz(-0.64223993) q[3];
sx q[3];
rz(3.0321966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1632094) q[2];
sx q[2];
rz(-2.3613112) q[2];
sx q[2];
rz(3.1261463) q[2];
rz(-0.77203006) q[3];
sx q[3];
rz(-0.49842146) q[3];
sx q[3];
rz(-1.8933403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60927272) q[0];
sx q[0];
rz(-0.63432223) q[0];
sx q[0];
rz(-2.3442205) q[0];
rz(-2.4728921) q[1];
sx q[1];
rz(-1.0182321) q[1];
sx q[1];
rz(-0.55577898) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41125339) q[0];
sx q[0];
rz(-2.3603435) q[0];
sx q[0];
rz(-2.8528105) q[0];
rz(0.90644353) q[2];
sx q[2];
rz(-1.4600233) q[2];
sx q[2];
rz(3.0701766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68986702) q[1];
sx q[1];
rz(-0.84459442) q[1];
sx q[1];
rz(1.5867342) q[1];
x q[2];
rz(1.1709514) q[3];
sx q[3];
rz(-2.2388864) q[3];
sx q[3];
rz(0.081351697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3681616) q[2];
sx q[2];
rz(-1.975446) q[2];
sx q[2];
rz(-1.2097166) q[2];
rz(-0.11956735) q[3];
sx q[3];
rz(-2.5366294) q[3];
sx q[3];
rz(2.1623478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46215737) q[0];
sx q[0];
rz(-1.4428416) q[0];
sx q[0];
rz(-0.47732842) q[0];
rz(3.0408472) q[1];
sx q[1];
rz(-2.0858177) q[1];
sx q[1];
rz(3.005952) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1208219) q[0];
sx q[0];
rz(-0.77540708) q[0];
sx q[0];
rz(2.1935619) q[0];
rz(-pi) q[1];
rz(1.2775947) q[2];
sx q[2];
rz(-2.0005156) q[2];
sx q[2];
rz(-1.7079084) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98235929) q[1];
sx q[1];
rz(-1.9894588) q[1];
sx q[1];
rz(-1.7348947) q[1];
rz(1.8036912) q[3];
sx q[3];
rz(-0.80162382) q[3];
sx q[3];
rz(0.1375141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.2123083) q[2];
sx q[2];
rz(-1.4886651) q[2];
sx q[2];
rz(0.51367122) q[2];
rz(-2.9786927) q[3];
sx q[3];
rz(-2.8859911) q[3];
sx q[3];
rz(-2.0388849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.61360079) q[0];
sx q[0];
rz(-0.010471424) q[0];
sx q[0];
rz(-2.7136059) q[0];
rz(-1.1780659) q[1];
sx q[1];
rz(-1.4617498) q[1];
sx q[1];
rz(-1.4129432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1278541) q[0];
sx q[0];
rz(-1.1098521) q[0];
sx q[0];
rz(-2.7147033) q[0];
rz(-pi) q[1];
rz(-2.2828083) q[2];
sx q[2];
rz(-3.0860215) q[2];
sx q[2];
rz(2.839599) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.481491) q[1];
sx q[1];
rz(-2.8347184) q[1];
sx q[1];
rz(0.18395384) q[1];
rz(2.574202) q[3];
sx q[3];
rz(-1.9741283) q[3];
sx q[3];
rz(1.698026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.596375) q[2];
sx q[2];
rz(-1.7866106) q[2];
sx q[2];
rz(2.6522816) q[2];
rz(-2.4423068) q[3];
sx q[3];
rz(-2.0977061) q[3];
sx q[3];
rz(-1.4780686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74772239) q[0];
sx q[0];
rz(-2.1339397) q[0];
sx q[0];
rz(1.7973768) q[0];
rz(-2.4463553) q[1];
sx q[1];
rz(-1.8724915) q[1];
sx q[1];
rz(0.42323798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76038928) q[0];
sx q[0];
rz(-1.7402599) q[0];
sx q[0];
rz(1.5080835) q[0];
rz(-pi) q[1];
rz(3.0401523) q[2];
sx q[2];
rz(-0.82221088) q[2];
sx q[2];
rz(1.492983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6072487) q[1];
sx q[1];
rz(-2.3869276) q[1];
sx q[1];
rz(0.21790803) q[1];
rz(-pi) q[2];
rz(0.18169348) q[3];
sx q[3];
rz(-1.6258929) q[3];
sx q[3];
rz(-2.4247935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77373475) q[2];
sx q[2];
rz(-2.9025142) q[2];
sx q[2];
rz(-1.2131005) q[2];
rz(3.0554981) q[3];
sx q[3];
rz(-1.2505069) q[3];
sx q[3];
rz(1.7566173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4940779) q[0];
sx q[0];
rz(-2.8120742) q[0];
sx q[0];
rz(2.9373017) q[0];
rz(-0.34945166) q[1];
sx q[1];
rz(-1.5710257) q[1];
sx q[1];
rz(-1.6250767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1384927) q[0];
sx q[0];
rz(-2.6647143) q[0];
sx q[0];
rz(1.80507) q[0];
rz(-pi) q[1];
rz(-1.1687241) q[2];
sx q[2];
rz(-1.8622145) q[2];
sx q[2];
rz(2.5411487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2353846) q[1];
sx q[1];
rz(-1.6170119) q[1];
sx q[1];
rz(-1.4704204) q[1];
rz(-pi) q[2];
x q[2];
rz(0.013850526) q[3];
sx q[3];
rz(-0.99273721) q[3];
sx q[3];
rz(-0.25307551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1427631) q[2];
sx q[2];
rz(-2.5122354) q[2];
sx q[2];
rz(-0.91267419) q[2];
rz(3.1257889) q[3];
sx q[3];
rz(-2.0533357) q[3];
sx q[3];
rz(-0.44939941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939746) q[0];
sx q[0];
rz(-2.0517218) q[0];
sx q[0];
rz(0.81899482) q[0];
rz(0.35797572) q[1];
sx q[1];
rz(-2.0134108) q[1];
sx q[1];
rz(-2.5111759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8837161) q[0];
sx q[0];
rz(-1.64332) q[0];
sx q[0];
rz(3.1255015) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2906456) q[2];
sx q[2];
rz(-1.2149879) q[2];
sx q[2];
rz(-1.8263888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9832657) q[1];
sx q[1];
rz(-0.7559146) q[1];
sx q[1];
rz(-2.1793038) q[1];
rz(-0.94576337) q[3];
sx q[3];
rz(-1.182036) q[3];
sx q[3];
rz(0.22076535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36998746) q[2];
sx q[2];
rz(-0.5094173) q[2];
sx q[2];
rz(-2.7484757) q[2];
rz(-2.3734132) q[3];
sx q[3];
rz(-0.79777515) q[3];
sx q[3];
rz(-1.2099077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8623675) q[0];
sx q[0];
rz(-2.6287862) q[0];
sx q[0];
rz(-2.7592036) q[0];
rz(2.6928316) q[1];
sx q[1];
rz(-0.02350137) q[1];
sx q[1];
rz(0.99865595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778393) q[0];
sx q[0];
rz(-2.4344517) q[0];
sx q[0];
rz(-0.36974837) q[0];
x q[1];
rz(1.9085541) q[2];
sx q[2];
rz(-1.6509368) q[2];
sx q[2];
rz(0.077827443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4289394) q[1];
sx q[1];
rz(-1.3725159) q[1];
sx q[1];
rz(-2.8444931) q[1];
x q[2];
rz(3.0581045) q[3];
sx q[3];
rz(-0.69056702) q[3];
sx q[3];
rz(2.4960828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0950332) q[2];
sx q[2];
rz(-0.79594374) q[2];
sx q[2];
rz(-0.98540068) q[2];
rz(-0.69917786) q[3];
sx q[3];
rz(-0.85088426) q[3];
sx q[3];
rz(-3.0881171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92626524) q[0];
sx q[0];
rz(-0.43020058) q[0];
sx q[0];
rz(-2.3660124) q[0];
rz(-2.8714478) q[1];
sx q[1];
rz(-1.6855449) q[1];
sx q[1];
rz(-0.57033479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1352976) q[0];
sx q[0];
rz(-2.0699602) q[0];
sx q[0];
rz(2.1111302) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3675466) q[2];
sx q[2];
rz(-0.48664868) q[2];
sx q[2];
rz(1.987285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8034035) q[1];
sx q[1];
rz(-1.3373151) q[1];
sx q[1];
rz(1.6031053) q[1];
x q[2];
rz(1.6472642) q[3];
sx q[3];
rz(-2.2588393) q[3];
sx q[3];
rz(2.9043759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21738416) q[2];
sx q[2];
rz(-1.9839828) q[2];
sx q[2];
rz(2.6317934) q[2];
rz(-2.9057251) q[3];
sx q[3];
rz(-0.1736975) q[3];
sx q[3];
rz(2.164446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8248642) q[0];
sx q[0];
rz(-1.4126128) q[0];
sx q[0];
rz(1.8870507) q[0];
rz(-2.6084689) q[1];
sx q[1];
rz(-1.6785379) q[1];
sx q[1];
rz(1.0482845) q[1];
rz(1.6021846) q[2];
sx q[2];
rz(-2.7161408) q[2];
sx q[2];
rz(-1.3795992) q[2];
rz(-1.6857311) q[3];
sx q[3];
rz(-2.0685932) q[3];
sx q[3];
rz(2.6017068) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
