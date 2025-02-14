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
rz(1.3156112) q[0];
sx q[0];
rz(-0.82807461) q[0];
sx q[0];
rz(-0.84870422) q[0];
rz(-1.6500213) q[1];
sx q[1];
rz(-2.4529011) q[1];
sx q[1];
rz(-0.42660776) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72683452) q[0];
sx q[0];
rz(-1.6462741) q[0];
sx q[0];
rz(3.0407136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83580221) q[2];
sx q[2];
rz(-2.9907132) q[2];
sx q[2];
rz(2.7250233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9200807) q[1];
sx q[1];
rz(-2.6903841) q[1];
sx q[1];
rz(0.16772224) q[1];
rz(-pi) q[2];
rz(2.3412104) q[3];
sx q[3];
rz(-1.1100148) q[3];
sx q[3];
rz(-0.15298259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2171057) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(-3.1021049) q[2];
rz(2.110179) q[3];
sx q[3];
rz(-2.1125427) q[3];
sx q[3];
rz(1.2379117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6295488) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(-1.7426096) q[0];
rz(1.5139187) q[1];
sx q[1];
rz(-0.84295034) q[1];
sx q[1];
rz(1.4167851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48081452) q[0];
sx q[0];
rz(-0.82624683) q[0];
sx q[0];
rz(1.8742524) q[0];
rz(-pi) q[1];
rz(-2.3624647) q[2];
sx q[2];
rz(-1.7567524) q[2];
sx q[2];
rz(0.27006876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4128215) q[1];
sx q[1];
rz(-1.9660453) q[1];
sx q[1];
rz(-1.7186307) q[1];
rz(-2.7339277) q[3];
sx q[3];
rz(-1.614794) q[3];
sx q[3];
rz(0.10213479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17025718) q[2];
sx q[2];
rz(-2.038326) q[2];
sx q[2];
rz(0.7106759) q[2];
rz(-0.014287861) q[3];
sx q[3];
rz(-0.24737869) q[3];
sx q[3];
rz(1.3361196) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9209442) q[0];
sx q[0];
rz(-2.1823688) q[0];
sx q[0];
rz(-0.4162108) q[0];
rz(1.2843708) q[1];
sx q[1];
rz(-1.4764079) q[1];
sx q[1];
rz(2.823337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98781768) q[0];
sx q[0];
rz(-1.9526228) q[0];
sx q[0];
rz(2.8433958) q[0];
x q[1];
rz(0.84789611) q[2];
sx q[2];
rz(-2.090914) q[2];
sx q[2];
rz(1.073357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8313731) q[1];
sx q[1];
rz(-1.7064269) q[1];
sx q[1];
rz(-2.4844205) q[1];
x q[2];
rz(1.2622358) q[3];
sx q[3];
rz(-1.3342382) q[3];
sx q[3];
rz(2.6738338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8754114) q[2];
sx q[2];
rz(-2.3616932) q[2];
sx q[2];
rz(-1.8511339) q[2];
rz(3.087888) q[3];
sx q[3];
rz(-1.3196245) q[3];
sx q[3];
rz(-1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5829492) q[0];
sx q[0];
rz(-0.68840331) q[0];
sx q[0];
rz(-2.131856) q[0];
rz(2.2263777) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-2.7161652) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4998326) q[0];
sx q[0];
rz(-1.7374037) q[0];
sx q[0];
rz(-0.84979041) q[0];
rz(2.3138626) q[2];
sx q[2];
rz(-1.6249738) q[2];
sx q[2];
rz(-2.3119213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.912251) q[1];
sx q[1];
rz(-2.3315363) q[1];
sx q[1];
rz(2.1476169) q[1];
rz(-pi) q[2];
rz(-1.6050379) q[3];
sx q[3];
rz(-2.9786907) q[3];
sx q[3];
rz(0.77252856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2414744) q[2];
sx q[2];
rz(-1.5510617) q[2];
sx q[2];
rz(0.11108622) q[2];
rz(-0.19242081) q[3];
sx q[3];
rz(-2.7399053) q[3];
sx q[3];
rz(-0.69453159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7531994) q[0];
sx q[0];
rz(-0.73047262) q[0];
sx q[0];
rz(2.2764192) q[0];
rz(-0.49258891) q[1];
sx q[1];
rz(-2.045423) q[1];
sx q[1];
rz(-1.7561087) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82768142) q[0];
sx q[0];
rz(-1.4340785) q[0];
sx q[0];
rz(0.78929957) q[0];
rz(2.3326068) q[2];
sx q[2];
rz(-1.8833767) q[2];
sx q[2];
rz(2.9650826) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4069479) q[1];
sx q[1];
rz(-2.4666767) q[1];
sx q[1];
rz(2.4823443) q[1];
rz(-pi) q[2];
rz(0.68180696) q[3];
sx q[3];
rz(-0.94577937) q[3];
sx q[3];
rz(2.8375219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8869141) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(2.4349507) q[2];
rz(0.32736579) q[3];
sx q[3];
rz(-1.2923765) q[3];
sx q[3];
rz(-0.65739003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3051598) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(-2.7850372) q[0];
rz(-0.56218475) q[1];
sx q[1];
rz(-1.1327344) q[1];
sx q[1];
rz(-0.98639375) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1323469) q[0];
sx q[0];
rz(-0.86902394) q[0];
sx q[0];
rz(-2.470507) q[0];
rz(-2.3190034) q[2];
sx q[2];
rz(-0.486136) q[2];
sx q[2];
rz(1.2631455) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79365208) q[1];
sx q[1];
rz(-0.75885443) q[1];
sx q[1];
rz(0.037686613) q[1];
x q[2];
rz(2.2067542) q[3];
sx q[3];
rz(-0.29948161) q[3];
sx q[3];
rz(0.071223473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3761313) q[2];
sx q[2];
rz(-2.4262846) q[2];
sx q[2];
rz(0.7005271) q[2];
rz(-0.96945196) q[3];
sx q[3];
rz(-1.0337831) q[3];
sx q[3];
rz(0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1932061) q[0];
sx q[0];
rz(-2.8128615) q[0];
sx q[0];
rz(2.9902003) q[0];
rz(-1.6311749) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(-1.4600533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3356381) q[0];
sx q[0];
rz(-0.80712576) q[0];
sx q[0];
rz(-2.4390079) q[0];
rz(-2.7786656) q[2];
sx q[2];
rz(-2.0690326) q[2];
sx q[2];
rz(-1.0403596) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5269176) q[1];
sx q[1];
rz(-1.9567031) q[1];
sx q[1];
rz(-0.98395394) q[1];
x q[2];
rz(2.9956508) q[3];
sx q[3];
rz(-1.3781007) q[3];
sx q[3];
rz(1.3789498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2363362) q[2];
sx q[2];
rz(-0.6826123) q[2];
sx q[2];
rz(-0.35161099) q[2];
rz(-2.6998399) q[3];
sx q[3];
rz(-2.2124002) q[3];
sx q[3];
rz(2.9898804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0996967) q[0];
sx q[0];
rz(-1.2956887) q[0];
sx q[0];
rz(1.1610485) q[0];
rz(2.4940122) q[1];
sx q[1];
rz(-1.7627962) q[1];
sx q[1];
rz(-2.3191998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0121545) q[0];
sx q[0];
rz(-2.4732051) q[0];
sx q[0];
rz(0.8968312) q[0];
rz(-pi) q[1];
rz(-2.6641409) q[2];
sx q[2];
rz(-2.3465354) q[2];
sx q[2];
rz(-0.24518046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31588512) q[1];
sx q[1];
rz(-1.2824821) q[1];
sx q[1];
rz(-0.79135367) q[1];
rz(-2.0045207) q[3];
sx q[3];
rz(-1.2081283) q[3];
sx q[3];
rz(0.39946242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9245727) q[2];
sx q[2];
rz(-1.9968888) q[2];
sx q[2];
rz(1.8355969) q[2];
rz(-0.71803391) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(-0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.43593916) q[0];
sx q[0];
rz(-2.0581364) q[0];
sx q[0];
rz(-3.1245681) q[0];
rz(1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(0.33504018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3075313) q[0];
sx q[0];
rz(-2.1387707) q[0];
sx q[0];
rz(-3.1295883) q[0];
rz(-2.4990988) q[2];
sx q[2];
rz(-1.9064925) q[2];
sx q[2];
rz(-0.030980008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2372053) q[1];
sx q[1];
rz(-0.6977302) q[1];
sx q[1];
rz(1.2751139) q[1];
rz(-pi) q[2];
rz(-2.0552956) q[3];
sx q[3];
rz(-1.1346176) q[3];
sx q[3];
rz(2.0791813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7177141) q[2];
sx q[2];
rz(-2.3953891) q[2];
sx q[2];
rz(0.27900532) q[2];
rz(1.772607) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(-2.2122808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170711) q[0];
sx q[0];
rz(-0.14685024) q[0];
sx q[0];
rz(2.3745234) q[0];
rz(2.7686139) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(0.17072089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0559346) q[0];
sx q[0];
rz(-1.5996337) q[0];
sx q[0];
rz(-0.7349073) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4095979) q[2];
sx q[2];
rz(-1.1606163) q[2];
sx q[2];
rz(0.31482492) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0788012) q[1];
sx q[1];
rz(-1.2578221) q[1];
sx q[1];
rz(-0.2614577) q[1];
rz(1.7307698) q[3];
sx q[3];
rz(-2.0334646) q[3];
sx q[3];
rz(1.3079974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79147044) q[2];
sx q[2];
rz(-2.5779724) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(0.9015829) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(-2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053454178) q[0];
sx q[0];
rz(-0.32129856) q[0];
sx q[0];
rz(-1.970485) q[0];
rz(0.83012719) q[1];
sx q[1];
rz(-1.8142038) q[1];
sx q[1];
rz(1.289191) q[1];
rz(-0.51132497) q[2];
sx q[2];
rz(-2.1187388) q[2];
sx q[2];
rz(-0.50730898) q[2];
rz(-2.4975111) q[3];
sx q[3];
rz(-0.91337503) q[3];
sx q[3];
rz(2.5144284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
