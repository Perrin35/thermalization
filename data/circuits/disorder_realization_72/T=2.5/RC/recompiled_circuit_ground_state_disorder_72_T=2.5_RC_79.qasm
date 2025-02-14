OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(-1.0847217) q[0];
rz(1.9864858) q[1];
sx q[1];
rz(-2.3218563) q[1];
sx q[1];
rz(-2.2302332) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93854967) q[0];
sx q[0];
rz(-1.3289551) q[0];
sx q[0];
rz(-2.3641402) q[0];
rz(-1.7167822) q[2];
sx q[2];
rz(-0.8643736) q[2];
sx q[2];
rz(-1.3024769) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45381308) q[1];
sx q[1];
rz(-1.1247083) q[1];
sx q[1];
rz(2.5519754) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0769071) q[3];
sx q[3];
rz(-2.0035335) q[3];
sx q[3];
rz(-2.329934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34065166) q[2];
sx q[2];
rz(-1.1110577) q[2];
sx q[2];
rz(-3.0621373) q[2];
rz(-2.5331412) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(2.2505545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4873753) q[0];
sx q[0];
rz(-2.2779164) q[0];
sx q[0];
rz(1.7864216) q[0];
rz(-1.0379418) q[1];
sx q[1];
rz(-2.4619921) q[1];
sx q[1];
rz(2.3341446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12804764) q[0];
sx q[0];
rz(-0.68251784) q[0];
sx q[0];
rz(-1.1004992) q[0];
rz(-pi) q[1];
rz(2.1218929) q[2];
sx q[2];
rz(-0.65097133) q[2];
sx q[2];
rz(2.5215182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5521621) q[1];
sx q[1];
rz(-2.8928601) q[1];
sx q[1];
rz(-2.8692607) q[1];
x q[2];
rz(1.7044675) q[3];
sx q[3];
rz(-1.6764056) q[3];
sx q[3];
rz(-1.228412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9281533) q[2];
sx q[2];
rz(-1.5636874) q[2];
sx q[2];
rz(-1.6820924) q[2];
rz(0.45808074) q[3];
sx q[3];
rz(-2.5638678) q[3];
sx q[3];
rz(-0.95988449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9816575) q[0];
sx q[0];
rz(-1.0502879) q[0];
sx q[0];
rz(-2.4666069) q[0];
rz(0.58810294) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(1.810422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0494306) q[0];
sx q[0];
rz(-1.9665008) q[0];
sx q[0];
rz(0.53073287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5187732) q[2];
sx q[2];
rz(-2.3732936) q[2];
sx q[2];
rz(1.7211421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4294376) q[1];
sx q[1];
rz(-0.78199103) q[1];
sx q[1];
rz(-1.9351134) q[1];
rz(-pi) q[2];
rz(2.939091) q[3];
sx q[3];
rz(-2.0612067) q[3];
sx q[3];
rz(-0.59668505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5171234) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(2.9208753) q[2];
rz(2.3366426) q[3];
sx q[3];
rz(-1.5996108) q[3];
sx q[3];
rz(1.8055003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6014366) q[0];
sx q[0];
rz(-2.8995081) q[0];
sx q[0];
rz(2.5228187) q[0];
rz(-0.082911804) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(-1.6212911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12099685) q[0];
sx q[0];
rz(-1.474664) q[0];
sx q[0];
rz(2.0970048) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5624406) q[2];
sx q[2];
rz(-1.7699827) q[2];
sx q[2];
rz(2.3483089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6449086) q[1];
sx q[1];
rz(-2.0808947) q[1];
sx q[1];
rz(0.094385191) q[1];
rz(-0.374745) q[3];
sx q[3];
rz(-1.7386769) q[3];
sx q[3];
rz(1.8207267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80359047) q[2];
sx q[2];
rz(-3.0323196) q[2];
sx q[2];
rz(2.9962311) q[2];
rz(-0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(-1.9947778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1266601) q[0];
sx q[0];
rz(-1.3260051) q[0];
sx q[0];
rz(0.10619157) q[0];
rz(-2.1821678) q[1];
sx q[1];
rz(-0.54490772) q[1];
sx q[1];
rz(-0.19439654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0036164408) q[0];
sx q[0];
rz(-2.4038393) q[0];
sx q[0];
rz(-3.0135148) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6856367) q[2];
sx q[2];
rz(-1.5756338) q[2];
sx q[2];
rz(-0.23813914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.186708) q[1];
sx q[1];
rz(-0.59046035) q[1];
sx q[1];
rz(0.073530274) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3931469) q[3];
sx q[3];
rz(-0.41209778) q[3];
sx q[3];
rz(0.92890152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83547366) q[2];
sx q[2];
rz(-1.4543616) q[2];
sx q[2];
rz(2.5731738) q[2];
rz(-3.0889555) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(-0.13681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1133872) q[0];
sx q[0];
rz(-2.5882692) q[0];
sx q[0];
rz(-0.18948874) q[0];
rz(2.1151309) q[1];
sx q[1];
rz(-2.2920513) q[1];
sx q[1];
rz(2.386327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24863923) q[0];
sx q[0];
rz(-1.0414855) q[0];
sx q[0];
rz(-2.3387699) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.009047) q[2];
sx q[2];
rz(-2.0239186) q[2];
sx q[2];
rz(2.8293186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89596302) q[1];
sx q[1];
rz(-1.4019483) q[1];
sx q[1];
rz(-2.9814475) q[1];
rz(-1.2047661) q[3];
sx q[3];
rz(-1.9037953) q[3];
sx q[3];
rz(2.7878417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75817627) q[2];
sx q[2];
rz(-1.5551609) q[2];
sx q[2];
rz(-1.4413393) q[2];
rz(1.7656743) q[3];
sx q[3];
rz(-0.91840363) q[3];
sx q[3];
rz(0.70022303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.43948424) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(-1.1442319) q[0];
rz(0.2598091) q[1];
sx q[1];
rz(-2.3016498) q[1];
sx q[1];
rz(0.98794404) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086083273) q[0];
sx q[0];
rz(-1.4713773) q[0];
sx q[0];
rz(-1.587633) q[0];
rz(-1.9960711) q[2];
sx q[2];
rz(-2.2605611) q[2];
sx q[2];
rz(2.9381616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0007947) q[1];
sx q[1];
rz(-2.4355445) q[1];
sx q[1];
rz(2.7435859) q[1];
x q[2];
rz(-2.0201398) q[3];
sx q[3];
rz(-0.62329037) q[3];
sx q[3];
rz(-1.2227071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7320431) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(-2.4260855) q[3];
sx q[3];
rz(-2.3698273) q[3];
sx q[3];
rz(-2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20329423) q[0];
sx q[0];
rz(-2.2977915) q[0];
sx q[0];
rz(2.8734558) q[0];
rz(-2.9109491) q[1];
sx q[1];
rz(-1.5430887) q[1];
sx q[1];
rz(3.1275829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.926353) q[0];
sx q[0];
rz(-0.82414675) q[0];
sx q[0];
rz(1.8570379) q[0];
rz(-1.1300541) q[2];
sx q[2];
rz(-1.4258372) q[2];
sx q[2];
rz(-0.10565378) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3611759) q[1];
sx q[1];
rz(-1.0339972) q[1];
sx q[1];
rz(-1.518834) q[1];
rz(-pi) q[2];
rz(-1.3669861) q[3];
sx q[3];
rz(-0.14575726) q[3];
sx q[3];
rz(0.83728079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9739428) q[2];
sx q[2];
rz(-1.3385945) q[2];
sx q[2];
rz(0.30926427) q[2];
rz(-2.9850128) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(2.1046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88368791) q[0];
sx q[0];
rz(-0.64249277) q[0];
sx q[0];
rz(-0.25074348) q[0];
rz(-0.58654395) q[1];
sx q[1];
rz(-1.5778912) q[1];
sx q[1];
rz(2.3768545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79478489) q[0];
sx q[0];
rz(-1.8612004) q[0];
sx q[0];
rz(1.8998763) q[0];
rz(0.20275764) q[2];
sx q[2];
rz(-2.360095) q[2];
sx q[2];
rz(-2.4818713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34652485) q[1];
sx q[1];
rz(-1.5800313) q[1];
sx q[1];
rz(3.0187876) q[1];
x q[2];
rz(-0.13651092) q[3];
sx q[3];
rz(-2.7595466) q[3];
sx q[3];
rz(0.34367022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6278729) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(2.9676843) q[2];
rz(-0.74226132) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-3.098587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0088418) q[0];
sx q[0];
rz(-2.371377) q[0];
sx q[0];
rz(0.26790628) q[0];
rz(1.335089) q[1];
sx q[1];
rz(-0.91064149) q[1];
sx q[1];
rz(-1.7693899) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55937476) q[0];
sx q[0];
rz(-1.1211044) q[0];
sx q[0];
rz(0.26828464) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29628654) q[2];
sx q[2];
rz(-1.6673267) q[2];
sx q[2];
rz(1.5636843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42555537) q[1];
sx q[1];
rz(-0.73409427) q[1];
sx q[1];
rz(-2.3922763) q[1];
rz(2.3389111) q[3];
sx q[3];
rz(-2.6978328) q[3];
sx q[3];
rz(-0.26324031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.09769663) q[2];
sx q[2];
rz(-2.1052269) q[2];
sx q[2];
rz(-1.3807266) q[2];
rz(-2.0128287) q[3];
sx q[3];
rz(-2.1916316) q[3];
sx q[3];
rz(-1.5560163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8783405) q[0];
sx q[0];
rz(-1.6077519) q[0];
sx q[0];
rz(2.0335249) q[0];
rz(0.68645984) q[1];
sx q[1];
rz(-1.7484799) q[1];
sx q[1];
rz(1.9304986) q[1];
rz(1.8615234) q[2];
sx q[2];
rz(-2.0474993) q[2];
sx q[2];
rz(-1.7119424) q[2];
rz(-0.93871212) q[3];
sx q[3];
rz(-0.6365255) q[3];
sx q[3];
rz(-0.33088007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
