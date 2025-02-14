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
rz(2.7316982) q[0];
sx q[0];
rz(-1.7562261) q[0];
sx q[0];
rz(2.07055) q[0];
rz(1.8154124) q[1];
sx q[1];
rz(-0.92547995) q[1];
sx q[1];
rz(2.7993536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54968327) q[0];
sx q[0];
rz(-0.91381493) q[0];
sx q[0];
rz(2.6651938) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5985122) q[2];
sx q[2];
rz(-2.2545071) q[2];
sx q[2];
rz(-1.3165084) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9896868) q[1];
sx q[1];
rz(-1.1380956) q[1];
sx q[1];
rz(-3.0244083) q[1];
rz(0.017090509) q[3];
sx q[3];
rz(-1.9668764) q[3];
sx q[3];
rz(-1.9310967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5388415) q[2];
sx q[2];
rz(-0.72776908) q[2];
sx q[2];
rz(1.7068498) q[2];
rz(-0.96628609) q[3];
sx q[3];
rz(-1.1937701) q[3];
sx q[3];
rz(0.57893354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12427881) q[0];
sx q[0];
rz(-0.74324981) q[0];
sx q[0];
rz(0.055835128) q[0];
rz(2.3827379) q[1];
sx q[1];
rz(-2.0003624) q[1];
sx q[1];
rz(2.8208044) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13981217) q[0];
sx q[0];
rz(-2.3235112) q[0];
sx q[0];
rz(-1.7745738) q[0];
rz(0.77938883) q[2];
sx q[2];
rz(-1.8223127) q[2];
sx q[2];
rz(-0.0096461065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6763416) q[1];
sx q[1];
rz(-0.23561978) q[1];
sx q[1];
rz(1.444677) q[1];
rz(2.4092968) q[3];
sx q[3];
rz(-1.2114015) q[3];
sx q[3];
rz(0.72903192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8212829) q[2];
sx q[2];
rz(-1.8822957) q[2];
sx q[2];
rz(0.28309509) q[2];
rz(2.8458332) q[3];
sx q[3];
rz(-1.7167973) q[3];
sx q[3];
rz(1.6368216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5931684) q[0];
sx q[0];
rz(-1.230509) q[0];
sx q[0];
rz(-3.044627) q[0];
rz(-1.4292258) q[1];
sx q[1];
rz(-1.2173419) q[1];
sx q[1];
rz(0.355535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9622346) q[0];
sx q[0];
rz(-1.0197687) q[0];
sx q[0];
rz(1.2310811) q[0];
rz(-1.468268) q[2];
sx q[2];
rz(-1.9838196) q[2];
sx q[2];
rz(-2.6669974) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83427338) q[1];
sx q[1];
rz(-1.6641698) q[1];
sx q[1];
rz(1.1591256) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.612298) q[3];
sx q[3];
rz(-1.8222162) q[3];
sx q[3];
rz(1.4700898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43704438) q[2];
sx q[2];
rz(-1.1532249) q[2];
sx q[2];
rz(-0.28771773) q[2];
rz(-0.25117609) q[3];
sx q[3];
rz(-2.0468678) q[3];
sx q[3];
rz(-1.8672966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1890202) q[0];
sx q[0];
rz(-2.548521) q[0];
sx q[0];
rz(-2.3841542) q[0];
rz(-2.3144552) q[1];
sx q[1];
rz(-0.65991455) q[1];
sx q[1];
rz(1.3077024) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73440512) q[0];
sx q[0];
rz(-1.4467753) q[0];
sx q[0];
rz(3.0888686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1266379) q[2];
sx q[2];
rz(-0.53947811) q[2];
sx q[2];
rz(-1.2459823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74119905) q[1];
sx q[1];
rz(-0.82655061) q[1];
sx q[1];
rz(-1.939102) q[1];
rz(1.3321628) q[3];
sx q[3];
rz(-0.25128579) q[3];
sx q[3];
rz(-0.3005614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4266251) q[2];
sx q[2];
rz(-1.7968618) q[2];
sx q[2];
rz(1.7460543) q[2];
rz(1.7914145) q[3];
sx q[3];
rz(-1.502864) q[3];
sx q[3];
rz(-3.1135528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.79513079) q[0];
sx q[0];
rz(-2.6030354) q[0];
sx q[0];
rz(1.7002456) q[0];
rz(0.92789188) q[1];
sx q[1];
rz(-0.83506942) q[1];
sx q[1];
rz(2.3477614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5466325) q[0];
sx q[0];
rz(-1.4671456) q[0];
sx q[0];
rz(-1.1569381) q[0];
x q[1];
rz(-0.16818856) q[2];
sx q[2];
rz(-1.6861746) q[2];
sx q[2];
rz(-0.66957973) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87463996) q[1];
sx q[1];
rz(-0.62332223) q[1];
sx q[1];
rz(-1.0250574) q[1];
rz(-pi) q[2];
rz(-2.6364274) q[3];
sx q[3];
rz(-1.9117711) q[3];
sx q[3];
rz(2.7974432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4980207) q[2];
sx q[2];
rz(-1.6732432) q[2];
sx q[2];
rz(-2.7872861) q[2];
rz(-1.0962567) q[3];
sx q[3];
rz(-2.8772964) q[3];
sx q[3];
rz(1.8335584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7588014) q[0];
sx q[0];
rz(-0.76171869) q[0];
sx q[0];
rz(0.85269165) q[0];
rz(-0.26607749) q[1];
sx q[1];
rz(-1.0238901) q[1];
sx q[1];
rz(-1.783225) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.81937) q[0];
sx q[0];
rz(-2.529105) q[0];
sx q[0];
rz(-1.9507381) q[0];
rz(1.9044962) q[2];
sx q[2];
rz(-0.70095567) q[2];
sx q[2];
rz(0.70065166) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7733506) q[1];
sx q[1];
rz(-2.1424286) q[1];
sx q[1];
rz(-1.9707457) q[1];
rz(-pi) q[2];
rz(0.90161277) q[3];
sx q[3];
rz(-1.7898149) q[3];
sx q[3];
rz(2.7215093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69409662) q[2];
sx q[2];
rz(-1.2022377) q[2];
sx q[2];
rz(0.48074943) q[2];
rz(2.2065744) q[3];
sx q[3];
rz(-0.98116773) q[3];
sx q[3];
rz(-0.41486129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9347436) q[0];
sx q[0];
rz(-2.594279) q[0];
sx q[0];
rz(-0.68688399) q[0];
rz(2.1737449) q[1];
sx q[1];
rz(-1.6136439) q[1];
sx q[1];
rz(0.17766775) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2813276) q[0];
sx q[0];
rz(-1.2277368) q[0];
sx q[0];
rz(2.706091) q[0];
rz(0.91089852) q[2];
sx q[2];
rz(-1.0086806) q[2];
sx q[2];
rz(-1.0719932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1473917) q[1];
sx q[1];
rz(-1.8462291) q[1];
sx q[1];
rz(-2.0250399) q[1];
x q[2];
rz(-1.4804041) q[3];
sx q[3];
rz(-1.5022583) q[3];
sx q[3];
rz(-1.7523505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.054242) q[2];
sx q[2];
rz(-0.44173104) q[2];
sx q[2];
rz(-0.98197118) q[2];
rz(-0.25518498) q[3];
sx q[3];
rz(-1.988215) q[3];
sx q[3];
rz(1.0758146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690689) q[0];
sx q[0];
rz(-0.60966063) q[0];
sx q[0];
rz(-0.12741086) q[0];
rz(-1.472507) q[1];
sx q[1];
rz(-1.1839048) q[1];
sx q[1];
rz(1.4471819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873647) q[0];
sx q[0];
rz(-1.7561098) q[0];
sx q[0];
rz(0.32295708) q[0];
rz(-0.54355343) q[2];
sx q[2];
rz(-1.8951891) q[2];
sx q[2];
rz(-0.75553644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.081320612) q[1];
sx q[1];
rz(-1.6439207) q[1];
sx q[1];
rz(-3.0036219) q[1];
rz(-pi) q[2];
rz(0.24932464) q[3];
sx q[3];
rz(-2.1405947) q[3];
sx q[3];
rz(-2.9227481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.24826) q[2];
sx q[2];
rz(-2.1620763) q[2];
sx q[2];
rz(0.40929201) q[2];
rz(2.4457757) q[3];
sx q[3];
rz(-2.4455363) q[3];
sx q[3];
rz(-2.5223562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5291587) q[0];
sx q[0];
rz(-0.22888628) q[0];
sx q[0];
rz(2.9869475) q[0];
rz(1.0908499) q[1];
sx q[1];
rz(-2.2583074) q[1];
sx q[1];
rz(-1.8152016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33783864) q[0];
sx q[0];
rz(-2.3804166) q[0];
sx q[0];
rz(1.9597783) q[0];
rz(0.048166231) q[2];
sx q[2];
rz(-2.0694975) q[2];
sx q[2];
rz(-1.1737385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0834956) q[1];
sx q[1];
rz(-0.82193437) q[1];
sx q[1];
rz(3.0739944) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6821413) q[3];
sx q[3];
rz(-1.3927476) q[3];
sx q[3];
rz(2.0083754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0632625) q[2];
sx q[2];
rz(-2.7276954) q[2];
sx q[2];
rz(2.4388893) q[2];
rz(-0.094376877) q[3];
sx q[3];
rz(-2.1636212) q[3];
sx q[3];
rz(2.2167061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5265882) q[0];
sx q[0];
rz(-2.1115392) q[0];
sx q[0];
rz(0.35453844) q[0];
rz(-1.4226557) q[1];
sx q[1];
rz(-1.4965897) q[1];
sx q[1];
rz(-2.6197701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1352393) q[0];
sx q[0];
rz(-1.420605) q[0];
sx q[0];
rz(1.2144258) q[0];
rz(-pi) q[1];
rz(-3.0031239) q[2];
sx q[2];
rz(-1.8731786) q[2];
sx q[2];
rz(-0.87463986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2592755) q[1];
sx q[1];
rz(-2.2188201) q[1];
sx q[1];
rz(-0.8671182) q[1];
rz(-0.75504889) q[3];
sx q[3];
rz(-0.88249373) q[3];
sx q[3];
rz(1.2699708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3598513) q[2];
sx q[2];
rz(-1.968911) q[2];
sx q[2];
rz(-0.68217984) q[2];
rz(1.3683246) q[3];
sx q[3];
rz(-1.5409639) q[3];
sx q[3];
rz(-1.0002331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2405887) q[0];
sx q[0];
rz(-2.1841342) q[0];
sx q[0];
rz(2.8484455) q[0];
rz(2.3121569) q[1];
sx q[1];
rz(-1.71143) q[1];
sx q[1];
rz(-1.5477187) q[1];
rz(-1.310036) q[2];
sx q[2];
rz(-1.3481067) q[2];
sx q[2];
rz(1.3630661) q[2];
rz(-0.33942038) q[3];
sx q[3];
rz(-3.0302553) q[3];
sx q[3];
rz(-0.29449022) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
