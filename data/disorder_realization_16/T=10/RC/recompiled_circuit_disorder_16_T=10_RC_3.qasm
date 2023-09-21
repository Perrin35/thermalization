OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(4.1339388) q[1];
sx q[1];
rz(9.0864656) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4727288) q[0];
sx q[0];
rz(-2.7358486) q[0];
sx q[0];
rz(-0.91234447) q[0];
rz(-0.89262427) q[2];
sx q[2];
rz(-1.2783588) q[2];
sx q[2];
rz(-2.6543648) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6602064) q[1];
sx q[1];
rz(-1.8567137) q[1];
sx q[1];
rz(1.5155161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5987414) q[3];
sx q[3];
rz(-0.51252796) q[3];
sx q[3];
rz(-0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9989495) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.1738698) q[2];
rz(-0.075803444) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-0.092806667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8006111) q[0];
sx q[0];
rz(-2.0759463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(-2.5669572) q[1];
sx q[1];
rz(-2.7119633) q[1];
sx q[1];
rz(-1.2423135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949462) q[0];
sx q[0];
rz(-1.6992237) q[0];
sx q[0];
rz(2.8917679) q[0];
rz(-0.93810268) q[2];
sx q[2];
rz(-1.4422851) q[2];
sx q[2];
rz(0.036966952) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17774432) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(0.53667712) q[1];
rz(-pi) q[2];
rz(1.0817238) q[3];
sx q[3];
rz(-0.92182577) q[3];
sx q[3];
rz(2.0826516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3339281) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(-0.58369613) q[2];
rz(-0.57404533) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(-0.13124245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211001) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-0.39847001) q[1];
sx q[1];
rz(-1.0167936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4749878) q[0];
sx q[0];
rz(-3.0525065) q[0];
sx q[0];
rz(0.41390093) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3967907) q[2];
sx q[2];
rz(-0.6664657) q[2];
sx q[2];
rz(-1.9772066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6779855) q[1];
sx q[1];
rz(-1.7079759) q[1];
sx q[1];
rz(-0.82171085) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0924762) q[3];
sx q[3];
rz(-1.8521063) q[3];
sx q[3];
rz(1.0055055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8016522) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-1.0495079) q[2];
rz(-2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.1857741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.6695492) q[0];
rz(2.4064348) q[1];
sx q[1];
rz(-2.3627294) q[1];
sx q[1];
rz(-0.24681117) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128162) q[0];
sx q[0];
rz(-1.719559) q[0];
sx q[0];
rz(2.2674198) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3601801) q[2];
sx q[2];
rz(-1.1014551) q[2];
sx q[2];
rz(2.5460555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3226763) q[1];
sx q[1];
rz(-2.7671742) q[1];
sx q[1];
rz(0.14426343) q[1];
rz(-pi) q[2];
x q[2];
rz(2.187192) q[3];
sx q[3];
rz(-1.8064926) q[3];
sx q[3];
rz(2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-2.0229979) q[2];
sx q[2];
rz(1.5412615) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-1.0975081) q[3];
sx q[3];
rz(0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8108869) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.5785626) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(-0.24838233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9527234) q[0];
sx q[0];
rz(-1.6048604) q[0];
sx q[0];
rz(-0.54038318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9551198) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(-1.1764256) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0866962) q[1];
sx q[1];
rz(-0.61453648) q[1];
sx q[1];
rz(-0.30026786) q[1];
x q[2];
rz(-0.6487209) q[3];
sx q[3];
rz(-1.7626764) q[3];
sx q[3];
rz(0.81024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(0.79745897) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1335063) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(1.7957934) q[0];
rz(1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(-3.016901) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9799177) q[0];
sx q[0];
rz(-1.7470164) q[0];
sx q[0];
rz(1.4833223) q[0];
rz(-pi) q[1];
rz(-0.40839809) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(2.1085395) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5387419) q[1];
sx q[1];
rz(-2.3067143) q[1];
sx q[1];
rz(1.49453) q[1];
rz(-pi) q[2];
rz(1.3521306) q[3];
sx q[3];
rz(-1.8018186) q[3];
sx q[3];
rz(-0.71393379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6283915) q[2];
sx q[2];
rz(-1.1867563) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597647) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(2.0986309) q[0];
rz(-0.46328059) q[1];
sx q[1];
rz(-1.1136585) q[1];
sx q[1];
rz(1.0707062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844855) q[0];
sx q[0];
rz(-2.6575408) q[0];
sx q[0];
rz(-0.0012782106) q[0];
x q[1];
rz(-1.7037017) q[2];
sx q[2];
rz(-1.226236) q[2];
sx q[2];
rz(2.8466356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5412022) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(-0.77460918) q[1];
rz(2.8093852) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36859194) q[2];
sx q[2];
rz(-2.4020782) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(2.0543082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535646) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.2063684) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-0.85314631) q[1];
sx q[1];
rz(-0.94747296) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065814106) q[0];
sx q[0];
rz(-1.4697207) q[0];
sx q[0];
rz(2.1094735) q[0];
rz(-pi) q[1];
rz(-0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(-1.8159602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87654963) q[1];
sx q[1];
rz(-1.8076841) q[1];
sx q[1];
rz(2.151728) q[1];
rz(1.2097589) q[3];
sx q[3];
rz(-2.3449538) q[3];
sx q[3];
rz(-2.3494997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(-0.46978152) q[2];
rz(1.3011159) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.3289733) q[0];
rz(-0.7912311) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(-0.20283094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79265362) q[0];
sx q[0];
rz(-1.5968423) q[0];
sx q[0];
rz(-3.1025725) q[0];
x q[1];
rz(-0.75195306) q[2];
sx q[2];
rz(-2.8402036) q[2];
sx q[2];
rz(1.6981268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5962332) q[1];
sx q[1];
rz(-1.6714449) q[1];
sx q[1];
rz(-0.049493162) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4167452) q[3];
sx q[3];
rz(-2.027958) q[3];
sx q[3];
rz(-2.9726213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41708502) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(0.99651304) q[2];
rz(-0.35342446) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.080169454) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(-2.9598575) q[0];
rz(3.0985447) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(2.8607686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1998394) q[0];
sx q[0];
rz(-2.0240677) q[0];
sx q[0];
rz(2.0487294) q[0];
rz(-pi) q[1];
rz(-2.4091987) q[2];
sx q[2];
rz(-2.8576982) q[2];
sx q[2];
rz(2.4925799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3481969) q[1];
sx q[1];
rz(-1.6284202) q[1];
sx q[1];
rz(0.036199526) q[1];
x q[2];
rz(0.19631581) q[3];
sx q[3];
rz(-1.2442949) q[3];
sx q[3];
rz(1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3165555) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(-1.0021707) q[3];
sx q[3];
rz(-1.802417) q[3];
sx q[3];
rz(-0.56308693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4939209) q[0];
sx q[0];
rz(-1.5734084) q[0];
sx q[0];
rz(-1.5403803) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.95332425) q[2];
sx q[2];
rz(-0.8669903) q[2];
sx q[2];
rz(0.16624761) q[2];
rz(-2.279083) q[3];
sx q[3];
rz(-2.6211092) q[3];
sx q[3];
rz(0.72343788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];