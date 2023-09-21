OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(5.5881349) q[1];
sx q[1];
rz(10.499428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63700919) q[0];
sx q[0];
rz(-1.5890117) q[0];
sx q[0];
rz(-1.5760742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4104112) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(-0.10057848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97779146) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(1.4896643) q[1];
x q[2];
rz(-0.53332897) q[3];
sx q[3];
rz(-1.9682069) q[3];
sx q[3];
rz(-1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.507507) q[0];
sx q[0];
rz(-1.5956143) q[0];
sx q[0];
rz(-0.12269845) q[0];
rz(-pi) q[1];
rz(-2.03962) q[2];
sx q[2];
rz(-1.8509682) q[2];
sx q[2];
rz(-2.3902992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4109107) q[1];
sx q[1];
rz(-0.69697471) q[1];
sx q[1];
rz(2.8721786) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1355751) q[3];
sx q[3];
rz(-1.9618417) q[3];
sx q[3];
rz(-1.0139549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30341879) q[0];
sx q[0];
rz(-1.6359328) q[0];
sx q[0];
rz(-1.4831878) q[0];
rz(-pi) q[1];
rz(2.1554699) q[2];
sx q[2];
rz(-1.976527) q[2];
sx q[2];
rz(3.0509146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9864195) q[1];
sx q[1];
rz(-1.6134312) q[1];
sx q[1];
rz(0.30100545) q[1];
x q[2];
rz(-0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(0.7437931) q[2];
rz(2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.1490885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3455032) q[0];
sx q[0];
rz(-0.26608276) q[0];
sx q[0];
rz(-0.43568133) q[0];
x q[1];
rz(0.97082246) q[2];
sx q[2];
rz(-2.4568825) q[2];
sx q[2];
rz(-0.1240571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.035605343) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(-1.1483907) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74242384) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(-2.1267668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(2.2573684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0852172) q[0];
sx q[0];
rz(-1.3803122) q[0];
sx q[0];
rz(-2.9604244) q[0];
rz(0.40712936) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(1.7083573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7985178) q[1];
sx q[1];
rz(-1.6371562) q[1];
sx q[1];
rz(0.14454486) q[1];
rz(-pi) q[2];
rz(0.36966427) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(-2.8408185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(-0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-0.96907369) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58020619) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(1.6798621) q[0];
rz(-pi) q[1];
rz(-3.063383) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(0.90925928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6094738) q[1];
sx q[1];
rz(-1.9874959) q[1];
sx q[1];
rz(2.6344224) q[1];
rz(-1.045741) q[3];
sx q[3];
rz(-1.9805555) q[3];
sx q[3];
rz(0.72011442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(-2.7039841) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.70599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3754417) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(0.96813162) q[0];
x q[1];
rz(-2.3195763) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(-0.4617304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85305271) q[1];
sx q[1];
rz(-1.9318046) q[1];
sx q[1];
rz(2.8193874) q[1];
rz(-pi) q[2];
rz(2.7123659) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(-2.6847117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4009565) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-2.9675972) q[2];
rz(-2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(2.984916) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(-2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93850219) q[0];
sx q[0];
rz(-1.4346968) q[0];
sx q[0];
rz(-0.076119856) q[0];
x q[1];
rz(-1.6985967) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(0.70805659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3268765) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(-0.53175064) q[1];
rz(-pi) q[2];
rz(1.385231) q[3];
sx q[3];
rz(-1.984333) q[3];
sx q[3];
rz(-2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(-0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16238427) q[0];
sx q[0];
rz(-0.93063336) q[0];
sx q[0];
rz(-1.0248653) q[0];
x q[1];
rz(1.5588435) q[2];
sx q[2];
rz(-1.0822312) q[2];
sx q[2];
rz(-1.919463) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3194771) q[1];
sx q[1];
rz(-1.9648106) q[1];
sx q[1];
rz(1.0072717) q[1];
rz(2.4273848) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(-0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.7473934) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(2.6031301) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(-1.753153) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5441226) q[0];
sx q[0];
rz(-0.96203502) q[0];
x q[1];
rz(1.795904) q[2];
sx q[2];
rz(-0.39160608) q[2];
sx q[2];
rz(0.12347808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0109509) q[1];
sx q[1];
rz(-0.37591939) q[1];
sx q[1];
rz(1.4569267) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5358701) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(-1.4981131) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(1.5126363) q[2];
sx q[2];
rz(-1.6806316) q[2];
sx q[2];
rz(-3.1113839) q[2];
rz(-0.31717832) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
