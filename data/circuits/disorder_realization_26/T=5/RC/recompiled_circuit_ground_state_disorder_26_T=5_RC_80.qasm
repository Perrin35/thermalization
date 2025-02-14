OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.77312624) q[0];
sx q[0];
rz(2.3946895) q[0];
sx q[0];
rz(11.725732) q[0];
rz(0.12159881) q[1];
sx q[1];
rz(-1.2727979) q[1];
sx q[1];
rz(-2.8425541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75343412) q[0];
sx q[0];
rz(-0.82601604) q[0];
sx q[0];
rz(-2.1623934) q[0];
rz(-pi) q[1];
rz(3.063721) q[2];
sx q[2];
rz(-1.0857333) q[2];
sx q[2];
rz(-2.9295706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2192397) q[1];
sx q[1];
rz(-1.427622) q[1];
sx q[1];
rz(0.24210614) q[1];
x q[2];
rz(2.998803) q[3];
sx q[3];
rz(-2.0215109) q[3];
sx q[3];
rz(0.55270178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3806939) q[2];
sx q[2];
rz(-1.0502522) q[2];
sx q[2];
rz(2.8139581) q[2];
rz(-1.7662175) q[3];
sx q[3];
rz(-1.4204493) q[3];
sx q[3];
rz(-0.49427858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7647917) q[0];
sx q[0];
rz(-1.6015653) q[0];
sx q[0];
rz(2.2862527) q[0];
rz(0.0080464706) q[1];
sx q[1];
rz(-1.2866373) q[1];
sx q[1];
rz(-2.6904552) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52569235) q[0];
sx q[0];
rz(-0.85291686) q[0];
sx q[0];
rz(-2.9327716) q[0];
x q[1];
rz(0.19719736) q[2];
sx q[2];
rz(-2.3596546) q[2];
sx q[2];
rz(0.3202084) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0135368) q[1];
sx q[1];
rz(-1.744413) q[1];
sx q[1];
rz(0.44349576) q[1];
rz(-1.1633918) q[3];
sx q[3];
rz(-2.0137069) q[3];
sx q[3];
rz(0.16678424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1796639) q[2];
sx q[2];
rz(-2.0018061) q[2];
sx q[2];
rz(0.12953225) q[2];
rz(0.1772964) q[3];
sx q[3];
rz(-0.57759053) q[3];
sx q[3];
rz(0.052848335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7496846) q[0];
sx q[0];
rz(-0.41202298) q[0];
sx q[0];
rz(-0.55927292) q[0];
rz(-0.094712146) q[1];
sx q[1];
rz(-1.6030703) q[1];
sx q[1];
rz(0.47725484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9608979) q[0];
sx q[0];
rz(-1.3154234) q[0];
sx q[0];
rz(-1.9810505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36156003) q[2];
sx q[2];
rz(-1.7107367) q[2];
sx q[2];
rz(-0.41658336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4871769) q[1];
sx q[1];
rz(-0.80779874) q[1];
sx q[1];
rz(-0.31262763) q[1];
x q[2];
rz(-0.84945143) q[3];
sx q[3];
rz(-1.9799383) q[3];
sx q[3];
rz(2.9144998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39530784) q[2];
sx q[2];
rz(-1.2664653) q[2];
sx q[2];
rz(2.2115808) q[2];
rz(-1.2578472) q[3];
sx q[3];
rz(-1.7141637) q[3];
sx q[3];
rz(0.045684489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607218) q[0];
sx q[0];
rz(-1.5683132) q[0];
sx q[0];
rz(-1.038653) q[0];
rz(-2.5301798) q[1];
sx q[1];
rz(-2.2544506) q[1];
sx q[1];
rz(-2.2183653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0221585) q[0];
sx q[0];
rz(-2.4537477) q[0];
sx q[0];
rz(-1.0115959) q[0];
rz(-pi) q[1];
rz(3.1298248) q[2];
sx q[2];
rz(-1.1530877) q[2];
sx q[2];
rz(-2.2209446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9070396) q[1];
sx q[1];
rz(-0.72826339) q[1];
sx q[1];
rz(1.6978463) q[1];
rz(-1.3096894) q[3];
sx q[3];
rz(-0.93702836) q[3];
sx q[3];
rz(-1.5696309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5103147) q[2];
sx q[2];
rz(-2.6738561) q[2];
sx q[2];
rz(-2.5941217) q[2];
rz(3.0771717) q[3];
sx q[3];
rz(-1.8378601) q[3];
sx q[3];
rz(2.2678383) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15544686) q[0];
sx q[0];
rz(-0.90279818) q[0];
sx q[0];
rz(-1.4601532) q[0];
rz(2.1414781) q[1];
sx q[1];
rz(-0.8668879) q[1];
sx q[1];
rz(0.11988457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45792199) q[0];
sx q[0];
rz(-2.409465) q[0];
sx q[0];
rz(1.9089655) q[0];
rz(-0.91953711) q[2];
sx q[2];
rz(-2.0555758) q[2];
sx q[2];
rz(0.19249053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5226871) q[1];
sx q[1];
rz(-2.9123016) q[1];
sx q[1];
rz(0.35405901) q[1];
rz(-0.67976953) q[3];
sx q[3];
rz(-2.418737) q[3];
sx q[3];
rz(-2.2732609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.039006058) q[2];
sx q[2];
rz(-0.90709364) q[2];
sx q[2];
rz(0.26068035) q[2];
rz(0.87812224) q[3];
sx q[3];
rz(-1.2295281) q[3];
sx q[3];
rz(-2.3584283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6108625) q[0];
sx q[0];
rz(-0.090228883) q[0];
sx q[0];
rz(-2.1395785) q[0];
rz(-0.37725457) q[1];
sx q[1];
rz(-2.1899624) q[1];
sx q[1];
rz(0.9446876) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8588036) q[0];
sx q[0];
rz(-1.7886046) q[0];
sx q[0];
rz(2.5522638) q[0];
rz(-pi) q[1];
rz(2.2425507) q[2];
sx q[2];
rz(-2.0625249) q[2];
sx q[2];
rz(2.0042302) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78996736) q[1];
sx q[1];
rz(-1.891544) q[1];
sx q[1];
rz(-1.7299537) q[1];
rz(-2.112859) q[3];
sx q[3];
rz(-1.6717864) q[3];
sx q[3];
rz(-2.2057057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79823309) q[2];
sx q[2];
rz(-1.6075906) q[2];
sx q[2];
rz(3.0976683) q[2];
rz(-1.6262866) q[3];
sx q[3];
rz(-2.01912) q[3];
sx q[3];
rz(-1.4656434) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2783022) q[0];
sx q[0];
rz(-2.589812) q[0];
sx q[0];
rz(-0.58468753) q[0];
rz(2.0901285) q[1];
sx q[1];
rz(-0.81948558) q[1];
sx q[1];
rz(-0.47031602) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7746587) q[0];
sx q[0];
rz(-2.1485188) q[0];
sx q[0];
rz(1.7455186) q[0];
x q[1];
rz(-1.9740747) q[2];
sx q[2];
rz(-1.9152616) q[2];
sx q[2];
rz(3.0967876) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0329758) q[1];
sx q[1];
rz(-1.931463) q[1];
sx q[1];
rz(2.3920139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5830273) q[3];
sx q[3];
rz(-1.4577971) q[3];
sx q[3];
rz(-1.5025653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26479244) q[2];
sx q[2];
rz(-0.53692836) q[2];
sx q[2];
rz(2.8301767) q[2];
rz(2.9733859) q[3];
sx q[3];
rz(-1.5752537) q[3];
sx q[3];
rz(0.28347191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6300221) q[0];
sx q[0];
rz(-3.1302852) q[0];
sx q[0];
rz(2.2054963) q[0];
rz(-0.49631897) q[1];
sx q[1];
rz(-0.66718188) q[1];
sx q[1];
rz(-0.47028968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0867827) q[0];
sx q[0];
rz(-2.5544689) q[0];
sx q[0];
rz(-0.55956383) q[0];
x q[1];
rz(1.3730714) q[2];
sx q[2];
rz(-1.573296) q[2];
sx q[2];
rz(-0.041415215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.072402231) q[1];
sx q[1];
rz(-1.6399929) q[1];
sx q[1];
rz(-1.2284159) q[1];
rz(1.7467198) q[3];
sx q[3];
rz(-1.298549) q[3];
sx q[3];
rz(1.7737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63142598) q[2];
sx q[2];
rz(-1.3672071) q[2];
sx q[2];
rz(-1.4754971) q[2];
rz(-0.79536074) q[3];
sx q[3];
rz(-0.16203351) q[3];
sx q[3];
rz(1.2696666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0610166) q[0];
sx q[0];
rz(-2.5503655) q[0];
sx q[0];
rz(-0.06037816) q[0];
rz(-0.16054842) q[1];
sx q[1];
rz(-1.5269273) q[1];
sx q[1];
rz(2.1626332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6290603) q[0];
sx q[0];
rz(-2.1950025) q[0];
sx q[0];
rz(-1.4275622) q[0];
x q[1];
rz(0.8884807) q[2];
sx q[2];
rz(-1.7856626) q[2];
sx q[2];
rz(2.9944098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1886602) q[1];
sx q[1];
rz(-2.1702936) q[1];
sx q[1];
rz(-0.99296928) q[1];
rz(1.5185202) q[3];
sx q[3];
rz(-2.5617122) q[3];
sx q[3];
rz(-1.1842022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0742566) q[2];
sx q[2];
rz(-2.8496075) q[2];
sx q[2];
rz(3.0687029) q[2];
rz(-0.59761754) q[3];
sx q[3];
rz(-1.3772734) q[3];
sx q[3];
rz(-3.1387175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.446796) q[0];
sx q[0];
rz(-1.0121166) q[0];
sx q[0];
rz(-0.52892518) q[0];
rz(-2.9534598) q[1];
sx q[1];
rz(-2.4362322) q[1];
sx q[1];
rz(-0.58427748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463294) q[0];
sx q[0];
rz(-1.3130616) q[0];
sx q[0];
rz(-1.3168174) q[0];
x q[1];
rz(2.4227117) q[2];
sx q[2];
rz(-2.5846057) q[2];
sx q[2];
rz(2.5124036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1297785) q[1];
sx q[1];
rz(-0.41167694) q[1];
sx q[1];
rz(-3.0342388) q[1];
rz(-pi) q[2];
rz(-0.27354555) q[3];
sx q[3];
rz(-0.97018948) q[3];
sx q[3];
rz(0.41714868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3820485) q[2];
sx q[2];
rz(-2.1376762) q[2];
sx q[2];
rz(-0.071852597) q[2];
rz(1.0233277) q[3];
sx q[3];
rz(-1.5724678) q[3];
sx q[3];
rz(-2.6527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857984) q[0];
sx q[0];
rz(-2.4837942) q[0];
sx q[0];
rz(2.876045) q[0];
rz(-2.4664948) q[1];
sx q[1];
rz(-1.594512) q[1];
sx q[1];
rz(-0.095269861) q[1];
rz(0.51495348) q[2];
sx q[2];
rz(-1.5025768) q[2];
sx q[2];
rz(1.2160355) q[2];
rz(-0.76473372) q[3];
sx q[3];
rz(-0.22976362) q[3];
sx q[3];
rz(0.094658628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
