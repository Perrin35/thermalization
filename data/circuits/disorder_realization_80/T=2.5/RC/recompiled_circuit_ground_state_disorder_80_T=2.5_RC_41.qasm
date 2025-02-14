OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(2.2194982) q[0];
sx q[0];
rz(10.194869) q[0];
rz(-2.4251179) q[1];
sx q[1];
rz(-0.94243503) q[1];
sx q[1];
rz(2.4997349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72015136) q[0];
sx q[0];
rz(-2.9099265) q[0];
sx q[0];
rz(-2.6257444) q[0];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5696758) q[2];
sx q[2];
rz(-1.6482774) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7388307) q[1];
sx q[1];
rz(-1.9301199) q[1];
sx q[1];
rz(-1.7372543) q[1];
x q[2];
rz(-0.83323102) q[3];
sx q[3];
rz(-1.1166443) q[3];
sx q[3];
rz(-1.3564701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98051071) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(-0.80992997) q[2];
rz(0.03446456) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(0.1035498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.506839) q[0];
sx q[0];
rz(-2.6320808) q[0];
sx q[0];
rz(0.65938812) q[0];
rz(-1.4393282) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(-2.6606681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5372807) q[0];
sx q[0];
rz(-1.3635646) q[0];
sx q[0];
rz(-2.0181433) q[0];
x q[1];
rz(-2.4056817) q[2];
sx q[2];
rz(-0.99756587) q[2];
sx q[2];
rz(2.5664751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2175423) q[1];
sx q[1];
rz(-0.21044193) q[1];
sx q[1];
rz(1.7173052) q[1];
rz(-pi) q[2];
rz(-2.5629077) q[3];
sx q[3];
rz(-1.1704786) q[3];
sx q[3];
rz(2.2861885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7646358) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(-0.62977201) q[2];
rz(-1.9624286) q[3];
sx q[3];
rz(-0.7032913) q[3];
sx q[3];
rz(-2.0298957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-0.010882219) q[0];
sx q[0];
rz(-0.96845281) q[0];
rz(3.0015266) q[1];
sx q[1];
rz(-1.3532956) q[1];
sx q[1];
rz(2.5586939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4980782) q[0];
sx q[0];
rz(-1.5257784) q[0];
sx q[0];
rz(2.974653) q[0];
rz(-pi) q[1];
rz(-1.8331362) q[2];
sx q[2];
rz(-0.749974) q[2];
sx q[2];
rz(-1.7469847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1000468) q[1];
sx q[1];
rz(-2.1963781) q[1];
sx q[1];
rz(2.2872674) q[1];
rz(-2.8606877) q[3];
sx q[3];
rz(-1.2634522) q[3];
sx q[3];
rz(-0.25562778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6110903) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(-0.21406847) q[2];
rz(2.8392082) q[3];
sx q[3];
rz(-2.3979135) q[3];
sx q[3];
rz(2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18836235) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(-2.0971712) q[0];
rz(-2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(0.16876076) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89655748) q[0];
sx q[0];
rz(-1.3447133) q[0];
sx q[0];
rz(-2.079179) q[0];
x q[1];
rz(0.20511583) q[2];
sx q[2];
rz(-1.2502708) q[2];
sx q[2];
rz(-1.4687302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16774878) q[1];
sx q[1];
rz(-2.1653321) q[1];
sx q[1];
rz(-1.3092625) q[1];
x q[2];
rz(1.6200284) q[3];
sx q[3];
rz(-2.4353046) q[3];
sx q[3];
rz(-0.8362706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7793444) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(-2.0812422) q[2];
rz(0.29252163) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5523858) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(-2.2733083) q[0];
rz(-0.6262511) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.7832322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5007965) q[0];
sx q[0];
rz(-1.2828886) q[0];
sx q[0];
rz(-1.4782216) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38154885) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(-0.82009456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4386126) q[1];
sx q[1];
rz(-1.3452521) q[1];
sx q[1];
rz(1.799052) q[1];
rz(-pi) q[2];
rz(2.503848) q[3];
sx q[3];
rz(-1.9336091) q[3];
sx q[3];
rz(-2.2408298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2061578) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(2.6182776) q[2];
rz(2.7715136) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(1.7962598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0321781) q[0];
sx q[0];
rz(-0.21571708) q[0];
sx q[0];
rz(2.7874462) q[0];
rz(-0.94193637) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(1.8249493) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7553332) q[0];
sx q[0];
rz(-1.0585896) q[0];
sx q[0];
rz(1.451158) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5923827) q[2];
sx q[2];
rz(-1.5592279) q[2];
sx q[2];
rz(-2.9988097) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5608974) q[1];
sx q[1];
rz(-1.3968916) q[1];
sx q[1];
rz(1.5087125) q[1];
rz(-pi) q[2];
rz(2.7735071) q[3];
sx q[3];
rz(-2.4637665) q[3];
sx q[3];
rz(-2.8742636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3256623) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-2.8032934) q[2];
rz(0.47652388) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(2.8009955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7020096) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(-1.4078377) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(0.62526155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940373) q[0];
sx q[0];
rz(-2.4522344) q[0];
sx q[0];
rz(2.5563142) q[0];
x q[1];
rz(-0.53886063) q[2];
sx q[2];
rz(-1.6340573) q[2];
sx q[2];
rz(-0.45743313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75692528) q[1];
sx q[1];
rz(-0.97434348) q[1];
sx q[1];
rz(-1.8310962) q[1];
rz(2.7688249) q[3];
sx q[3];
rz(-0.81171747) q[3];
sx q[3];
rz(1.852664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(-1.3803049) q[2];
rz(2.7050833) q[3];
sx q[3];
rz(-1.0218388) q[3];
sx q[3];
rz(0.71353394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9647144) q[0];
sx q[0];
rz(-0.62129337) q[0];
sx q[0];
rz(2.6253413) q[0];
rz(-0.77975726) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-1.0468743) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488432) q[0];
sx q[0];
rz(-1.7147281) q[0];
sx q[0];
rz(-2.8400665) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5604232) q[2];
sx q[2];
rz(-2.0072674) q[2];
sx q[2];
rz(1.3429067) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88491625) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(-1.0156471) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1740999) q[3];
sx q[3];
rz(-1.124212) q[3];
sx q[3];
rz(0.74029884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20389916) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(-2.1843809) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-2.011994) q[3];
sx q[3];
rz(-2.8637776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99943632) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(0.18705046) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(1.4923219) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87875116) q[0];
sx q[0];
rz(-0.10058144) q[0];
sx q[0];
rz(-0.79128964) q[0];
x q[1];
rz(0.19047462) q[2];
sx q[2];
rz(-1.8996793) q[2];
sx q[2];
rz(-2.637907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8111251) q[1];
sx q[1];
rz(-0.44422517) q[1];
sx q[1];
rz(-1.1796239) q[1];
rz(-0.093676693) q[3];
sx q[3];
rz(-1.4149349) q[3];
sx q[3];
rz(-0.81934281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6734068) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(-0.51472384) q[3];
sx q[3];
rz(-2.616021) q[3];
sx q[3];
rz(-1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.89727) q[0];
sx q[0];
rz(-1.6615302) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(1.9482535) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(1.4512482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.25989) q[0];
sx q[0];
rz(-2.8326747) q[0];
sx q[0];
rz(-2.3007459) q[0];
x q[1];
rz(-1.7733943) q[2];
sx q[2];
rz(-2.9349064) q[2];
sx q[2];
rz(-2.8967173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2047233) q[1];
sx q[1];
rz(-1.544049) q[1];
sx q[1];
rz(0.95519798) q[1];
rz(-1.162446) q[3];
sx q[3];
rz(-1.8682042) q[3];
sx q[3];
rz(1.4434467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(-0.80121458) q[2];
rz(-2.2090705) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(0.41509375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5781317) q[0];
sx q[0];
rz(-1.3787855) q[0];
sx q[0];
rz(-0.61510573) q[0];
rz(-2.8201132) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(-2.9762815) q[2];
sx q[2];
rz(-1.4480235) q[2];
sx q[2];
rz(1.5245246) q[2];
rz(2.2533909) q[3];
sx q[3];
rz(-2.785688) q[3];
sx q[3];
rz(1.9733081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
