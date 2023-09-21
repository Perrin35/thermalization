OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(-1.0441138) q[1];
sx q[1];
rz(1.27966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5582433) q[0];
sx q[0];
rz(-1.3238751) q[0];
sx q[0];
rz(0.14351828) q[0];
rz(-2.6718164) q[2];
sx q[2];
rz(-2.6758286) q[2];
sx q[2];
rz(-2.2848406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76525926) q[1];
sx q[1];
rz(-1.1632803) q[1];
sx q[1];
rz(-2.5024662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7198556) q[3];
sx q[3];
rz(-0.30656439) q[3];
sx q[3];
rz(-0.43678624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(0.65650666) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1332557) q[0];
sx q[0];
rz(-0.79611859) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(2.6541236) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66747626) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(1.4897896) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0080645) q[2];
sx q[2];
rz(-1.184706) q[2];
sx q[2];
rz(1.3646477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65709719) q[1];
sx q[1];
rz(-0.24689455) q[1];
sx q[1];
rz(0.8685649) q[1];
x q[2];
rz(-1.8295733) q[3];
sx q[3];
rz(-2.8142455) q[3];
sx q[3];
rz(0.60206383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-1.0428492) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(-0.22274676) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6762786) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(0.42155427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1026956) q[2];
sx q[2];
rz(-0.42429081) q[2];
sx q[2];
rz(1.588856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1258771) q[1];
sx q[1];
rz(-1.6385957) q[1];
sx q[1];
rz(2.5281639) q[1];
x q[2];
rz(-1.2336897) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(-3.0984578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.2732482) q[2];
sx q[2];
rz(-1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(0.96364337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(-0.9202252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5266787) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(-0.99292361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.780519) q[2];
sx q[2];
rz(-1.4898584) q[2];
sx q[2];
rz(-2.0780448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(2.4469417) q[1];
x q[2];
rz(-1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(0.40311381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(-2.3222893) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.3764494) q[0];
sx q[0];
rz(-1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-2.4914425) q[1];
sx q[1];
rz(-2.6884902) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74005175) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(3.1042276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3327417) q[2];
sx q[2];
rz(-0.14698725) q[2];
sx q[2];
rz(0.57839314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15935005) q[1];
sx q[1];
rz(-0.23540007) q[1];
sx q[1];
rz(-1.2118641) q[1];
x q[2];
rz(2.681053) q[3];
sx q[3];
rz(-1.0214877) q[3];
sx q[3];
rz(-0.82563803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(1.5117234) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081651) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(0.23813716) q[0];
rz(0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.5135117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7662402) q[0];
sx q[0];
rz(-2.0413114) q[0];
sx q[0];
rz(-1.9348295) q[0];
x q[1];
rz(-3.0797708) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(-2.2668554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39735079) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(2.1703297) q[1];
x q[2];
rz(2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(0.57487956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(-1.1435821) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(-2.9706764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(1.827084) q[1];
sx q[1];
rz(-2.59771) q[1];
sx q[1];
rz(1.1869173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-2.5573686) q[0];
sx q[0];
rz(-2.0969735) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6830964) q[2];
sx q[2];
rz(-3.0756604) q[2];
sx q[2];
rz(-2.3022848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0901511) q[1];
sx q[1];
rz(-1.4370343) q[1];
sx q[1];
rz(2.306288) q[1];
rz(-pi) q[2];
rz(-1.8625453) q[3];
sx q[3];
rz(-0.96127931) q[3];
sx q[3];
rz(2.9677344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(1.0303248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6894585) q[0];
sx q[0];
rz(-1.2297213) q[0];
sx q[0];
rz(0.80490168) q[0];
x q[1];
rz(1.6542997) q[2];
sx q[2];
rz(-0.52831542) q[2];
sx q[2];
rz(-0.63937843) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.8712908) q[1];
x q[2];
rz(2.6183073) q[3];
sx q[3];
rz(-1.0417632) q[3];
sx q[3];
rz(-2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29256233) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(-2.0896185) q[1];
sx q[1];
rz(-2.4119222) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90056706) q[0];
sx q[0];
rz(-1.7591488) q[0];
sx q[0];
rz(-2.7450949) q[0];
rz(0.029149292) q[2];
sx q[2];
rz(-1.2505194) q[2];
sx q[2];
rz(2.0707891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4371722) q[1];
sx q[1];
rz(-1.51209) q[1];
sx q[1];
rz(1.0874332) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.023852392) q[3];
sx q[3];
rz(-2.1937222) q[3];
sx q[3];
rz(-0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51745522) q[2];
sx q[2];
rz(-2.3949261) q[2];
sx q[2];
rz(-1.5967782) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88403945) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(0.89417017) q[0];
x q[1];
rz(-1.7061383) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(-2.4382255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0386104) q[1];
sx q[1];
rz(-0.83676941) q[1];
sx q[1];
rz(-3.1208913) q[1];
rz(-pi) q[2];
rz(0.65532834) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(-2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(3.0604559) q[2];
rz(-1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1186196) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.8267869) q[1];
sx q[1];
rz(-2.343315) q[1];
sx q[1];
rz(1.5409484) q[1];
rz(1.5224456) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
rz(0.046394596) q[3];
sx q[3];
rz(-1.8002602) q[3];
sx q[3];
rz(0.044063448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
