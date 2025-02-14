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
rz(1.7477859) q[0];
sx q[0];
rz(-2.2414247) q[0];
sx q[0];
rz(1.409344) q[0];
rz(-2.4318168) q[1];
sx q[1];
rz(-0.33325279) q[1];
sx q[1];
rz(-0.38212734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2173806) q[0];
sx q[0];
rz(-1.7048302) q[0];
sx q[0];
rz(2.4045375) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0067031) q[2];
sx q[2];
rz(-1.8230048) q[2];
sx q[2];
rz(-1.2812966) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5772878) q[1];
sx q[1];
rz(-1.9365942) q[1];
sx q[1];
rz(-1.9655966) q[1];
x q[2];
rz(-3.0218533) q[3];
sx q[3];
rz(-2.2594707) q[3];
sx q[3];
rz(-2.7514091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0673151) q[2];
sx q[2];
rz(-1.6396319) q[2];
sx q[2];
rz(0.074946694) q[2];
rz(-1.8986757) q[3];
sx q[3];
rz(-2.8800745) q[3];
sx q[3];
rz(0.3956795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0089834) q[0];
sx q[0];
rz(-2.7217396) q[0];
sx q[0];
rz(-2.5123151) q[0];
rz(2.3043326) q[1];
sx q[1];
rz(-0.75051296) q[1];
sx q[1];
rz(0.57964051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0934187) q[0];
sx q[0];
rz(-2.3701128) q[0];
sx q[0];
rz(2.8475965) q[0];
x q[1];
rz(-1.8331892) q[2];
sx q[2];
rz(-1.1875404) q[2];
sx q[2];
rz(0.33363909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9847199) q[1];
sx q[1];
rz(-0.52238363) q[1];
sx q[1];
rz(0.70981437) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5509255) q[3];
sx q[3];
rz(-1.7345034) q[3];
sx q[3];
rz(-0.57396561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2418182) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(-2.3972798) q[2];
rz(-0.36738473) q[3];
sx q[3];
rz(-1.8495879) q[3];
sx q[3];
rz(-1.415409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431417) q[0];
sx q[0];
rz(-0.95040584) q[0];
sx q[0];
rz(2.1344192) q[0];
rz(3.1305283) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(-0.99753582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580306) q[0];
sx q[0];
rz(-3.0016388) q[0];
sx q[0];
rz(1.3849694) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80834016) q[2];
sx q[2];
rz(-1.488121) q[2];
sx q[2];
rz(0.17553628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1147194) q[1];
sx q[1];
rz(-0.51843868) q[1];
sx q[1];
rz(0.015109574) q[1];
rz(-pi) q[2];
rz(0.51181958) q[3];
sx q[3];
rz(-1.5580172) q[3];
sx q[3];
rz(0.46841533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.02115383) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(0.8417449) q[2];
rz(-0.65819955) q[3];
sx q[3];
rz(-0.87116146) q[3];
sx q[3];
rz(3.1200718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.72092527) q[0];
sx q[0];
rz(-0.1249211) q[0];
sx q[0];
rz(-0.7290054) q[0];
rz(-2.3782102) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(0.36300945) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.027452) q[0];
sx q[0];
rz(-1.318256) q[0];
sx q[0];
rz(2.8414498) q[0];
x q[1];
rz(-0.82421397) q[2];
sx q[2];
rz(-1.1477648) q[2];
sx q[2];
rz(-1.3049478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22717366) q[1];
sx q[1];
rz(-2.8110831) q[1];
sx q[1];
rz(0.4497437) q[1];
x q[2];
rz(-1.5176124) q[3];
sx q[3];
rz(-1.6973572) q[3];
sx q[3];
rz(0.80985957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9558941) q[2];
sx q[2];
rz(-2.318435) q[2];
sx q[2];
rz(-1.4996747) q[2];
rz(-2.5540292) q[3];
sx q[3];
rz(-0.9891808) q[3];
sx q[3];
rz(-0.51830083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8374306) q[0];
sx q[0];
rz(-2.4346133) q[0];
sx q[0];
rz(-0.28513232) q[0];
rz(2.888491) q[1];
sx q[1];
rz(-1.0292091) q[1];
sx q[1];
rz(2.0957799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6308002) q[0];
sx q[0];
rz(-1.0655154) q[0];
sx q[0];
rz(2.0204161) q[0];
rz(2.3836993) q[2];
sx q[2];
rz(-0.57595384) q[2];
sx q[2];
rz(-0.26612258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1872028) q[1];
sx q[1];
rz(-2.7437401) q[1];
sx q[1];
rz(0.14601645) q[1];
x q[2];
rz(-0.24744125) q[3];
sx q[3];
rz(-2.4486922) q[3];
sx q[3];
rz(-3.122356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73695856) q[2];
sx q[2];
rz(-1.0598695) q[2];
sx q[2];
rz(-0.57445478) q[2];
rz(-0.61070853) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059998) q[0];
sx q[0];
rz(-1.3845504) q[0];
sx q[0];
rz(2.7912676) q[0];
rz(-0.2925182) q[1];
sx q[1];
rz(-0.12648335) q[1];
sx q[1];
rz(-0.77936053) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.821601) q[0];
sx q[0];
rz(-1.6896833) q[0];
sx q[0];
rz(3.0827808) q[0];
rz(-1.3931403) q[2];
sx q[2];
rz(-1.0972995) q[2];
sx q[2];
rz(2.908978) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8694589) q[1];
sx q[1];
rz(-1.773352) q[1];
sx q[1];
rz(2.814568) q[1];
rz(-pi) q[2];
rz(-0.27542122) q[3];
sx q[3];
rz(-0.53715992) q[3];
sx q[3];
rz(0.28322434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.902035) q[2];
sx q[2];
rz(-1.5152479) q[2];
sx q[2];
rz(0.96417344) q[2];
rz(0.17298175) q[3];
sx q[3];
rz(-2.1292584) q[3];
sx q[3];
rz(3.0103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5445589) q[0];
sx q[0];
rz(-1.9439161) q[0];
sx q[0];
rz(0.069393754) q[0];
rz(1.767905) q[1];
sx q[1];
rz(-1.8927788) q[1];
sx q[1];
rz(2.6232041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98250472) q[0];
sx q[0];
rz(-2.2347274) q[0];
sx q[0];
rz(0.21592617) q[0];
rz(-pi) q[1];
rz(2.0608299) q[2];
sx q[2];
rz(-1.9980717) q[2];
sx q[2];
rz(0.13908041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6024149) q[1];
sx q[1];
rz(-1.5966192) q[1];
sx q[1];
rz(0.13823814) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9767738) q[3];
sx q[3];
rz(-1.6281152) q[3];
sx q[3];
rz(0.27637339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18846506) q[2];
sx q[2];
rz(-1.943482) q[2];
sx q[2];
rz(-2.8971064) q[2];
rz(1.9237579) q[3];
sx q[3];
rz(-3.0471424) q[3];
sx q[3];
rz(-0.34734669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090488) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(0.69945139) q[0];
rz(-2.4199016) q[1];
sx q[1];
rz(-1.7803918) q[1];
sx q[1];
rz(-1.1915421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31909865) q[0];
sx q[0];
rz(-2.8945277) q[0];
sx q[0];
rz(-2.6074227) q[0];
rz(1.4100685) q[2];
sx q[2];
rz(-1.3597915) q[2];
sx q[2];
rz(-2.4461022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6859795) q[1];
sx q[1];
rz(-2.1081104) q[1];
sx q[1];
rz(-3.0431391) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6939387) q[3];
sx q[3];
rz(-1.8869683) q[3];
sx q[3];
rz(2.3929736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2193853) q[2];
sx q[2];
rz(-0.81874138) q[2];
sx q[2];
rz(1.7259664) q[2];
rz(-2.7042232) q[3];
sx q[3];
rz(-2.8919817) q[3];
sx q[3];
rz(-2.6387446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9472083) q[0];
sx q[0];
rz(-1.9358862) q[0];
sx q[0];
rz(-2.6956287) q[0];
rz(1.3392316) q[1];
sx q[1];
rz(-0.65473348) q[1];
sx q[1];
rz(1.9816678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9438749) q[0];
sx q[0];
rz(-2.3811584) q[0];
sx q[0];
rz(2.6307153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8322938) q[2];
sx q[2];
rz(-1.6435677) q[2];
sx q[2];
rz(2.6852599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36148188) q[1];
sx q[1];
rz(-2.105753) q[1];
sx q[1];
rz(-1.8423716) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2307617) q[3];
sx q[3];
rz(-1.0928705) q[3];
sx q[3];
rz(-2.7895989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7356073) q[2];
sx q[2];
rz(-1.0286101) q[2];
sx q[2];
rz(-2.410991) q[2];
rz(-2.2405911) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(-0.034339529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63448298) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(0.46257567) q[0];
rz(-2.0787461) q[1];
sx q[1];
rz(-1.3674419) q[1];
sx q[1];
rz(-0.06632334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49100599) q[0];
sx q[0];
rz(-0.9977451) q[0];
sx q[0];
rz(-1.2038403) q[0];
rz(-1.2356051) q[2];
sx q[2];
rz(-1.7370213) q[2];
sx q[2];
rz(-0.34632296) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6219225) q[1];
sx q[1];
rz(-1.038045) q[1];
sx q[1];
rz(-0.68343917) q[1];
rz(-pi) q[2];
rz(2.9154882) q[3];
sx q[3];
rz(-1.9852828) q[3];
sx q[3];
rz(-1.7630486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.037420951) q[2];
sx q[2];
rz(-0.53356844) q[2];
sx q[2];
rz(-1.9083692) q[2];
rz(-2.6836266) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(2.7428194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11252277) q[0];
sx q[0];
rz(-1.7081013) q[0];
sx q[0];
rz(-1.3992455) q[0];
rz(2.9872672) q[1];
sx q[1];
rz(-1.4230774) q[1];
sx q[1];
rz(1.9565061) q[1];
rz(-0.90171705) q[2];
sx q[2];
rz(-2.6606885) q[2];
sx q[2];
rz(-0.056405141) q[2];
rz(2.1460228) q[3];
sx q[3];
rz(-1.155326) q[3];
sx q[3];
rz(-0.028459878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
