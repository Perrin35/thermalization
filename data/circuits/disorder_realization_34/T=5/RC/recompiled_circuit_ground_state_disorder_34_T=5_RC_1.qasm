OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(0.85564268) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(0.56245437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320975) q[0];
sx q[0];
rz(-1.2065071) q[0];
sx q[0];
rz(-3.0845736) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15563528) q[2];
sx q[2];
rz(-2.0726207) q[2];
sx q[2];
rz(1.6574455) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7592647) q[1];
sx q[1];
rz(-0.60324962) q[1];
sx q[1];
rz(-2.3690577) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9770369) q[3];
sx q[3];
rz(-1.9349758) q[3];
sx q[3];
rz(-1.4573569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2354551) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(2.3423024) q[2];
rz(-0.47131395) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(-0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.8164182) q[0];
sx q[0];
rz(1.7934196) q[0];
rz(-1.7680291) q[1];
sx q[1];
rz(-1.9919688) q[1];
sx q[1];
rz(1.1522393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2179759) q[0];
sx q[0];
rz(-2.2647175) q[0];
sx q[0];
rz(-0.38461916) q[0];
rz(-pi) q[1];
x q[1];
rz(2.466408) q[2];
sx q[2];
rz(-2.084888) q[2];
sx q[2];
rz(-2.6509283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79346953) q[1];
sx q[1];
rz(-1.7382227) q[1];
sx q[1];
rz(0.69116418) q[1];
rz(-pi) q[2];
rz(-1.0188607) q[3];
sx q[3];
rz(-2.0235007) q[3];
sx q[3];
rz(-0.34358968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79933244) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(-0.93117923) q[2];
rz(-0.10654199) q[3];
sx q[3];
rz(-2.0276766) q[3];
sx q[3];
rz(0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057864144) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(2.6237543) q[0];
rz(-2.2528516) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(2.6944366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155375) q[0];
sx q[0];
rz(-1.8319905) q[0];
sx q[0];
rz(-1.5047856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4141623) q[2];
sx q[2];
rz(-1.1377678) q[2];
sx q[2];
rz(-0.43011452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4616926) q[1];
sx q[1];
rz(-2.5352806) q[1];
sx q[1];
rz(1.3316505) q[1];
rz(0.52618653) q[3];
sx q[3];
rz(-1.8136214) q[3];
sx q[3];
rz(-0.015004166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7670224) q[2];
sx q[2];
rz(-0.76408237) q[2];
sx q[2];
rz(-2.6089597) q[2];
rz(0.24752188) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240876) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(0.10661539) q[0];
rz(2.7952349) q[1];
sx q[1];
rz(-0.84500161) q[1];
sx q[1];
rz(1.9812298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713523) q[0];
sx q[0];
rz(-1.3246857) q[0];
sx q[0];
rz(-1.7497803) q[0];
x q[1];
rz(-0.71765064) q[2];
sx q[2];
rz(-0.96390488) q[2];
sx q[2];
rz(-0.90536149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47487709) q[1];
sx q[1];
rz(-1.8711539) q[1];
sx q[1];
rz(-1.1984065) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0456234) q[3];
sx q[3];
rz(-1.9672638) q[3];
sx q[3];
rz(3.0399655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0059263) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(0.076210991) q[2];
rz(2.5557319) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(-1.3425672) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12014408) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(-0.35010499) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(-1.2299445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0098388) q[0];
sx q[0];
rz(-1.3309877) q[0];
sx q[0];
rz(-1.9480223) q[0];
rz(-1.7361264) q[2];
sx q[2];
rz(-1.606719) q[2];
sx q[2];
rz(-1.3816116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1175244) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(2.9419521) q[1];
rz(-0.1996207) q[3];
sx q[3];
rz(-0.99310447) q[3];
sx q[3];
rz(-2.9132089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(-1.492307) q[2];
rz(0.40488511) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(-0.58240044) q[0];
rz(-2.4549585) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(-1.2581717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63623896) q[0];
sx q[0];
rz(-0.4551783) q[0];
sx q[0];
rz(1.2824461) q[0];
rz(-pi) q[1];
rz(-2.1297087) q[2];
sx q[2];
rz(-0.61043533) q[2];
sx q[2];
rz(-1.8545146) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88731474) q[1];
sx q[1];
rz(-1.5149024) q[1];
sx q[1];
rz(0.069737597) q[1];
rz(-pi) q[2];
rz(-2.9223676) q[3];
sx q[3];
rz(-0.98837438) q[3];
sx q[3];
rz(1.6265914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76576343) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(1.0180391) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(-1.6233981) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4356284) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(0.14353453) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-0.20763436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589446) q[0];
sx q[0];
rz(-2.2981854) q[0];
sx q[0];
rz(1.101053) q[0];
rz(-pi) q[1];
rz(-2.1499499) q[2];
sx q[2];
rz(-1.8657547) q[2];
sx q[2];
rz(0.3504914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27397284) q[1];
sx q[1];
rz(-1.311353) q[1];
sx q[1];
rz(0.8670437) q[1];
rz(-0.36334857) q[3];
sx q[3];
rz(-1.4814113) q[3];
sx q[3];
rz(-2.9575728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(0.60619727) q[2];
rz(-2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(-2.4351951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48924482) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(-3.0319013) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.441997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2588651) q[0];
sx q[0];
rz(-1.8804714) q[0];
sx q[0];
rz(0.41559269) q[0];
rz(-pi) q[1];
rz(2.6410854) q[2];
sx q[2];
rz(-1.8138514) q[2];
sx q[2];
rz(-0.99832035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7322526) q[1];
sx q[1];
rz(-0.85459083) q[1];
sx q[1];
rz(0.90998896) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23542102) q[3];
sx q[3];
rz(-1.306433) q[3];
sx q[3];
rz(2.2478769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66120061) q[2];
sx q[2];
rz(-2.9471687) q[2];
sx q[2];
rz(-1.9332168) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(-1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97063589) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(3.048625) q[0];
rz(-1.2804821) q[1];
sx q[1];
rz(-0.72851506) q[1];
sx q[1];
rz(-0.13519898) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057892) q[0];
sx q[0];
rz(-0.069878526) q[0];
sx q[0];
rz(-2.4326434) q[0];
rz(0.5298631) q[2];
sx q[2];
rz(-0.29680291) q[2];
sx q[2];
rz(-0.57020818) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0307903) q[1];
sx q[1];
rz(-2.8867278) q[1];
sx q[1];
rz(-0.084666208) q[1];
rz(-pi) q[2];
rz(1.1208833) q[3];
sx q[3];
rz(-2.9723047) q[3];
sx q[3];
rz(-1.5811046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70904237) q[2];
sx q[2];
rz(-1.9229527) q[2];
sx q[2];
rz(2.3700628) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.8621657) q[3];
sx q[3];
rz(0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5563357) q[0];
sx q[0];
rz(-2.519643) q[0];
sx q[0];
rz(-2.0167895) q[0];
rz(1.2987761) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(-2.3769456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9287531) q[0];
sx q[0];
rz(-0.1360341) q[0];
sx q[0];
rz(-0.68599756) q[0];
rz(-pi) q[1];
rz(-2.1642288) q[2];
sx q[2];
rz(-0.77319169) q[2];
sx q[2];
rz(0.14466454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.56988) q[1];
sx q[1];
rz(-2.9015744) q[1];
sx q[1];
rz(-2.1313138) q[1];
rz(-2.5705757) q[3];
sx q[3];
rz(-0.63435508) q[3];
sx q[3];
rz(2.3650124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(-0.20720227) q[2];
rz(2.1616705) q[3];
sx q[3];
rz(-2.0082974) q[3];
sx q[3];
rz(0.19237147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(0.021066396) q[0];
sx q[0];
rz(-1.6854032) q[0];
sx q[0];
rz(2.2717463) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(1.6899213) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(-2.5475827) q[3];
sx q[3];
rz(-0.28492155) q[3];
sx q[3];
rz(-0.84656885) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
