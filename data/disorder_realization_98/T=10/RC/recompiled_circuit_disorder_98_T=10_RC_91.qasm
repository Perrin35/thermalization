OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5538841) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(0.61650886) q[0];
rz(-1.0295463) q[2];
sx q[2];
rz(-1.1272578) q[2];
sx q[2];
rz(0.29543791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40587273) q[1];
sx q[1];
rz(-0.18438965) q[1];
sx q[1];
rz(-1.7806446) q[1];
rz(2.9505694) q[3];
sx q[3];
rz(-1.028562) q[3];
sx q[3];
rz(0.99381002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(1.4037508) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(-2.4893563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61440496) q[0];
sx q[0];
rz(-1.5831278) q[0];
sx q[0];
rz(3.1129818) q[0];
rz(-pi) q[1];
rz(2.2488238) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(-0.54112753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10416874) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(1.8038521) q[1];
x q[2];
rz(-2.5370876) q[3];
sx q[3];
rz(-1.0009871) q[3];
sx q[3];
rz(-1.2165716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.3495548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916546) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(1.2664938) q[0];
x q[1];
rz(-1.0565874) q[2];
sx q[2];
rz(-3.0658256) q[2];
sx q[2];
rz(0.55759341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50074358) q[1];
sx q[1];
rz(-1.8060246) q[1];
sx q[1];
rz(-0.51847036) q[1];
x q[2];
rz(-1.549167) q[3];
sx q[3];
rz(-0.4261741) q[3];
sx q[3];
rz(2.222995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0321908) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0050126652) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(1.2269292) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47469791) q[2];
sx q[2];
rz(-0.63650741) q[2];
sx q[2];
rz(-2.0091025) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8424884) q[1];
sx q[1];
rz(-2.4310388) q[1];
sx q[1];
rz(-1.3683661) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23457228) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(-2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(2.9072445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34412947) q[0];
sx q[0];
rz(-1.2327694) q[0];
sx q[0];
rz(-0.56030886) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0460065) q[2];
sx q[2];
rz(-1.0786973) q[2];
sx q[2];
rz(1.0620067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.105141) q[1];
sx q[1];
rz(-0.76117939) q[1];
sx q[1];
rz(-0.21703227) q[1];
x q[2];
rz(-1.8946394) q[3];
sx q[3];
rz(-1.4959335) q[3];
sx q[3];
rz(0.74136855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(-2.9210572) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.9979427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8457984) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(-0.11008115) q[0];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(-0.93271819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3559349) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(1.7788586) q[1];
x q[2];
rz(-0.46623047) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(-3.1371491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75366655) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.9082665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(-1.5412488) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(-0.90460888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.308097) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(-1.584946) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2193905) q[3];
sx q[3];
rz(-2.0316342) q[3];
sx q[3];
rz(1.2724862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.258237) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90342605) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(-1.1358791) q[0];
x q[1];
rz(0.7779185) q[2];
sx q[2];
rz(-2.3292543) q[2];
sx q[2];
rz(1.9922436) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5453323) q[1];
sx q[1];
rz(-0.55570554) q[1];
sx q[1];
rz(-0.19821367) q[1];
x q[2];
rz(0.97790896) q[3];
sx q[3];
rz(-0.80855723) q[3];
sx q[3];
rz(-2.9582994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(0.2510221) q[0];
rz(0.42731467) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-0.02773157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160307) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-0.82657878) q[0];
rz(-0.8717732) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(-2.2733462) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30222826) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(2.539413) q[1];
x q[2];
rz(0.025891993) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47782183) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.1429943) q[0];
rz(-0.54458877) q[2];
sx q[2];
rz(-1.1368183) q[2];
sx q[2];
rz(1.7314272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.232302) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(0.69625744) q[1];
rz(1.0912283) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(-1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(0.67129927) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(-2.0444617) q[3];
sx q[3];
rz(-2.6087425) q[3];
sx q[3];
rz(-2.9125924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
