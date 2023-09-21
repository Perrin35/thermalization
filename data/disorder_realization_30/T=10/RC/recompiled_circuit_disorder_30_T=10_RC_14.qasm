OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8631247) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
rz(1.1236905) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(-2.6927039) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(2.0631454) q[1];
rz(3.0311534) q[3];
sx q[3];
rz(-2.7613598) q[3];
sx q[3];
rz(1.5649753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(0.079205714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4488892) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(-3.0485857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5835186) q[1];
sx q[1];
rz(-1.1221702) q[1];
sx q[1];
rz(-3.0706057) q[1];
x q[2];
rz(0.077640688) q[3];
sx q[3];
rz(-1.75378) q[3];
sx q[3];
rz(1.3444572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(1.0650939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247892) q[0];
sx q[0];
rz(-1.5771126) q[0];
sx q[0];
rz(1.9301231) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58653411) q[2];
sx q[2];
rz(-0.9622935) q[2];
sx q[2];
rz(1.8412794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0703735) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25032708) q[3];
sx q[3];
rz(-2.2776789) q[3];
sx q[3];
rz(1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(0.68850368) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(-0.011118523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.714689) q[0];
sx q[0];
rz(-0.62781292) q[0];
sx q[0];
rz(-2.3054302) q[0];
rz(2.501802) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(-1.421979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74777714) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(-1.0541037) q[1];
rz(-pi) q[2];
rz(2.6075881) q[3];
sx q[3];
rz(-2.7760193) q[3];
sx q[3];
rz(0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-0.28309506) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(0.33183137) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(1.3269075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834284) q[0];
sx q[0];
rz(-2.313662) q[0];
sx q[0];
rz(1.3616256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1660216) q[2];
sx q[2];
rz(-0.86704463) q[2];
sx q[2];
rz(0.032035839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9304282) q[1];
sx q[1];
rz(-1.2439562) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(-pi) q[2];
rz(-1.9782449) q[3];
sx q[3];
rz(-1.3031928) q[3];
sx q[3];
rz(-2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(-2.025827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003446) q[0];
sx q[0];
rz(-1.9755409) q[0];
sx q[0];
rz(-0.83165283) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0954779) q[2];
sx q[2];
rz(-2.4002053) q[2];
sx q[2];
rz(-1.807883) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(0.58847217) q[1];
rz(-2.9168105) q[3];
sx q[3];
rz(-0.68192476) q[3];
sx q[3];
rz(-3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(0.3195233) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4903957) q[0];
sx q[0];
rz(-2.2389452) q[0];
sx q[0];
rz(0.018880318) q[0];
rz(-pi) q[1];
x q[1];
rz(2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(-0.069318511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9971804) q[1];
sx q[1];
rz(-1.715853) q[1];
sx q[1];
rz(-0.1952862) q[1];
rz(-pi) q[2];
rz(-0.087248487) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(-1.5213535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053145807) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(-1.501207) q[0];
rz(-pi) q[1];
rz(-2.0422158) q[2];
sx q[2];
rz(-1.2087012) q[2];
sx q[2];
rz(2.9277756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13008598) q[1];
sx q[1];
rz(-0.79400051) q[1];
sx q[1];
rz(1.6542875) q[1];
rz(-pi) q[2];
rz(2.2500854) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(-0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-2.14595) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(2.8709581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7662738) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(-0.71592285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.052202) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(-2.7871728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.7120413) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(-0.045750381) q[1];
x q[2];
rz(0.6587894) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(-0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(0.45544004) q[2];
rz(-2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(2.4023138) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56787581) q[0];
sx q[0];
rz(-2.0567237) q[0];
sx q[0];
rz(3.0987415) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32785373) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(1.9890131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.478185) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(0.138476) q[1];
rz(-pi) q[2];
rz(-3.0462618) q[3];
sx q[3];
rz(-0.87626002) q[3];
sx q[3];
rz(-1.8180083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(2.949529) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(1.0275502) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(-2.8752747) q[3];
sx q[3];
rz(-1.8177633) q[3];
sx q[3];
rz(0.37533356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
