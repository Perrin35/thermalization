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
rz(-1.3574358) q[0];
sx q[0];
rz(3.7665851) q[0];
sx q[0];
rz(7.9023043) q[0];
rz(-1.9982665) q[1];
sx q[1];
rz(-0.63653094) q[1];
sx q[1];
rz(-1.7826537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5464294) q[0];
sx q[0];
rz(-0.63906416) q[0];
sx q[0];
rz(1.0733814) q[0];
rz(1.2011365) q[2];
sx q[2];
rz(-1.2944752) q[2];
sx q[2];
rz(-1.2075961) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5660043) q[1];
sx q[1];
rz(-1.8270565) q[1];
sx q[1];
rz(-0.23988597) q[1];
rz(-0.10349689) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(-0.58799839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4437359) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(-0.13992986) q[3];
sx q[3];
rz(-1.6401451) q[3];
sx q[3];
rz(-3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.3530465) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(2.7781558) q[0];
rz(2.4711171) q[1];
sx q[1];
rz(-1.8164219) q[1];
sx q[1];
rz(2.061981) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.712942) q[0];
sx q[0];
rz(-2.0976557) q[0];
sx q[0];
rz(0.27609472) q[0];
rz(-pi) q[1];
x q[1];
rz(3.011376) q[2];
sx q[2];
rz(-2.0998345) q[2];
sx q[2];
rz(2.0646273) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1324001) q[1];
sx q[1];
rz(-1.8120488) q[1];
sx q[1];
rz(3.0962837) q[1];
x q[2];
rz(0.3877181) q[3];
sx q[3];
rz(-2.5250375) q[3];
sx q[3];
rz(1.316837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4662027) q[2];
sx q[2];
rz(-2.4066996) q[2];
sx q[2];
rz(-0.43851635) q[2];
rz(2.130326) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(-1.4810168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587875) q[0];
sx q[0];
rz(-0.52849448) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(1.8939182) q[1];
sx q[1];
rz(-1.1670466) q[1];
sx q[1];
rz(-0.25965986) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6557216) q[0];
sx q[0];
rz(-1.8410793) q[0];
sx q[0];
rz(-0.87422687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38085085) q[2];
sx q[2];
rz(-1.2253739) q[2];
sx q[2];
rz(-2.7598515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5071755) q[1];
sx q[1];
rz(-0.89752642) q[1];
sx q[1];
rz(-0.61051621) q[1];
x q[2];
rz(-0.91783701) q[3];
sx q[3];
rz(-1.3363095) q[3];
sx q[3];
rz(-0.10374903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20716771) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(-2.3210607) q[2];
rz(-2.3151243) q[3];
sx q[3];
rz(-2.1491437) q[3];
sx q[3];
rz(1.8291399) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8823223) q[0];
sx q[0];
rz(-0.81040183) q[0];
sx q[0];
rz(-0.93233863) q[0];
rz(1.8238292) q[1];
sx q[1];
rz(-1.1612786) q[1];
sx q[1];
rz(-3.0188149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6500191) q[0];
sx q[0];
rz(-0.54560018) q[0];
sx q[0];
rz(0.16262098) q[0];
rz(1.2447692) q[2];
sx q[2];
rz(-2.3609128) q[2];
sx q[2];
rz(-0.29034607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9726561) q[1];
sx q[1];
rz(-2.3984809) q[1];
sx q[1];
rz(1.5692706) q[1];
x q[2];
rz(0.44466059) q[3];
sx q[3];
rz(-2.529699) q[3];
sx q[3];
rz(1.9652308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5968898) q[2];
sx q[2];
rz(-2.2201316) q[2];
sx q[2];
rz(0.60453647) q[2];
rz(-2.4321411) q[3];
sx q[3];
rz(-2.3716898) q[3];
sx q[3];
rz(2.3179222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17968793) q[0];
sx q[0];
rz(-0.12374319) q[0];
sx q[0];
rz(1.7748348) q[0];
rz(2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(-2.6008115) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3872469) q[0];
sx q[0];
rz(-1.7779113) q[0];
sx q[0];
rz(-1.4282754) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86414637) q[2];
sx q[2];
rz(-0.96071488) q[2];
sx q[2];
rz(2.6928201) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1652226) q[1];
sx q[1];
rz(-1.3787377) q[1];
sx q[1];
rz(-2.8403175) q[1];
rz(0.28022061) q[3];
sx q[3];
rz(-2.5570405) q[3];
sx q[3];
rz(-1.9909797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6775386) q[2];
sx q[2];
rz(-0.81391922) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(1.4520491) q[3];
sx q[3];
rz(-1.6477081) q[3];
sx q[3];
rz(-2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0238817) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(1.1718132) q[0];
rz(0.68012971) q[1];
sx q[1];
rz(-1.2260194) q[1];
sx q[1];
rz(-0.5562869) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6345917) q[0];
sx q[0];
rz(-1.1959574) q[0];
sx q[0];
rz(-2.8681173) q[0];
rz(-pi) q[1];
x q[1];
rz(1.589619) q[2];
sx q[2];
rz(-0.14655098) q[2];
sx q[2];
rz(2.3702247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7032675) q[1];
sx q[1];
rz(-2.5653953) q[1];
sx q[1];
rz(0.11771113) q[1];
x q[2];
rz(-2.049796) q[3];
sx q[3];
rz(-2.1508039) q[3];
sx q[3];
rz(0.73778462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80643001) q[2];
sx q[2];
rz(-1.8589636) q[2];
sx q[2];
rz(0.20509091) q[2];
rz(-2.3388376) q[3];
sx q[3];
rz(-1.1057248) q[3];
sx q[3];
rz(0.55454379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270808) q[0];
sx q[0];
rz(-2.0682122) q[0];
sx q[0];
rz(-2.9141973) q[0];
rz(-2.3511476) q[1];
sx q[1];
rz(-2.3643654) q[1];
sx q[1];
rz(0.8367742) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28115434) q[0];
sx q[0];
rz(-1.7895164) q[0];
sx q[0];
rz(0.52978306) q[0];
rz(-pi) q[1];
rz(-0.58800943) q[2];
sx q[2];
rz(-1.2704029) q[2];
sx q[2];
rz(-0.38640768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5143699) q[1];
sx q[1];
rz(-2.6191776) q[1];
sx q[1];
rz(0.66607968) q[1];
rz(-pi) q[2];
rz(-0.76117875) q[3];
sx q[3];
rz(-1.2952779) q[3];
sx q[3];
rz(-0.31106424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5922015) q[2];
sx q[2];
rz(-1.6624007) q[2];
sx q[2];
rz(-2.5679892) q[2];
rz(0.32522374) q[3];
sx q[3];
rz(-2.2461788) q[3];
sx q[3];
rz(-2.7371244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(-0.67894116) q[0];
rz(3.0067054) q[1];
sx q[1];
rz(-1.0604246) q[1];
sx q[1];
rz(1.4237684) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2436062) q[0];
sx q[0];
rz(-0.73315128) q[0];
sx q[0];
rz(0.85444684) q[0];
rz(-pi) q[1];
rz(-2.277454) q[2];
sx q[2];
rz(-1.1306964) q[2];
sx q[2];
rz(-2.3789483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9805124) q[1];
sx q[1];
rz(-0.80787611) q[1];
sx q[1];
rz(-0.041461583) q[1];
rz(-pi) q[2];
rz(2.1072949) q[3];
sx q[3];
rz(-0.52742672) q[3];
sx q[3];
rz(-2.3303383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5474825) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(-0.2317079) q[2];
rz(1.0101275) q[3];
sx q[3];
rz(-1.908327) q[3];
sx q[3];
rz(2.2505984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6758839) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(-2.6278507) q[0];
rz(1.4328009) q[1];
sx q[1];
rz(-1.2756196) q[1];
sx q[1];
rz(1.6339711) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5861478) q[0];
sx q[0];
rz(-2.9990104) q[0];
sx q[0];
rz(0.34839387) q[0];
x q[1];
rz(0.25357004) q[2];
sx q[2];
rz(-1.7268983) q[2];
sx q[2];
rz(-1.9954322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6027939) q[1];
sx q[1];
rz(-1.4414296) q[1];
sx q[1];
rz(-1.6954697) q[1];
x q[2];
rz(1.474103) q[3];
sx q[3];
rz(-2.8887199) q[3];
sx q[3];
rz(-2.0677572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12019176) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(-2.8821442) q[2];
rz(-1.1546968) q[3];
sx q[3];
rz(-2.3754933) q[3];
sx q[3];
rz(1.8416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(1.4890626) q[0];
rz(0.64000714) q[1];
sx q[1];
rz(-2.044544) q[1];
sx q[1];
rz(-2.1515813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3418808) q[0];
sx q[0];
rz(-2.9721301) q[0];
sx q[0];
rz(0.54903309) q[0];
x q[1];
rz(2.3732164) q[2];
sx q[2];
rz(-1.2207165) q[2];
sx q[2];
rz(-0.96736747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9089936) q[1];
sx q[1];
rz(-2.6910344) q[1];
sx q[1];
rz(-1.1934936) q[1];
rz(-3.1166273) q[3];
sx q[3];
rz(-1.7933266) q[3];
sx q[3];
rz(-2.3937267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2085569) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(-2.9162858) q[2];
rz(0.92132583) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(-3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0113572) q[0];
sx q[0];
rz(-2.7738032) q[0];
sx q[0];
rz(-2.2875447) q[0];
rz(-1.3182974) q[1];
sx q[1];
rz(-1.6164936) q[1];
sx q[1];
rz(2.2093538) q[1];
rz(1.83056) q[2];
sx q[2];
rz(-1.4767892) q[2];
sx q[2];
rz(-1.3108419) q[2];
rz(2.9574838) q[3];
sx q[3];
rz(-1.0297903) q[3];
sx q[3];
rz(-0.31576706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
