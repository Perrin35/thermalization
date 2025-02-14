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
rz(1.3266069) q[0];
sx q[0];
rz(-2.9171483) q[0];
sx q[0];
rz(1.8087968) q[0];
rz(-0.83238554) q[1];
sx q[1];
rz(4.8424911) q[1];
sx q[1];
rz(10.535156) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8827986) q[0];
sx q[0];
rz(-1.4075532) q[0];
sx q[0];
rz(1.4397653) q[0];
rz(-0.85904636) q[2];
sx q[2];
rz(-1.8084322) q[2];
sx q[2];
rz(0.29881921) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.367358) q[1];
sx q[1];
rz(-1.5176766) q[1];
sx q[1];
rz(-1.5426643) q[1];
rz(-2.7361462) q[3];
sx q[3];
rz(-0.67908248) q[3];
sx q[3];
rz(0.17857851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9139468) q[2];
sx q[2];
rz(-3.1326742) q[2];
sx q[2];
rz(1.439636) q[2];
rz(1.7281744) q[3];
sx q[3];
rz(-3.1296802) q[3];
sx q[3];
rz(0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.444376) q[0];
sx q[0];
rz(-1.603729) q[0];
sx q[0];
rz(-2.5268396) q[0];
rz(-0.53601021) q[1];
sx q[1];
rz(-0.025608048) q[1];
sx q[1];
rz(0.33686179) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0885411) q[0];
sx q[0];
rz(-1.4362925) q[0];
sx q[0];
rz(0.10945871) q[0];
rz(0.031177715) q[2];
sx q[2];
rz(-2.665069) q[2];
sx q[2];
rz(-0.88821238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.58648169) q[1];
sx q[1];
rz(-0.046849061) q[1];
sx q[1];
rz(1.2270801) q[1];
rz(-pi) q[2];
rz(2.0982127) q[3];
sx q[3];
rz(-2.7763397) q[3];
sx q[3];
rz(3.0556847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2577995) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(-2.4419355) q[2];
rz(2.9720225) q[3];
sx q[3];
rz(-3.1278059) q[3];
sx q[3];
rz(-0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70497847) q[0];
sx q[0];
rz(-0.545937) q[0];
sx q[0];
rz(-0.28526947) q[0];
rz(0.26372313) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(2.347351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6503822) q[0];
sx q[0];
rz(-0.69612078) q[0];
sx q[0];
rz(0.18265035) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2474485) q[2];
sx q[2];
rz(-1.4686606) q[2];
sx q[2];
rz(-2.7030526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31595597) q[1];
sx q[1];
rz(-1.5772305) q[1];
sx q[1];
rz(3.107162) q[1];
rz(-pi) q[2];
rz(0.66624347) q[3];
sx q[3];
rz(-1.0292205) q[3];
sx q[3];
rz(-0.2592087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11702015) q[2];
sx q[2];
rz(-0.046155013) q[2];
sx q[2];
rz(-2.2771007) q[2];
rz(0.74508673) q[3];
sx q[3];
rz(-0.86747187) q[3];
sx q[3];
rz(-3.0985221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2660148) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(2.2752046) q[0];
rz(0.59151793) q[1];
sx q[1];
rz(-2.1951127) q[1];
sx q[1];
rz(-2.0902858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54615462) q[0];
sx q[0];
rz(-0.0020448472) q[0];
sx q[0];
rz(-0.05495222) q[0];
rz(-3.1411416) q[2];
sx q[2];
rz(-1.573083) q[2];
sx q[2];
rz(-0.0046943517) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2308553) q[1];
sx q[1];
rz(-1.1713109) q[1];
sx q[1];
rz(-2.3302156) q[1];
rz(-pi) q[2];
rz(2.0962001) q[3];
sx q[3];
rz(-2.0936493) q[3];
sx q[3];
rz(-2.9749572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5279348) q[2];
sx q[2];
rz(-0.064585678) q[2];
sx q[2];
rz(-0.85395542) q[2];
rz(2.4438786) q[3];
sx q[3];
rz(-1.8055975) q[3];
sx q[3];
rz(-0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39603221) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(-2.7295617) q[0];
rz(-1.283006) q[1];
sx q[1];
rz(-0.7737776) q[1];
sx q[1];
rz(2.3903019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728516) q[0];
sx q[0];
rz(-1.6882467) q[0];
sx q[0];
rz(-0.23953481) q[0];
x q[1];
rz(2.2504248) q[2];
sx q[2];
rz(-1.749265) q[2];
sx q[2];
rz(2.8986069) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48482286) q[1];
sx q[1];
rz(-2.1570766) q[1];
sx q[1];
rz(2.1692679) q[1];
rz(-pi) q[2];
rz(0.86814084) q[3];
sx q[3];
rz(-2.9538547) q[3];
sx q[3];
rz(-0.38892285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46386197) q[2];
sx q[2];
rz(-3.1158713) q[2];
sx q[2];
rz(-2.1417997) q[2];
rz(2.540588) q[3];
sx q[3];
rz(-0.060636245) q[3];
sx q[3];
rz(2.291688) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8782225) q[0];
sx q[0];
rz(-0.090494089) q[0];
sx q[0];
rz(0.40060842) q[0];
rz(1.7349617) q[1];
sx q[1];
rz(-0.35146439) q[1];
sx q[1];
rz(1.7469143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0907441) q[0];
sx q[0];
rz(-1.5904477) q[0];
sx q[0];
rz(-3.1369563) q[0];
rz(-0.80446316) q[2];
sx q[2];
rz(-1.7906252) q[2];
sx q[2];
rz(1.1457844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38401022) q[1];
sx q[1];
rz(-1.3802802) q[1];
sx q[1];
rz(-1.5662929) q[1];
rz(-1.3974779) q[3];
sx q[3];
rz(-2.299274) q[3];
sx q[3];
rz(2.4761379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3957735) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(-2.0664717) q[2];
rz(0.51946688) q[3];
sx q[3];
rz(-0.1778917) q[3];
sx q[3];
rz(1.547267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2125242) q[0];
sx q[0];
rz(-1.9269749) q[0];
sx q[0];
rz(1.6125096) q[0];
rz(-0.56135881) q[1];
sx q[1];
rz(-3.1293271) q[1];
sx q[1];
rz(0.57048172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9425526) q[0];
sx q[0];
rz(-1.4496452) q[0];
sx q[0];
rz(3.0899307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9134741) q[2];
sx q[2];
rz(-1.2608796) q[2];
sx q[2];
rz(0.44911929) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.884932) q[1];
sx q[1];
rz(-0.015555582) q[1];
sx q[1];
rz(1.6473098) q[1];
rz(-0.99659427) q[3];
sx q[3];
rz(-2.3086267) q[3];
sx q[3];
rz(-1.0723237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.059171112) q[2];
sx q[2];
rz(-2.8736281) q[2];
sx q[2];
rz(-1.9783665) q[2];
rz(1.5393114) q[3];
sx q[3];
rz(-3.0395165) q[3];
sx q[3];
rz(0.99360895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7546643) q[0];
sx q[0];
rz(-0.29351497) q[0];
sx q[0];
rz(-2.5013404) q[0];
rz(-0.86296588) q[1];
sx q[1];
rz(-0.14185618) q[1];
sx q[1];
rz(-0.77756768) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1877768) q[0];
sx q[0];
rz(-1.1788771) q[0];
sx q[0];
rz(-1.1338393) q[0];
rz(-pi) q[1];
rz(0.24224327) q[2];
sx q[2];
rz(-2.3628334) q[2];
sx q[2];
rz(-0.81204108) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.150084) q[1];
sx q[1];
rz(-1.5822295) q[1];
sx q[1];
rz(1.6437794) q[1];
x q[2];
rz(-0.74109091) q[3];
sx q[3];
rz(-2.5020385) q[3];
sx q[3];
rz(-3.0921788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61206478) q[2];
sx q[2];
rz(-0.42403388) q[2];
sx q[2];
rz(2.6112134) q[2];
rz(-0.027675962) q[3];
sx q[3];
rz(-3.0961302) q[3];
sx q[3];
rz(1.8224705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46488047) q[0];
sx q[0];
rz(-0.16093971) q[0];
sx q[0];
rz(0.93223923) q[0];
rz(1.5981916) q[1];
sx q[1];
rz(-0.80359572) q[1];
sx q[1];
rz(-2.7995321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8300104) q[0];
sx q[0];
rz(-1.5556635) q[0];
sx q[0];
rz(-0.084447817) q[0];
x q[1];
rz(1.507296) q[2];
sx q[2];
rz(-0.8339774) q[2];
sx q[2];
rz(-2.8066016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2310357) q[1];
sx q[1];
rz(-2.0879832) q[1];
sx q[1];
rz(1.789235) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4985648) q[3];
sx q[3];
rz(-0.22966239) q[3];
sx q[3];
rz(0.1230474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0590234) q[2];
sx q[2];
rz(-3.1408568) q[2];
sx q[2];
rz(-1.9334582) q[2];
rz(2.1900603) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(2.4212196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7893938) q[0];
sx q[0];
rz(-2.3645526) q[0];
sx q[0];
rz(0.081789516) q[0];
rz(0.31546047) q[1];
sx q[1];
rz(-3.0888562) q[1];
sx q[1];
rz(-1.2299406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953461) q[0];
sx q[0];
rz(-2.917756) q[0];
sx q[0];
rz(2.361103) q[0];
x q[1];
rz(-0.17773262) q[2];
sx q[2];
rz(-1.4338253) q[2];
sx q[2];
rz(-2.9229119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9435159) q[1];
sx q[1];
rz(-1.480989) q[1];
sx q[1];
rz(-1.4526558) q[1];
rz(1.0395524) q[3];
sx q[3];
rz(-1.6907516) q[3];
sx q[3];
rz(-2.1841072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6007467) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(2.02796) q[2];
rz(-0.039637808) q[3];
sx q[3];
rz(-0.0097291917) q[3];
sx q[3];
rz(2.5447194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46350805) q[0];
sx q[0];
rz(-1.2152553) q[0];
sx q[0];
rz(-1.8334462) q[0];
rz(-2.5254163) q[1];
sx q[1];
rz(-1.0721075) q[1];
sx q[1];
rz(0.20872605) q[1];
rz(-1.8088593) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.0669591) q[2];
rz(1.3499089) q[3];
sx q[3];
rz(-1.532423) q[3];
sx q[3];
rz(1.5929089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
