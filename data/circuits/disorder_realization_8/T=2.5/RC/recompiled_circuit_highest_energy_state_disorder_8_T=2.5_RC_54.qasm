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
rz(-1.8149858) q[0];
sx q[0];
rz(-0.22444434) q[0];
sx q[0];
rz(1.3327959) q[0];
rz(-0.83238554) q[1];
sx q[1];
rz(4.8424911) q[1];
sx q[1];
rz(10.535156) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334167) q[0];
sx q[0];
rz(-1.7000755) q[0];
sx q[0];
rz(0.16462932) q[0];
rz(-0.85904636) q[2];
sx q[2];
rz(-1.8084322) q[2];
sx q[2];
rz(-2.8427734) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8800115) q[1];
sx q[1];
rz(-0.060102988) q[1];
sx q[1];
rz(2.6549935) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63817231) q[3];
sx q[3];
rz(-1.8211357) q[3];
sx q[3];
rz(1.7145715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9139468) q[2];
sx q[2];
rz(-0.0089184428) q[2];
sx q[2];
rz(1.7019567) q[2];
rz(1.4134183) q[3];
sx q[3];
rz(-0.011912502) q[3];
sx q[3];
rz(0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.444376) q[0];
sx q[0];
rz(-1.603729) q[0];
sx q[0];
rz(2.5268396) q[0];
rz(0.53601021) q[1];
sx q[1];
rz(-0.025608048) q[1];
sx q[1];
rz(-0.33686179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63348544) q[0];
sx q[0];
rz(-2.9683873) q[0];
sx q[0];
rz(2.2499535) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1104149) q[2];
sx q[2];
rz(-2.665069) q[2];
sx q[2];
rz(-2.2533803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24241703) q[1];
sx q[1];
rz(-1.6149033) q[1];
sx q[1];
rz(0.015797829) q[1];
rz(2.9514489) q[3];
sx q[3];
rz(-1.8845673) q[3];
sx q[3];
rz(0.47167512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88379318) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(-0.69965714) q[2];
rz(-0.16957016) q[3];
sx q[3];
rz(-0.01378672) q[3];
sx q[3];
rz(-2.5794896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4366142) q[0];
sx q[0];
rz(-2.5956557) q[0];
sx q[0];
rz(2.8563232) q[0];
rz(-0.26372313) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(-2.347351) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4141615) q[0];
sx q[0];
rz(-0.88849706) q[0];
sx q[0];
rz(-1.4201565) q[0];
rz(-pi) q[1];
rz(-1.8828234) q[2];
sx q[2];
rz(-0.3385545) q[2];
sx q[2];
rz(0.83690137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4395098) q[1];
sx q[1];
rz(-3.1065662) q[1];
sx q[1];
rz(0.18478025) q[1];
x q[2];
rz(-0.91762791) q[3];
sx q[3];
rz(-1.0126739) q[3];
sx q[3];
rz(2.2150326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11702015) q[2];
sx q[2];
rz(-0.046155013) q[2];
sx q[2];
rz(2.2771007) q[2];
rz(-0.74508673) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(0.043070506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87557781) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(0.86638802) q[0];
rz(2.5500747) q[1];
sx q[1];
rz(-0.94647995) q[1];
sx q[1];
rz(-2.0902858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54615462) q[0];
sx q[0];
rz(-3.1395478) q[0];
sx q[0];
rz(3.0866404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5685097) q[2];
sx q[2];
rz(-1.5703452) q[2];
sx q[2];
rz(1.5754896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4128139) q[1];
sx q[1];
rz(-2.3024913) q[1];
sx q[1];
rz(-1.0207291) q[1];
rz(-pi) q[2];
rz(1.0453926) q[3];
sx q[3];
rz(-2.0936493) q[3];
sx q[3];
rz(2.9749572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5279348) q[2];
sx q[2];
rz(-3.077007) q[2];
sx q[2];
rz(-2.2876372) q[2];
rz(-2.4438786) q[3];
sx q[3];
rz(-1.3359952) q[3];
sx q[3];
rz(-0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39603221) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(0.41203099) q[0];
rz(-1.8585867) q[1];
sx q[1];
rz(-0.7737776) q[1];
sx q[1];
rz(0.7512908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17344483) q[0];
sx q[0];
rz(-1.8086489) q[0];
sx q[0];
rz(1.449927) q[0];
rz(0.89116781) q[2];
sx q[2];
rz(-1.749265) q[2];
sx q[2];
rz(-2.8986069) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4163782) q[1];
sx q[1];
rz(-2.0592923) q[1];
sx q[1];
rz(0.67711551) q[1];
rz(3.0194332) q[3];
sx q[3];
rz(-1.7137104) q[3];
sx q[3];
rz(1.1003332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6777307) q[2];
sx q[2];
rz(-3.1158713) q[2];
sx q[2];
rz(-2.1417997) q[2];
rz(0.60100466) q[3];
sx q[3];
rz(-3.0809564) q[3];
sx q[3];
rz(2.291688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18086223) q[0];
sx q[0];
rz(-3.1214018) q[0];
sx q[0];
rz(-1.8024615) q[0];
rz(-0.30072923) q[2];
sx q[2];
rz(-0.82737714) q[2];
sx q[2];
rz(-2.9234795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38401022) q[1];
sx q[1];
rz(-1.7613125) q[1];
sx q[1];
rz(1.5752998) q[1];
x q[2];
rz(-2.4056099) q[3];
sx q[3];
rz(-1.6998359) q[3];
sx q[3];
rz(2.1202212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74581915) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(2.0664717) q[2];
rz(-0.51946688) q[3];
sx q[3];
rz(-0.1778917) q[3];
sx q[3];
rz(1.5943257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2125242) q[0];
sx q[0];
rz(-1.9269749) q[0];
sx q[0];
rz(-1.6125096) q[0];
rz(-0.56135881) q[1];
sx q[1];
rz(-3.1293271) q[1];
sx q[1];
rz(0.57048172) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9425526) q[0];
sx q[0];
rz(-1.4496452) q[0];
sx q[0];
rz(0.051661927) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2531742) q[2];
sx q[2];
rz(-1.7878685) q[2];
sx q[2];
rz(-2.0906015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25666061) q[1];
sx q[1];
rz(-3.1260371) q[1];
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
rz(3.0824215) q[2];
sx q[2];
rz(-0.26796451) q[2];
sx q[2];
rz(-1.9783665) q[2];
rz(-1.6022812) q[3];
sx q[3];
rz(-3.0395165) q[3];
sx q[3];
rz(0.99360895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38692835) q[0];
sx q[0];
rz(-0.29351497) q[0];
sx q[0];
rz(0.64025229) q[0];
rz(2.2786268) q[1];
sx q[1];
rz(-2.9997365) q[1];
sx q[1];
rz(0.77756768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8391188) q[0];
sx q[0];
rz(-0.57838744) q[0];
sx q[0];
rz(0.79721398) q[0];
x q[1];
rz(1.338358) q[2];
sx q[2];
rz(-0.82045499) q[2];
sx q[2];
rz(-2.6636843) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41987644) q[1];
sx q[1];
rz(-1.6437746) q[1];
sx q[1];
rz(-3.130129) q[1];
rz(0.50189314) q[3];
sx q[3];
rz(-1.9855026) q[3];
sx q[3];
rz(2.1547174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61206478) q[2];
sx q[2];
rz(-2.7175588) q[2];
sx q[2];
rz(-0.5303793) q[2];
rz(0.027675962) q[3];
sx q[3];
rz(-0.045462463) q[3];
sx q[3];
rz(-1.3191222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767122) q[0];
sx q[0];
rz(-0.16093971) q[0];
sx q[0];
rz(2.2093534) q[0];
rz(1.5981916) q[1];
sx q[1];
rz(-0.80359572) q[1];
sx q[1];
rz(0.34206051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0592744) q[0];
sx q[0];
rz(-0.085789811) q[0];
sx q[0];
rz(0.17753521) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4037696) q[2];
sx q[2];
rz(-1.5237817) q[2];
sx q[2];
rz(1.948487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4885713) q[1];
sx q[1];
rz(-2.584051) q[1];
sx q[1];
rz(-0.36399059) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1247221) q[3];
sx q[3];
rz(-1.3417435) q[3];
sx q[3];
rz(3.0927174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.082569294) q[2];
sx q[2];
rz(-0.00073585357) q[2];
sx q[2];
rz(1.2081344) q[2];
rz(0.9515323) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(-2.4212196) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35219881) q[0];
sx q[0];
rz(-0.77704) q[0];
sx q[0];
rz(-0.081789516) q[0];
rz(-2.8261322) q[1];
sx q[1];
rz(-3.0888562) q[1];
sx q[1];
rz(1.9116521) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0462466) q[0];
sx q[0];
rz(-2.917756) q[0];
sx q[0];
rz(-2.361103) q[0];
x q[1];
rz(0.66218485) q[2];
sx q[2];
rz(-2.9176468) q[2];
sx q[2];
rz(1.1395154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3620748) q[1];
sx q[1];
rz(-1.6884585) q[1];
sx q[1];
rz(-3.0511583) q[1];
rz(-pi) q[2];
x q[2];
rz(1.337215) q[3];
sx q[3];
rz(-0.54334917) q[3];
sx q[3];
rz(0.41239967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54084593) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(2.02796) q[2];
rz(0.039637808) q[3];
sx q[3];
rz(-0.0097291917) q[3];
sx q[3];
rz(-2.5447194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780846) q[0];
sx q[0];
rz(-1.2152553) q[0];
sx q[0];
rz(-1.8334462) q[0];
rz(-2.5254163) q[1];
sx q[1];
rz(-1.0721075) q[1];
sx q[1];
rz(0.20872605) q[1];
rz(-2.9504572) q[2];
sx q[2];
rz(-1.8046817) q[2];
sx q[2];
rz(2.6006151) q[2];
rz(-1.7916837) q[3];
sx q[3];
rz(-1.532423) q[3];
sx q[3];
rz(1.5929089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
