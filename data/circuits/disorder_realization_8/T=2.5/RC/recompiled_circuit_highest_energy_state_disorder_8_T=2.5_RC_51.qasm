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
rz(2.3092071) q[1];
sx q[1];
rz(-1.7008984) q[1];
sx q[1];
rz(2.031215) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8081759) q[0];
sx q[0];
rz(-1.7000755) q[0];
sx q[0];
rz(0.16462932) q[0];
rz(-pi) q[1];
rz(-2.832007) q[2];
sx q[2];
rz(-0.88304115) q[2];
sx q[2];
rz(-1.0716455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.367358) q[1];
sx q[1];
rz(-1.5176766) q[1];
sx q[1];
rz(-1.5989283) q[1];
x q[2];
rz(-1.8790122) q[3];
sx q[3];
rz(-2.1860414) q[3];
sx q[3];
rz(0.32549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2276459) q[2];
sx q[2];
rz(-3.1326742) q[2];
sx q[2];
rz(-1.7019567) q[2];
rz(1.4134183) q[3];
sx q[3];
rz(-3.1296802) q[3];
sx q[3];
rz(-0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.444376) q[0];
sx q[0];
rz(-1.603729) q[0];
sx q[0];
rz(0.61475301) q[0];
rz(0.53601021) q[1];
sx q[1];
rz(-3.1159846) q[1];
sx q[1];
rz(0.33686179) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5081072) q[0];
sx q[0];
rz(-0.17320536) q[0];
sx q[0];
rz(-2.2499535) q[0];
rz(-pi) q[1];
rz(1.5868863) q[2];
sx q[2];
rz(-2.0470691) q[2];
sx q[2];
rz(-2.2884638) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.81391) q[1];
sx q[1];
rz(-1.5865788) q[1];
sx q[1];
rz(-1.6149088) q[1];
x q[2];
rz(-0.19014374) q[3];
sx q[3];
rz(-1.8845673) q[3];
sx q[3];
rz(0.47167512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88379318) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(0.69965714) q[2];
rz(2.9720225) q[3];
sx q[3];
rz(-0.01378672) q[3];
sx q[3];
rz(0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70497847) q[0];
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
rz(-3.0803826) q[0];
sx q[0];
rz(-1.4540577) q[0];
sx q[0];
rz(2.453713) q[0];
x q[1];
rz(3.0339166) q[2];
sx q[2];
rz(-1.2491944) q[2];
sx q[2];
rz(1.9751825) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8869739) q[1];
sx q[1];
rz(-1.5363664) q[1];
sx q[1];
rz(1.5772343) q[1];
rz(-pi) q[2];
rz(2.3696926) q[3];
sx q[3];
rz(-2.3099314) q[3];
sx q[3];
rz(-1.2496304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0245725) q[2];
sx q[2];
rz(-3.0954376) q[2];
sx q[2];
rz(2.2771007) q[2];
rz(-2.3965059) q[3];
sx q[3];
rz(-2.2741208) q[3];
sx q[3];
rz(3.0985221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2660148) q[0];
sx q[0];
rz(-3.0953396) q[0];
sx q[0];
rz(-2.2752046) q[0];
rz(0.59151793) q[1];
sx q[1];
rz(-0.94647995) q[1];
sx q[1];
rz(-1.0513069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619988) q[0];
sx q[0];
rz(-1.5709086) q[0];
sx q[0];
rz(-3.1395509) q[0];
x q[1];
rz(1.5685097) q[2];
sx q[2];
rz(-1.5703452) q[2];
sx q[2];
rz(1.5754896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72877872) q[1];
sx q[1];
rz(-0.83910131) q[1];
sx q[1];
rz(-1.0207291) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0962001) q[3];
sx q[3];
rz(-2.0936493) q[3];
sx q[3];
rz(0.16663545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6136578) q[2];
sx q[2];
rz(-3.077007) q[2];
sx q[2];
rz(-2.2876372) q[2];
rz(2.4438786) q[3];
sx q[3];
rz(-1.3359952) q[3];
sx q[3];
rz(0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39603221) q[0];
sx q[0];
rz(-3.0389391) q[0];
sx q[0];
rz(-0.41203099) q[0];
rz(-1.283006) q[1];
sx q[1];
rz(-0.7737776) q[1];
sx q[1];
rz(-0.7512908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9681478) q[0];
sx q[0];
rz(-1.3329437) q[0];
sx q[0];
rz(1.6916656) q[0];
rz(-pi) q[1];
rz(-2.9137026) q[2];
sx q[2];
rz(-0.90393674) q[2];
sx q[2];
rz(-1.6713072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4163782) q[1];
sx q[1];
rz(-2.0592923) q[1];
sx q[1];
rz(0.67711551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2734518) q[3];
sx q[3];
rz(-2.9538547) q[3];
sx q[3];
rz(2.7526698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6777307) q[2];
sx q[2];
rz(-3.1158713) q[2];
sx q[2];
rz(-0.99979293) q[2];
rz(2.540588) q[3];
sx q[3];
rz(-3.0809564) q[3];
sx q[3];
rz(0.84990466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26337013) q[0];
sx q[0];
rz(-3.0510986) q[0];
sx q[0];
rz(0.40060842) q[0];
rz(1.406631) q[1];
sx q[1];
rz(-2.7901283) q[1];
sx q[1];
rz(1.7469143) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9607304) q[0];
sx q[0];
rz(-3.1214018) q[0];
sx q[0];
rz(-1.8024615) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30072923) q[2];
sx q[2];
rz(-2.3142155) q[2];
sx q[2];
rz(0.21811315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36023271) q[1];
sx q[1];
rz(-0.1905687) q[1];
sx q[1];
rz(3.1182454) q[1];
rz(-pi) q[2];
rz(2.9506545) q[3];
sx q[3];
rz(-0.74511601) q[3];
sx q[3];
rz(-2.7332954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3957735) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(1.0751209) q[2];
rz(2.6221258) q[3];
sx q[3];
rz(-2.963701) q[3];
sx q[3];
rz(-1.5943257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290685) q[0];
sx q[0];
rz(-1.2146177) q[0];
sx q[0];
rz(1.6125096) q[0];
rz(2.5802338) q[1];
sx q[1];
rz(-0.012265597) q[1];
sx q[1];
rz(-0.57048172) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5382696) q[0];
sx q[0];
rz(-3.0099359) q[0];
sx q[0];
rz(-1.1696474) q[0];
x q[1];
rz(-0.95593234) q[2];
sx q[2];
rz(-2.758965) q[2];
sx q[2];
rz(-1.0997694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25666061) q[1];
sx q[1];
rz(-3.1260371) q[1];
sx q[1];
rz(1.4942829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6030406) q[3];
sx q[3];
rz(-2.2411514) q[3];
sx q[3];
rz(-2.8351401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.059171112) q[2];
sx q[2];
rz(-0.26796451) q[2];
sx q[2];
rz(1.1632261) q[2];
rz(1.6022812) q[3];
sx q[3];
rz(-0.10207615) q[3];
sx q[3];
rz(-2.1479837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7546643) q[0];
sx q[0];
rz(-0.29351497) q[0];
sx q[0];
rz(2.5013404) q[0];
rz(0.86296588) q[1];
sx q[1];
rz(-0.14185618) q[1];
sx q[1];
rz(0.77756768) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7011576) q[0];
sx q[0];
rz(-1.1689725) q[0];
sx q[0];
rz(-0.42796414) q[0];
rz(-pi) q[1];
rz(-0.76394947) q[2];
sx q[2];
rz(-1.4014971) q[2];
sx q[2];
rz(0.93285376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7217162) q[1];
sx q[1];
rz(-1.4978181) q[1];
sx q[1];
rz(3.130129) q[1];
rz(-pi) q[2];
rz(-0.50189314) q[3];
sx q[3];
rz(-1.9855026) q[3];
sx q[3];
rz(0.98687526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5295279) q[2];
sx q[2];
rz(-0.42403388) q[2];
sx q[2];
rz(0.5303793) q[2];
rz(3.1139167) q[3];
sx q[3];
rz(-0.045462463) q[3];
sx q[3];
rz(1.3191222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.2093534) q[0];
rz(-1.543401) q[1];
sx q[1];
rz(-2.3379969) q[1];
sx q[1];
rz(-0.34206051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0823183) q[0];
sx q[0];
rz(-0.085789811) q[0];
sx q[0];
rz(0.17753521) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069829536) q[2];
sx q[2];
rz(-2.402555) q[2];
sx q[2];
rz(-0.42934092) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6530214) q[1];
sx q[1];
rz(-2.584051) q[1];
sx q[1];
rz(2.7776021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7998806) q[3];
sx q[3];
rz(-1.5872262) q[3];
sx q[3];
rz(-1.5180902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.082569294) q[2];
sx q[2];
rz(-0.00073585357) q[2];
sx q[2];
rz(1.2081344) q[2];
rz(2.1900603) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(2.4212196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7893938) q[0];
sx q[0];
rz(-2.3645526) q[0];
sx q[0];
rz(-3.0598031) q[0];
rz(2.8261322) q[1];
sx q[1];
rz(-3.0888562) q[1];
sx q[1];
rz(1.2299406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30222505) q[0];
sx q[0];
rz(-1.4124083) q[0];
sx q[0];
rz(1.7296289) q[0];
x q[1];
rz(0.66218485) q[2];
sx q[2];
rz(-0.22394584) q[2];
sx q[2];
rz(2.0020773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1217864) q[1];
sx q[1];
rz(-0.14827327) q[1];
sx q[1];
rz(2.2231977) q[1];
rz(1.8043777) q[3];
sx q[3];
rz(-0.54334917) q[3];
sx q[3];
rz(-0.41239967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54084593) q[2];
sx q[2];
rz(-0.019286152) q[2];
sx q[2];
rz(1.1136327) q[2];
rz(-0.039637808) q[3];
sx q[3];
rz(-3.1318635) q[3];
sx q[3];
rz(-2.5447194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.5254163) q[1];
sx q[1];
rz(-2.0694852) q[1];
sx q[1];
rz(-2.9328666) q[1];
rz(1.3327333) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.0669591) q[2];
rz(0.039327903) q[3];
sx q[3];
rz(-1.3500742) q[3];
sx q[3];
rz(-3.1108656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
