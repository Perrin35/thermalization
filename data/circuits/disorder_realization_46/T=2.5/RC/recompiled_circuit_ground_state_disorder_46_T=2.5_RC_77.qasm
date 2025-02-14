OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1281066) q[0];
sx q[0];
rz(3.9689316) q[0];
sx q[0];
rz(7.6710424) q[0];
rz(0.85343051) q[1];
sx q[1];
rz(-0.4018468) q[1];
sx q[1];
rz(2.0274577) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83682293) q[0];
sx q[0];
rz(-2.3802505) q[0];
sx q[0];
rz(-1.9403275) q[0];
rz(-pi) q[1];
rz(-2.4365455) q[2];
sx q[2];
rz(-1.4234666) q[2];
sx q[2];
rz(1.6697869) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88190313) q[1];
sx q[1];
rz(-1.5756541) q[1];
sx q[1];
rz(1.5905981) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2306289) q[3];
sx q[3];
rz(-2.0475884) q[3];
sx q[3];
rz(-1.7309675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16093291) q[2];
sx q[2];
rz(-2.7200343) q[2];
sx q[2];
rz(-1.6932311) q[2];
rz(2.8990922) q[3];
sx q[3];
rz(-2.226053) q[3];
sx q[3];
rz(-1.8814794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777957) q[0];
sx q[0];
rz(-1.3160492) q[0];
sx q[0];
rz(-1.0003566) q[0];
rz(2.4099804) q[1];
sx q[1];
rz(-1.4294521) q[1];
sx q[1];
rz(-1.1801205) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6574616) q[0];
sx q[0];
rz(-3.0041782) q[0];
sx q[0];
rz(0.71533393) q[0];
rz(-pi) q[1];
rz(1.2928455) q[2];
sx q[2];
rz(-1.6033844) q[2];
sx q[2];
rz(-1.1162835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0457256) q[1];
sx q[1];
rz(-1.5689227) q[1];
sx q[1];
rz(1.7049403) q[1];
x q[2];
rz(2.3445156) q[3];
sx q[3];
rz(-0.38603544) q[3];
sx q[3];
rz(-0.71938709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0483027) q[2];
sx q[2];
rz(-2.2403109) q[2];
sx q[2];
rz(-1.0789336) q[2];
rz(-0.087873936) q[3];
sx q[3];
rz(-1.352997) q[3];
sx q[3];
rz(-2.8972076) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3946149) q[0];
sx q[0];
rz(-1.9786388) q[0];
sx q[0];
rz(1.8633457) q[0];
rz(-0.46319115) q[1];
sx q[1];
rz(-1.7852424) q[1];
sx q[1];
rz(2.3203704) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44959497) q[0];
sx q[0];
rz(-2.2173081) q[0];
sx q[0];
rz(-2.9086935) q[0];
x q[1];
rz(1.6326007) q[2];
sx q[2];
rz(-1.9414177) q[2];
sx q[2];
rz(0.48423094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6602409) q[1];
sx q[1];
rz(-1.7180182) q[1];
sx q[1];
rz(0.46105701) q[1];
rz(-pi) q[2];
rz(2.4850003) q[3];
sx q[3];
rz(-0.87513211) q[3];
sx q[3];
rz(-0.10564239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2440196) q[2];
sx q[2];
rz(-1.5622746) q[2];
sx q[2];
rz(1.0025586) q[2];
rz(1.9133441) q[3];
sx q[3];
rz(-1.4017665) q[3];
sx q[3];
rz(-1.7206515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9035852) q[0];
sx q[0];
rz(-2.0915732) q[0];
sx q[0];
rz(-0.28582698) q[0];
rz(-2.8803275) q[1];
sx q[1];
rz(-1.3328054) q[1];
sx q[1];
rz(3.1108943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54346976) q[0];
sx q[0];
rz(-1.3075519) q[0];
sx q[0];
rz(-1.4953161) q[0];
rz(-pi) q[1];
rz(-0.89392406) q[2];
sx q[2];
rz(-1.6325762) q[2];
sx q[2];
rz(0.0047574818) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0847807) q[1];
sx q[1];
rz(-1.1263945) q[1];
sx q[1];
rz(-0.2803352) q[1];
rz(-pi) q[2];
rz(-1.2571583) q[3];
sx q[3];
rz(-1.4935857) q[3];
sx q[3];
rz(1.9208637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8135445) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(3.1385341) q[2];
rz(-2.3004153) q[3];
sx q[3];
rz(-0.58924651) q[3];
sx q[3];
rz(0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641814) q[0];
sx q[0];
rz(-2.2941636) q[0];
sx q[0];
rz(-1.6743976) q[0];
rz(-1.6293619) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(-2.1893952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6124355) q[0];
sx q[0];
rz(-1.5780492) q[0];
sx q[0];
rz(-3.0713505) q[0];
rz(-pi) q[1];
rz(-1.9604076) q[2];
sx q[2];
rz(-0.91105538) q[2];
sx q[2];
rz(2.8174741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84633652) q[1];
sx q[1];
rz(-1.6115115) q[1];
sx q[1];
rz(1.4671007) q[1];
x q[2];
rz(0.85785474) q[3];
sx q[3];
rz(-2.3393173) q[3];
sx q[3];
rz(-1.0209393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3516922) q[2];
sx q[2];
rz(-1.60195) q[2];
sx q[2];
rz(-2.6969625) q[2];
rz(0.45186684) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(-3.0795081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43279466) q[0];
sx q[0];
rz(-2.4169156) q[0];
sx q[0];
rz(0.92025796) q[0];
rz(0.96011773) q[1];
sx q[1];
rz(-1.7476387) q[1];
sx q[1];
rz(0.49930176) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647674) q[0];
sx q[0];
rz(-1.9070325) q[0];
sx q[0];
rz(-1.7306843) q[0];
rz(2.9174706) q[2];
sx q[2];
rz(-2.4680228) q[2];
sx q[2];
rz(2.3526255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0681794) q[1];
sx q[1];
rz(-1.2559895) q[1];
sx q[1];
rz(0.59609658) q[1];
rz(-pi) q[2];
rz(2.4003827) q[3];
sx q[3];
rz(-1.831372) q[3];
sx q[3];
rz(-1.2737361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7509191) q[2];
sx q[2];
rz(-1.5385188) q[2];
sx q[2];
rz(2.1092559) q[2];
rz(-2.27683) q[3];
sx q[3];
rz(-1.4356177) q[3];
sx q[3];
rz(0.56149703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.432935) q[0];
sx q[0];
rz(-1.9517887) q[0];
sx q[0];
rz(1.7565961) q[0];
rz(0.9224836) q[1];
sx q[1];
rz(-1.9360767) q[1];
sx q[1];
rz(-2.9965957) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3475048) q[0];
sx q[0];
rz(-1.4055168) q[0];
sx q[0];
rz(1.7519622) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.112193) q[2];
sx q[2];
rz(-2.2839632) q[2];
sx q[2];
rz(0.61427639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.662949) q[1];
sx q[1];
rz(-1.4976209) q[1];
sx q[1];
rz(-3.1107145) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0078405) q[3];
sx q[3];
rz(-2.2641695) q[3];
sx q[3];
rz(-1.0348606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21236913) q[2];
sx q[2];
rz(-0.23304686) q[2];
sx q[2];
rz(1.93335) q[2];
rz(2.3182747) q[3];
sx q[3];
rz(-1.1813141) q[3];
sx q[3];
rz(1.3594782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061148297) q[0];
sx q[0];
rz(-2.3539982) q[0];
sx q[0];
rz(-2.0620692) q[0];
rz(0.98006717) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(-3.0316839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11342587) q[0];
sx q[0];
rz(-0.47385212) q[0];
sx q[0];
rz(-2.4555311) q[0];
x q[1];
rz(1.5797516) q[2];
sx q[2];
rz(-1.2924465) q[2];
sx q[2];
rz(3.0305733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3802783) q[1];
sx q[1];
rz(-1.2732693) q[1];
sx q[1];
rz(0.44400906) q[1];
x q[2];
rz(-1.1989459) q[3];
sx q[3];
rz(-0.56104413) q[3];
sx q[3];
rz(-2.4435465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93426934) q[2];
sx q[2];
rz(-1.0730275) q[2];
sx q[2];
rz(-2.4369924) q[2];
rz(2.8568824) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(-2.708191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339742) q[0];
sx q[0];
rz(-1.3173137) q[0];
sx q[0];
rz(-0.91868573) q[0];
rz(0.48565117) q[1];
sx q[1];
rz(-0.44356569) q[1];
sx q[1];
rz(2.0408911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1087295) q[0];
sx q[0];
rz(-2.9447365) q[0];
sx q[0];
rz(0.78702505) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2186949) q[2];
sx q[2];
rz(-0.37247286) q[2];
sx q[2];
rz(3.0650919) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6599258) q[1];
sx q[1];
rz(-1.1705913) q[1];
sx q[1];
rz(3.0277962) q[1];
x q[2];
rz(2.5060008) q[3];
sx q[3];
rz(-1.6219536) q[3];
sx q[3];
rz(3.1380461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47621581) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(0.7594792) q[2];
rz(2.1935943) q[3];
sx q[3];
rz(-1.4576603) q[3];
sx q[3];
rz(-1.52004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043902472) q[0];
sx q[0];
rz(-2.6444785) q[0];
sx q[0];
rz(-0.33101606) q[0];
rz(-2.379592) q[1];
sx q[1];
rz(-0.29251978) q[1];
sx q[1];
rz(-2.5133572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8549919) q[0];
sx q[0];
rz(-1.0074662) q[0];
sx q[0];
rz(1.2391633) q[0];
rz(-pi) q[1];
rz(0.64543311) q[2];
sx q[2];
rz(-2.7629768) q[2];
sx q[2];
rz(0.23832527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8208295) q[1];
sx q[1];
rz(-2.4270128) q[1];
sx q[1];
rz(-2.5933867) q[1];
rz(-pi) q[2];
rz(-0.68898983) q[3];
sx q[3];
rz(-0.51577079) q[3];
sx q[3];
rz(-2.6067579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.021775333) q[2];
sx q[2];
rz(-2.6675318) q[2];
sx q[2];
rz(2.9277756) q[2];
rz(2.1553433) q[3];
sx q[3];
rz(-1.8503559) q[3];
sx q[3];
rz(3.0333062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4719791) q[0];
sx q[0];
rz(-1.4431974) q[0];
sx q[0];
rz(1.6786014) q[0];
rz(1.9646473) q[1];
sx q[1];
rz(-1.4851478) q[1];
sx q[1];
rz(-3.1335395) q[1];
rz(-0.17433737) q[2];
sx q[2];
rz(-1.871289) q[2];
sx q[2];
rz(1.540779) q[2];
rz(-0.10569345) q[3];
sx q[3];
rz(-0.65423818) q[3];
sx q[3];
rz(-0.26442179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
