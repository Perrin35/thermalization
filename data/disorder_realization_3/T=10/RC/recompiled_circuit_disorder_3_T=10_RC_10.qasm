OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(2.3303633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1404999) q[2];
sx q[2];
rz(-1.6506519) q[2];
sx q[2];
rz(-0.18730883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2990883) q[1];
sx q[1];
rz(-0.43049225) q[1];
sx q[1];
rz(1.7249785) q[1];
rz(-pi) q[2];
rz(-0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-2.6020715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182506) q[0];
sx q[0];
rz(-1.6174416) q[0];
sx q[0];
rz(-1.1211066) q[0];
rz(-2.4194854) q[2];
sx q[2];
rz(-1.2566084) q[2];
sx q[2];
rz(2.9177641) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19903781) q[1];
sx q[1];
rz(-2.6032762) q[1];
sx q[1];
rz(1.946279) q[1];
rz(-pi) q[2];
rz(0.55595056) q[3];
sx q[3];
rz(-2.2009938) q[3];
sx q[3];
rz(2.6549908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(-0.2562491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218131) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(0.42827423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0923549) q[2];
sx q[2];
rz(-1.1698327) q[2];
sx q[2];
rz(1.608633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2910034) q[1];
sx q[1];
rz(-2.2271419) q[1];
sx q[1];
rz(-2.7235051) q[1];
x q[2];
rz(2.933421) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(0.65073035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(2.5668872) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(-2.8821049) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-0.73192275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703092) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(2.1156103) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4243083) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(0.35536534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4983474) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(0.64756067) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22051208) q[3];
sx q[3];
rz(-1.6933428) q[3];
sx q[3];
rz(1.3361271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(-0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5916409) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(-0.79777065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3482462) q[0];
sx q[0];
rz(-0.12013809) q[0];
sx q[0];
rz(-1.5500463) q[0];
rz(-pi) q[1];
rz(-1.351379) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(2.9079633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.55802) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-0.11548345) q[1];
x q[2];
rz(1.9005152) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72567155) q[0];
sx q[0];
rz(-1.3999181) q[0];
sx q[0];
rz(-0.15539774) q[0];
x q[1];
rz(-2.1224535) q[2];
sx q[2];
rz(-0.86691228) q[2];
sx q[2];
rz(0.64955074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1706714) q[1];
sx q[1];
rz(-2.1340003) q[1];
sx q[1];
rz(-1.3794273) q[1];
rz(-pi) q[2];
rz(3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(-1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(0.0079356114) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82848362) q[0];
sx q[0];
rz(-0.61921739) q[0];
sx q[0];
rz(0.45066582) q[0];
x q[1];
rz(-0.43086149) q[2];
sx q[2];
rz(-2.2164946) q[2];
sx q[2];
rz(-0.25260392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.36564) q[1];
sx q[1];
rz(-2.0740293) q[1];
sx q[1];
rz(-0.788049) q[1];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(0.41637862) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84425981) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(2.8914333) q[0];
rz(-pi) q[1];
rz(-2.084923) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(-2.4664997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2545401) q[1];
sx q[1];
rz(-0.92003838) q[1];
sx q[1];
rz(2.410143) q[1];
rz(-1.7765462) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1736261) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(0.047441479) q[0];
rz(-pi) q[1];
rz(2.4531104) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(0.85306963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3560564) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(-2.9836125) q[1];
x q[2];
rz(0.74420332) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(-0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.8189925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(1.9508719) q[0];
rz(1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.10510124) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(2.3829616) q[1];
x q[2];
rz(-2.9880452) q[3];
sx q[3];
rz(-0.921075) q[3];
sx q[3];
rz(-2.6447907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(1.580711) q[2];
sx q[2];
rz(-1.0615988) q[2];
sx q[2];
rz(1.4524729) q[2];
rz(1.9382523) q[3];
sx q[3];
rz(-2.8946579) q[3];
sx q[3];
rz(-2.2263087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
