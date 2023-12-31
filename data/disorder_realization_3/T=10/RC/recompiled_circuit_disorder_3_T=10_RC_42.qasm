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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(2.3303633) q[0];
rz(-pi) q[1];
rz(1.7180874) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(1.6342083) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.011885) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(3.0711864) q[1];
rz(-pi) q[2];
rz(-0.36177735) q[3];
sx q[3];
rz(-2.8594115) q[3];
sx q[3];
rz(2.2892945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(0.52662915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3027705) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-2.6020715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(3.0898068) q[0];
x q[1];
rz(-2.6846634) q[2];
sx q[2];
rz(-2.3655342) q[2];
sx q[2];
rz(-1.6844144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0454228) q[1];
sx q[1];
rz(-1.3816557) q[1];
sx q[1];
rz(1.063709) q[1];
rz(-pi) q[2];
rz(0.86124729) q[3];
sx q[3];
rz(-2.0112787) q[3];
sx q[3];
rz(-1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-0.2562491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(1.8947253) q[0];
x q[1];
rz(-2.696064) q[2];
sx q[2];
rz(-2.0085213) q[2];
sx q[2];
rz(2.90403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5970522) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(-0.87045963) q[1];
x q[2];
rz(2.933421) q[3];
sx q[3];
rz(-0.18067193) q[3];
sx q[3];
rz(2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8024575) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-0.25948778) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(2.4096699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5607802) q[0];
sx q[0];
rz(-2.4111351) q[0];
sx q[0];
rz(0.74303445) q[0];
x q[1];
rz(3.0078997) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(-1.114811) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1720097) q[1];
sx q[1];
rz(-1.3130377) q[1];
sx q[1];
rz(1.7715363) q[1];
rz(-0.22051208) q[3];
sx q[3];
rz(-1.4482499) q[3];
sx q[3];
rz(-1.8054655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-2.990492) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20194963) q[0];
sx q[0];
rz(-1.5683096) q[0];
sx q[0];
rz(-1.4506838) q[0];
x q[1];
rz(-1.351379) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(2.9079633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3623912) q[1];
sx q[1];
rz(-0.63450846) q[1];
sx q[1];
rz(-1.7290551) q[1];
rz(-pi) q[2];
rz(-0.62545125) q[3];
sx q[3];
rz(-0.52895412) q[3];
sx q[3];
rz(1.7367425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4159211) q[0];
sx q[0];
rz(-1.3999181) q[0];
sx q[0];
rz(2.9861949) q[0];
rz(-pi) q[1];
rz(2.5885133) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(-1.4097708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-pi/2) q[1];
sx q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.4403884) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(-0.0079356114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7769653) q[0];
sx q[0];
rz(-2.1205175) q[0];
sx q[0];
rz(1.2697898) q[0];
rz(-pi) q[1];
rz(0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(1.0483339) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4850033) q[1];
sx q[1];
rz(-0.90066972) q[1];
sx q[1];
rz(0.90799241) q[1];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.5023124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72909268) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(-1.5604707) q[0];
x q[1];
rz(2.6697568) q[2];
sx q[2];
rz(-2.2484357) q[2];
sx q[2];
rz(0.0080646947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9598436) q[1];
sx q[1];
rz(-2.1310924) q[1];
sx q[1];
rz(-2.3676141) q[1];
x q[2];
rz(-3.0495166) q[3];
sx q[3];
rz(-1.3658938) q[3];
sx q[3];
rz(-2.3931707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(0.41697821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0957085) q[0];
sx q[0];
rz(-0.65512148) q[0];
sx q[0];
rz(1.5089633) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99879361) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(-1.7172161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13874395) q[1];
sx q[1];
rz(-1.7094304) q[1];
sx q[1];
rz(-2.0744051) q[1];
rz(-pi) q[2];
rz(-3.0005089) q[3];
sx q[3];
rz(-1.4420296) q[3];
sx q[3];
rz(1.5537378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.3226002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.246829) q[0];
sx q[0];
rz(-1.6920977) q[0];
sx q[0];
rz(1.9508719) q[0];
rz(0.25090353) q[2];
sx q[2];
rz(-1.13399) q[2];
sx q[2];
rz(-1.8032339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1111856) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(0.91099693) q[1];
rz(-pi) q[2];
rz(2.9880452) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-2.6447907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14324698) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.3636419) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-0.017756391) q[2];
sx q[2];
rz(-0.50928558) q[2];
sx q[2];
rz(-1.7094564) q[2];
rz(-1.3397459) q[3];
sx q[3];
rz(-1.6587202) q[3];
sx q[3];
rz(-0.29826577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
