OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(-2.6383658) q[0];
sx q[0];
rz(0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0137579) q[0];
sx q[0];
rz(-0.36359596) q[0];
sx q[0];
rz(-0.6291997) q[0];
rz(-pi) q[1];
rz(2.7230524) q[2];
sx q[2];
rz(-1.6633908) q[2];
sx q[2];
rz(-1.6469524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9468294) q[1];
sx q[1];
rz(-1.0461055) q[1];
sx q[1];
rz(-2.8835433) q[1];
rz(-pi) q[2];
rz(3.1078994) q[3];
sx q[3];
rz(-1.1564621) q[3];
sx q[3];
rz(1.3106443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.2215866) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-0.63175732) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97665635) q[0];
sx q[0];
rz(-1.5979842) q[0];
sx q[0];
rz(1.8261766) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9594876) q[2];
sx q[2];
rz(-1.9027862) q[2];
sx q[2];
rz(1.7322844) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.020097453) q[1];
sx q[1];
rz(-1.1250245) q[1];
sx q[1];
rz(-0.14264588) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6173108) q[3];
sx q[3];
rz(-0.97913137) q[3];
sx q[3];
rz(-0.072629645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0624861) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(0.24965723) q[2];
rz(-2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.867884) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0262336) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(-0.94985234) q[0];
rz(-pi) q[1];
rz(1.6689698) q[2];
sx q[2];
rz(-1.86491) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.102036) q[1];
sx q[1];
rz(-1.9196871) q[1];
sx q[1];
rz(2.4400913) q[1];
rz(-pi) q[2];
rz(-0.7234296) q[3];
sx q[3];
rz(-1.1838786) q[3];
sx q[3];
rz(-2.5748411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(2.188142) q[2];
rz(-3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-3.0674556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1124681) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(2.8549457) q[0];
rz(0.3481625) q[2];
sx q[2];
rz(-2.2194038) q[2];
sx q[2];
rz(-1.5319038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6756235) q[1];
sx q[1];
rz(-0.33756653) q[1];
sx q[1];
rz(-1.9338495) q[1];
rz(1.3965963) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(-3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-0.038643535) q[2];
rz(0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(1.3624396) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.1626676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029318854) q[0];
sx q[0];
rz(-0.11538878) q[0];
sx q[0];
rz(1.3089887) q[0];
rz(-pi) q[1];
rz(2.9722399) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(-2.7404075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7971024) q[1];
sx q[1];
rz(-1.7624197) q[1];
sx q[1];
rz(-1.9413129) q[1];
rz(-pi) q[2];
rz(2.2546222) q[3];
sx q[3];
rz(-0.63347048) q[3];
sx q[3];
rz(2.5630643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5148619) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(2.999372) q[2];
rz(-2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(2.8335559) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7361569) q[0];
sx q[0];
rz(-1.476256) q[0];
sx q[0];
rz(-0.12673881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61200895) q[2];
sx q[2];
rz(-1.0527805) q[2];
sx q[2];
rz(0.81105622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6412515) q[1];
sx q[1];
rz(-0.89741035) q[1];
sx q[1];
rz(-0.25232368) q[1];
x q[2];
rz(-2.8940053) q[3];
sx q[3];
rz(-2.7831804) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797453) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(0.67108363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2790047) q[0];
sx q[0];
rz(-3.0775078) q[0];
sx q[0];
rz(-2.5628753) q[0];
x q[1];
rz(1.3865878) q[2];
sx q[2];
rz(-1.5407908) q[2];
sx q[2];
rz(1.9629994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96324254) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(-0.11022719) q[1];
x q[2];
rz(2.2460849) q[3];
sx q[3];
rz(-1.5338147) q[3];
sx q[3];
rz(1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.122763) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(3.1130062) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.2602497) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(1.9746045) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056571) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(-0.94353326) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86161676) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(2.784563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1046909) q[1];
sx q[1];
rz(-2.9505886) q[1];
sx q[1];
rz(0.16049338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(-0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.6519201) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4091858) q[0];
sx q[0];
rz(-0.86940765) q[0];
sx q[0];
rz(2.5138445) q[0];
x q[1];
rz(-1.0987368) q[2];
sx q[2];
rz(-1.5257116) q[2];
sx q[2];
rz(-0.91143196) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97396321) q[1];
sx q[1];
rz(-1.3988252) q[1];
sx q[1];
rz(-0.19584882) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1065815) q[3];
sx q[3];
rz(-2.0097369) q[3];
sx q[3];
rz(-0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.757471) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(2.8531895) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-2.9945701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4924016) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(-2.6834821) q[0];
x q[1];
rz(-1.64098) q[2];
sx q[2];
rz(-1.3314221) q[2];
sx q[2];
rz(1.0603051) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.920267) q[1];
sx q[1];
rz(-0.95278554) q[1];
sx q[1];
rz(3.0926535) q[1];
rz(0.23707323) q[3];
sx q[3];
rz(-0.58963886) q[3];
sx q[3];
rz(-1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(2.384281) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(3.0030737) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(2.1327303) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
