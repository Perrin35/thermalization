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
rz(-1.8259814) q[0];
sx q[0];
rz(-2.313518) q[0];
sx q[0];
rz(-2.2928884) q[0];
rz(1.4915713) q[1];
sx q[1];
rz(-0.68869156) q[1];
sx q[1];
rz(-2.7149849) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899982) q[0];
sx q[0];
rz(-1.6713872) q[0];
sx q[0];
rz(-1.4949343) q[0];
rz(-pi) q[1];
rz(1.6831065) q[2];
sx q[2];
rz(-1.6717615) q[2];
sx q[2];
rz(-1.2580521) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22151193) q[1];
sx q[1];
rz(-0.4512085) q[1];
sx q[1];
rz(2.9738704) q[1];
x q[2];
rz(-2.1900666) q[3];
sx q[3];
rz(-0.87267704) q[3];
sx q[3];
rz(-1.8474471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9244869) q[2];
sx q[2];
rz(-1.7542398) q[2];
sx q[2];
rz(0.039487751) q[2];
rz(-2.110179) q[3];
sx q[3];
rz(-1.02905) q[3];
sx q[3];
rz(-1.903681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6295488) q[0];
sx q[0];
rz(-1.5271674) q[0];
sx q[0];
rz(1.7426096) q[0];
rz(1.6276739) q[1];
sx q[1];
rz(-2.2986423) q[1];
sx q[1];
rz(-1.7248076) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91367224) q[0];
sx q[0];
rz(-2.3487072) q[0];
sx q[0];
rz(2.8280054) q[0];
rz(-pi) q[1];
rz(-2.879989) q[2];
sx q[2];
rz(-0.79643476) q[2];
sx q[2];
rz(1.4858044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.09782) q[1];
sx q[1];
rz(-2.7209681) q[1];
sx q[1];
rz(2.8021977) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6187158) q[3];
sx q[3];
rz(-1.1635492) q[3];
sx q[3];
rz(-1.6919235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17025718) q[2];
sx q[2];
rz(-1.1032666) q[2];
sx q[2];
rz(0.7106759) q[2];
rz(0.014287861) q[3];
sx q[3];
rz(-2.894214) q[3];
sx q[3];
rz(-1.8054731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22064848) q[0];
sx q[0];
rz(-2.1823688) q[0];
sx q[0];
rz(-2.7253819) q[0];
rz(-1.2843708) q[1];
sx q[1];
rz(-1.6651848) q[1];
sx q[1];
rz(-0.31825569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46894761) q[0];
sx q[0];
rz(-1.8469317) q[0];
sx q[0];
rz(1.1731082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2936965) q[2];
sx q[2];
rz(-1.0506786) q[2];
sx q[2];
rz(2.0682356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(15/(11*pi)) q[1];
sx q[1];
rz(-0.66898443) q[1];
sx q[1];
rz(0.21978746) q[1];
rz(-pi) q[2];
rz(2.2417547) q[3];
sx q[3];
rz(-2.7551015) q[3];
sx q[3];
rz(0.46903601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2661813) q[2];
sx q[2];
rz(-2.3616932) q[2];
sx q[2];
rz(1.2904588) q[2];
rz(-0.053704638) q[3];
sx q[3];
rz(-1.8219681) q[3];
sx q[3];
rz(1.3857589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5829492) q[0];
sx q[0];
rz(-2.4531893) q[0];
sx q[0];
rz(-1.0097367) q[0];
rz(2.2263777) q[1];
sx q[1];
rz(-1.5292294) q[1];
sx q[1];
rz(-2.7161652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4998326) q[0];
sx q[0];
rz(-1.7374037) q[0];
sx q[0];
rz(0.84979041) q[0];
rz(-pi) q[1];
rz(3.0680823) q[2];
sx q[2];
rz(-2.3125154) q[2];
sx q[2];
rz(2.4501767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.912251) q[1];
sx q[1];
rz(-2.3315363) q[1];
sx q[1];
rz(-0.99397578) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5365547) q[3];
sx q[3];
rz(-2.9786907) q[3];
sx q[3];
rz(2.3690641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90011826) q[2];
sx q[2];
rz(-1.5510617) q[2];
sx q[2];
rz(0.11108622) q[2];
rz(0.19242081) q[3];
sx q[3];
rz(-2.7399053) q[3];
sx q[3];
rz(-2.4470611) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38839328) q[0];
sx q[0];
rz(-2.41112) q[0];
sx q[0];
rz(-2.2764192) q[0];
rz(0.49258891) q[1];
sx q[1];
rz(-2.045423) q[1];
sx q[1];
rz(1.7561087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139112) q[0];
sx q[0];
rz(-1.7075141) q[0];
sx q[0];
rz(-0.78929957) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4200468) q[2];
sx q[2];
rz(-2.2873452) q[2];
sx q[2];
rz(-1.6796215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.046987586) q[1];
sx q[1];
rz(-1.054227) q[1];
sx q[1];
rz(2.0265685) q[1];
rz(-pi) q[2];
rz(2.4597857) q[3];
sx q[3];
rz(-2.1958133) q[3];
sx q[3];
rz(2.8375219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2546786) q[2];
sx q[2];
rz(-2.6723537) q[2];
sx q[2];
rz(-2.4349507) q[2];
rz(0.32736579) q[3];
sx q[3];
rz(-1.8492161) q[3];
sx q[3];
rz(-2.4842026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3051598) q[0];
sx q[0];
rz(-1.3494455) q[0];
sx q[0];
rz(2.7850372) q[0];
rz(2.5794079) q[1];
sx q[1];
rz(-2.0088582) q[1];
sx q[1];
rz(-2.1551989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1323469) q[0];
sx q[0];
rz(-2.2725687) q[0];
sx q[0];
rz(2.470507) q[0];
x q[1];
rz(-2.3190034) q[2];
sx q[2];
rz(-2.6554567) q[2];
sx q[2];
rz(-1.2631455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3918009) q[1];
sx q[1];
rz(-1.5967249) q[1];
sx q[1];
rz(0.75849979) q[1];
x q[2];
rz(2.2067542) q[3];
sx q[3];
rz(-2.842111) q[3];
sx q[3];
rz(-0.071223473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76546136) q[2];
sx q[2];
rz(-0.71530801) q[2];
sx q[2];
rz(-0.7005271) q[2];
rz(0.96945196) q[3];
sx q[3];
rz(-2.1078096) q[3];
sx q[3];
rz(0.9303003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9483865) q[0];
sx q[0];
rz(-0.32873118) q[0];
sx q[0];
rz(2.9902003) q[0];
rz(-1.6311749) q[1];
sx q[1];
rz(-2.0163586) q[1];
sx q[1];
rz(1.6815394) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4497283) q[0];
sx q[0];
rz(-0.98694333) q[0];
sx q[0];
rz(-2.1644781) q[0];
rz(-2.1489819) q[2];
sx q[2];
rz(-2.5343203) q[2];
sx q[2];
rz(-0.36888514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55893873) q[1];
sx q[1];
rz(-2.451921) q[1];
sx q[1];
rz(2.2037872) q[1];
x q[2];
rz(-2.2112956) q[3];
sx q[3];
rz(-2.9004164) q[3];
sx q[3];
rz(-1.1080139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90525642) q[2];
sx q[2];
rz(-0.6826123) q[2];
sx q[2];
rz(-0.35161099) q[2];
rz(-2.6998399) q[3];
sx q[3];
rz(-2.2124002) q[3];
sx q[3];
rz(-0.15171224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.041895954) q[0];
sx q[0];
rz(-1.2956887) q[0];
sx q[0];
rz(-1.1610485) q[0];
rz(2.4940122) q[1];
sx q[1];
rz(-1.3787965) q[1];
sx q[1];
rz(2.3191998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3352584) q[0];
sx q[0];
rz(-1.0653235) q[0];
sx q[0];
rz(2.6837206) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0089125) q[2];
sx q[2];
rz(-0.88399502) q[2];
sx q[2];
rz(-2.7504454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1668561) q[1];
sx q[1];
rz(-2.3212163) q[1];
sx q[1];
rz(-1.9700762) q[1];
rz(-0.83637442) q[3];
sx q[3];
rz(-0.55787239) q[3];
sx q[3];
rz(1.8250067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9245727) q[2];
sx q[2];
rz(-1.1447039) q[2];
sx q[2];
rz(1.8355969) q[2];
rz(0.71803391) q[3];
sx q[3];
rz(-0.82457232) q[3];
sx q[3];
rz(0.048862783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.43593916) q[0];
sx q[0];
rz(-1.0834563) q[0];
sx q[0];
rz(-3.1245681) q[0];
rz(1.1964993) q[1];
sx q[1];
rz(-2.7462609) q[1];
sx q[1];
rz(-2.8065525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3075313) q[0];
sx q[0];
rz(-1.002822) q[0];
sx q[0];
rz(-0.012004367) q[0];
rz(-0.52729221) q[2];
sx q[2];
rz(-0.71374536) q[2];
sx q[2];
rz(1.1251768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10427231) q[1];
sx q[1];
rz(-1.7591216) q[1];
sx q[1];
rz(-0.89481797) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7850168) q[3];
sx q[3];
rz(-0.64006539) q[3];
sx q[3];
rz(2.9736116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7177141) q[2];
sx q[2];
rz(-0.7462036) q[2];
sx q[2];
rz(-2.8625873) q[2];
rz(-1.3689857) q[3];
sx q[3];
rz(-1.5149346) q[3];
sx q[3];
rz(0.92931187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245215) q[0];
sx q[0];
rz(-0.14685024) q[0];
sx q[0];
rz(-0.76706925) q[0];
rz(2.7686139) q[1];
sx q[1];
rz(-1.6423128) q[1];
sx q[1];
rz(0.17072089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883461) q[0];
sx q[0];
rz(-2.4062254) q[0];
sx q[0];
rz(3.0985996) q[0];
rz(1.7319948) q[2];
sx q[2];
rz(-1.1606163) q[2];
sx q[2];
rz(-2.8267677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3633299) q[1];
sx q[1];
rz(-0.40503854) q[1];
sx q[1];
rz(0.89680775) q[1];
rz(-2.6737718) q[3];
sx q[3];
rz(-1.7138283) q[3];
sx q[3];
rz(0.19090548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3501222) q[2];
sx q[2];
rz(-0.56362027) q[2];
sx q[2];
rz(-0.0027837022) q[2];
rz(-0.9015829) q[3];
sx q[3];
rz(-1.0222579) q[3];
sx q[3];
rz(2.3367052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881385) q[0];
sx q[0];
rz(-2.8202941) q[0];
sx q[0];
rz(1.1711076) q[0];
rz(2.3114655) q[1];
sx q[1];
rz(-1.3273888) q[1];
sx q[1];
rz(-1.8524016) q[1];
rz(-2.6302677) q[2];
sx q[2];
rz(-1.0228538) q[2];
sx q[2];
rz(2.6342837) q[2];
rz(-0.8030025) q[3];
sx q[3];
rz(-2.0661125) q[3];
sx q[3];
rz(-2.6281602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
