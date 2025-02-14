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
rz(-2.6900238) q[0];
sx q[0];
rz(-1.4181925) q[0];
sx q[0];
rz(-0.42887846) q[0];
rz(2.9984527) q[1];
sx q[1];
rz(-1.0909456) q[1];
sx q[1];
rz(-0.080088869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1096261) q[0];
sx q[0];
rz(-1.6499777) q[0];
sx q[0];
rz(-1.5840785) q[0];
x q[1];
rz(1.8378788) q[2];
sx q[2];
rz(-0.3729698) q[2];
sx q[2];
rz(-2.3139062) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3428719) q[1];
sx q[1];
rz(-1.8340543) q[1];
sx q[1];
rz(-0.85204551) q[1];
x q[2];
rz(0.57755263) q[3];
sx q[3];
rz(-1.7087382) q[3];
sx q[3];
rz(0.70226125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.074778883) q[2];
sx q[2];
rz(-2.0902233) q[2];
sx q[2];
rz(1.7195513) q[2];
rz(1.7285796) q[3];
sx q[3];
rz(-1.6001817) q[3];
sx q[3];
rz(1.448267) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220575) q[0];
sx q[0];
rz(-1.7449361) q[0];
sx q[0];
rz(-2.8080217) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(-2.3894892) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81588868) q[0];
sx q[0];
rz(-2.3518343) q[0];
sx q[0];
rz(2.2237097) q[0];
rz(-pi) q[1];
rz(-3.0235748) q[2];
sx q[2];
rz(-2.9213597) q[2];
sx q[2];
rz(1.3853467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57422728) q[1];
sx q[1];
rz(-0.53814473) q[1];
sx q[1];
rz(-1.3567011) q[1];
x q[2];
rz(2.8798298) q[3];
sx q[3];
rz(-1.7883375) q[3];
sx q[3];
rz(2.1037773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0636474) q[2];
sx q[2];
rz(-0.91219488) q[2];
sx q[2];
rz(2.4845541) q[2];
rz(0.71197236) q[3];
sx q[3];
rz(-0.65198055) q[3];
sx q[3];
rz(2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.7065358) q[0];
sx q[0];
rz(-0.61921316) q[0];
sx q[0];
rz(-1.0414498) q[0];
rz(-0.374818) q[1];
sx q[1];
rz(-1.638214) q[1];
sx q[1];
rz(-2.5846438) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14675635) q[0];
sx q[0];
rz(-2.2731011) q[0];
sx q[0];
rz(1.102314) q[0];
x q[1];
rz(-1.416549) q[2];
sx q[2];
rz(-2.6852312) q[2];
sx q[2];
rz(0.92044324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8136602) q[1];
sx q[1];
rz(-2.2371462) q[1];
sx q[1];
rz(-1.7009237) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65105533) q[3];
sx q[3];
rz(-1.3811033) q[3];
sx q[3];
rz(2.9168791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(2.7064145) q[2];
rz(-2.6168881) q[3];
sx q[3];
rz(-0.33354959) q[3];
sx q[3];
rz(0.13478002) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(-1.6271628) q[0];
rz(1.0487522) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.4357766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.393261) q[0];
sx q[0];
rz(-0.18042856) q[0];
sx q[0];
rz(0.24112757) q[0];
rz(-2.5982791) q[2];
sx q[2];
rz(-1.2676953) q[2];
sx q[2];
rz(1.4599279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.857261) q[1];
sx q[1];
rz(-1.5067717) q[1];
sx q[1];
rz(-1.048719) q[1];
rz(2.0346927) q[3];
sx q[3];
rz(-1.8827065) q[3];
sx q[3];
rz(-0.31181617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.783739) q[2];
sx q[2];
rz(0.025156585) q[2];
rz(-1.6545506) q[3];
sx q[3];
rz(-0.39792037) q[3];
sx q[3];
rz(0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0461034) q[0];
sx q[0];
rz(-2.4920721) q[0];
sx q[0];
rz(2.982614) q[0];
rz(-2.036463) q[1];
sx q[1];
rz(-2.4159894) q[1];
sx q[1];
rz(-2.9172858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858356) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(2.0465536) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0105466) q[2];
sx q[2];
rz(-1.5429826) q[2];
sx q[2];
rz(1.5956439) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4558002) q[1];
sx q[1];
rz(-0.53736254) q[1];
sx q[1];
rz(-1.3440029) q[1];
rz(-2.5644825) q[3];
sx q[3];
rz(-1.464586) q[3];
sx q[3];
rz(2.4662227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9362681) q[2];
sx q[2];
rz(-1.779413) q[2];
sx q[2];
rz(0.74860191) q[2];
rz(-0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(0.41142472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1614138) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(2.4112716) q[0];
rz(-1.6397363) q[1];
sx q[1];
rz(-1.390099) q[1];
sx q[1];
rz(-2.9072442) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47539179) q[0];
sx q[0];
rz(-2.0115543) q[0];
sx q[0];
rz(2.4077364) q[0];
x q[1];
rz(-2.5426834) q[2];
sx q[2];
rz(-1.0157758) q[2];
sx q[2];
rz(2.2428494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1533656) q[1];
sx q[1];
rz(-1.4082137) q[1];
sx q[1];
rz(-1.52262) q[1];
x q[2];
rz(-1.9803011) q[3];
sx q[3];
rz(-0.1529158) q[3];
sx q[3];
rz(-0.59329882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.948287) q[2];
sx q[2];
rz(-2.327658) q[2];
sx q[2];
rz(1.8602271) q[2];
rz(2.5365601) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(-1.9090778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.989711) q[0];
sx q[0];
rz(-1.2622156) q[0];
sx q[0];
rz(2.5285517) q[0];
rz(2.4609861) q[1];
sx q[1];
rz(-1.1111518) q[1];
sx q[1];
rz(-1.8162762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8794358) q[0];
sx q[0];
rz(-0.61816321) q[0];
sx q[0];
rz(-2.6690041) q[0];
rz(-2.9284186) q[2];
sx q[2];
rz(-1.6521378) q[2];
sx q[2];
rz(0.48254642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5449574) q[1];
sx q[1];
rz(-1.6973293) q[1];
sx q[1];
rz(2.5759349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.16991) q[3];
sx q[3];
rz(-2.1155778) q[3];
sx q[3];
rz(-1.137103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1943872) q[2];
sx q[2];
rz(-1.2014061) q[2];
sx q[2];
rz(-2.3335333) q[2];
rz(-2.4957073) q[3];
sx q[3];
rz(-2.8736727) q[3];
sx q[3];
rz(-2.8382235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5757669) q[0];
sx q[0];
rz(-2.4149826) q[0];
sx q[0];
rz(-0.44276825) q[0];
rz(2.5495461) q[1];
sx q[1];
rz(-0.7213842) q[1];
sx q[1];
rz(0.17568849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21766996) q[0];
sx q[0];
rz(-0.45638535) q[0];
sx q[0];
rz(-0.66500591) q[0];
rz(-pi) q[1];
rz(0.4515516) q[2];
sx q[2];
rz(-1.6066666) q[2];
sx q[2];
rz(0.6605688) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86398839) q[1];
sx q[1];
rz(-2.2114807) q[1];
sx q[1];
rz(-1.5169657) q[1];
rz(-pi) q[2];
rz(2.1814303) q[3];
sx q[3];
rz(-1.9879793) q[3];
sx q[3];
rz(-0.3822592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(2.234484) q[2];
rz(2.234327) q[3];
sx q[3];
rz(-1.7472569) q[3];
sx q[3];
rz(3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-1.1842781) q[0];
sx q[0];
rz(-2.2397581) q[0];
sx q[0];
rz(-2.9515475) q[0];
rz(-1.5826781) q[1];
sx q[1];
rz(-1.5946486) q[1];
sx q[1];
rz(1.83439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060576749) q[0];
sx q[0];
rz(-1.2071484) q[0];
sx q[0];
rz(0.15369291) q[0];
x q[1];
rz(-2.2707269) q[2];
sx q[2];
rz(-1.8339089) q[2];
sx q[2];
rz(1.7990636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9865454) q[1];
sx q[1];
rz(-1.6035756) q[1];
sx q[1];
rz(-1.678534) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43453026) q[3];
sx q[3];
rz(-1.6756246) q[3];
sx q[3];
rz(-1.2140345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.068453161) q[2];
sx q[2];
rz(-1.0871202) q[2];
sx q[2];
rz(-2.9404409) q[2];
rz(1.0226095) q[3];
sx q[3];
rz(-2.4590731) q[3];
sx q[3];
rz(1.539591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1744743) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(0.87646595) q[0];
rz(-0.62249741) q[1];
sx q[1];
rz(-2.9298156) q[1];
sx q[1];
rz(2.7988787) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22607813) q[0];
sx q[0];
rz(-2.5097339) q[0];
sx q[0];
rz(-1.4591239) q[0];
rz(-pi) q[1];
rz(-1.2340401) q[2];
sx q[2];
rz(-1.6861262) q[2];
sx q[2];
rz(2.8591836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73574793) q[1];
sx q[1];
rz(-2.1830507) q[1];
sx q[1];
rz(-1.014018) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38767918) q[3];
sx q[3];
rz(-1.3362543) q[3];
sx q[3];
rz(2.925433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1891302) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-2.6978037) q[2];
rz(1.1187547) q[3];
sx q[3];
rz(-1.5951364) q[3];
sx q[3];
rz(2.5877623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583869) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(-0.063477909) q[1];
sx q[1];
rz(-1.7964446) q[1];
sx q[1];
rz(-1.7367015) q[1];
rz(1.9012801) q[2];
sx q[2];
rz(-1.5122432) q[2];
sx q[2];
rz(1.6846364) q[2];
rz(2.8982481) q[3];
sx q[3];
rz(-1.018033) q[3];
sx q[3];
rz(1.3924934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
