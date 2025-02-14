OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6168851) q[0];
sx q[0];
rz(5.5445639) q[0];
sx q[0];
rz(9.5417547) q[0];
rz(0.95247954) q[1];
sx q[1];
rz(-1.6713961) q[1];
sx q[1];
rz(0.87364668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40634333) q[0];
sx q[0];
rz(-0.85052089) q[0];
sx q[0];
rz(-0.27168897) q[0];
rz(1.123465) q[2];
sx q[2];
rz(-1.8904933) q[2];
sx q[2];
rz(1.5699707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80109875) q[1];
sx q[1];
rz(-1.5135161) q[1];
sx q[1];
rz(2.1199473) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1659032) q[3];
sx q[3];
rz(-1.3171853) q[3];
sx q[3];
rz(-0.028279956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17254193) q[2];
sx q[2];
rz(-1.3323063) q[2];
sx q[2];
rz(1.8633899) q[2];
rz(1.5365907) q[3];
sx q[3];
rz(-1.2383818) q[3];
sx q[3];
rz(2.8240375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4441967) q[0];
sx q[0];
rz(-2.9228656) q[0];
sx q[0];
rz(-0.87795192) q[0];
rz(-0.76389337) q[1];
sx q[1];
rz(-1.0822252) q[1];
sx q[1];
rz(2.0108932) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5754841) q[0];
sx q[0];
rz(-2.974925) q[0];
sx q[0];
rz(-1.1995633) q[0];
x q[1];
rz(1.2487407) q[2];
sx q[2];
rz(-2.6832504) q[2];
sx q[2];
rz(-2.2001668) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69399) q[1];
sx q[1];
rz(-0.75335767) q[1];
sx q[1];
rz(-1.5894458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5155792) q[3];
sx q[3];
rz(-1.6580868) q[3];
sx q[3];
rz(-1.5777301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6905602) q[2];
sx q[2];
rz(-1.6310383) q[2];
sx q[2];
rz(-2.4859599) q[2];
rz(1.2304652) q[3];
sx q[3];
rz(-0.96223193) q[3];
sx q[3];
rz(-0.83278304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7448298) q[0];
sx q[0];
rz(-1.9704882) q[0];
sx q[0];
rz(0.50286621) q[0];
rz(-0.14547959) q[1];
sx q[1];
rz(-1.611404) q[1];
sx q[1];
rz(1.1605638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6737719) q[0];
sx q[0];
rz(-1.2611817) q[0];
sx q[0];
rz(1.2542115) q[0];
rz(-1.7021178) q[2];
sx q[2];
rz(-1.9177556) q[2];
sx q[2];
rz(-1.0441213) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8050011) q[1];
sx q[1];
rz(-2.1865926) q[1];
sx q[1];
rz(2.6686882) q[1];
rz(-pi) q[2];
x q[2];
rz(1.840402) q[3];
sx q[3];
rz(-1.2343711) q[3];
sx q[3];
rz(0.92576448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9460556) q[2];
sx q[2];
rz(-1.035752) q[2];
sx q[2];
rz(-0.76988402) q[2];
rz(2.5464673) q[3];
sx q[3];
rz(-2.6502521) q[3];
sx q[3];
rz(2.6814521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6339517) q[0];
sx q[0];
rz(-1.0469629) q[0];
sx q[0];
rz(1.3557583) q[0];
rz(-0.62375623) q[1];
sx q[1];
rz(-1.535894) q[1];
sx q[1];
rz(2.6502868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2441933) q[0];
sx q[0];
rz(-1.4174875) q[0];
sx q[0];
rz(-0.53863948) q[0];
x q[1];
rz(2.86209) q[2];
sx q[2];
rz(-1.3692999) q[2];
sx q[2];
rz(-3.1147) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92190336) q[1];
sx q[1];
rz(-1.2619201) q[1];
sx q[1];
rz(-2.3050453) q[1];
x q[2];
rz(-2.4748865) q[3];
sx q[3];
rz(-1.1475862) q[3];
sx q[3];
rz(0.70021399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9451311) q[2];
sx q[2];
rz(-1.6170231) q[2];
sx q[2];
rz(0.54836908) q[2];
rz(-1.8114926) q[3];
sx q[3];
rz(-0.29812223) q[3];
sx q[3];
rz(-2.8032081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86692989) q[0];
sx q[0];
rz(-0.45329705) q[0];
sx q[0];
rz(-1.0484265) q[0];
rz(0.10920814) q[1];
sx q[1];
rz(-1.5914773) q[1];
sx q[1];
rz(-1.106326) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7503742) q[0];
sx q[0];
rz(-0.24570172) q[0];
sx q[0];
rz(1.1121763) q[0];
rz(-0.84594251) q[2];
sx q[2];
rz(-1.8258313) q[2];
sx q[2];
rz(-0.77470335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8804923) q[1];
sx q[1];
rz(-1.893781) q[1];
sx q[1];
rz(2.9862981) q[1];
x q[2];
rz(0.23979964) q[3];
sx q[3];
rz(-1.8604401) q[3];
sx q[3];
rz(-2.0568796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3722374) q[2];
sx q[2];
rz(-1.8753588) q[2];
sx q[2];
rz(-0.20277578) q[2];
rz(2.3619704) q[3];
sx q[3];
rz(-1.0896261) q[3];
sx q[3];
rz(0.37105086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5146273) q[0];
sx q[0];
rz(-0.045128673) q[0];
sx q[0];
rz(0.87920642) q[0];
rz(2.5458287) q[1];
sx q[1];
rz(-2.2674982) q[1];
sx q[1];
rz(1.8858058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0107675) q[0];
sx q[0];
rz(-1.6629476) q[0];
sx q[0];
rz(-1.5177726) q[0];
rz(2.6379616) q[2];
sx q[2];
rz(-2.5066895) q[2];
sx q[2];
rz(-1.1213999) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3869898) q[1];
sx q[1];
rz(-2.5856254) q[1];
sx q[1];
rz(-0.43773361) q[1];
x q[2];
rz(2.5257931) q[3];
sx q[3];
rz(-1.7753505) q[3];
sx q[3];
rz(1.2990862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51861989) q[2];
sx q[2];
rz(-2.8039248) q[2];
sx q[2];
rz(-1.6597623) q[2];
rz(1.4431813) q[3];
sx q[3];
rz(-1.7549763) q[3];
sx q[3];
rz(2.0229009) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5959394) q[0];
sx q[0];
rz(-0.82939363) q[0];
sx q[0];
rz(-0.19350061) q[0];
rz(0.045348383) q[1];
sx q[1];
rz(-1.9769316) q[1];
sx q[1];
rz(-2.4605816) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15475965) q[0];
sx q[0];
rz(-1.4082419) q[0];
sx q[0];
rz(1.6587371) q[0];
rz(-pi) q[1];
rz(0.25410507) q[2];
sx q[2];
rz(-2.0850344) q[2];
sx q[2];
rz(0.50287535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4894708) q[1];
sx q[1];
rz(-0.014111405) q[1];
sx q[1];
rz(-1.1443787) q[1];
rz(-pi) q[2];
rz(-2.8425334) q[3];
sx q[3];
rz(-2.7388696) q[3];
sx q[3];
rz(-0.85351588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27651522) q[2];
sx q[2];
rz(-2.4962208) q[2];
sx q[2];
rz(-1.3481677) q[2];
rz(-1.3639785) q[3];
sx q[3];
rz(-1.8230702) q[3];
sx q[3];
rz(-1.9510423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50787038) q[0];
sx q[0];
rz(-0.81955925) q[0];
sx q[0];
rz(0.67935294) q[0];
rz(-0.23718111) q[1];
sx q[1];
rz(-0.91865426) q[1];
sx q[1];
rz(1.0666581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5954067) q[0];
sx q[0];
rz(-2.0810938) q[0];
sx q[0];
rz(1.575281) q[0];
rz(-0.96644281) q[2];
sx q[2];
rz(-1.6487412) q[2];
sx q[2];
rz(-1.9029531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.28344) q[1];
sx q[1];
rz(-1.4464708) q[1];
sx q[1];
rz(1.0470301) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2150463) q[3];
sx q[3];
rz(-0.51623453) q[3];
sx q[3];
rz(-0.16027361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25881585) q[2];
sx q[2];
rz(-1.2315742) q[2];
sx q[2];
rz(-2.7302177) q[2];
rz(1.8607633) q[3];
sx q[3];
rz(-2.4954093) q[3];
sx q[3];
rz(-2.083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48788747) q[0];
sx q[0];
rz(-1.3478841) q[0];
sx q[0];
rz(2.2109798) q[0];
rz(2.6764684) q[1];
sx q[1];
rz(-2.5587176) q[1];
sx q[1];
rz(1.4177657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2965558) q[0];
sx q[0];
rz(-0.64528685) q[0];
sx q[0];
rz(2.7373903) q[0];
rz(0.88820712) q[2];
sx q[2];
rz(-2.0964551) q[2];
sx q[2];
rz(-1.5789248) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.666856) q[1];
sx q[1];
rz(-2.6856093) q[1];
sx q[1];
rz(0.70959301) q[1];
x q[2];
rz(-2.1499436) q[3];
sx q[3];
rz(-1.4738899) q[3];
sx q[3];
rz(-2.4976418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.174939) q[2];
sx q[2];
rz(-2.311309) q[2];
sx q[2];
rz(1.2154328) q[2];
rz(-2.4231966) q[3];
sx q[3];
rz(-1.6738439) q[3];
sx q[3];
rz(-0.64232701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4833118) q[0];
sx q[0];
rz(-0.42228666) q[0];
sx q[0];
rz(1.8102113) q[0];
rz(-2.4347958) q[1];
sx q[1];
rz(-1.7168609) q[1];
sx q[1];
rz(-1.7165064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9658222) q[0];
sx q[0];
rz(-0.81908161) q[0];
sx q[0];
rz(0.62721647) q[0];
rz(-0.064878929) q[2];
sx q[2];
rz(-2.2264495) q[2];
sx q[2];
rz(0.82198373) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7692868) q[1];
sx q[1];
rz(-2.0897749) q[1];
sx q[1];
rz(1.0610046) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44300191) q[3];
sx q[3];
rz(-2.2878929) q[3];
sx q[3];
rz(2.9505961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66296545) q[2];
sx q[2];
rz(-2.4586283) q[2];
sx q[2];
rz(1.1753725) q[2];
rz(3.1183682) q[3];
sx q[3];
rz(-2.569779) q[3];
sx q[3];
rz(0.2636675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5671253) q[0];
sx q[0];
rz(-0.99526417) q[0];
sx q[0];
rz(-1.9010726) q[0];
rz(2.4868838) q[1];
sx q[1];
rz(-1.2798825) q[1];
sx q[1];
rz(1.9725694) q[1];
rz(-1.7083009) q[2];
sx q[2];
rz(-1.6824097) q[2];
sx q[2];
rz(-2.8472441) q[2];
rz(0.8620048) q[3];
sx q[3];
rz(-1.0222407) q[3];
sx q[3];
rz(0.23289451) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
