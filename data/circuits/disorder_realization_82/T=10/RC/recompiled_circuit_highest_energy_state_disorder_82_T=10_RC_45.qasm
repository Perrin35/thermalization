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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(-1.9544344) q[0];
rz(4.123426) q[1];
sx q[1];
rz(0.47458664) q[1];
sx q[1];
rz(8.4835806) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1008074) q[0];
sx q[0];
rz(-1.5724764) q[0];
sx q[0];
rz(-2.6462862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0082939) q[2];
sx q[2];
rz(-2.9037711) q[2];
sx q[2];
rz(2.651525) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57454771) q[1];
sx q[1];
rz(-0.68463782) q[1];
sx q[1];
rz(2.5570832) q[1];
x q[2];
rz(1.4210798) q[3];
sx q[3];
rz(-2.6805373) q[3];
sx q[3];
rz(-2.113101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5473951) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(-1.0198062) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(-1.3297133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52784598) q[0];
sx q[0];
rz(-2.3044523) q[0];
sx q[0];
rz(-1.6380731) q[0];
rz(-1.9095406) q[1];
sx q[1];
rz(-2.3161395) q[1];
sx q[1];
rz(-1.253461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.764956) q[0];
sx q[0];
rz(-2.6692163) q[0];
sx q[0];
rz(2.0852023) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49053662) q[2];
sx q[2];
rz(-2.5744458) q[2];
sx q[2];
rz(2.4012152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4644891) q[1];
sx q[1];
rz(-1.3687412) q[1];
sx q[1];
rz(1.4239763) q[1];
x q[2];
rz(2.0157865) q[3];
sx q[3];
rz(-0.87712327) q[3];
sx q[3];
rz(-1.7852269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85679179) q[2];
sx q[2];
rz(-2.5992726) q[2];
sx q[2];
rz(-2.8112603) q[2];
rz(3.0387943) q[3];
sx q[3];
rz(-1.8967352) q[3];
sx q[3];
rz(-2.5225294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5264346) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(-1.0636299) q[0];
rz(3.0377153) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(2.0361384) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5828667) q[0];
sx q[0];
rz(-1.5387282) q[0];
sx q[0];
rz(1.523287) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8145829) q[2];
sx q[2];
rz(-2.3578224) q[2];
sx q[2];
rz(-1.6683589) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3893095) q[1];
sx q[1];
rz(-1.0573309) q[1];
sx q[1];
rz(-3.0016243) q[1];
x q[2];
rz(-1.6126339) q[3];
sx q[3];
rz(-1.5825019) q[3];
sx q[3];
rz(-0.042543471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1198472) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(0.22551192) q[2];
rz(1.8146023) q[3];
sx q[3];
rz(-2.3046389) q[3];
sx q[3];
rz(0.64083159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.91318146) q[0];
sx q[0];
rz(-1.5629733) q[0];
sx q[0];
rz(-2.5033503) q[0];
rz(-2.0300716) q[1];
sx q[1];
rz(-2.1941954) q[1];
sx q[1];
rz(-0.18724719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603724) q[0];
sx q[0];
rz(-0.6938254) q[0];
sx q[0];
rz(1.1129473) q[0];
rz(-1.6333605) q[2];
sx q[2];
rz(-2.3507032) q[2];
sx q[2];
rz(-3.084201) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79689246) q[1];
sx q[1];
rz(-2.6019367) q[1];
sx q[1];
rz(2.9873288) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10113206) q[3];
sx q[3];
rz(-0.71590483) q[3];
sx q[3];
rz(-2.9734395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1049261) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(-0.46910134) q[2];
rz(0.30096287) q[3];
sx q[3];
rz(-0.260869) q[3];
sx q[3];
rz(-1.4891967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8982573) q[0];
sx q[0];
rz(-1.6247592) q[0];
sx q[0];
rz(-2.6014056) q[0];
rz(0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(-0.25925055) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2592372) q[0];
sx q[0];
rz(-1.8412388) q[0];
sx q[0];
rz(-2.8195802) q[0];
rz(0.44633003) q[2];
sx q[2];
rz(-2.0392373) q[2];
sx q[2];
rz(-0.51821741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4673678) q[1];
sx q[1];
rz(-1.4150054) q[1];
sx q[1];
rz(2.1464661) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8794925) q[3];
sx q[3];
rz(-1.8912695) q[3];
sx q[3];
rz(2.4985034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79869142) q[2];
sx q[2];
rz(-0.68621245) q[2];
sx q[2];
rz(0.9101103) q[2];
rz(-0.56509194) q[3];
sx q[3];
rz(-1.7730954) q[3];
sx q[3];
rz(-0.71649396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703576) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(1.2275335) q[0];
rz(-0.50499376) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(1.203677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6617463) q[0];
sx q[0];
rz(-1.8176314) q[0];
sx q[0];
rz(1.4479475) q[0];
rz(-2.2199322) q[2];
sx q[2];
rz(-1.1716915) q[2];
sx q[2];
rz(-1.468935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5745637) q[1];
sx q[1];
rz(-2.1334799) q[1];
sx q[1];
rz(-0.74498727) q[1];
x q[2];
rz(-0.26546462) q[3];
sx q[3];
rz(-2.0164312) q[3];
sx q[3];
rz(-1.7441986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(2.1825979) q[2];
rz(0.2995019) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(1.449077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9129979) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(-0.49384299) q[0];
rz(-2.2834868) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(1.3546622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49789594) q[0];
sx q[0];
rz(-2.5312811) q[0];
sx q[0];
rz(1.9998788) q[0];
x q[1];
rz(-1.3184403) q[2];
sx q[2];
rz(-0.50559536) q[2];
sx q[2];
rz(-2.8211968) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5112016) q[1];
sx q[1];
rz(-1.5005209) q[1];
sx q[1];
rz(-0.75842885) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89833791) q[3];
sx q[3];
rz(-1.9914989) q[3];
sx q[3];
rz(1.6572233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4570423) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(0.54614145) q[2];
rz(2.9438733) q[3];
sx q[3];
rz(-2.1106014) q[3];
sx q[3];
rz(-0.24699591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9207183) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(-1.1065296) q[0];
rz(2.5571892) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(-1.0027142) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4179392) q[0];
sx q[0];
rz(-1.1375715) q[0];
sx q[0];
rz(2.9648925) q[0];
x q[1];
rz(1.2757589) q[2];
sx q[2];
rz(-1.8031839) q[2];
sx q[2];
rz(-0.7600998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.20368491) q[1];
sx q[1];
rz(-1.6069429) q[1];
sx q[1];
rz(0.51236492) q[1];
x q[2];
rz(2.5224696) q[3];
sx q[3];
rz(-1.8612218) q[3];
sx q[3];
rz(2.8324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5661085) q[2];
sx q[2];
rz(-0.9245975) q[2];
sx q[2];
rz(1.7522579) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(-2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0480807) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(2.804948) q[0];
rz(3.0632784) q[1];
sx q[1];
rz(-1.9322194) q[1];
sx q[1];
rz(1.1557109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9243568) q[0];
sx q[0];
rz(-1.9072394) q[0];
sx q[0];
rz(2.2948059) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9843654) q[2];
sx q[2];
rz(-1.509827) q[2];
sx q[2];
rz(-2.3496036) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3382321) q[1];
sx q[1];
rz(-2.6322717) q[1];
sx q[1];
rz(-1.2912911) q[1];
rz(-pi) q[2];
rz(-1.6011681) q[3];
sx q[3];
rz(-1.2541459) q[3];
sx q[3];
rz(0.79398549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0597421) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(2.3183863) q[2];
rz(1.0169704) q[3];
sx q[3];
rz(-0.40139324) q[3];
sx q[3];
rz(0.038173525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391649) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(-1.7145351) q[0];
rz(1.6629705) q[1];
sx q[1];
rz(-1.0717816) q[1];
sx q[1];
rz(-0.85085416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9768391) q[0];
sx q[0];
rz(-2.0434336) q[0];
sx q[0];
rz(-0.12906277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97423221) q[2];
sx q[2];
rz(-2.9880004) q[2];
sx q[2];
rz(-3.0418015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58356279) q[1];
sx q[1];
rz(-2.0940015) q[1];
sx q[1];
rz(-1.5091404) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6353278) q[3];
sx q[3];
rz(-1.4307662) q[3];
sx q[3];
rz(-2.9794995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6829546) q[2];
sx q[2];
rz(-1.7913603) q[2];
sx q[2];
rz(2.3214935) q[2];
rz(-1.5470777) q[3];
sx q[3];
rz(-1.7087414) q[3];
sx q[3];
rz(1.1944176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9783258) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(-2.3691879) q[1];
sx q[1];
rz(-0.62320566) q[1];
sx q[1];
rz(-2.1355245) q[1];
rz(1.7987235) q[2];
sx q[2];
rz(-2.2493636) q[2];
sx q[2];
rz(2.0668088) q[2];
rz(2.3631572) q[3];
sx q[3];
rz(-2.5166471) q[3];
sx q[3];
rz(-2.2770721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
