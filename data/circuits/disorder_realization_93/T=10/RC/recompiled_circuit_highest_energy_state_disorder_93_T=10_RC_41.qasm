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
rz(-2.2095069) q[0];
sx q[0];
rz(-1.0340438) q[0];
sx q[0];
rz(0.75420585) q[0];
rz(-1.3940613) q[1];
sx q[1];
rz(-2.5480707) q[1];
sx q[1];
rz(1.4374179) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93255723) q[0];
sx q[0];
rz(-0.13929312) q[0];
sx q[0];
rz(-0.4362696) q[0];
x q[1];
rz(1.4778334) q[2];
sx q[2];
rz(-2.7507901) q[2];
sx q[2];
rz(-1.2338232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5240548) q[1];
sx q[1];
rz(-2.1270942) q[1];
sx q[1];
rz(1.9281045) q[1];
x q[2];
rz(-1.2487605) q[3];
sx q[3];
rz(-0.40169558) q[3];
sx q[3];
rz(-2.5046949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6393911) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(0.26546738) q[2];
rz(-1.4207077) q[3];
sx q[3];
rz(-1.9749494) q[3];
sx q[3];
rz(-1.407912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.079064) q[0];
sx q[0];
rz(-1.9363576) q[0];
sx q[0];
rz(-2.089654) q[0];
rz(2.1417292) q[1];
sx q[1];
rz(-1.5971767) q[1];
sx q[1];
rz(-0.76905191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6175943) q[0];
sx q[0];
rz(-2.4193561) q[0];
sx q[0];
rz(-2.2666032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.090754358) q[2];
sx q[2];
rz(-0.57951515) q[2];
sx q[2];
rz(1.3578292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5271142) q[1];
sx q[1];
rz(-2.0172146) q[1];
sx q[1];
rz(2.9698257) q[1];
rz(2.6399986) q[3];
sx q[3];
rz(-0.47712773) q[3];
sx q[3];
rz(0.17233822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96180463) q[2];
sx q[2];
rz(-0.804681) q[2];
sx q[2];
rz(-1.5847607) q[2];
rz(-2.4368317) q[3];
sx q[3];
rz(-0.2125936) q[3];
sx q[3];
rz(2.052665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253118) q[0];
sx q[0];
rz(-2.3820057) q[0];
sx q[0];
rz(2.3749206) q[0];
rz(0.39241544) q[1];
sx q[1];
rz(-2.2800443) q[1];
sx q[1];
rz(-0.79889417) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93358675) q[0];
sx q[0];
rz(-1.0867991) q[0];
sx q[0];
rz(-3.1307277) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33943265) q[2];
sx q[2];
rz(-2.8763874) q[2];
sx q[2];
rz(0.79193927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70452141) q[1];
sx q[1];
rz(-2.066631) q[1];
sx q[1];
rz(2.7896672) q[1];
rz(0.53839012) q[3];
sx q[3];
rz(-0.205919) q[3];
sx q[3];
rz(-1.8507322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7400292) q[2];
sx q[2];
rz(-2.8758958) q[2];
sx q[2];
rz(-0.024988739) q[2];
rz(-1.3966903) q[3];
sx q[3];
rz(-1.9091505) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(1.5806737) q[0];
sx q[0];
rz(-0.54045254) q[0];
sx q[0];
rz(2.8448291) q[0];
rz(-2.1538323) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(0.74772778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4241705) q[0];
sx q[0];
rz(-1.7892196) q[0];
sx q[0];
rz(0.79346795) q[0];
x q[1];
rz(2.426366) q[2];
sx q[2];
rz(-1.919121) q[2];
sx q[2];
rz(-2.0722318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8750413) q[1];
sx q[1];
rz(-1.5822486) q[1];
sx q[1];
rz(-0.40240605) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70246299) q[3];
sx q[3];
rz(-0.80804658) q[3];
sx q[3];
rz(-1.74287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1838386) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(2.6217065) q[2];
rz(0.94684354) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(-3.0900743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87676048) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(1.9448036) q[0];
rz(1.4747249) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(-2.8428452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7525879) q[0];
sx q[0];
rz(-2.2225131) q[0];
sx q[0];
rz(2.0155922) q[0];
rz(2.6514154) q[2];
sx q[2];
rz(-0.85974795) q[2];
sx q[2];
rz(2.4012964) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50089624) q[1];
sx q[1];
rz(-1.6846034) q[1];
sx q[1];
rz(-2.3747185) q[1];
rz(-pi) q[2];
rz(-2.7389333) q[3];
sx q[3];
rz(-2.23051) q[3];
sx q[3];
rz(0.92056489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8220736) q[2];
sx q[2];
rz(-0.75054979) q[2];
sx q[2];
rz(0.015241148) q[2];
rz(-0.12766734) q[3];
sx q[3];
rz(-0.36104194) q[3];
sx q[3];
rz(-0.84967363) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3540038) q[0];
sx q[0];
rz(-2.6241527) q[0];
sx q[0];
rz(-2.2678243) q[0];
rz(-1.9173701) q[1];
sx q[1];
rz(-2.1189549) q[1];
sx q[1];
rz(-1.9220985) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8003094) q[0];
sx q[0];
rz(-2.5195846) q[0];
sx q[0];
rz(-0.84026297) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0245514) q[2];
sx q[2];
rz(-0.92140475) q[2];
sx q[2];
rz(-1.3017728) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4095597) q[1];
sx q[1];
rz(-0.28342208) q[1];
sx q[1];
rz(2.3054512) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1732916) q[3];
sx q[3];
rz(-2.1263206) q[3];
sx q[3];
rz(1.3501957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.445861) q[2];
sx q[2];
rz(-1.4706688) q[2];
sx q[2];
rz(-2.640558) q[2];
rz(-1.3387574) q[3];
sx q[3];
rz(-2.7300291) q[3];
sx q[3];
rz(1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.5434791) q[0];
sx q[0];
rz(-1.1971594) q[0];
sx q[0];
rz(-0.57394779) q[0];
rz(-2.5698938) q[1];
sx q[1];
rz(-1.0400925) q[1];
sx q[1];
rz(-1.5843102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1645419) q[0];
sx q[0];
rz(-2.2254311) q[0];
sx q[0];
rz(-0.47619168) q[0];
rz(0.90960501) q[2];
sx q[2];
rz(-1.8665458) q[2];
sx q[2];
rz(0.32509229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6819463) q[1];
sx q[1];
rz(-0.78725029) q[1];
sx q[1];
rz(-0.86904591) q[1];
rz(-1.9628223) q[3];
sx q[3];
rz(-1.2268664) q[3];
sx q[3];
rz(-2.1363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68465662) q[2];
sx q[2];
rz(-1.1494278) q[2];
sx q[2];
rz(-1.9880902) q[2];
rz(1.614511) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1077147) q[0];
sx q[0];
rz(-0.77924538) q[0];
sx q[0];
rz(-2.9039134) q[0];
rz(0.73860812) q[1];
sx q[1];
rz(-1.2146726) q[1];
sx q[1];
rz(-0.019066378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23526351) q[0];
sx q[0];
rz(-1.5625347) q[0];
sx q[0];
rz(-1.4981734) q[0];
x q[1];
rz(0.76808651) q[2];
sx q[2];
rz(-1.0407018) q[2];
sx q[2];
rz(1.7045316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.993666) q[1];
sx q[1];
rz(-1.6523517) q[1];
sx q[1];
rz(3.102735) q[1];
rz(-pi) q[2];
rz(2.158659) q[3];
sx q[3];
rz(-2.1275828) q[3];
sx q[3];
rz(-2.9012321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1261403) q[2];
sx q[2];
rz(-1.3498053) q[2];
sx q[2];
rz(0.17507412) q[2];
rz(-1.763688) q[3];
sx q[3];
rz(-0.62052369) q[3];
sx q[3];
rz(-0.11848816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058067583) q[0];
sx q[0];
rz(-1.7836934) q[0];
sx q[0];
rz(2.626626) q[0];
rz(-2.1452451) q[1];
sx q[1];
rz(-2.0774272) q[1];
sx q[1];
rz(1.4774342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9436059) q[0];
sx q[0];
rz(-2.1353525) q[0];
sx q[0];
rz(-2.415231) q[0];
x q[1];
rz(-1.5375877) q[2];
sx q[2];
rz(-2.1589212) q[2];
sx q[2];
rz(2.3606481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96326376) q[1];
sx q[1];
rz(-0.70581268) q[1];
sx q[1];
rz(2.3030445) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9733125) q[3];
sx q[3];
rz(-0.63197836) q[3];
sx q[3];
rz(-3.0807854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.065206334) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(0.15929407) q[2];
rz(1.6431036) q[3];
sx q[3];
rz(-0.67125932) q[3];
sx q[3];
rz(-1.338965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945187) q[0];
sx q[0];
rz(-0.044059489) q[0];
sx q[0];
rz(-0.86167589) q[0];
rz(-1.2791951) q[1];
sx q[1];
rz(-0.28293124) q[1];
sx q[1];
rz(-0.18712015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9368071) q[0];
sx q[0];
rz(-1.7996368) q[0];
sx q[0];
rz(-0.28099999) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8420024) q[2];
sx q[2];
rz(-2.0496297) q[2];
sx q[2];
rz(-1.2287272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0963073) q[1];
sx q[1];
rz(-0.86353179) q[1];
sx q[1];
rz(1.2169669) q[1];
rz(-pi) q[2];
rz(0.39463233) q[3];
sx q[3];
rz(-1.2773716) q[3];
sx q[3];
rz(0.59673264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44237915) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(1.3628192) q[2];
rz(1.5022701) q[3];
sx q[3];
rz(-0.5539186) q[3];
sx q[3];
rz(1.4208008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61158553) q[0];
sx q[0];
rz(-1.8066318) q[0];
sx q[0];
rz(-3.1351177) q[0];
rz(-0.50637983) q[1];
sx q[1];
rz(-2.9492999) q[1];
sx q[1];
rz(0.86023387) q[1];
rz(-1.1477825) q[2];
sx q[2];
rz(-1.7034875) q[2];
sx q[2];
rz(1.7166058) q[2];
rz(-1.7888482) q[3];
sx q[3];
rz(-1.3094123) q[3];
sx q[3];
rz(-1.9121758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
