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
rz(3.7351146) q[1];
sx q[1];
rz(10.862196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0707998) q[0];
sx q[0];
rz(-1.6294998) q[0];
sx q[0];
rz(-3.015201) q[0];
x q[1];
rz(0.038226323) q[2];
sx q[2];
rz(-1.9598205) q[2];
sx q[2];
rz(-1.3343177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9086409) q[1];
sx q[1];
rz(-2.490762) q[1];
sx q[1];
rz(2.6292168) q[1];
rz(-pi) q[2];
rz(-1.8928321) q[3];
sx q[3];
rz(-2.7398971) q[3];
sx q[3];
rz(-2.5046949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5022016) q[2];
sx q[2];
rz(-0.23362564) q[2];
sx q[2];
rz(0.26546738) q[2];
rz(-1.4207077) q[3];
sx q[3];
rz(-1.9749494) q[3];
sx q[3];
rz(1.7336806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.0625286) q[0];
sx q[0];
rz(-1.205235) q[0];
sx q[0];
rz(-1.0519387) q[0];
rz(-2.1417292) q[1];
sx q[1];
rz(-1.5971767) q[1];
sx q[1];
rz(-2.3725407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3628415) q[0];
sx q[0];
rz(-1.0386416) q[0];
sx q[0];
rz(-0.51409419) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0508383) q[2];
sx q[2];
rz(-0.57951515) q[2];
sx q[2];
rz(1.7837634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96892541) q[1];
sx q[1];
rz(-1.7255867) q[1];
sx q[1];
rz(-1.118578) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6399986) q[3];
sx q[3];
rz(-0.47712773) q[3];
sx q[3];
rz(2.9692544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.179788) q[2];
sx q[2];
rz(-0.804681) q[2];
sx q[2];
rz(-1.556832) q[2];
rz(2.4368317) q[3];
sx q[3];
rz(-2.9289991) q[3];
sx q[3];
rz(2.052665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253118) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(-2.3749206) q[0];
rz(2.7491772) q[1];
sx q[1];
rz(-2.2800443) q[1];
sx q[1];
rz(0.79889417) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2080059) q[0];
sx q[0];
rz(-2.0547935) q[0];
sx q[0];
rz(3.1307277) q[0];
rz(-0.33943265) q[2];
sx q[2];
rz(-2.8763874) q[2];
sx q[2];
rz(-0.79193927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.047192725) q[1];
sx q[1];
rz(-2.5421738) q[1];
sx q[1];
rz(-2.1381738) q[1];
x q[2];
rz(-1.464099) q[3];
sx q[3];
rz(-1.7472526) q[3];
sx q[3];
rz(0.7430232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40156349) q[2];
sx q[2];
rz(-2.8758958) q[2];
sx q[2];
rz(3.1166039) q[2];
rz(-1.3966903) q[3];
sx q[3];
rz(-1.9091505) q[3];
sx q[3];
rz(3.0241844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5609189) q[0];
sx q[0];
rz(-2.6011401) q[0];
sx q[0];
rz(0.29676357) q[0];
rz(0.98776039) q[1];
sx q[1];
rz(-1.2514273) q[1];
sx q[1];
rz(2.3938649) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0714617) q[0];
sx q[0];
rz(-0.80118766) q[0];
sx q[0];
rz(-1.2642994) q[0];
rz(2.6358887) q[2];
sx q[2];
rz(-2.3597368) q[2];
sx q[2];
rz(-0.1270369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8104483) q[1];
sx q[1];
rz(-2.7390326) q[1];
sx q[1];
rz(3.1123575) q[1];
rz(-0.70246299) q[3];
sx q[3];
rz(-2.3335461) q[3];
sx q[3];
rz(-1.3987227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1838386) q[2];
sx q[2];
rz(-2.2695393) q[2];
sx q[2];
rz(0.5198861) q[2];
rz(-0.94684354) q[3];
sx q[3];
rz(-1.9498473) q[3];
sx q[3];
rz(-3.0900743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2648322) q[0];
sx q[0];
rz(-0.89986372) q[0];
sx q[0];
rz(1.9448036) q[0];
rz(1.6668677) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(2.8428452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0419436) q[0];
sx q[0];
rz(-1.22166) q[0];
sx q[0];
rz(2.4399202) q[0];
x q[1];
rz(-2.3442106) q[2];
sx q[2];
rz(-1.2060617) q[2];
sx q[2];
rz(-1.1656176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6406964) q[1];
sx q[1];
rz(-1.6846034) q[1];
sx q[1];
rz(-0.76687419) q[1];
rz(-2.7389333) q[3];
sx q[3];
rz(-0.91108262) q[3];
sx q[3];
rz(-0.92056489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8220736) q[2];
sx q[2];
rz(-2.3910429) q[2];
sx q[2];
rz(-0.015241148) q[2];
rz(0.12766734) q[3];
sx q[3];
rz(-0.36104194) q[3];
sx q[3];
rz(-2.291919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3540038) q[0];
sx q[0];
rz(-0.51743999) q[0];
sx q[0];
rz(0.87376839) q[0];
rz(1.2242225) q[1];
sx q[1];
rz(-2.1189549) q[1];
sx q[1];
rz(-1.9220985) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6343459) q[0];
sx q[0];
rz(-1.1218881) q[0];
sx q[0];
rz(0.44621356) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60004931) q[2];
sx q[2];
rz(-2.3193137) q[2];
sx q[2];
rz(-2.6278969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.022910206) q[1];
sx q[1];
rz(-1.7798276) q[1];
sx q[1];
rz(2.9487756) q[1];
rz(-pi) q[2];
rz(-1.9683011) q[3];
sx q[3];
rz(-1.0152721) q[3];
sx q[3];
rz(1.3501957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.445861) q[2];
sx q[2];
rz(-1.4706688) q[2];
sx q[2];
rz(0.50103465) q[2];
rz(-1.3387574) q[3];
sx q[3];
rz(-0.41156358) q[3];
sx q[3];
rz(-1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5981136) q[0];
sx q[0];
rz(-1.9444332) q[0];
sx q[0];
rz(0.57394779) q[0];
rz(-0.57169882) q[1];
sx q[1];
rz(-2.1015002) q[1];
sx q[1];
rz(1.5572825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2435574) q[0];
sx q[0];
rz(-1.1986309) q[0];
sx q[0];
rz(2.283147) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90960501) q[2];
sx q[2];
rz(-1.2750468) q[2];
sx q[2];
rz(-2.8165004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3347335) q[1];
sx q[1];
rz(-0.99913952) q[1];
sx q[1];
rz(0.57493361) q[1];
rz(-pi) q[2];
rz(-0.81767003) q[3];
sx q[3];
rz(-2.6260321) q[3];
sx q[3];
rz(1.8918693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68465662) q[2];
sx q[2];
rz(-1.1494278) q[2];
sx q[2];
rz(-1.9880902) q[2];
rz(-1.5270816) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1077147) q[0];
sx q[0];
rz(-2.3623473) q[0];
sx q[0];
rz(2.9039134) q[0];
rz(-0.73860812) q[1];
sx q[1];
rz(-1.9269201) q[1];
sx q[1];
rz(3.1225263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9063291) q[0];
sx q[0];
rz(-1.5625347) q[0];
sx q[0];
rz(-1.6434192) q[0];
x q[1];
rz(0.88709082) q[2];
sx q[2];
rz(-0.92803144) q[2];
sx q[2];
rz(-0.32059344) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7155559) q[1];
sx q[1];
rz(-1.6095248) q[1];
sx q[1];
rz(1.6524131) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.158659) q[3];
sx q[3];
rz(-2.1275828) q[3];
sx q[3];
rz(-0.24036053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0154524) q[2];
sx q[2];
rz(-1.7917874) q[2];
sx q[2];
rz(2.9665185) q[2];
rz(-1.3779047) q[3];
sx q[3];
rz(-0.62052369) q[3];
sx q[3];
rz(0.11848816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0835251) q[0];
sx q[0];
rz(-1.3578992) q[0];
sx q[0];
rz(0.51496664) q[0];
rz(2.1452451) q[1];
sx q[1];
rz(-1.0641655) q[1];
sx q[1];
rz(1.4774342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1979868) q[0];
sx q[0];
rz(-1.0062402) q[0];
sx q[0];
rz(0.72636162) q[0];
x q[1];
rz(-1.6040049) q[2];
sx q[2];
rz(-2.1589212) q[2];
sx q[2];
rz(-2.3606481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96326376) q[1];
sx q[1];
rz(-2.43578) q[1];
sx q[1];
rz(-0.8385482) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9733125) q[3];
sx q[3];
rz(-2.5096143) q[3];
sx q[3];
rz(-0.06080725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.065206334) q[2];
sx q[2];
rz(-1.6519158) q[2];
sx q[2];
rz(-0.15929407) q[2];
rz(-1.6431036) q[3];
sx q[3];
rz(-2.4703333) q[3];
sx q[3];
rz(-1.338965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047073929) q[0];
sx q[0];
rz(-0.044059489) q[0];
sx q[0];
rz(-2.2799168) q[0];
rz(1.2791951) q[1];
sx q[1];
rz(-0.28293124) q[1];
sx q[1];
rz(0.18712015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1094799) q[0];
sx q[0];
rz(-2.7811227) q[0];
sx q[0];
rz(-0.69860639) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8420024) q[2];
sx q[2];
rz(-1.091963) q[2];
sx q[2];
rz(-1.9128654) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0963073) q[1];
sx q[1];
rz(-2.2780609) q[1];
sx q[1];
rz(-1.2169669) q[1];
x q[2];
rz(-2.7469603) q[3];
sx q[3];
rz(-1.2773716) q[3];
sx q[3];
rz(-2.54486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44237915) q[2];
sx q[2];
rz(-0.97457653) q[2];
sx q[2];
rz(-1.3628192) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(-1.9938102) q[2];
sx q[2];
rz(-1.4381051) q[2];
sx q[2];
rz(-1.4249868) q[2];
rz(1.7888482) q[3];
sx q[3];
rz(-1.8321804) q[3];
sx q[3];
rz(1.2294168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
