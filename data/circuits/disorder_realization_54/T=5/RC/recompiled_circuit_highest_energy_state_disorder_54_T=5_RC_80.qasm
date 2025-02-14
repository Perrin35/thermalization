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
rz(0.74066585) q[0];
sx q[0];
rz(-1.292115) q[0];
sx q[0];
rz(2.3722755) q[0];
rz(1.5953335) q[1];
sx q[1];
rz(-2.9268664) q[1];
sx q[1];
rz(-0.62155849) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21827605) q[0];
sx q[0];
rz(-1.0729323) q[0];
sx q[0];
rz(2.455869) q[0];
rz(-pi) q[1];
rz(-1.1530959) q[2];
sx q[2];
rz(-1.5868895) q[2];
sx q[2];
rz(2.180445) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5468419) q[1];
sx q[1];
rz(-1.2870645) q[1];
sx q[1];
rz(0.8304412) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6042799) q[3];
sx q[3];
rz(-2.9859324) q[3];
sx q[3];
rz(-0.90720219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2810716) q[2];
sx q[2];
rz(-0.29025429) q[2];
sx q[2];
rz(-2.2104635) q[2];
rz(-2.9562601) q[3];
sx q[3];
rz(-1.0328181) q[3];
sx q[3];
rz(-1.4144271) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0952048) q[0];
sx q[0];
rz(-0.33400184) q[0];
sx q[0];
rz(-0.97292501) q[0];
rz(0.75353777) q[1];
sx q[1];
rz(-0.85644186) q[1];
sx q[1];
rz(-0.10261745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55089968) q[0];
sx q[0];
rz(-1.4859938) q[0];
sx q[0];
rz(1.4728611) q[0];
rz(-0.63541143) q[2];
sx q[2];
rz(-2.9757068) q[2];
sx q[2];
rz(2.7775641) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8637789) q[1];
sx q[1];
rz(-1.4881388) q[1];
sx q[1];
rz(-1.4138926) q[1];
rz(-0.13861175) q[3];
sx q[3];
rz(-1.9034121) q[3];
sx q[3];
rz(-0.0083855586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40833452) q[2];
sx q[2];
rz(-1.4380941) q[2];
sx q[2];
rz(-2.5737838) q[2];
rz(-2.7693977) q[3];
sx q[3];
rz(-1.8122383) q[3];
sx q[3];
rz(-1.6224434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55783015) q[0];
sx q[0];
rz(-2.3985641) q[0];
sx q[0];
rz(1.2416191) q[0];
rz(-0.81635967) q[1];
sx q[1];
rz(-0.38874778) q[1];
sx q[1];
rz(0.98966086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0123204) q[0];
sx q[0];
rz(-2.7033157) q[0];
sx q[0];
rz(0.52198903) q[0];
rz(0.092421345) q[2];
sx q[2];
rz(-1.7310609) q[2];
sx q[2];
rz(-1.607464) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1898797) q[1];
sx q[1];
rz(-0.44766176) q[1];
sx q[1];
rz(-2.7192804) q[1];
x q[2];
rz(-1.8258404) q[3];
sx q[3];
rz(-0.38557426) q[3];
sx q[3];
rz(0.60690597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.376754) q[2];
sx q[2];
rz(-1.2412485) q[2];
sx q[2];
rz(-1.494701) q[2];
rz(0.74690789) q[3];
sx q[3];
rz(-2.5266095) q[3];
sx q[3];
rz(-0.99228215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.295739) q[0];
sx q[0];
rz(-0.60765147) q[0];
sx q[0];
rz(-2.1920152) q[0];
rz(-1.1664561) q[1];
sx q[1];
rz(-1.8828853) q[1];
sx q[1];
rz(-0.20225254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324449) q[0];
sx q[0];
rz(-1.5209001) q[0];
sx q[0];
rz(1.6014535) q[0];
rz(-pi) q[1];
rz(2.9390088) q[2];
sx q[2];
rz(-0.69949818) q[2];
sx q[2];
rz(-1.0765178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3551605) q[1];
sx q[1];
rz(-2.0789685) q[1];
sx q[1];
rz(0.90759214) q[1];
rz(-pi) q[2];
rz(0.99503354) q[3];
sx q[3];
rz(-1.0594133) q[3];
sx q[3];
rz(0.56727876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4919081) q[2];
sx q[2];
rz(-1.7077571) q[2];
sx q[2];
rz(1.8053619) q[2];
rz(-0.46184552) q[3];
sx q[3];
rz(-1.0577842) q[3];
sx q[3];
rz(-2.113078) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.860054) q[0];
sx q[0];
rz(-3.0335732) q[0];
sx q[0];
rz(-0.39529133) q[0];
rz(-2.1832502) q[1];
sx q[1];
rz(-1.7610565) q[1];
sx q[1];
rz(2.7470632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.677934) q[0];
sx q[0];
rz(-1.5531085) q[0];
sx q[0];
rz(1.836548) q[0];
rz(2.7005211) q[2];
sx q[2];
rz(-1.8243378) q[2];
sx q[2];
rz(2.2833786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1198152) q[1];
sx q[1];
rz(-1.2477326) q[1];
sx q[1];
rz(0.82746877) q[1];
x q[2];
rz(1.35704) q[3];
sx q[3];
rz(-2.1780664) q[3];
sx q[3];
rz(0.34363817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40256527) q[2];
sx q[2];
rz(-1.6141011) q[2];
sx q[2];
rz(-1.019545) q[2];
rz(-1.839365) q[3];
sx q[3];
rz(-0.052509382) q[3];
sx q[3];
rz(-2.4359865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5604945) q[0];
sx q[0];
rz(-0.56466931) q[0];
sx q[0];
rz(-2.1288921) q[0];
rz(0.46214354) q[1];
sx q[1];
rz(-1.7299165) q[1];
sx q[1];
rz(0.046646811) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7830392) q[0];
sx q[0];
rz(-1.5601242) q[0];
sx q[0];
rz(1.693876) q[0];
x q[1];
rz(2.279206) q[2];
sx q[2];
rz(-2.7593331) q[2];
sx q[2];
rz(-1.6406989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0477284) q[1];
sx q[1];
rz(-1.0055491) q[1];
sx q[1];
rz(0.39611343) q[1];
rz(1.9479475) q[3];
sx q[3];
rz(-1.5483678) q[3];
sx q[3];
rz(-0.42322657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8320273) q[2];
sx q[2];
rz(-3.053061) q[2];
sx q[2];
rz(2.1868165) q[2];
rz(-0.66983062) q[3];
sx q[3];
rz(-1.9570743) q[3];
sx q[3];
rz(-0.36439103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1466115) q[0];
sx q[0];
rz(-2.8611188) q[0];
sx q[0];
rz(1.9914419) q[0];
rz(3.0448044) q[1];
sx q[1];
rz(-1.140118) q[1];
sx q[1];
rz(-1.4801625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.262275) q[0];
sx q[0];
rz(-1.9497383) q[0];
sx q[0];
rz(2.3728601) q[0];
x q[1];
rz(-1.8100914) q[2];
sx q[2];
rz(-1.8836438) q[2];
sx q[2];
rz(-2.8205736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5642359) q[1];
sx q[1];
rz(-0.32315394) q[1];
sx q[1];
rz(-1.1250457) q[1];
x q[2];
rz(1.893793) q[3];
sx q[3];
rz(-0.45489254) q[3];
sx q[3];
rz(-0.98857075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73411554) q[2];
sx q[2];
rz(-0.87782562) q[2];
sx q[2];
rz(0.62848282) q[2];
rz(-2.8001522) q[3];
sx q[3];
rz(-0.98613685) q[3];
sx q[3];
rz(-0.83610523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8516561) q[0];
sx q[0];
rz(-2.4735232) q[0];
sx q[0];
rz(3.0776899) q[0];
rz(-0.54197657) q[1];
sx q[1];
rz(-2.0109476) q[1];
sx q[1];
rz(-0.37905395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0444894) q[0];
sx q[0];
rz(-1.7428723) q[0];
sx q[0];
rz(-1.340036) q[0];
rz(2.3325868) q[2];
sx q[2];
rz(-1.6651285) q[2];
sx q[2];
rz(1.091711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8833911) q[1];
sx q[1];
rz(-0.54377006) q[1];
sx q[1];
rz(-3.1115502) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7160801) q[3];
sx q[3];
rz(-2.1149016) q[3];
sx q[3];
rz(1.5489006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1698251) q[2];
sx q[2];
rz(-1.1173893) q[2];
sx q[2];
rz(2.8525412) q[2];
rz(2.1257832) q[3];
sx q[3];
rz(-0.93520516) q[3];
sx q[3];
rz(-2.9142006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.10580258) q[0];
sx q[0];
rz(-2.7547024) q[0];
sx q[0];
rz(-2.8111358) q[0];
rz(2.6608048) q[1];
sx q[1];
rz(-1.5590706) q[1];
sx q[1];
rz(2.6343583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1523154) q[0];
sx q[0];
rz(-2.0880654) q[0];
sx q[0];
rz(-0.74801561) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45926913) q[2];
sx q[2];
rz(-1.9063564) q[2];
sx q[2];
rz(-0.84268779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1363234) q[1];
sx q[1];
rz(-1.988896) q[1];
sx q[1];
rz(0.24009248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60023579) q[3];
sx q[3];
rz(-1.6199937) q[3];
sx q[3];
rz(-0.95358301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5759739) q[2];
sx q[2];
rz(-2.2810292) q[2];
sx q[2];
rz(-0.38798517) q[2];
rz(0.062813736) q[3];
sx q[3];
rz(-0.7395491) q[3];
sx q[3];
rz(1.9529823) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0906618) q[0];
sx q[0];
rz(-3.1212786) q[0];
sx q[0];
rz(-0.055572979) q[0];
rz(0.22124258) q[1];
sx q[1];
rz(-1.0098927) q[1];
sx q[1];
rz(1.6411068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243603) q[0];
sx q[0];
rz(-1.8460197) q[0];
sx q[0];
rz(2.0936398) q[0];
rz(-pi) q[1];
rz(0.13672853) q[2];
sx q[2];
rz(-1.5216344) q[2];
sx q[2];
rz(-1.3280185) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55602395) q[1];
sx q[1];
rz(-2.635286) q[1];
sx q[1];
rz(2.271818) q[1];
rz(-0.9110578) q[3];
sx q[3];
rz(-0.83759901) q[3];
sx q[3];
rz(-0.78105951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1859493) q[2];
sx q[2];
rz(-0.187749) q[2];
sx q[2];
rz(2.032568) q[2];
rz(3.1254613) q[3];
sx q[3];
rz(-1.4343836) q[3];
sx q[3];
rz(-0.46512887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7598509) q[0];
sx q[0];
rz(-2.1325337) q[0];
sx q[0];
rz(2.7251563) q[0];
rz(-0.79553678) q[1];
sx q[1];
rz(-1.7557314) q[1];
sx q[1];
rz(3.0605127) q[1];
rz(-0.14079413) q[2];
sx q[2];
rz(-1.9119605) q[2];
sx q[2];
rz(-0.20902363) q[2];
rz(0.53989177) q[3];
sx q[3];
rz(-0.85689982) q[3];
sx q[3];
rz(0.11107437) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
