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
rz(0.93208575) q[0];
sx q[0];
rz(-2.1075489) q[0];
sx q[0];
rz(2.3873868) q[0];
rz(-1.3940613) q[1];
sx q[1];
rz(3.7351146) q[1];
sx q[1];
rz(10.862196) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49254848) q[0];
sx q[0];
rz(-1.6969691) q[0];
sx q[0];
rz(-1.5116219) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.038226323) q[2];
sx q[2];
rz(-1.1817721) q[2];
sx q[2];
rz(1.807275) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2329518) q[1];
sx q[1];
rz(-0.65083069) q[1];
sx q[1];
rz(0.51237583) q[1];
x q[2];
rz(1.1877443) q[3];
sx q[3];
rz(-1.6948587) q[3];
sx q[3];
rz(-1.2318486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6393911) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(0.26546738) q[2];
rz(1.7208849) q[3];
sx q[3];
rz(-1.1666433) q[3];
sx q[3];
rz(-1.7336806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(1.079064) q[0];
sx q[0];
rz(-1.205235) q[0];
sx q[0];
rz(2.089654) q[0];
rz(2.1417292) q[1];
sx q[1];
rz(-1.544416) q[1];
sx q[1];
rz(-2.3725407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3628415) q[0];
sx q[0];
rz(-2.102951) q[0];
sx q[0];
rz(0.51409419) q[0];
x q[1];
rz(3.0508383) q[2];
sx q[2];
rz(-2.5620775) q[2];
sx q[2];
rz(1.3578292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2324128) q[1];
sx q[1];
rz(-0.47623834) q[1];
sx q[1];
rz(1.9137812) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6399986) q[3];
sx q[3];
rz(-2.6644649) q[3];
sx q[3];
rz(-0.17233822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.179788) q[2];
sx q[2];
rz(-2.3369117) q[2];
sx q[2];
rz(1.5847607) q[2];
rz(-0.70476091) q[3];
sx q[3];
rz(-0.2125936) q[3];
sx q[3];
rz(-2.052665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8162808) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(-0.76667205) q[0];
rz(0.39241544) q[1];
sx q[1];
rz(-0.86154834) q[1];
sx q[1];
rz(-2.3426985) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1846599) q[0];
sx q[0];
rz(-2.6574832) q[0];
sx q[0];
rz(1.5501322) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33943265) q[2];
sx q[2];
rz(-0.26520525) q[2];
sx q[2];
rz(-2.3496534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4370712) q[1];
sx q[1];
rz(-1.0749617) q[1];
sx q[1];
rz(-0.35192546) q[1];
x q[2];
rz(-2.9641482) q[3];
sx q[3];
rz(-1.6758306) q[3];
sx q[3];
rz(-2.3326186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40156349) q[2];
sx q[2];
rz(-2.8758958) q[2];
sx q[2];
rz(3.1166039) q[2];
rz(-1.7449024) q[3];
sx q[3];
rz(-1.2324421) q[3];
sx q[3];
rz(3.0241844) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5806737) q[0];
sx q[0];
rz(-2.6011401) q[0];
sx q[0];
rz(-0.29676357) q[0];
rz(-2.1538323) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(-2.3938649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0714617) q[0];
sx q[0];
rz(-2.340405) q[0];
sx q[0];
rz(1.2642994) q[0];
rz(-pi) q[1];
rz(1.1224611) q[2];
sx q[2];
rz(-2.2349662) q[2];
sx q[2];
rz(0.78967204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2993707) q[1];
sx q[1];
rz(-1.9731745) q[1];
sx q[1];
rz(-1.55835) q[1];
rz(-2.4391297) q[3];
sx q[3];
rz(-2.3335461) q[3];
sx q[3];
rz(-1.74287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.957754) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(2.6217065) q[2];
rz(0.94684354) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(0.051518353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2648322) q[0];
sx q[0];
rz(-0.89986372) q[0];
sx q[0];
rz(-1.9448036) q[0];
rz(1.4747249) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(-2.8428452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099649053) q[0];
sx q[0];
rz(-1.22166) q[0];
sx q[0];
rz(2.4399202) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3442106) q[2];
sx q[2];
rz(-1.935531) q[2];
sx q[2];
rz(-1.9759751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9626959) q[1];
sx q[1];
rz(-0.81013727) q[1];
sx q[1];
rz(-1.728216) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87033724) q[3];
sx q[3];
rz(-1.8856109) q[3];
sx q[3];
rz(-0.39488068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3195191) q[2];
sx q[2];
rz(-0.75054979) q[2];
sx q[2];
rz(3.1263515) q[2];
rz(3.0139253) q[3];
sx q[3];
rz(-0.36104194) q[3];
sx q[3];
rz(2.291919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.78758883) q[0];
sx q[0];
rz(-2.6241527) q[0];
sx q[0];
rz(0.87376839) q[0];
rz(1.9173701) q[1];
sx q[1];
rz(-2.1189549) q[1];
sx q[1];
rz(-1.2194941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50724678) q[0];
sx q[0];
rz(-2.0197045) q[0];
sx q[0];
rz(0.44621356) q[0];
x q[1];
rz(2.4151426) q[2];
sx q[2];
rz(-1.9973576) q[2];
sx q[2];
rz(-2.5202765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5883792) q[1];
sx q[1];
rz(-1.7593652) q[1];
sx q[1];
rz(-1.3579353) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59252809) q[3];
sx q[3];
rz(-1.9059407) q[3];
sx q[3];
rz(0.0026897653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.445861) q[2];
sx q[2];
rz(-1.6709238) q[2];
sx q[2];
rz(-2.640558) q[2];
rz(-1.8028353) q[3];
sx q[3];
rz(-2.7300291) q[3];
sx q[3];
rz(-1.1494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5434791) q[0];
sx q[0];
rz(-1.9444332) q[0];
sx q[0];
rz(-0.57394779) q[0];
rz(2.5698938) q[1];
sx q[1];
rz(-1.0400925) q[1];
sx q[1];
rz(-1.5572825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97705072) q[0];
sx q[0];
rz(-0.91616154) q[0];
sx q[0];
rz(2.665401) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90960501) q[2];
sx q[2];
rz(-1.2750468) q[2];
sx q[2];
rz(0.32509229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3347335) q[1];
sx q[1];
rz(-0.99913952) q[1];
sx q[1];
rz(2.566659) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36974233) q[3];
sx q[3];
rz(-1.9387285) q[3];
sx q[3];
rz(-2.7145601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.456936) q[2];
sx q[2];
rz(-1.9921649) q[2];
sx q[2];
rz(-1.1535025) q[2];
rz(-1.614511) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.8089186) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1077147) q[0];
sx q[0];
rz(-0.77924538) q[0];
sx q[0];
rz(-0.23767924) q[0];
rz(-2.4029845) q[1];
sx q[1];
rz(-1.2146726) q[1];
sx q[1];
rz(-0.019066378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23526351) q[0];
sx q[0];
rz(-1.5790579) q[0];
sx q[0];
rz(-1.4981734) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76808651) q[2];
sx q[2];
rz(-1.0407018) q[2];
sx q[2];
rz(-1.7045316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.42603675) q[1];
sx q[1];
rz(-1.5320679) q[1];
sx q[1];
rz(-1.4891796) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4138149) q[3];
sx q[3];
rz(-0.78633868) q[3];
sx q[3];
rz(-1.1407408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0154524) q[2];
sx q[2];
rz(-1.3498053) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0835251) q[0];
sx q[0];
rz(-1.3578992) q[0];
sx q[0];
rz(-0.51496664) q[0];
rz(-2.1452451) q[1];
sx q[1];
rz(-2.0774272) q[1];
sx q[1];
rz(1.4774342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2124951) q[0];
sx q[0];
rz(-0.97508865) q[0];
sx q[0];
rz(2.2736249) q[0];
x q[1];
rz(2.5532132) q[2];
sx q[2];
rz(-1.5431689) q[2];
sx q[2];
rz(2.333312) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9340201) q[1];
sx q[1];
rz(-1.1222543) q[1];
sx q[1];
rz(-2.1357082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9733125) q[3];
sx q[3];
rz(-0.63197836) q[3];
sx q[3];
rz(0.06080725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0763863) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(-0.15929407) q[2];
rz(-1.4984891) q[3];
sx q[3];
rz(-0.67125932) q[3];
sx q[3];
rz(1.8026277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047073929) q[0];
sx q[0];
rz(-3.0975332) q[0];
sx q[0];
rz(-2.2799168) q[0];
rz(-1.2791951) q[1];
sx q[1];
rz(-2.8586614) q[1];
sx q[1];
rz(0.18712015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0321127) q[0];
sx q[0];
rz(-2.7811227) q[0];
sx q[0];
rz(2.4429863) q[0];
rz(-0.47640064) q[2];
sx q[2];
rz(-0.54506732) q[2];
sx q[2];
rz(-2.4558512) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0963073) q[1];
sx q[1];
rz(-2.2780609) q[1];
sx q[1];
rz(1.2169669) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.475521) q[3];
sx q[3];
rz(-2.6544673) q[3];
sx q[3];
rz(0.36706007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6992135) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(1.7787735) q[2];
rz(-1.5022701) q[3];
sx q[3];
rz(-0.5539186) q[3];
sx q[3];
rz(1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5300071) q[0];
sx q[0];
rz(-1.3349608) q[0];
sx q[0];
rz(0.0064749574) q[0];
rz(-2.6352128) q[1];
sx q[1];
rz(-0.19229278) q[1];
sx q[1];
rz(-2.2813588) q[1];
rz(1.2564332) q[2];
sx q[2];
rz(-0.4421352) q[2];
sx q[2];
rz(3.0015702) q[2];
rz(-2.8741638) q[3];
sx q[3];
rz(-1.3602644) q[3];
sx q[3];
rz(-0.2841831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
