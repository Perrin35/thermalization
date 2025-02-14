OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.75093961) q[0];
sx q[0];
rz(3.7151389) q[0];
sx q[0];
rz(9.4555785) q[0];
rz(-2.7544694) q[1];
sx q[1];
rz(-0.20063278) q[1];
sx q[1];
rz(0.90955847) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5542339) q[0];
sx q[0];
rz(-2.2116714) q[0];
sx q[0];
rz(2.4184722) q[0];
rz(-pi) q[1];
rz(2.4417905) q[2];
sx q[2];
rz(-1.095311) q[2];
sx q[2];
rz(2.3043927) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2594257) q[1];
sx q[1];
rz(-1.7710238) q[1];
sx q[1];
rz(3.0012896) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2991139) q[3];
sx q[3];
rz(-1.4702904) q[3];
sx q[3];
rz(2.1372014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6261141) q[2];
sx q[2];
rz(-1.6823744) q[2];
sx q[2];
rz(0.98575753) q[2];
rz(-1.9043026) q[3];
sx q[3];
rz(-2.3111549) q[3];
sx q[3];
rz(-1.5498836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80376959) q[0];
sx q[0];
rz(-2.0812806) q[0];
sx q[0];
rz(0.34326237) q[0];
rz(0.23974165) q[1];
sx q[1];
rz(-2.7570351) q[1];
sx q[1];
rz(-3.1125617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.507258) q[0];
sx q[0];
rz(-1.0936948) q[0];
sx q[0];
rz(1.7768589) q[0];
rz(-pi) q[1];
rz(3.0986039) q[2];
sx q[2];
rz(-2.413297) q[2];
sx q[2];
rz(-2.9790619) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5933725) q[1];
sx q[1];
rz(-2.1637329) q[1];
sx q[1];
rz(1.5895459) q[1];
rz(-pi) q[2];
rz(0.47201158) q[3];
sx q[3];
rz(-2.4084512) q[3];
sx q[3];
rz(2.5705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79641882) q[2];
sx q[2];
rz(-2.5721305) q[2];
sx q[2];
rz(2.2929537) q[2];
rz(2.7560915) q[3];
sx q[3];
rz(-1.4717088) q[3];
sx q[3];
rz(2.8482385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.327453) q[0];
sx q[0];
rz(-1.7610981) q[0];
sx q[0];
rz(3.0019548) q[0];
rz(-2.8756554) q[1];
sx q[1];
rz(-1.1775002) q[1];
sx q[1];
rz(-1.0136484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3861198) q[0];
sx q[0];
rz(-1.9546224) q[0];
sx q[0];
rz(2.1912133) q[0];
x q[1];
rz(-0.75516197) q[2];
sx q[2];
rz(-0.56729672) q[2];
sx q[2];
rz(-3.0586363) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2964281) q[1];
sx q[1];
rz(-1.0966874) q[1];
sx q[1];
rz(-1.7557322) q[1];
rz(1.6740213) q[3];
sx q[3];
rz(-1.6397392) q[3];
sx q[3];
rz(-0.063223039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93918148) q[2];
sx q[2];
rz(-2.1911669) q[2];
sx q[2];
rz(-2.6964296) q[2];
rz(-2.4116481) q[3];
sx q[3];
rz(-1.4533307) q[3];
sx q[3];
rz(-1.9385519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79778033) q[0];
sx q[0];
rz(-0.36425632) q[0];
sx q[0];
rz(-1.9601747) q[0];
rz(1.5161318) q[1];
sx q[1];
rz(-1.5189891) q[1];
sx q[1];
rz(-2.6536062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0868091) q[0];
sx q[0];
rz(-1.8300608) q[0];
sx q[0];
rz(-1.7340585) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4242742) q[2];
sx q[2];
rz(-1.0861168) q[2];
sx q[2];
rz(-0.48004638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5637602) q[1];
sx q[1];
rz(-1.8298733) q[1];
sx q[1];
rz(2.1772846) q[1];
rz(1.2156297) q[3];
sx q[3];
rz(-1.9653335) q[3];
sx q[3];
rz(1.153314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92859149) q[2];
sx q[2];
rz(-1.4139621) q[2];
sx q[2];
rz(-2.6334527) q[2];
rz(-0.43934923) q[3];
sx q[3];
rz(-0.63826799) q[3];
sx q[3];
rz(-2.5590069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4947434) q[0];
sx q[0];
rz(-1.2127533) q[0];
sx q[0];
rz(-2.5272722) q[0];
rz(0.99614227) q[1];
sx q[1];
rz(-2.3042945) q[1];
sx q[1];
rz(-0.070929758) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1114601) q[0];
sx q[0];
rz(-1.8778525) q[0];
sx q[0];
rz(-2.2024406) q[0];
x q[1];
rz(1.0164169) q[2];
sx q[2];
rz(-0.51018894) q[2];
sx q[2];
rz(-1.4025127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48551895) q[1];
sx q[1];
rz(-0.89886256) q[1];
sx q[1];
rz(0.80826064) q[1];
x q[2];
rz(-0.67885375) q[3];
sx q[3];
rz(-2.1311789) q[3];
sx q[3];
rz(-0.43401516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7190242) q[2];
sx q[2];
rz(-1.7302128) q[2];
sx q[2];
rz(-0.076449797) q[2];
rz(0.092175305) q[3];
sx q[3];
rz(-1.1276827) q[3];
sx q[3];
rz(2.5341212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344264) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(-0.46519753) q[0];
rz(3.0472962) q[1];
sx q[1];
rz(-2.1262719) q[1];
sx q[1];
rz(1.2634855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2282361) q[0];
sx q[0];
rz(-0.37274088) q[0];
sx q[0];
rz(1.2380225) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1390148) q[2];
sx q[2];
rz(-2.5093304) q[2];
sx q[2];
rz(-0.31429502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3249619) q[1];
sx q[1];
rz(-2.0010608) q[1];
sx q[1];
rz(0.63741405) q[1];
x q[2];
rz(-2.6776836) q[3];
sx q[3];
rz(-1.2281443) q[3];
sx q[3];
rz(-0.040657834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5491526) q[2];
sx q[2];
rz(-2.0024039) q[2];
sx q[2];
rz(-0.69063866) q[2];
rz(-1.45951) q[3];
sx q[3];
rz(-0.032460902) q[3];
sx q[3];
rz(2.8800817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9948147) q[0];
sx q[0];
rz(-1.0228461) q[0];
sx q[0];
rz(0.47635517) q[0];
rz(-1.9626544) q[1];
sx q[1];
rz(-0.6911239) q[1];
sx q[1];
rz(-1.8144511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0020719) q[0];
sx q[0];
rz(-2.5304573) q[0];
sx q[0];
rz(0.72211653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05169241) q[2];
sx q[2];
rz(-1.0982656) q[2];
sx q[2];
rz(-2.56531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33846291) q[1];
sx q[1];
rz(-1.6410562) q[1];
sx q[1];
rz(-0.52609013) q[1];
rz(-pi) q[2];
x q[2];
rz(0.031257625) q[3];
sx q[3];
rz(-0.61504902) q[3];
sx q[3];
rz(0.94155967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2771161) q[2];
sx q[2];
rz(-1.6552552) q[2];
sx q[2];
rz(-2.3790996) q[2];
rz(1.7755194) q[3];
sx q[3];
rz(-1.450918) q[3];
sx q[3];
rz(-2.8097927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6890126) q[0];
sx q[0];
rz(-3.0735425) q[0];
sx q[0];
rz(-2.3633603) q[0];
rz(1.4100086) q[1];
sx q[1];
rz(-0.86763132) q[1];
sx q[1];
rz(-0.57077879) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3616791) q[0];
sx q[0];
rz(-0.81244367) q[0];
sx q[0];
rz(-0.98093428) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2879175) q[2];
sx q[2];
rz(-1.520051) q[2];
sx q[2];
rz(3.1078448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5321884) q[1];
sx q[1];
rz(-1.6488607) q[1];
sx q[1];
rz(2.0933269) q[1];
rz(0.024939288) q[3];
sx q[3];
rz(-2.20633) q[3];
sx q[3];
rz(2.3221515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0290231) q[2];
sx q[2];
rz(-1.3444291) q[2];
sx q[2];
rz(2.6018108) q[2];
rz(-2.6953186) q[3];
sx q[3];
rz(-2.4036784) q[3];
sx q[3];
rz(-2.3257183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096263252) q[0];
sx q[0];
rz(-2.6427866) q[0];
sx q[0];
rz(2.2742284) q[0];
rz(1.438104) q[1];
sx q[1];
rz(-2.4080364) q[1];
sx q[1];
rz(-2.6677456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60199262) q[0];
sx q[0];
rz(-1.3253544) q[0];
sx q[0];
rz(-0.28781366) q[0];
x q[1];
rz(0.3963741) q[2];
sx q[2];
rz(-2.1389942) q[2];
sx q[2];
rz(3.0432354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0666536) q[1];
sx q[1];
rz(-2.3015296) q[1];
sx q[1];
rz(1.380062) q[1];
rz(-pi) q[2];
rz(2.2603277) q[3];
sx q[3];
rz(-1.9816958) q[3];
sx q[3];
rz(-1.2311862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69502407) q[2];
sx q[2];
rz(-2.5769672) q[2];
sx q[2];
rz(1.1308283) q[2];
rz(1.1787777) q[3];
sx q[3];
rz(-1.3935573) q[3];
sx q[3];
rz(2.668837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9216565) q[0];
sx q[0];
rz(-2.642785) q[0];
sx q[0];
rz(-0.98861277) q[0];
rz(-1.0722718) q[1];
sx q[1];
rz(-1.5360473) q[1];
sx q[1];
rz(-1.706121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0611183) q[0];
sx q[0];
rz(-2.1466564) q[0];
sx q[0];
rz(0.030035069) q[0];
rz(1.5631235) q[2];
sx q[2];
rz(-1.3425217) q[2];
sx q[2];
rz(1.829819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2769673) q[1];
sx q[1];
rz(-1.5709779) q[1];
sx q[1];
rz(1.5110043) q[1];
rz(1.9063453) q[3];
sx q[3];
rz(-0.90651262) q[3];
sx q[3];
rz(-1.2227525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84206218) q[2];
sx q[2];
rz(-2.5008423) q[2];
sx q[2];
rz(-0.93471849) q[2];
rz(2.0958021) q[3];
sx q[3];
rz(-0.17242923) q[3];
sx q[3];
rz(0.24833965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1187779) q[0];
sx q[0];
rz(-1.7775443) q[0];
sx q[0];
rz(-1.4015629) q[0];
rz(-0.10764311) q[1];
sx q[1];
rz(-1.5168774) q[1];
sx q[1];
rz(0.37227896) q[1];
rz(-2.8714928) q[2];
sx q[2];
rz(-0.80672376) q[2];
sx q[2];
rz(-3.0835019) q[2];
rz(2.6519898) q[3];
sx q[3];
rz(-2.1402749) q[3];
sx q[3];
rz(-0.28771852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
