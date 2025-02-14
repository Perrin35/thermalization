OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.14804949) q[0];
sx q[0];
rz(-0.43819675) q[0];
sx q[0];
rz(0.73102695) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(-2.5875523) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06023887) q[0];
sx q[0];
rz(-2.0610272) q[0];
sx q[0];
rz(-2.1054563) q[0];
x q[1];
rz(1.117716) q[2];
sx q[2];
rz(-2.6181698) q[2];
sx q[2];
rz(-2.7414309) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82137992) q[1];
sx q[1];
rz(-2.2578635) q[1];
sx q[1];
rz(0.043619556) q[1];
x q[2];
rz(-1.0022121) q[3];
sx q[3];
rz(-0.7665638) q[3];
sx q[3];
rz(-2.1994906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76250166) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(-0.16538922) q[2];
rz(2.9108544) q[3];
sx q[3];
rz(-0.24300353) q[3];
sx q[3];
rz(0.86133426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.47925258) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(1.3586556) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(-0.84017909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5206487) q[0];
sx q[0];
rz(-2.0387702) q[0];
sx q[0];
rz(-2.8415989) q[0];
rz(-pi) q[1];
rz(-1.1112059) q[2];
sx q[2];
rz(-0.75558096) q[2];
sx q[2];
rz(-0.43967512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11586861) q[1];
sx q[1];
rz(-2.3626973) q[1];
sx q[1];
rz(0.42515691) q[1];
rz(1.718037) q[3];
sx q[3];
rz(-2.3309532) q[3];
sx q[3];
rz(2.7642706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.75300616) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(2.2071655) q[2];
rz(1.8560483) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(-2.2540895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282106) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(-1.5895948) q[0];
rz(-1.7734843) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(-2.0210463) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89884886) q[0];
sx q[0];
rz(-1.1087497) q[0];
sx q[0];
rz(0.36926134) q[0];
rz(1.1324544) q[2];
sx q[2];
rz(-1.8016644) q[2];
sx q[2];
rz(-1.9491652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4996329) q[1];
sx q[1];
rz(-1.4781535) q[1];
sx q[1];
rz(-0.53643061) q[1];
rz(-1.2854101) q[3];
sx q[3];
rz(-2.719398) q[3];
sx q[3];
rz(-0.58870047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0997194) q[2];
sx q[2];
rz(-2.9283044) q[2];
sx q[2];
rz(-0.29207692) q[2];
rz(1.9000351) q[3];
sx q[3];
rz(-1.7009267) q[3];
sx q[3];
rz(0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9860155) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(2.7098932) q[0];
rz(-0.15402928) q[1];
sx q[1];
rz(-1.5700211) q[1];
sx q[1];
rz(-1.0607176) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3531279) q[0];
sx q[0];
rz(-1.2305088) q[0];
sx q[0];
rz(2.8022604) q[0];
x q[1];
rz(-0.94580146) q[2];
sx q[2];
rz(-0.29871395) q[2];
sx q[2];
rz(-2.309989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4894197) q[1];
sx q[1];
rz(-1.2448346) q[1];
sx q[1];
rz(0.23613813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4522375) q[3];
sx q[3];
rz(-0.24447799) q[3];
sx q[3];
rz(2.1357128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3920307) q[2];
sx q[2];
rz(-1.6101086) q[2];
sx q[2];
rz(-0.57382601) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.4797689) q[3];
sx q[3];
rz(-2.3315232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60657984) q[0];
sx q[0];
rz(-0.25595328) q[0];
sx q[0];
rz(0.25273299) q[0];
rz(-2.9622954) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.4089233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64302639) q[0];
sx q[0];
rz(-1.6448978) q[0];
sx q[0];
rz(1.3331494) q[0];
x q[1];
rz(-0.73136605) q[2];
sx q[2];
rz(-2.593267) q[2];
sx q[2];
rz(0.62125118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.682413) q[1];
sx q[1];
rz(-2.3485567) q[1];
sx q[1];
rz(-0.17211087) q[1];
rz(2.5645676) q[3];
sx q[3];
rz(-2.5379556) q[3];
sx q[3];
rz(-1.1567504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1285642) q[2];
sx q[2];
rz(-1.3821673) q[2];
sx q[2];
rz(1.2010835) q[2];
rz(2.3342093) q[3];
sx q[3];
rz(-1.7816593) q[3];
sx q[3];
rz(-2.0090296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0627237) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(-0.8412745) q[0];
rz(1.2048644) q[1];
sx q[1];
rz(-1.3179904) q[1];
sx q[1];
rz(1.0296317) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1372927) q[0];
sx q[0];
rz(-0.51057112) q[0];
sx q[0];
rz(-1.592203) q[0];
rz(0.44544051) q[2];
sx q[2];
rz(-2.5800319) q[2];
sx q[2];
rz(-0.5024006) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8334044) q[1];
sx q[1];
rz(-0.17886111) q[1];
sx q[1];
rz(-0.043828242) q[1];
x q[2];
rz(-1.383841) q[3];
sx q[3];
rz(-1.2724981) q[3];
sx q[3];
rz(-1.3307856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.131375) q[2];
sx q[2];
rz(-2.3848332) q[2];
sx q[2];
rz(-2.6436464) q[2];
rz(-2.61854) q[3];
sx q[3];
rz(-1.7224576) q[3];
sx q[3];
rz(-0.12652346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2987357) q[0];
sx q[0];
rz(-0.13164483) q[0];
sx q[0];
rz(-2.329622) q[0];
rz(0.89093351) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-2.1422211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1645714) q[0];
sx q[0];
rz(-1.1081105) q[0];
sx q[0];
rz(-2.0624731) q[0];
rz(-pi) q[1];
rz(-0.81951253) q[2];
sx q[2];
rz(-3.0009821) q[2];
sx q[2];
rz(-2.484172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68518406) q[1];
sx q[1];
rz(-1.0899836) q[1];
sx q[1];
rz(2.494704) q[1];
rz(0.90535594) q[3];
sx q[3];
rz(-1.8017343) q[3];
sx q[3];
rz(-3.0311716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.823395) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(-1.8434175) q[2];
rz(3.1031109) q[3];
sx q[3];
rz(-2.0352071) q[3];
sx q[3];
rz(1.5255671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0549952) q[0];
sx q[0];
rz(-2.0107858) q[0];
sx q[0];
rz(2.8259592) q[0];
rz(1.9075958) q[1];
sx q[1];
rz(-2.4256568) q[1];
sx q[1];
rz(1.8826145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9568898) q[0];
sx q[0];
rz(-1.8543921) q[0];
sx q[0];
rz(-2.3882939) q[0];
rz(2.9614212) q[2];
sx q[2];
rz(-0.77356662) q[2];
sx q[2];
rz(-0.13828466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8371305) q[1];
sx q[1];
rz(-1.9726957) q[1];
sx q[1];
rz(-2.7086651) q[1];
rz(-pi) q[2];
rz(-3.0837667) q[3];
sx q[3];
rz(-3.0396121) q[3];
sx q[3];
rz(-2.0174055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2235609) q[2];
sx q[2];
rz(-0.74350244) q[2];
sx q[2];
rz(-0.027916748) q[2];
rz(2.3935086) q[3];
sx q[3];
rz(-1.4192162) q[3];
sx q[3];
rz(2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7276723) q[0];
sx q[0];
rz(-0.99069178) q[0];
sx q[0];
rz(0.45898166) q[0];
rz(-2.0052295) q[1];
sx q[1];
rz(-1.7106979) q[1];
sx q[1];
rz(2.9232025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9680351) q[0];
sx q[0];
rz(-1.5941248) q[0];
sx q[0];
rz(1.5541881) q[0];
x q[1];
rz(-1.893546) q[2];
sx q[2];
rz(-0.40229978) q[2];
sx q[2];
rz(-0.93076555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23109197) q[1];
sx q[1];
rz(-1.0808966) q[1];
sx q[1];
rz(1.4040158) q[1];
rz(-pi) q[2];
rz(-1.9588542) q[3];
sx q[3];
rz(-1.4493296) q[3];
sx q[3];
rz(1.1282008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2961262) q[2];
sx q[2];
rz(-0.25456905) q[2];
sx q[2];
rz(1.0477585) q[2];
rz(-2.3962077) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(1.6489702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02136852) q[0];
sx q[0];
rz(-2.5701018) q[0];
sx q[0];
rz(2.7987203) q[0];
rz(0.19035569) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(2.6284133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7452538) q[0];
sx q[0];
rz(-1.7293735) q[0];
sx q[0];
rz(-2.2157085) q[0];
rz(-2.1634899) q[2];
sx q[2];
rz(-2.5524271) q[2];
sx q[2];
rz(1.9798673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3575705) q[1];
sx q[1];
rz(-2.5139132) q[1];
sx q[1];
rz(-0.72369544) q[1];
x q[2];
rz(-1.9916198) q[3];
sx q[3];
rz(-1.4412035) q[3];
sx q[3];
rz(2.0146927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43789431) q[2];
sx q[2];
rz(-2.5056705) q[2];
sx q[2];
rz(-0.61100125) q[2];
rz(2.8827187) q[3];
sx q[3];
rz(-1.2114108) q[3];
sx q[3];
rz(-1.3795616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3872869) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(-2.1517131) q[1];
sx q[1];
rz(-1.2073333) q[1];
sx q[1];
rz(-2.5055199) q[1];
rz(0.0089916747) q[2];
sx q[2];
rz(-1.7594382) q[2];
sx q[2];
rz(0.3450763) q[2];
rz(-1.5045139) q[3];
sx q[3];
rz(-0.87285973) q[3];
sx q[3];
rz(-1.9535337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
