OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5759597) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(2.5984882) q[0];
rz(-1.7898516) q[2];
sx q[2];
rz(-0.54974216) q[2];
sx q[2];
rz(1.4866536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5719205) q[1];
sx q[1];
rz(-1.306428) q[1];
sx q[1];
rz(0.7260679) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8957542) q[3];
sx q[3];
rz(-0.80245362) q[3];
sx q[3];
rz(-1.6368395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(2.6057459) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(-2.9931086) q[0];
rz(2.8500154) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(0.72327327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9745001) q[1];
sx q[1];
rz(-1.8247316) q[1];
sx q[1];
rz(2.7212935) q[1];
x q[2];
rz(-1.7137394) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(0.95834857) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(2.6843605) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31608554) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(-0.020629701) q[0];
rz(-pi) q[1];
rz(-0.089695887) q[2];
sx q[2];
rz(-1.2604453) q[2];
sx q[2];
rz(-2.2657564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64297134) q[1];
sx q[1];
rz(-1.1585304) q[1];
sx q[1];
rz(0.69795124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2529536) q[3];
sx q[3];
rz(-2.7405973) q[3];
sx q[3];
rz(-0.81438118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(2.5562111) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90081763) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6648383) q[0];
sx q[0];
rz(-1.1676482) q[0];
sx q[0];
rz(2.320015) q[0];
rz(-pi) q[1];
rz(1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(-1.1255217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.199898) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(-0.56337507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67646497) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(1.3456618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-1.0513603) q[0];
rz(1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2252786) q[0];
sx q[0];
rz(-1.3818704) q[0];
sx q[0];
rz(-0.96154763) q[0];
rz(-pi) q[1];
rz(-2.8565065) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(-1.8269055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6592308) q[1];
sx q[1];
rz(-1.8419957) q[1];
sx q[1];
rz(-1.2960474) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.352036) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.4279799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(-2.7894003) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(1.0587143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4689894) q[0];
sx q[0];
rz(-2.6575343) q[0];
sx q[0];
rz(2.6812535) q[0];
rz(0.58745678) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(-2.6641252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9204273) q[1];
sx q[1];
rz(-1.7943008) q[1];
sx q[1];
rz(1.8799053) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22535725) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.9865215) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(-2.2999433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8182897) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(-1.3099758) q[0];
rz(-pi) q[1];
rz(-0.26288962) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(2.7123244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0591653) q[1];
sx q[1];
rz(-1.6470243) q[1];
sx q[1];
rz(0.04870292) q[1];
x q[2];
rz(-2.3119885) q[3];
sx q[3];
rz(-1.7332352) q[3];
sx q[3];
rz(0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-2.9522827) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(2.7239674) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4293489) q[0];
sx q[0];
rz(-1.1375543) q[0];
sx q[0];
rz(0.54746763) q[0];
rz(-pi) q[1];
rz(-1.8629486) q[2];
sx q[2];
rz(-0.38049618) q[2];
sx q[2];
rz(-0.7064864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8199181) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(0.8701156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91248625) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(-2.130079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(1.0428838) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-2.1946857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.2456018) q[0];
sx q[0];
rz(-1.8878216) q[0];
rz(-1.5309179) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(-3.1181042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0526035) q[1];
sx q[1];
rz(-0.34480428) q[1];
sx q[1];
rz(-0.25119541) q[1];
rz(-pi) q[2];
rz(-1.9584719) q[3];
sx q[3];
rz(-0.47260731) q[3];
sx q[3];
rz(-1.2848867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6039186) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(0.3084348) q[0];
x q[1];
rz(-1.3078493) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(-0.76542379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88313738) q[1];
sx q[1];
rz(-1.9069792) q[1];
sx q[1];
rz(-1.2698445) q[1];
rz(2.6976922) q[3];
sx q[3];
rz(-1.1392986) q[3];
sx q[3];
rz(0.75507009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(-0.59990668) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.29466378) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(1.7462126) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(0.89623981) q[3];
sx q[3];
rz(-1.2590209) q[3];
sx q[3];
rz(-0.38929064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
