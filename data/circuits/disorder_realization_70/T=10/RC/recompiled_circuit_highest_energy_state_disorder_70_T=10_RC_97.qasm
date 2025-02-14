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
rz(1.7477859) q[0];
sx q[0];
rz(4.0417606) q[0];
sx q[0];
rz(10.834122) q[0];
rz(0.70977587) q[1];
sx q[1];
rz(-2.8083399) q[1];
sx q[1];
rz(0.38212734) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2173806) q[0];
sx q[0];
rz(-1.7048302) q[0];
sx q[0];
rz(-2.4045375) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3163739) q[2];
sx q[2];
rz(-1.7013936) q[2];
sx q[2];
rz(-0.25564627) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.00013) q[1];
sx q[1];
rz(-1.9381818) q[1];
sx q[1];
rz(-2.7482583) q[1];
rz(3.0218533) q[3];
sx q[3];
rz(-2.2594707) q[3];
sx q[3];
rz(-0.39018351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0742775) q[2];
sx q[2];
rz(-1.5019608) q[2];
sx q[2];
rz(0.074946694) q[2];
rz(1.8986757) q[3];
sx q[3];
rz(-2.8800745) q[3];
sx q[3];
rz(-0.3956795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260929) q[0];
sx q[0];
rz(-0.41985303) q[0];
sx q[0];
rz(-2.5123151) q[0];
rz(-0.83726007) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(-0.57964051) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.048174) q[0];
sx q[0];
rz(-0.77147986) q[0];
sx q[0];
rz(2.8475965) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7461065) q[2];
sx q[2];
rz(-1.3278393) q[2];
sx q[2];
rz(-2.0045351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2036248) q[1];
sx q[1];
rz(-1.9589099) q[1];
sx q[1];
rz(-1.2118503) q[1];
rz(-2.9778538) q[3];
sx q[3];
rz(-1.5511912) q[3];
sx q[3];
rz(2.1480008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2418182) q[2];
sx q[2];
rz(-0.16336975) q[2];
sx q[2];
rz(-2.3972798) q[2];
rz(-2.7742079) q[3];
sx q[3];
rz(-1.8495879) q[3];
sx q[3];
rz(-1.7261837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431417) q[0];
sx q[0];
rz(-2.1911868) q[0];
sx q[0];
rz(-2.1344192) q[0];
rz(-0.011064359) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(-0.99753582) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384093) q[0];
sx q[0];
rz(-1.5965726) q[0];
sx q[0];
rz(-1.4332214) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4513955) q[2];
sx q[2];
rz(-2.3755671) q[2];
sx q[2];
rz(-1.8325782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1147194) q[1];
sx q[1];
rz(-2.623154) q[1];
sx q[1];
rz(0.015109574) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51181958) q[3];
sx q[3];
rz(-1.5835754) q[3];
sx q[3];
rz(-0.46841533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.02115383) q[2];
sx q[2];
rz(-0.28205559) q[2];
sx q[2];
rz(0.8417449) q[2];
rz(-2.4833931) q[3];
sx q[3];
rz(-0.87116146) q[3];
sx q[3];
rz(0.021520821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4206674) q[0];
sx q[0];
rz(-0.1249211) q[0];
sx q[0];
rz(-0.7290054) q[0];
rz(-0.76338243) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(-0.36300945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0057392) q[0];
sx q[0];
rz(-2.7518138) q[0];
sx q[0];
rz(-0.71758349) q[0];
rz(-pi) q[1];
rz(2.1562026) q[2];
sx q[2];
rz(-2.3040758) q[2];
sx q[2];
rz(-0.68345165) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91528401) q[1];
sx q[1];
rz(-1.7123509) q[1];
sx q[1];
rz(-0.29967163) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7458526) q[3];
sx q[3];
rz(-3.0043663) q[3];
sx q[3];
rz(-1.2089704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9558941) q[2];
sx q[2];
rz(-0.82315767) q[2];
sx q[2];
rz(1.4996747) q[2];
rz(2.5540292) q[3];
sx q[3];
rz(-0.9891808) q[3];
sx q[3];
rz(-2.6232918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8374306) q[0];
sx q[0];
rz(-0.70697933) q[0];
sx q[0];
rz(-2.8564603) q[0];
rz(0.25310165) q[1];
sx q[1];
rz(-2.1123835) q[1];
sx q[1];
rz(-1.0458127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8305215) q[0];
sx q[0];
rz(-1.1806628) q[0];
sx q[0];
rz(-2.5907787) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1509456) q[2];
sx q[2];
rz(-1.9774745) q[2];
sx q[2];
rz(0.57957725) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95438984) q[1];
sx q[1];
rz(-2.7437401) q[1];
sx q[1];
rz(2.9955762) q[1];
rz(-2.8941514) q[3];
sx q[3];
rz(-0.69290042) q[3];
sx q[3];
rz(0.019236658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4046341) q[2];
sx q[2];
rz(-1.0598695) q[2];
sx q[2];
rz(-2.5671379) q[2];
rz(-0.61070853) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33559281) q[0];
sx q[0];
rz(-1.7570423) q[0];
sx q[0];
rz(2.7912676) q[0];
rz(0.2925182) q[1];
sx q[1];
rz(-0.12648335) q[1];
sx q[1];
rz(-2.3622321) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85920632) q[0];
sx q[0];
rz(-3.0090158) q[0];
sx q[0];
rz(-2.0280806) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6616277) q[2];
sx q[2];
rz(-1.4128608) q[2];
sx q[2];
rz(-1.2564893) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36679572) q[1];
sx q[1];
rz(-1.8908943) q[1];
sx q[1];
rz(1.3572378) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.73137) q[3];
sx q[3];
rz(-1.055937) q[3];
sx q[3];
rz(2.5405879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23955762) q[2];
sx q[2];
rz(-1.6263447) q[2];
sx q[2];
rz(0.96417344) q[2];
rz(2.9686109) q[3];
sx q[3];
rz(-2.1292584) q[3];
sx q[3];
rz(-3.0103736) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59703374) q[0];
sx q[0];
rz(-1.1976765) q[0];
sx q[0];
rz(0.069393754) q[0];
rz(1.3736877) q[1];
sx q[1];
rz(-1.2488139) q[1];
sx q[1];
rz(-0.51838851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1590879) q[0];
sx q[0];
rz(-2.2347274) q[0];
sx q[0];
rz(0.21592617) q[0];
rz(2.0608299) q[2];
sx q[2];
rz(-1.9980717) q[2];
sx q[2];
rz(-3.0025122) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1063819) q[1];
sx q[1];
rz(-1.7089881) q[1];
sx q[1];
rz(-1.5968678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9767738) q[3];
sx q[3];
rz(-1.6281152) q[3];
sx q[3];
rz(0.27637339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18846506) q[2];
sx q[2];
rz(-1.1981107) q[2];
sx q[2];
rz(2.8971064) q[2];
rz(1.9237579) q[3];
sx q[3];
rz(-0.094450258) q[3];
sx q[3];
rz(0.34734669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090488) q[0];
sx q[0];
rz(-2.8489887) q[0];
sx q[0];
rz(-0.69945139) q[0];
rz(-2.4199016) q[1];
sx q[1];
rz(-1.3612008) q[1];
sx q[1];
rz(-1.9500505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9129097) q[0];
sx q[0];
rz(-1.7828724) q[0];
sx q[0];
rz(1.4430844) q[0];
x q[1];
rz(-1.4100685) q[2];
sx q[2];
rz(-1.3597915) q[2];
sx q[2];
rz(2.4461022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4556132) q[1];
sx q[1];
rz(-1.0334823) q[1];
sx q[1];
rz(0.098453589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6939387) q[3];
sx q[3];
rz(-1.2546243) q[3];
sx q[3];
rz(-2.3929736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9222074) q[2];
sx q[2];
rz(-2.3228513) q[2];
sx q[2];
rz(-1.7259664) q[2];
rz(2.7042232) q[3];
sx q[3];
rz(-0.24961095) q[3];
sx q[3];
rz(0.50284809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9472083) q[0];
sx q[0];
rz(-1.2057065) q[0];
sx q[0];
rz(0.44596392) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-2.4868592) q[1];
sx q[1];
rz(1.9816678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98731536) q[0];
sx q[0];
rz(-1.9145218) q[0];
sx q[0];
rz(2.4489342) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.075322876) q[2];
sx q[2];
rz(-1.8315856) q[2];
sx q[2];
rz(-1.0950077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7912631) q[1];
sx q[1];
rz(-1.3379249) q[1];
sx q[1];
rz(-2.5901152) q[1];
rz(-pi) q[2];
rz(2.5612513) q[3];
sx q[3];
rz(-0.99513061) q[3];
sx q[3];
rz(-2.2656296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40598536) q[2];
sx q[2];
rz(-1.0286101) q[2];
sx q[2];
rz(2.410991) q[2];
rz(2.2405911) q[3];
sx q[3];
rz(-2.7121219) q[3];
sx q[3];
rz(-0.034339529) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5071097) q[0];
sx q[0];
rz(-0.49840885) q[0];
sx q[0];
rz(-2.679017) q[0];
rz(2.0787461) q[1];
sx q[1];
rz(-1.3674419) q[1];
sx q[1];
rz(0.06632334) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0339113) q[0];
sx q[0];
rz(-2.4723791) q[0];
sx q[0];
rz(2.6341556) q[0];
x q[1];
rz(-2.0424329) q[2];
sx q[2];
rz(-2.7688469) q[2];
sx q[2];
rz(2.3605704) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6219225) q[1];
sx q[1];
rz(-1.038045) q[1];
sx q[1];
rz(0.68343917) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22610449) q[3];
sx q[3];
rz(-1.1563099) q[3];
sx q[3];
rz(1.3785441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1041717) q[2];
sx q[2];
rz(-2.6080242) q[2];
sx q[2];
rz(1.9083692) q[2];
rz(-2.6836266) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(-0.39877322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.11252277) q[0];
sx q[0];
rz(-1.7081013) q[0];
sx q[0];
rz(-1.3992455) q[0];
rz(-0.15432547) q[1];
sx q[1];
rz(-1.4230774) q[1];
sx q[1];
rz(1.9565061) q[1];
rz(-1.1823282) q[2];
sx q[2];
rz(-1.8618089) q[2];
sx q[2];
rz(2.1258327) q[2];
rz(-2.1460228) q[3];
sx q[3];
rz(-1.9862666) q[3];
sx q[3];
rz(3.1131328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
