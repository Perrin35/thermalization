OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29785922) q[0];
sx q[0];
rz(3.7552667) q[0];
sx q[0];
rz(11.847191) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2253217) q[0];
sx q[0];
rz(-2.7164408) q[0];
sx q[0];
rz(-1.7122373) q[0];
x q[1];
rz(-2.1626986) q[2];
sx q[2];
rz(-1.7040633) q[2];
sx q[2];
rz(2.8507581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72585427) q[1];
sx q[1];
rz(-2.038318) q[1];
sx q[1];
rz(2.095925) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8139008) q[3];
sx q[3];
rz(-1.5339359) q[3];
sx q[3];
rz(2.778361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44895479) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(-0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(-2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.73873591) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(0.36303315) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(-1.8889069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371212) q[0];
sx q[0];
rz(-1.7303109) q[0];
sx q[0];
rz(-1.479854) q[0];
x q[1];
rz(0.91141869) q[2];
sx q[2];
rz(-2.0479925) q[2];
sx q[2];
rz(-2.0267682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0536086) q[1];
sx q[1];
rz(-1.0454185) q[1];
sx q[1];
rz(-2.5768075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1329123) q[3];
sx q[3];
rz(-1.4012194) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0210555) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(-3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(2.4287756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8757652) q[0];
sx q[0];
rz(-1.7808502) q[0];
sx q[0];
rz(-1.7957627) q[0];
x q[1];
rz(-1.1586645) q[2];
sx q[2];
rz(-0.40908989) q[2];
sx q[2];
rz(1.2661753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7506999) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(-0.23551029) q[1];
rz(-0.080428877) q[3];
sx q[3];
rz(-2.0518528) q[3];
sx q[3];
rz(-1.8773432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5793844) q[2];
sx q[2];
rz(2.4776069) q[2];
rz(-1.9021696) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.5749213) q[0];
sx q[0];
rz(-3.1118588) q[0];
sx q[0];
rz(-0.57408875) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(-2.2132197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95604491) q[0];
sx q[0];
rz(-1.5717686) q[0];
sx q[0];
rz(-2.2768343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8027014) q[2];
sx q[2];
rz(-2.3232984) q[2];
sx q[2];
rz(-1.3202867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91765431) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(2.0018105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5023605) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(2.3140698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(-0.14285764) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(-0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-0.074247867) q[0];
rz(-1.6429365) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(0.9517076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18626801) q[0];
sx q[0];
rz(-0.82337472) q[0];
sx q[0];
rz(3.0427409) q[0];
rz(-pi) q[1];
rz(-0.56065083) q[2];
sx q[2];
rz(-0.89386212) q[2];
sx q[2];
rz(-1.4844984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73478991) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(-1.2924679) q[1];
x q[2];
rz(-2.1577155) q[3];
sx q[3];
rz(-2.4873599) q[3];
sx q[3];
rz(-0.65140843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.1523694) q[2];
rz(2.0955829) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(-3.1402821) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(0.48962012) q[0];
rz(-1.024225) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.6061868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55035704) q[0];
sx q[0];
rz(-0.11002692) q[0];
sx q[0];
rz(-0.45022924) q[0];
rz(-pi) q[1];
rz(2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(-0.51598179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9163497) q[1];
sx q[1];
rz(-2.8147329) q[1];
sx q[1];
rz(-1.0570231) q[1];
rz(-pi) q[2];
rz(-2.9165886) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9635222) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(2.9658588) q[2];
rz(1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734633) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(-1.3659182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1195927) q[0];
sx q[0];
rz(-2.3455142) q[0];
sx q[0];
rz(2.6636366) q[0];
rz(-pi) q[1];
rz(1.8132509) q[2];
sx q[2];
rz(-1.889466) q[2];
sx q[2];
rz(-2.6048425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3599638) q[1];
sx q[1];
rz(-2.8287323) q[1];
sx q[1];
rz(-0.60731213) q[1];
rz(-pi) q[2];
rz(-0.52929438) q[3];
sx q[3];
rz(-1.4636453) q[3];
sx q[3];
rz(-1.2308434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(-0.89795566) q[2];
rz(0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.9550069) q[0];
sx q[0];
rz(-2.2877559) q[0];
sx q[0];
rz(0.32178497) q[0];
rz(-2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(0.61703533) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962284) q[0];
sx q[0];
rz(-1.3188688) q[0];
sx q[0];
rz(-0.20586254) q[0];
x q[1];
rz(1.4204942) q[2];
sx q[2];
rz(-0.21810025) q[2];
sx q[2];
rz(2.4739166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.033213381) q[1];
sx q[1];
rz(-1.2011659) q[1];
sx q[1];
rz(-2.0481678) q[1];
rz(-pi) q[2];
rz(3.0934107) q[3];
sx q[3];
rz(-2.391624) q[3];
sx q[3];
rz(0.90660209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5499605) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(1.6581416) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9634919) q[0];
sx q[0];
rz(-3.0577457) q[0];
sx q[0];
rz(-1.447178) q[0];
x q[1];
rz(-1.4610897) q[2];
sx q[2];
rz(-0.98329558) q[2];
sx q[2];
rz(2.7320931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2918313) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(0.15177365) q[1];
rz(1.9547192) q[3];
sx q[3];
rz(-1.0537638) q[3];
sx q[3];
rz(-2.8843055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(1.8971987) q[2];
rz(-0.42516285) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(-1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72220951) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(1.8006905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140708) q[0];
sx q[0];
rz(-0.43826575) q[0];
sx q[0];
rz(1.0036432) q[0];
rz(-1.4341899) q[2];
sx q[2];
rz(-0.86106664) q[2];
sx q[2];
rz(3.1106069) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(-0.96149573) q[1];
rz(-pi) q[2];
rz(-3.020346) q[3];
sx q[3];
rz(-1.5910501) q[3];
sx q[3];
rz(-2.1840087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(-0.86087888) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027949) q[0];
sx q[0];
rz(-1.4807985) q[0];
sx q[0];
rz(2.280622) q[0];
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.6265097) q[2];
sx q[2];
rz(-2.734676) q[2];
sx q[2];
rz(2.717072) q[2];
rz(0.38701804) q[3];
sx q[3];
rz(-2.0123464) q[3];
sx q[3];
rz(0.46586566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];