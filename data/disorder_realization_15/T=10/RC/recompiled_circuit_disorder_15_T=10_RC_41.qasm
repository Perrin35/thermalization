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
rz(-1.7742046) q[1];
sx q[1];
rz(-2.8957638) q[1];
sx q[1];
rz(-2.153102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52553015) q[0];
sx q[0];
rz(-1.5126192) q[0];
sx q[0];
rz(-1.9921897) q[0];
rz(-pi) q[1];
rz(-1.3349873) q[2];
sx q[2];
rz(-0.60496444) q[2];
sx q[2];
rz(-1.6665104) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72585427) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(1.0456677) q[1];
rz(-pi) q[2];
rz(3.1036166) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(-1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6926379) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4028567) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(-2.1733213) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(-1.2526858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38329268) q[0];
sx q[0];
rz(-0.18342605) q[0];
sx q[0];
rz(0.51390506) q[0];
x q[1];
rz(-2.230174) q[2];
sx q[2];
rz(-2.0479925) q[2];
sx q[2];
rz(1.1148244) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15236552) q[1];
sx q[1];
rz(-2.3902635) q[1];
sx q[1];
rz(-2.3163124) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1329123) q[3];
sx q[3];
rz(-1.7403733) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0210555) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-0.34376124) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(0.37653157) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(2.4287756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3526488) q[0];
sx q[0];
rz(-1.3508571) q[0];
sx q[0];
rz(0.21531944) q[0];
rz(-1.9829282) q[2];
sx q[2];
rz(-0.40908989) q[2];
sx q[2];
rz(-1.2661753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3908927) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(-0.23551029) q[1];
x q[2];
rz(1.4180693) q[3];
sx q[3];
rz(-2.6543791) q[3];
sx q[3];
rz(-1.0917851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-0.66398579) q[2];
rz(1.9021696) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(2.5675039) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(0.92837292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5279816) q[0];
sx q[0];
rz(-2.4355542) q[0];
sx q[0];
rz(-1.5722949) q[0];
rz(-2.8027014) q[2];
sx q[2];
rz(-0.8182943) q[2];
sx q[2];
rz(-1.8213059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2239383) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(2.0018105) q[1];
rz(2.6392322) q[3];
sx q[3];
rz(-0.97548786) q[3];
sx q[3];
rz(-0.82752284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(-0.14285764) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083754152) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(-3.0673448) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(2.1898851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4518406) q[0];
sx q[0];
rz(-1.6432439) q[0];
sx q[0];
rz(-2.3206582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81163002) q[2];
sx q[2];
rz(-1.9981761) q[2];
sx q[2];
rz(0.4610093) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4068027) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(1.8491247) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1577155) q[3];
sx q[3];
rz(-0.65423274) q[3];
sx q[3];
rz(2.4901842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(1.1523694) q[2];
rz(-1.0460098) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(-3.1402821) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9933269) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(0.48962012) q[0];
rz(-1.024225) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(1.5354059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682966) q[0];
sx q[0];
rz(-1.6185986) q[0];
sx q[0];
rz(3.0424546) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9027009) q[2];
sx q[2];
rz(-1.0039582) q[2];
sx q[2];
rz(-1.6859646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8363758) q[1];
sx q[1];
rz(-1.7292542) q[1];
sx q[1];
rz(1.8578908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.51227) q[3];
sx q[3];
rz(-1.704708) q[3];
sx q[3];
rz(-0.32979345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(-0.17573389) q[2];
rz(-2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068129383) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(-2.9636256) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.7756745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38458347) q[0];
sx q[0];
rz(-0.88338822) q[0];
sx q[0];
rz(-2.0100726) q[0];
x q[1];
rz(-0.62909158) q[2];
sx q[2];
rz(-0.39789879) q[2];
sx q[2];
rz(1.9366022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9907896) q[1];
sx q[1];
rz(-1.3152796) q[1];
sx q[1];
rz(-1.3882511) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6947721) q[3];
sx q[3];
rz(-1.0448536) q[3];
sx q[3];
rz(-2.7391609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(-2.9979624) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(0.34415054) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550069) q[0];
sx q[0];
rz(-2.2877559) q[0];
sx q[0];
rz(-0.32178497) q[0];
rz(-2.2161662) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(-2.5245573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2935534) q[0];
sx q[0];
rz(-2.8176421) q[0];
sx q[0];
rz(2.241961) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7865137) q[2];
sx q[2];
rz(-1.6032013) q[2];
sx q[2];
rz(1.0499133) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.9279234) q[1];
sx q[1];
rz(-0.59487768) q[1];
sx q[1];
rz(0.87022123) q[1];
x q[2];
rz(-1.6156322) q[3];
sx q[3];
rz(-2.3196844) q[3];
sx q[3];
rz(-0.84079784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5499605) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(3.0275596) q[2];
rz(-0.36241254) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(-2.3214031) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(1.483451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950715) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(1.6540065) q[0];
rz(2.9786417) q[2];
sx q[2];
rz(-2.5451247) q[2];
sx q[2];
rz(-2.5359254) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8497613) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(2.989819) q[1];
rz(1.9547192) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(-0.25728713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(1.2443939) q[2];
rz(-2.7164298) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(-1.6311133) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72220951) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(-0.37049946) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(-1.8006905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271034) q[0];
sx q[0];
rz(-1.9368441) q[0];
sx q[0];
rz(-2.8949379) q[0];
rz(-2.9843763) q[2];
sx q[2];
rz(-2.4210857) q[2];
sx q[2];
rz(-0.23888982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(0.96149573) q[1];
rz(3.020346) q[3];
sx q[3];
rz(-1.5505425) q[3];
sx q[3];
rz(0.95758394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.6293388) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(-1.6188999) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(-1.6265097) q[2];
sx q[2];
rz(-0.40691661) q[2];
sx q[2];
rz(-0.42452068) q[2];
rz(2.7545746) q[3];
sx q[3];
rz(-1.1292463) q[3];
sx q[3];
rz(-2.675727) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
