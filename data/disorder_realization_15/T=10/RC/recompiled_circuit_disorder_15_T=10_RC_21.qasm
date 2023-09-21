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
rz(-2.5279186) q[0];
sx q[0];
rz(2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(2.8957638) q[1];
sx q[1];
rz(10.413269) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0713232) q[0];
sx q[0];
rz(-1.9914314) q[0];
sx q[0];
rz(3.0778528) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97889401) q[2];
sx q[2];
rz(-1.7040633) q[2];
sx q[2];
rz(-2.8507581) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4157384) q[1];
sx q[1];
rz(-2.038318) q[1];
sx q[1];
rz(-1.0456677) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.037976102) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(1.2167041) q[3];
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
rz(-2.5732178) q[2];
rz(1.9521936) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-2.7785595) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(-1.8889069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371212) q[0];
sx q[0];
rz(-1.7303109) q[0];
sx q[0];
rz(-1.479854) q[0];
rz(-pi) q[1];
rz(2.5622796) q[2];
sx q[2];
rz(-0.99537731) q[2];
sx q[2];
rz(2.343611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0536086) q[1];
sx q[1];
rz(-1.0454185) q[1];
sx q[1];
rz(2.5768075) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1870175) q[3];
sx q[3];
rz(-0.46758258) q[3];
sx q[3];
rz(-0.13320696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0210555) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6702328) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(-0.37653157) q[1];
sx q[1];
rz(-0.9286325) q[1];
sx q[1];
rz(0.71281707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56601554) q[0];
sx q[0];
rz(-0.30656719) q[0];
sx q[0];
rz(2.3335639) q[0];
rz(1.9829282) q[2];
sx q[2];
rz(-0.40908989) q[2];
sx q[2];
rz(-1.8754174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3908927) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(2.9060824) q[1];
x q[2];
rz(-1.0884104) q[3];
sx q[3];
rz(-1.4995121) q[3];
sx q[3];
rz(-2.7977668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3537447) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-2.4776069) q[2];
rz(-1.239423) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(2.2263039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(-0.92837292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61361109) q[0];
sx q[0];
rz(-2.4355542) q[0];
sx q[0];
rz(1.5692977) q[0];
x q[1];
rz(-2.8027014) q[2];
sx q[2];
rz(-2.3232984) q[2];
sx q[2];
rz(1.8213059) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1410364) q[1];
sx q[1];
rz(-0.49266854) q[1];
sx q[1];
rz(2.1129184) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5023605) q[3];
sx q[3];
rz(-0.97548786) q[3];
sx q[3];
rz(2.3140698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(-2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0578385) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(-0.074247867) q[0];
rz(-1.6429365) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18626801) q[0];
sx q[0];
rz(-2.3182179) q[0];
sx q[0];
rz(0.098851725) q[0];
rz(-2.5809418) q[2];
sx q[2];
rz(-0.89386212) q[2];
sx q[2];
rz(1.4844984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4068027) q[1];
sx q[1];
rz(-0.64195913) q[1];
sx q[1];
rz(1.8491247) q[1];
rz(-pi) q[2];
rz(2.1577155) q[3];
sx q[3];
rz(-0.65423274) q[3];
sx q[3];
rz(2.4901842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64530659) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.1523694) q[2];
rz(2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9933269) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(-0.48962012) q[0];
rz(-2.1173677) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.5354059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55035704) q[0];
sx q[0];
rz(-0.11002692) q[0];
sx q[0];
rz(2.6913634) q[0];
x q[1];
rz(2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(-0.51598179) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8363758) q[1];
sx q[1];
rz(-1.7292542) q[1];
sx q[1];
rz(1.8578908) q[1];
rz(2.9165886) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(-1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(2.9658588) q[2];
rz(-1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(-0.1779671) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(1.7756745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89643909) q[0];
sx q[0];
rz(-1.23587) q[0];
sx q[0];
rz(0.73672898) q[0];
rz(-pi) q[1];
rz(-1.3283417) q[2];
sx q[2];
rz(-1.889466) q[2];
sx q[2];
rz(-2.6048425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7682225) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(-0.25964398) q[1];
rz(-pi) q[2];
rz(1.6947721) q[3];
sx q[3];
rz(-2.0967391) q[3];
sx q[3];
rz(2.7391609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0718677) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(-0.89795566) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(0.32178497) q[0];
rz(-0.92542648) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641545) q[0];
sx q[0];
rz(-1.7700717) q[0];
sx q[0];
rz(-1.827924) q[0];
x q[1];
rz(-1.355079) q[2];
sx q[2];
rz(-1.5383913) q[2];
sx q[2];
rz(2.0916794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1083793) q[1];
sx q[1];
rz(-1.2011659) q[1];
sx q[1];
rz(2.0481678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3922032) q[3];
sx q[3];
rz(-1.5379616) q[3];
sx q[3];
rz(2.4421305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(0.36241254) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.9053649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(-1.6581416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950715) q[0];
sx q[0];
rz(-1.5811231) q[0];
sx q[0];
rz(-1.6540065) q[0];
rz(-pi) q[1];
rz(-0.59028583) q[2];
sx q[2];
rz(-1.4795408) q[2];
sx q[2];
rz(1.1003189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22944731) q[1];
sx q[1];
rz(-1.427269) q[1];
sx q[1];
rz(-1.2372644) q[1];
rz(-pi) q[2];
rz(-2.5591194) q[3];
sx q[3];
rz(-2.5081722) q[3];
sx q[3];
rz(2.7137091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.2443939) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(1.6311133) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(1.1100948) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(1.8006905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140708) q[0];
sx q[0];
rz(-0.43826575) q[0];
sx q[0];
rz(1.0036432) q[0];
rz(-pi) q[1];
rz(1.7074028) q[2];
sx q[2];
rz(-0.86106664) q[2];
sx q[2];
rz(-0.030985706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1827613) q[1];
sx q[1];
rz(-1.3211831) q[1];
sx q[1];
rz(1.9447437) q[1];
rz(-pi) q[2];
rz(2.9756536) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(-0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(-0.86087888) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(-1.6265097) q[2];
sx q[2];
rz(-0.40691661) q[2];
sx q[2];
rz(-0.42452068) q[2];
rz(-2.2446185) q[3];
sx q[3];
rz(-2.5629808) q[3];
sx q[3];
rz(-1.9140011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];