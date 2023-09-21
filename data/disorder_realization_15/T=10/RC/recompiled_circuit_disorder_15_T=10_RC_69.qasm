OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(0.71917978) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(-0.9884907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0713232) q[0];
sx q[0];
rz(-1.1501612) q[0];
sx q[0];
rz(-3.0778528) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3349873) q[2];
sx q[2];
rz(-0.60496444) q[2];
sx q[2];
rz(-1.6665104) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5520652) q[1];
sx q[1];
rz(-2.0347934) q[1];
sx q[1];
rz(-0.5281756) q[1];
rz(-pi) q[2];
x q[2];
rz(1.418781) q[3];
sx q[3];
rz(-0.24582836) q[3];
sx q[3];
rz(-2.0815401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44895479) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(-2.5732178) q[2];
rz(1.9521936) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4028567) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.2526858) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68080901) q[0];
sx q[0];
rz(-1.660581) q[0];
sx q[0];
rz(-0.16016527) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91141869) q[2];
sx q[2];
rz(-1.0936001) q[2];
sx q[2];
rz(2.0267682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0536086) q[1];
sx q[1];
rz(-2.0961742) q[1];
sx q[1];
rz(2.5768075) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9545752) q[3];
sx q[3];
rz(-0.46758258) q[3];
sx q[3];
rz(0.13320696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(-3.1345471) q[0];
rz(-0.37653157) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-0.71281707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7889439) q[0];
sx q[0];
rz(-1.7907356) q[0];
sx q[0];
rz(-2.9262732) q[0];
rz(-pi) q[1];
rz(1.9489281) q[2];
sx q[2];
rz(-1.7308123) q[2];
sx q[2];
rz(0.68607054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3908927) q[1];
sx q[1];
rz(-1.0415959) q[1];
sx q[1];
rz(-0.23551029) q[1];
x q[2];
rz(-0.080428877) q[3];
sx q[3];
rz(-2.0518528) q[3];
sx q[3];
rz(1.2642494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(0.66398579) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(0.57408875) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-0.92837292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95604491) q[0];
sx q[0];
rz(-1.5717686) q[0];
sx q[0];
rz(-2.2768343) q[0];
x q[1];
rz(2.8027014) q[2];
sx q[2];
rz(-2.3232984) q[2];
sx q[2];
rz(1.3202867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6002027) q[1];
sx q[1];
rz(-1.9879513) q[1];
sx q[1];
rz(-0.27020176) q[1];
rz(-0.91289642) q[3];
sx q[3];
rz(-1.1607338) q[3];
sx q[3];
rz(1.0421154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(1.1553923) q[2];
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
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083754152) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(-3.0673448) q[0];
rz(1.4986562) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-0.9517076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9553246) q[0];
sx q[0];
rz(-0.82337472) q[0];
sx q[0];
rz(3.0427409) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81163002) q[2];
sx q[2];
rz(-1.9981761) q[2];
sx q[2];
rz(-0.4610093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0641099) q[1];
sx q[1];
rz(-2.1842842) q[1];
sx q[1];
rz(0.20258278) q[1];
x q[2];
rz(0.9838772) q[3];
sx q[3];
rz(-2.4873599) q[3];
sx q[3];
rz(-0.65140843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.9933269) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(0.48962012) q[0];
rz(1.024225) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(-1.5354059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55035704) q[0];
sx q[0];
rz(-3.0315657) q[0];
sx q[0];
rz(-2.6913634) q[0];
rz(-0.9027009) q[2];
sx q[2];
rz(-1.0039582) q[2];
sx q[2];
rz(-1.6859646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3052169) q[1];
sx q[1];
rz(-1.7292542) q[1];
sx q[1];
rz(-1.2837019) q[1];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(2.0819506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(2.9658588) q[2];
rz(-1.0774311) q[3];
sx q[3];
rz(-1.9947937) q[3];
sx q[3];
rz(3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0734633) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(1.8716795) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.3659182) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89643909) q[0];
sx q[0];
rz(-1.23587) q[0];
sx q[0];
rz(-2.4048637) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32760746) q[2];
sx q[2];
rz(-1.3407747) q[2];
sx q[2];
rz(2.1848781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3599638) q[1];
sx q[1];
rz(-2.8287323) q[1];
sx q[1];
rz(0.60731213) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4468206) q[3];
sx q[3];
rz(-2.0967391) q[3];
sx q[3];
rz(-0.40243173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0697249) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-2.8198077) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8480393) q[0];
sx q[0];
rz(-0.32395054) q[0];
sx q[0];
rz(2.241961) q[0];
rz(-pi) q[1];
rz(3.1084193) q[2];
sx q[2];
rz(-1.355194) q[2];
sx q[2];
rz(0.51378358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7223282) q[1];
sx q[1];
rz(-1.1280578) q[1];
sx q[1];
rz(-2.7302713) q[1];
rz(-2.3922032) q[3];
sx q[3];
rz(-1.6036311) q[3];
sx q[3];
rz(0.69946214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(0.11403306) q[2];
rz(-0.36241254) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645901) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.483451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720855) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(1.6540065) q[0];
rz(-pi) q[1];
rz(2.9786417) q[2];
sx q[2];
rz(-2.5451247) q[2];
sx q[2];
rz(0.60566723) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2918313) q[1];
sx q[1];
rz(-1.2408248) q[1];
sx q[1];
rz(-0.15177365) q[1];
rz(1.9547192) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(2.8843055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(1.8971987) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(-1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.3409021) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46643022) q[0];
sx q[0];
rz(-1.8008045) q[0];
sx q[0];
rz(-1.1943597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9843763) q[2];
sx q[2];
rz(-2.4210857) q[2];
sx q[2];
rz(0.23888982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(2.1800969) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5503928) q[3];
sx q[3];
rz(-1.692018) q[3];
sx q[3];
rz(-0.61074475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(-1.6188999) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(-0.32342708) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(1.1644438) q[2];
sx q[2];
rz(-1.5487557) q[2];
sx q[2];
rz(1.1974481) q[2];
rz(2.042751) q[3];
sx q[3];
rz(-1.2225716) q[3];
sx q[3];
rz(-0.93248758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
