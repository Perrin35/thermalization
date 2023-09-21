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
rz(-2.0702695) q[0];
sx q[0];
rz(-1.9914314) q[0];
sx q[0];
rz(3.0778528) q[0];
x q[1];
rz(1.8066054) q[2];
sx q[2];
rz(-2.5366282) q[2];
sx q[2];
rz(-1.4750823) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72585427) q[1];
sx q[1];
rz(-2.038318) q[1];
sx q[1];
rz(-2.095925) q[1];
rz(-pi) q[2];
rz(-1.7228116) q[3];
sx q[3];
rz(-2.8957643) q[3];
sx q[3];
rz(2.0815401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44895479) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(-2.5732178) q[2];
rz(1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4028567) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(2.7785595) q[0];
rz(-0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.2526858) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90447146) q[0];
sx q[0];
rz(-1.7303109) q[0];
sx q[0];
rz(1.6617387) q[0];
x q[1];
rz(-2.230174) q[2];
sx q[2];
rz(-1.0936001) q[2];
sx q[2];
rz(2.0267682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9892271) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(-2.3163124) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0086803) q[3];
sx q[3];
rz(-1.4012194) q[3];
sx q[3];
rz(-1.7835497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-2.7978314) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(0.0070455889) q[0];
rz(2.7650611) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-0.71281707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8757652) q[0];
sx q[0];
rz(-1.3607425) q[0];
sx q[0];
rz(-1.34583) q[0];
rz(-pi) q[1];
rz(-1.9829282) q[2];
sx q[2];
rz(-0.40908989) q[2];
sx q[2];
rz(-1.2661753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94718139) q[1];
sx q[1];
rz(-0.57465034) q[1];
sx q[1];
rz(1.9504207) q[1];
rz(-pi) q[2];
rz(-2.0531822) q[3];
sx q[3];
rz(-1.4995121) q[3];
sx q[3];
rz(-0.34382581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(2.4776069) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(2.2263039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5749213) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-0.57408875) q[0];
rz(2.8494049) q[1];
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
rz(2.5260122) q[0];
sx q[0];
rz(-2.2768339) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(-pi) q[1];
rz(0.33889126) q[2];
sx q[2];
rz(-0.8182943) q[2];
sx q[2];
rz(-1.8213059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1410364) q[1];
sx q[1];
rz(-0.49266854) q[1];
sx q[1];
rz(1.0286742) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2286962) q[3];
sx q[3];
rz(-1.9808589) q[3];
sx q[3];
rz(-1.0421154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.092768) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(1.9862004) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(-2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578385) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-0.074247867) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18626801) q[0];
sx q[0];
rz(-0.82337472) q[0];
sx q[0];
rz(-0.098851725) q[0];
rz(0.81163002) q[2];
sx q[2];
rz(-1.1434165) q[2];
sx q[2];
rz(-0.4610093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0774827) q[1];
sx q[1];
rz(-0.9573084) q[1];
sx q[1];
rz(2.9390099) q[1];
x q[2];
rz(2.7399671) q[3];
sx q[3];
rz(-2.1021608) q[3];
sx q[3];
rz(1.348996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(1.9892233) q[2];
rz(1.0460098) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(-0.48962012) q[0];
rz(-1.024225) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.6061868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55035704) q[0];
sx q[0];
rz(-0.11002692) q[0];
sx q[0];
rz(0.45022924) q[0];
rz(-2.3697482) q[2];
sx q[2];
rz(-2.2945885) q[2];
sx q[2];
rz(0.48230241) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3121351) q[1];
sx q[1];
rz(-1.8541938) q[1];
sx q[1];
rz(2.9764922) q[1];
x q[2];
rz(1.7359211) q[3];
sx q[3];
rz(-0.94797665) q[3];
sx q[3];
rz(-1.8036872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(0.17573389) q[2];
rz(-2.0641616) q[3];
sx q[3];
rz(-1.9947937) q[3];
sx q[3];
rz(0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(1.3659182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458347) q[0];
sx q[0];
rz(-0.88338822) q[0];
sx q[0];
rz(2.0100726) q[0];
rz(1.3283417) q[2];
sx q[2];
rz(-1.889466) q[2];
sx q[2];
rz(2.6048425) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7682225) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(0.25964398) q[1];
x q[2];
rz(-1.4468206) q[3];
sx q[3];
rz(-1.0448536) q[3];
sx q[3];
rz(0.40243173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0697249) q[2];
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
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550069) q[0];
sx q[0];
rz(-2.2877559) q[0];
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
rz(-2.8176421) q[0];
sx q[0];
rz(0.89963161) q[0];
x q[1];
rz(0.033173325) q[2];
sx q[2];
rz(-1.355194) q[2];
sx q[2];
rz(-0.51378358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2136692) q[1];
sx q[1];
rz(-0.59487768) q[1];
sx q[1];
rz(-2.2713714) q[1];
x q[2];
rz(1.5259605) q[3];
sx q[3];
rz(-0.82190824) q[3];
sx q[3];
rz(0.84079784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5499605) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(0.11403306) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(-1.9053649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(1.3569008) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.483451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26950715) q[0];
sx q[0];
rz(-1.5811231) q[0];
sx q[0];
rz(-1.6540065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9786417) q[2];
sx q[2];
rz(-2.5451247) q[2];
sx q[2];
rz(-0.60566723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8497613) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(-2.989819) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1868735) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(-0.25728713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-0.65594643) q[2];
sx q[2];
rz(-1.8971987) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(-0.37049946) q[0];
rz(2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(-1.3409021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46643022) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(1.1943597) q[0];
rz(-pi) q[1];
rz(0.15721639) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(-2.9027028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1738051) q[1];
sx q[1];
rz(-2.6952744) q[1];
sx q[1];
rz(2.1800969) q[1];
rz(-pi) q[2];
rz(-1.5911999) q[3];
sx q[3];
rz(-1.4495747) q[3];
sx q[3];
rz(0.61074475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(-1.6188999) q[2];
rz(0.86087888) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-3.1057152) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027949) q[0];
sx q[0];
rz(-1.4807985) q[0];
sx q[0];
rz(2.280622) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(-1.515083) q[2];
sx q[2];
rz(-2.734676) q[2];
sx q[2];
rz(2.717072) q[2];
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
