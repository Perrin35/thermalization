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
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0713232) q[0];
sx q[0];
rz(-1.1501612) q[0];
sx q[0];
rz(3.0778528) q[0];
rz(-pi) q[1];
rz(0.97889401) q[2];
sx q[2];
rz(-1.7040633) q[2];
sx q[2];
rz(2.8507581) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4157384) q[1];
sx q[1];
rz(-2.038318) q[1];
sx q[1];
rz(-2.095925) q[1];
x q[2];
rz(-1.418781) q[3];
sx q[3];
rz(-2.8957643) q[3];
sx q[3];
rz(-2.0815401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6926379) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(-0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4028567) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(0.36303315) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.8889069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607836) q[0];
sx q[0];
rz(-1.660581) q[0];
sx q[0];
rz(-0.16016527) q[0];
x q[1];
rz(0.57931309) q[2];
sx q[2];
rz(-2.1462153) q[2];
sx q[2];
rz(2.343611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9892271) q[1];
sx q[1];
rz(-2.3902635) q[1];
sx q[1];
rz(2.3163124) q[1];
x q[2];
rz(1.9545752) q[3];
sx q[3];
rz(-2.6740101) q[3];
sx q[3];
rz(-0.13320696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-2.7978314) q[2];
rz(-2.8072642) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702328) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(0.0070455889) q[0];
rz(2.7650611) q[1];
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
rz(-1.3508571) q[0];
sx q[0];
rz(2.9262732) q[0];
rz(-1.9489281) q[2];
sx q[2];
rz(-1.7308123) q[2];
sx q[2];
rz(-0.68607054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3908927) q[1];
sx q[1];
rz(-2.0999968) q[1];
sx q[1];
rz(-0.23551029) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0611638) q[3];
sx q[3];
rz(-1.0897398) q[3];
sx q[3];
rz(-1.8773432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3537447) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(2.4776069) q[2];
rz(1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(-2.8494049) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(0.92837292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61558047) q[0];
sx q[0];
rz(-2.2768339) q[0];
sx q[0];
rz(3.1403149) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3525535) q[2];
sx q[2];
rz(-1.3256729) q[2];
sx q[2];
rz(0.48691985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0005562) q[1];
sx q[1];
rz(-2.6489241) q[1];
sx q[1];
rz(2.1129184) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5023605) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(-2.3140698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0488247) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(0.14285764) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(0.75240451) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(3.0673448) q[0];
rz(-1.4986562) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518406) q[0];
sx q[0];
rz(-1.6432439) q[0];
sx q[0];
rz(2.3206582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3299626) q[2];
sx q[2];
rz(-1.9981761) q[2];
sx q[2];
rz(-0.4610093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61102066) q[1];
sx q[1];
rz(-1.7360577) q[1];
sx q[1];
rz(0.9475488) q[1];
rz(-2.7399671) q[3];
sx q[3];
rz(-2.1021608) q[3];
sx q[3];
rz(1.7925966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(1.9892233) q[2];
rz(-1.0460098) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(3.1402821) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(2.6519725) q[0];
rz(1.024225) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(-1.6061868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4682966) q[0];
sx q[0];
rz(-1.522994) q[0];
sx q[0];
rz(-0.099138069) q[0];
rz(-pi) q[1];
rz(-2.4602731) q[2];
sx q[2];
rz(-2.1207003) q[2];
sx q[2];
rz(2.6256109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.225243) q[1];
sx q[1];
rz(-0.32685977) q[1];
sx q[1];
rz(-1.0570231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4056715) q[3];
sx q[3];
rz(-2.193616) q[3];
sx q[3];
rz(1.3379054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(-2.0641616) q[3];
sx q[3];
rz(-1.9947937) q[3];
sx q[3];
rz(-3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068129383) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(1.2699132) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-2.3033419) q[1];
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
rz(-pi) q[1];
x q[1];
rz(0.32760746) q[2];
sx q[2];
rz(-1.8008179) q[2];
sx q[2];
rz(-2.1848781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7682225) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(0.25964398) q[1];
x q[2];
rz(1.6947721) q[3];
sx q[3];
rz(-1.0448536) q[3];
sx q[3];
rz(-2.7391609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(-2.9979624) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(-2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(2.8198077) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5962284) q[0];
sx q[0];
rz(-1.3188688) q[0];
sx q[0];
rz(-2.9357301) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4204942) q[2];
sx q[2];
rz(-2.9234924) q[2];
sx q[2];
rz(0.66767603) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2136692) q[1];
sx q[1];
rz(-0.59487768) q[1];
sx q[1];
rz(-0.87022123) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3922032) q[3];
sx q[3];
rz(-1.5379616) q[3];
sx q[3];
rz(-0.69946214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5499605) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(3.0275596) q[2];
rz(0.36241254) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.9053649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(-0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.483451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950715) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(-1.6540065) q[0];
x q[1];
rz(2.5513068) q[2];
sx q[2];
rz(-1.4795408) q[2];
sx q[2];
rz(1.1003189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2918313) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(0.15177365) q[1];
x q[2];
rz(1.9547192) q[3];
sx q[3];
rz(-1.0537638) q[3];
sx q[3];
rz(0.25728713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.8971987) q[2];
rz(-0.42516285) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.6311133) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72220951) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(-2.7710932) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(1.3409021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46643022) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(1.1943597) q[0];
rz(-pi) q[1];
rz(0.71435931) q[2];
sx q[2];
rz(-1.4673125) q[2];
sx q[2];
rz(1.6911182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1738051) q[1];
sx q[1];
rz(-2.6952744) q[1];
sx q[1];
rz(-2.1800969) q[1];
x q[2];
rz(0.16593905) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5122539) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(0.86087888) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2387977) q[0];
sx q[0];
rz(-1.4807985) q[0];
sx q[0];
rz(2.280622) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(-1.1644438) q[2];
sx q[2];
rz(-1.5928369) q[2];
sx q[2];
rz(-1.9441446) q[2];
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