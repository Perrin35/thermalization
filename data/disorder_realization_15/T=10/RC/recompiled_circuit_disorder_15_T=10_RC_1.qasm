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
rz(2.8957638) q[1];
sx q[1];
rz(10.413269) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52553015) q[0];
sx q[0];
rz(-1.5126192) q[0];
sx q[0];
rz(1.149403) q[0];
rz(-pi) q[1];
rz(0.97889401) q[2];
sx q[2];
rz(-1.4375293) q[2];
sx q[2];
rz(-2.8507581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6355675) q[1];
sx q[1];
rz(-2.4534181) q[1];
sx q[1];
rz(-2.3597005) q[1];
rz(-1.8139008) q[3];
sx q[3];
rz(-1.5339359) q[3];
sx q[3];
rz(-0.36323162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44895479) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(2.1733213) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(-1.2526858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7583) q[0];
sx q[0];
rz(-2.9581666) q[0];
sx q[0];
rz(-0.51390506) q[0];
rz(0.57931309) q[2];
sx q[2];
rz(-0.99537731) q[2];
sx q[2];
rz(-2.343611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15236552) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(2.3163124) q[1];
rz(-pi) q[2];
rz(-1.1329123) q[3];
sx q[3];
rz(-1.4012194) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0210555) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(-0.3343285) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702328) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(-3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-2.4287756) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26582742) q[0];
sx q[0];
rz(-1.7808502) q[0];
sx q[0];
rz(-1.34583) q[0];
x q[1];
rz(1.1926646) q[2];
sx q[2];
rz(-1.4107804) q[2];
sx q[2];
rz(0.68607054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94718139) q[1];
sx q[1];
rz(-0.57465034) q[1];
sx q[1];
rz(1.9504207) q[1];
x q[2];
rz(0.080428877) q[3];
sx q[3];
rz(-1.0897398) q[3];
sx q[3];
rz(1.2642494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-2.4776069) q[2];
rz(-1.9021696) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(-0.91528875) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5666714) q[0];
sx q[0];
rz(-3.1118588) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(0.92837292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61558047) q[0];
sx q[0];
rz(-0.86475879) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(-1.2296154) q[2];
sx q[2];
rz(-0.8114292) q[2];
sx q[2];
rz(-0.84412837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1410364) q[1];
sx q[1];
rz(-0.49266854) q[1];
sx q[1];
rz(-1.0286742) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9527693) q[3];
sx q[3];
rz(-0.7587772) q[3];
sx q[3];
rz(-3.0892059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0488247) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(1.1553923) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(-0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-3.0673448) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.689752) q[0];
sx q[0];
rz(-1.4983488) q[0];
sx q[0];
rz(-2.3206582) q[0];
x q[1];
rz(2.1553467) q[2];
sx q[2];
rz(-2.2918321) q[2];
sx q[2];
rz(-2.4436488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0641099) q[1];
sx q[1];
rz(-0.9573084) q[1];
sx q[1];
rz(0.20258278) q[1];
rz(-0.40162556) q[3];
sx q[3];
rz(-2.1021608) q[3];
sx q[3];
rz(1.348996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(1.1523694) q[2];
rz(-2.0955829) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(0.48962012) q[0];
rz(2.1173677) q[1];
sx q[1];
rz(-0.63413292) q[1];
sx q[1];
rz(1.5354059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438457) q[0];
sx q[0];
rz(-1.6698208) q[0];
sx q[0];
rz(-1.5227585) q[0];
x q[1];
rz(2.3697482) q[2];
sx q[2];
rz(-0.8470042) q[2];
sx q[2];
rz(-2.6592902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3052169) q[1];
sx q[1];
rz(-1.7292542) q[1];
sx q[1];
rz(-1.2837019) q[1];
x q[2];
rz(-1.7359211) q[3];
sx q[3];
rz(-2.193616) q[3];
sx q[3];
rz(1.3379054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9635222) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734633) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.3659182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.022) q[0];
sx q[0];
rz(-0.79607841) q[0];
sx q[0];
rz(0.47795602) q[0];
rz(-pi) q[1];
rz(-1.3283417) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(-0.53675011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9907896) q[1];
sx q[1];
rz(-1.3152796) q[1];
sx q[1];
rz(1.7533416) q[1];
x q[2];
rz(-2.9317022) q[3];
sx q[3];
rz(-2.6025747) q[3];
sx q[3];
rz(-2.982466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0697249) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(-2.7974421) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(0.32178497) q[0];
rz(0.92542648) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(0.61703533) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962284) q[0];
sx q[0];
rz(-1.8227238) q[0];
sx q[0];
rz(0.20586254) q[0];
x q[1];
rz(1.7210984) q[2];
sx q[2];
rz(-2.9234924) q[2];
sx q[2];
rz(-0.66767603) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2136692) q[1];
sx q[1];
rz(-2.546715) q[1];
sx q[1];
rz(2.2713714) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0934107) q[3];
sx q[3];
rz(-2.391624) q[3];
sx q[3];
rz(2.2349906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(3.0275596) q[2];
rz(-0.36241254) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(1.9053649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(-1.7846918) q[0];
rz(0.82018954) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(1.6581416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1781007) q[0];
sx q[0];
rz(-0.083847001) q[0];
sx q[0];
rz(1.6944147) q[0];
rz(1.4610897) q[2];
sx q[2];
rz(-2.1582971) q[2];
sx q[2];
rz(2.7320931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4087275) q[1];
sx q[1];
rz(-0.36204007) q[1];
sx q[1];
rz(-1.9865127) q[1];
rz(-pi) q[2];
rz(1.1868735) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(0.25728713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1411529) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.2443939) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72220951) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(-2.7710932) q[0];
rz(-1.1100948) q[1];
sx q[1];
rz(-2.4659174) q[1];
sx q[1];
rz(1.8006905) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46643022) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(1.947233) q[0];
x q[1];
rz(-2.9843763) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(0.23888982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1827613) q[1];
sx q[1];
rz(-1.3211831) q[1];
sx q[1];
rz(1.9447437) q[1];
x q[2];
rz(1.5503928) q[3];
sx q[3];
rz(-1.692018) q[3];
sx q[3];
rz(-0.61074475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5122539) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(-1.5226927) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387977) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(-3.1175989) q[2];
sx q[2];
rz(-1.1645483) q[2];
sx q[2];
rz(2.7777274) q[2];
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
