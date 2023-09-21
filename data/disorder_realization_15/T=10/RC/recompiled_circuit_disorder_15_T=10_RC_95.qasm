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
rz(-1.7742046) q[1];
sx q[1];
rz(-2.8957638) q[1];
sx q[1];
rz(-2.153102) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6160625) q[0];
sx q[0];
rz(-1.5126192) q[0];
sx q[0];
rz(-1.149403) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9814331) q[2];
sx q[2];
rz(-0.9848435) q[2];
sx q[2];
rz(1.1908659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5520652) q[1];
sx q[1];
rz(-2.0347934) q[1];
sx q[1];
rz(-2.6134171) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3276919) q[3];
sx q[3];
rz(-1.6076568) q[3];
sx q[3];
rz(-0.36323162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73873591) q[0];
sx q[0];
rz(-2.0500654) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(1.8889069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90447146) q[0];
sx q[0];
rz(-1.4112817) q[0];
sx q[0];
rz(1.479854) q[0];
rz(-pi) q[1];
rz(2.5622796) q[2];
sx q[2];
rz(-2.1462153) q[2];
sx q[2];
rz(0.79798165) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0536086) q[1];
sx q[1];
rz(-2.0961742) q[1];
sx q[1];
rz(-2.5768075) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9545752) q[3];
sx q[3];
rz(-2.6740101) q[3];
sx q[3];
rz(-3.0083857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0210555) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(2.7978314) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-0.25046644) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-2.4287756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3526488) q[0];
sx q[0];
rz(-1.3508571) q[0];
sx q[0];
rz(-0.21531944) q[0];
rz(-pi) q[1];
rz(-1.9829282) q[2];
sx q[2];
rz(-2.7325028) q[2];
sx q[2];
rz(-1.8754174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1944113) q[1];
sx q[1];
rz(-0.57465034) q[1];
sx q[1];
rz(-1.191172) q[1];
x q[2];
rz(2.0531822) q[3];
sx q[3];
rz(-1.4995121) q[3];
sx q[3];
rz(-2.7977668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(-2.4776069) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5749213) q[0];
sx q[0];
rz(-3.1118588) q[0];
sx q[0];
rz(2.5675039) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-2.2132197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1855477) q[0];
sx q[0];
rz(-1.5717686) q[0];
sx q[0];
rz(0.86475839) q[0];
rz(-pi) q[1];
rz(-2.3525535) q[2];
sx q[2];
rz(-1.3256729) q[2];
sx q[2];
rz(-0.48691985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91765431) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(2.0018105) q[1];
rz(-0.5023605) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(0.82752284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0488247) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(1.9862004) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(-0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578385) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(3.0673448) q[0];
rz(1.6429365) q[1];
sx q[1];
rz(-1.628412) q[1];
sx q[1];
rz(0.9517076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041391011) q[0];
sx q[0];
rz(-0.75267422) q[0];
sx q[0];
rz(1.676883) q[0];
rz(0.98624595) q[2];
sx q[2];
rz(-2.2918321) q[2];
sx q[2];
rz(-0.69794387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73478991) q[1];
sx q[1];
rz(-0.64195913) q[1];
sx q[1];
rz(-1.8491247) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0025025) q[3];
sx q[3];
rz(-1.9145402) q[3];
sx q[3];
rz(-2.7078201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9933269) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(1.024225) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(-1.6061868) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55035704) q[0];
sx q[0];
rz(-3.0315657) q[0];
sx q[0];
rz(-0.45022924) q[0];
rz(-pi) q[1];
rz(2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(2.6256109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3121351) q[1];
sx q[1];
rz(-1.8541938) q[1];
sx q[1];
rz(0.16510041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62932265) q[3];
sx q[3];
rz(-1.4368847) q[3];
sx q[3];
rz(-2.8117992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9635222) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(-0.17573389) q[2];
rz(2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734633) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.8716795) q[0];
rz(2.9636256) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(1.3659182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1195927) q[0];
sx q[0];
rz(-2.3455142) q[0];
sx q[0];
rz(-0.47795602) q[0];
x q[1];
rz(0.32760746) q[2];
sx q[2];
rz(-1.8008179) q[2];
sx q[2];
rz(0.95671456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7682225) q[1];
sx q[1];
rz(-1.3942413) q[1];
sx q[1];
rz(-0.25964398) q[1];
rz(2.6122983) q[3];
sx q[3];
rz(-1.6779473) q[3];
sx q[3];
rz(-1.9107493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(2.243637) q[2];
rz(2.9979624) q[3];
sx q[3];
rz(-1.8840021) q[3];
sx q[3];
rz(0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550069) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(-2.8198077) q[0];
rz(2.2161662) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(-2.5245573) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5453642) q[0];
sx q[0];
rz(-1.8227238) q[0];
sx q[0];
rz(-2.9357301) q[0];
rz(-1.7865137) q[2];
sx q[2];
rz(-1.5383913) q[2];
sx q[2];
rz(-2.0916794) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.033213381) q[1];
sx q[1];
rz(-1.9404267) q[1];
sx q[1];
rz(1.0934248) q[1];
x q[2];
rz(-0.048181941) q[3];
sx q[3];
rz(-0.74996862) q[3];
sx q[3];
rz(-0.90660209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5499605) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(-3.0275596) q[2];
rz(-0.36241254) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(1.9053649) q[3];
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
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(-1.3569008) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(1.483451) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9634919) q[0];
sx q[0];
rz(-3.0577457) q[0];
sx q[0];
rz(-1.6944147) q[0];
rz(0.59028583) q[2];
sx q[2];
rz(-1.4795408) q[2];
sx q[2];
rz(-1.1003189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7328651) q[1];
sx q[1];
rz(-2.7795526) q[1];
sx q[1];
rz(-1.15508) q[1];
rz(-0.58247329) q[3];
sx q[3];
rz(-2.5081722) q[3];
sx q[3];
rz(-2.7137091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(1.5104793) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4193831) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(0.37049946) q[0];
rz(2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.8006905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6751624) q[0];
sx q[0];
rz(-1.3407882) q[0];
sx q[0];
rz(-1.947233) q[0];
rz(-pi) q[1];
rz(-2.4272333) q[2];
sx q[2];
rz(-1.4673125) q[2];
sx q[2];
rz(-1.4504745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(-0.96149573) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.020346) q[3];
sx q[3];
rz(-1.5910501) q[3];
sx q[3];
rz(-2.1840087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5122539) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(-1.5226927) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387977) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(0.32342708) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(-0.023993775) q[2];
sx q[2];
rz(-1.9770443) q[2];
sx q[2];
rz(-0.36386522) q[2];
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