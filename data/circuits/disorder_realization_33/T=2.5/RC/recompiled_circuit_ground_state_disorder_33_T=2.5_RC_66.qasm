OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7140952) q[0];
sx q[0];
rz(3.7035898) q[0];
sx q[0];
rz(9.193767) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481349) q[0];
sx q[0];
rz(-2.0525816) q[0];
sx q[0];
rz(-2.8346377) q[0];
rz(0.76164772) q[2];
sx q[2];
rz(-2.8084694) q[2];
sx q[2];
rz(1.8322577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2778863) q[1];
sx q[1];
rz(-1.4046298) q[1];
sx q[1];
rz(-1.6707129) q[1];
rz(-pi) q[2];
rz(3.0662905) q[3];
sx q[3];
rz(-1.3960394) q[3];
sx q[3];
rz(-0.14698262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(-1.1245842) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(-1.2688961) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26038134) q[0];
sx q[0];
rz(-2.9660048) q[0];
sx q[0];
rz(-2.4847109) q[0];
rz(-0.33292133) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(-0.93516707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774076) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(0.10459374) q[0];
rz(-0.33719535) q[2];
sx q[2];
rz(-2.8697578) q[2];
sx q[2];
rz(1.4120917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.024927372) q[1];
sx q[1];
rz(-1.6102734) q[1];
sx q[1];
rz(0.14703207) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35547361) q[3];
sx q[3];
rz(-1.7661816) q[3];
sx q[3];
rz(3.1110142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3071345) q[2];
sx q[2];
rz(-1.5156526) q[2];
sx q[2];
rz(-1.2163986) q[2];
rz(2.6136716) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.886604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6127748) q[0];
sx q[0];
rz(-0.82614326) q[0];
sx q[0];
rz(2.5566027) q[0];
rz(-3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(1.4488719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9609279) q[0];
sx q[0];
rz(-1.573473) q[0];
sx q[0];
rz(-1.5737783) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18529202) q[2];
sx q[2];
rz(-3.0681562) q[2];
sx q[2];
rz(-0.18723182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6022012) q[1];
sx q[1];
rz(-2.3219206) q[1];
sx q[1];
rz(-2.17893) q[1];
x q[2];
rz(2.4873729) q[3];
sx q[3];
rz(-0.4163792) q[3];
sx q[3];
rz(2.6423796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2584194) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(-1.4460571) q[2];
rz(-2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676753) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(2.0539334) q[0];
rz(-2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(2.1813724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4488569) q[0];
sx q[0];
rz(-1.6105284) q[0];
sx q[0];
rz(1.3055152) q[0];
x q[1];
rz(-1.4713418) q[2];
sx q[2];
rz(-1.5410454) q[2];
sx q[2];
rz(-0.1386252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31611495) q[1];
sx q[1];
rz(-1.1255472) q[1];
sx q[1];
rz(-2.6632505) q[1];
rz(2.8720565) q[3];
sx q[3];
rz(-1.663066) q[3];
sx q[3];
rz(0.59449457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.20597657) q[2];
sx q[2];
rz(-1.874186) q[2];
sx q[2];
rz(-1.4975632) q[2];
rz(-2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(-0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(0.60428756) q[0];
rz(-1.9970278) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(0.42246517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6910962) q[0];
sx q[0];
rz(-1.4022777) q[0];
sx q[0];
rz(-3.0590579) q[0];
x q[1];
rz(-2.1429135) q[2];
sx q[2];
rz(-2.0586176) q[2];
sx q[2];
rz(1.3582548) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8026655) q[1];
sx q[1];
rz(-0.69310729) q[1];
sx q[1];
rz(0.1130123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38099184) q[3];
sx q[3];
rz(-1.2514827) q[3];
sx q[3];
rz(-1.3407941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8233238) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(2.4853415) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(-2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4055279) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(-2.3896353) q[1];
sx q[1];
rz(-1.0686921) q[1];
sx q[1];
rz(-1.6065074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4055734) q[0];
sx q[0];
rz(-1.6200119) q[0];
sx q[0];
rz(0.28926138) q[0];
rz(-2.4322832) q[2];
sx q[2];
rz(-1.0934208) q[2];
sx q[2];
rz(0.8698744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62404406) q[1];
sx q[1];
rz(-0.66429783) q[1];
sx q[1];
rz(0.37599998) q[1];
rz(0.84523369) q[3];
sx q[3];
rz(-0.74016011) q[3];
sx q[3];
rz(-2.6773767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9525166) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(2.2124186) q[2];
rz(-0.82516986) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(-2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.3744542) q[0];
rz(2.6311686) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(0.62320954) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685112) q[0];
sx q[0];
rz(-0.36751908) q[0];
sx q[0];
rz(2.0042242) q[0];
rz(-pi) q[1];
rz(1.7665777) q[2];
sx q[2];
rz(-2.0718241) q[2];
sx q[2];
rz(0.50925469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0556866) q[1];
sx q[1];
rz(-1.4594363) q[1];
sx q[1];
rz(1.423905) q[1];
x q[2];
rz(0.74759746) q[3];
sx q[3];
rz(-1.4211402) q[3];
sx q[3];
rz(-0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(0.022631571) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(-2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8959494) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(1.0173215) q[0];
rz(-1.3941049) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(2.5140433) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857972) q[0];
sx q[0];
rz(-1.9667407) q[0];
sx q[0];
rz(2.2084479) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4862138) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(1.6135474) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1422954) q[1];
sx q[1];
rz(-2.8311756) q[1];
sx q[1];
rz(2.4503972) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9110319) q[3];
sx q[3];
rz(-1.1200532) q[3];
sx q[3];
rz(2.4204202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5158186) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(0.41075692) q[2];
rz(-0.050203236) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(0.075411782) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.582616) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(-0.26891747) q[0];
rz(-0.31002632) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(-0.1300098) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01154218) q[0];
sx q[0];
rz(-0.39549144) q[0];
sx q[0];
rz(0.96590913) q[0];
rz(0.35138826) q[2];
sx q[2];
rz(-1.2341502) q[2];
sx q[2];
rz(1.2996246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7183154) q[1];
sx q[1];
rz(-2.1287103) q[1];
sx q[1];
rz(-1.7608791) q[1];
x q[2];
rz(-0.9762398) q[3];
sx q[3];
rz(-2.3701982) q[3];
sx q[3];
rz(0.4679799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93398634) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-3.0450191) q[2];
rz(-1.5038331) q[3];
sx q[3];
rz(-1.3373172) q[3];
sx q[3];
rz(-0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-2.9533563) q[0];
sx q[0];
rz(2.6085594) q[0];
rz(3.10532) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(-1.689555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17981681) q[0];
sx q[0];
rz(-0.70588934) q[0];
sx q[0];
rz(-2.174586) q[0];
rz(-0.31923652) q[2];
sx q[2];
rz(-1.6743273) q[2];
sx q[2];
rz(1.5872623) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.328985) q[1];
sx q[1];
rz(-1.9300864) q[1];
sx q[1];
rz(1.6027662) q[1];
rz(1.6598654) q[3];
sx q[3];
rz(-1.1067821) q[3];
sx q[3];
rz(2.9651412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4286246) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(0.0326322) q[2];
rz(-2.2964358) q[3];
sx q[3];
rz(-1.2266351) q[3];
sx q[3];
rz(-2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(-2.3445917) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-1.8664411) q[2];
sx q[2];
rz(-1.0545803) q[2];
sx q[2];
rz(2.7068356) q[2];
rz(-3.0211021) q[3];
sx q[3];
rz(-2.3269666) q[3];
sx q[3];
rz(2.5757488) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
