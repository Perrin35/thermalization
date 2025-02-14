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
rz(-0.40060488) q[0];
sx q[0];
rz(-2.7317943) q[0];
sx q[0];
rz(2.0706489) q[0];
rz(-2.5228956) q[1];
sx q[1];
rz(-0.56801152) q[1];
sx q[1];
rz(-2.2979589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6370359) q[0];
sx q[0];
rz(-1.3063653) q[0];
sx q[0];
rz(-3.1110071) q[0];
rz(-pi) q[1];
rz(0.29542342) q[2];
sx q[2];
rz(-1.7082126) q[2];
sx q[2];
rz(-0.58770056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93701332) q[1];
sx q[1];
rz(-2.3790303) q[1];
sx q[1];
rz(-2.049905) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0759962) q[3];
sx q[3];
rz(-1.9952979) q[3];
sx q[3];
rz(-0.44874661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7817276) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(3.1245933) q[2];
rz(1.7403691) q[3];
sx q[3];
rz(-2.5798116) q[3];
sx q[3];
rz(1.2167654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.568999) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(-2.3542985) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.2063824) q[1];
sx q[1];
rz(-1.9040727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71522349) q[0];
sx q[0];
rz(-2.5551717) q[0];
sx q[0];
rz(-3.0297902) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95713116) q[2];
sx q[2];
rz(-1.8120013) q[2];
sx q[2];
rz(1.7231736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54531714) q[1];
sx q[1];
rz(-1.8870237) q[1];
sx q[1];
rz(-1.113722) q[1];
rz(-pi) q[2];
rz(-2.5982598) q[3];
sx q[3];
rz(-0.50341922) q[3];
sx q[3];
rz(1.0375298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42147192) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(-3.0778911) q[2];
rz(1.8675768) q[3];
sx q[3];
rz(-0.84309045) q[3];
sx q[3];
rz(1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0434697) q[0];
sx q[0];
rz(-1.1935357) q[0];
sx q[0];
rz(-2.3980339) q[0];
rz(-2.6877563) q[1];
sx q[1];
rz(-1.3498787) q[1];
sx q[1];
rz(-0.41890621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48506698) q[0];
sx q[0];
rz(-2.65184) q[0];
sx q[0];
rz(1.9974677) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9100283) q[2];
sx q[2];
rz(-0.76398173) q[2];
sx q[2];
rz(0.52896777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6082904) q[1];
sx q[1];
rz(-0.69606298) q[1];
sx q[1];
rz(-2.0230369) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7396853) q[3];
sx q[3];
rz(-2.4554376) q[3];
sx q[3];
rz(-0.14631937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0539187) q[2];
sx q[2];
rz(-0.81834617) q[2];
sx q[2];
rz(1.9107001) q[2];
rz(-2.6026717) q[3];
sx q[3];
rz(-1.6659707) q[3];
sx q[3];
rz(-1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8777799) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(0.60212773) q[0];
rz(1.4471794) q[1];
sx q[1];
rz(-0.68232957) q[1];
sx q[1];
rz(-2.2785861) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569538) q[0];
sx q[0];
rz(-2.2881329) q[0];
sx q[0];
rz(1.7216645) q[0];
rz(2.237488) q[2];
sx q[2];
rz(-2.6764884) q[2];
sx q[2];
rz(2.1423774) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3366886) q[1];
sx q[1];
rz(-0.64283481) q[1];
sx q[1];
rz(1.1240608) q[1];
rz(-2.8252271) q[3];
sx q[3];
rz(-2.2834407) q[3];
sx q[3];
rz(-1.0567088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0756695) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(-0.90265957) q[2];
rz(2.3265808) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0937423) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(2.2064741) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-2.8159499) q[1];
sx q[1];
rz(-2.8757222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96052051) q[0];
sx q[0];
rz(-1.4645394) q[0];
sx q[0];
rz(-1.8828859) q[0];
x q[1];
rz(-2.8433617) q[2];
sx q[2];
rz(-2.3687393) q[2];
sx q[2];
rz(1.7333584) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8740377) q[1];
sx q[1];
rz(-1.4103762) q[1];
sx q[1];
rz(-0.35355132) q[1];
rz(-pi) q[2];
rz(1.83559) q[3];
sx q[3];
rz(-1.4715428) q[3];
sx q[3];
rz(0.81391993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68044668) q[2];
sx q[2];
rz(-1.6411883) q[2];
sx q[2];
rz(0.45073304) q[2];
rz(-1.9510673) q[3];
sx q[3];
rz(-2.4400986) q[3];
sx q[3];
rz(2.7019971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180069) q[0];
sx q[0];
rz(-0.30338767) q[0];
sx q[0];
rz(0.044483749) q[0];
rz(-2.8728409) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(-0.43089795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6857032) q[0];
sx q[0];
rz(-0.95818943) q[0];
sx q[0];
rz(-1.8926748) q[0];
x q[1];
rz(2.2399432) q[2];
sx q[2];
rz(-1.2001654) q[2];
sx q[2];
rz(1.8805875) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.84316501) q[1];
sx q[1];
rz(-0.80363217) q[1];
sx q[1];
rz(0.9982361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.481302) q[3];
sx q[3];
rz(-1.5753645) q[3];
sx q[3];
rz(0.82483236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1162794) q[2];
sx q[2];
rz(-0.37386027) q[2];
sx q[2];
rz(2.5171793) q[2];
rz(2.0404909) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-0.98744923) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91200149) q[0];
sx q[0];
rz(-2.0822552) q[0];
sx q[0];
rz(-2.4578995) q[0];
rz(-1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(-2.9396465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366906) q[0];
sx q[0];
rz(-0.32819191) q[0];
sx q[0];
rz(1.5668014) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5187289) q[2];
sx q[2];
rz(-0.74801842) q[2];
sx q[2];
rz(-1.6287273) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4836639) q[1];
sx q[1];
rz(-0.52573181) q[1];
sx q[1];
rz(-0.13327285) q[1];
x q[2];
rz(1.6138329) q[3];
sx q[3];
rz(-1.1133872) q[3];
sx q[3];
rz(2.4105173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0553637) q[2];
sx q[2];
rz(-1.4198885) q[2];
sx q[2];
rz(1.281338) q[2];
rz(-2.4751439) q[3];
sx q[3];
rz(-1.0969578) q[3];
sx q[3];
rz(-2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(2.0228588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7549158) q[0];
sx q[0];
rz(-0.43114063) q[0];
sx q[0];
rz(0.19622959) q[0];
x q[1];
rz(2.7114026) q[2];
sx q[2];
rz(-2.3025844) q[2];
sx q[2];
rz(-0.48626394) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9671849) q[1];
sx q[1];
rz(-2.6749415) q[1];
sx q[1];
rz(1.0789167) q[1];
rz(-pi) q[2];
rz(0.54688485) q[3];
sx q[3];
rz(-0.62084475) q[3];
sx q[3];
rz(-1.8016694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.755456) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(0.90292162) q[2];
rz(-1.6821945) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(0.82227796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4626386) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(1.1902887) q[0];
rz(-2.8207488) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(-1.920059) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4325334) q[0];
sx q[0];
rz(-1.4343059) q[0];
sx q[0];
rz(0.20705072) q[0];
x q[1];
rz(-2.8655445) q[2];
sx q[2];
rz(-0.85535565) q[2];
sx q[2];
rz(0.1705585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1561191) q[1];
sx q[1];
rz(-1.7745943) q[1];
sx q[1];
rz(0.76038313) q[1];
x q[2];
rz(1.3876347) q[3];
sx q[3];
rz(-1.3792999) q[3];
sx q[3];
rz(-0.93246704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(0.3332738) q[2];
rz(-2.1244369) q[3];
sx q[3];
rz(-0.95484304) q[3];
sx q[3];
rz(-0.3199544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79621133) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(-2.2228125) q[0];
rz(-0.19110075) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(2.3474615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988172) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(0.88943435) q[0];
rz(1.3599858) q[2];
sx q[2];
rz(-1.9242058) q[2];
sx q[2];
rz(1.4251874) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0550784) q[1];
sx q[1];
rz(-0.76755133) q[1];
sx q[1];
rz(-1.961019) q[1];
rz(2.8201032) q[3];
sx q[3];
rz(-1.7712542) q[3];
sx q[3];
rz(3.0908302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11067757) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(0.87016726) q[3];
sx q[3];
rz(-0.78799677) q[3];
sx q[3];
rz(-1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(-0.43791804) q[1];
sx q[1];
rz(-1.7057849) q[1];
sx q[1];
rz(1.5954856) q[1];
rz(0.6952473) q[2];
sx q[2];
rz(-2.5378479) q[2];
sx q[2];
rz(-2.9369864) q[2];
rz(-2.6043456) q[3];
sx q[3];
rz(-0.8575079) q[3];
sx q[3];
rz(-1.433123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
