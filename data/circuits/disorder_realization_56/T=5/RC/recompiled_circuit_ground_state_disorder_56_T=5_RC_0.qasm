OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(0.78328744) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(4.6286197) q[1];
sx q[1];
rz(11.774465) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5675674) q[0];
sx q[0];
rz(-1.3114531) q[0];
sx q[0];
rz(-1.3784842) q[0];
rz(-pi) q[1];
rz(1.804084) q[2];
sx q[2];
rz(-2.983494) q[2];
sx q[2];
rz(0.563941) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79719964) q[1];
sx q[1];
rz(-0.54423344) q[1];
sx q[1];
rz(-0.060677008) q[1];
rz(-0.57609419) q[3];
sx q[3];
rz(-0.8818501) q[3];
sx q[3];
rz(2.9722633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(0.099451065) q[2];
rz(0.24672306) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-0.39768404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1219015) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(-1.6909201) q[1];
sx q[1];
rz(-1.2932777) q[1];
sx q[1];
rz(-0.64847747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4546616) q[0];
sx q[0];
rz(-1.0333916) q[0];
sx q[0];
rz(-0.65673687) q[0];
rz(-pi) q[1];
rz(-1.8869867) q[2];
sx q[2];
rz(-2.11497) q[2];
sx q[2];
rz(-2.3426825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8295171) q[1];
sx q[1];
rz(-1.5863998) q[1];
sx q[1];
rz(-0.14194686) q[1];
rz(-3.0233211) q[3];
sx q[3];
rz(-1.857548) q[3];
sx q[3];
rz(-0.82055446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6917981) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(2.2890384) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6344675) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(1.0852098) q[0];
rz(-2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(-1.6403713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417898) q[0];
sx q[0];
rz(-2.9405624) q[0];
sx q[0];
rz(-0.19636671) q[0];
x q[1];
rz(1.3365715) q[2];
sx q[2];
rz(-1.7570436) q[2];
sx q[2];
rz(-3.0795003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13255626) q[1];
sx q[1];
rz(-1.471613) q[1];
sx q[1];
rz(2.1445091) q[1];
rz(-pi) q[2];
rz(1.3543868) q[3];
sx q[3];
rz(-1.4167656) q[3];
sx q[3];
rz(1.1104402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40470716) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(2.2331623) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(-2.0250208) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92723769) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(-1.9450564) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(1.1845142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9036983) q[0];
sx q[0];
rz(-1.333295) q[0];
sx q[0];
rz(-0.71126513) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3886613) q[2];
sx q[2];
rz(-1.7106461) q[2];
sx q[2];
rz(1.1876196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8912011) q[1];
sx q[1];
rz(-1.3445092) q[1];
sx q[1];
rz(-1.0550642) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5585737) q[3];
sx q[3];
rz(-0.32430092) q[3];
sx q[3];
rz(0.7079269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2148332) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(-0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(-3.1121837) q[0];
rz(0.75621653) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(-0.75278935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469893) q[0];
sx q[0];
rz(-1.5651817) q[0];
sx q[0];
rz(3.1412197) q[0];
rz(-pi) q[1];
rz(1.3026313) q[2];
sx q[2];
rz(-1.6577621) q[2];
sx q[2];
rz(3.027107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1323586) q[1];
sx q[1];
rz(-2.1387324) q[1];
sx q[1];
rz(-1.4465989) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54075586) q[3];
sx q[3];
rz(-1.1782559) q[3];
sx q[3];
rz(0.3730216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14043643) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(-2.746554) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(-0.37676677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-2.4212615) q[0];
sx q[0];
rz(2.1790867) q[0];
rz(-0.51482254) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(0.33822507) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047236) q[0];
sx q[0];
rz(-2.3677808) q[0];
sx q[0];
rz(1.85352) q[0];
rz(0.2372431) q[2];
sx q[2];
rz(-1.8117935) q[2];
sx q[2];
rz(-2.4324696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4938062) q[1];
sx q[1];
rz(-1.567651) q[1];
sx q[1];
rz(2.6814744e-05) q[1];
x q[2];
rz(-1.148726) q[3];
sx q[3];
rz(-1.5079807) q[3];
sx q[3];
rz(2.738225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6414791) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(2.2902299) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(2.5135777) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250658) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(0.75041962) q[0];
rz(-0.57506192) q[1];
sx q[1];
rz(-1.3651747) q[1];
sx q[1];
rz(-2.4837928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019264) q[0];
sx q[0];
rz(-2.3453379) q[0];
sx q[0];
rz(-2.5523561) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70722981) q[2];
sx q[2];
rz(-1.6042738) q[2];
sx q[2];
rz(2.1113124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9283674) q[1];
sx q[1];
rz(-1.8568653) q[1];
sx q[1];
rz(0.017315344) q[1];
rz(-pi) q[2];
rz(1.1893473) q[3];
sx q[3];
rz(-1.7717517) q[3];
sx q[3];
rz(-2.902346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49360069) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.3841265) q[3];
sx q[3];
rz(2.306126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(-1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.15377741) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7936121) q[0];
sx q[0];
rz(-1.6456933) q[0];
sx q[0];
rz(-2.8872284) q[0];
x q[1];
rz(2.2258899) q[2];
sx q[2];
rz(-2.2580552) q[2];
sx q[2];
rz(-2.5926431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1250026) q[1];
sx q[1];
rz(-1.224813) q[1];
sx q[1];
rz(-1.6660652) q[1];
x q[2];
rz(1.1438391) q[3];
sx q[3];
rz(-1.6817963) q[3];
sx q[3];
rz(-0.75253651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.095857233) q[2];
sx q[2];
rz(-2.6726674) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(2.5921953) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-2.8909491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9662387) q[0];
sx q[0];
rz(-2.4099775) q[0];
sx q[0];
rz(1.0989972) q[0];
rz(-0.05392404) q[2];
sx q[2];
rz(-2.8802875) q[2];
sx q[2];
rz(0.037029412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5329001) q[1];
sx q[1];
rz(-1.6402555) q[1];
sx q[1];
rz(0.77381247) q[1];
rz(-pi) q[2];
rz(2.9467877) q[3];
sx q[3];
rz(-1.4470295) q[3];
sx q[3];
rz(-0.74461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15829076) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(-0.41900751) q[0];
rz(-2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(-0.73582617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9402855) q[0];
sx q[0];
rz(-0.92865151) q[0];
sx q[0];
rz(-1.8117732) q[0];
rz(-pi) q[1];
rz(0.67086733) q[2];
sx q[2];
rz(-0.76422526) q[2];
sx q[2];
rz(-2.4799181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5246626) q[1];
sx q[1];
rz(-1.3549651) q[1];
sx q[1];
rz(0.642435) q[1];
x q[2];
rz(-0.60634585) q[3];
sx q[3];
rz(-0.56713533) q[3];
sx q[3];
rz(-2.4727269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(-0.62270069) q[2];
rz(2.4032118) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(-2.8625989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57711346) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(-2.4979757) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(-1.0254597) q[2];
sx q[2];
rz(-2.7164216) q[2];
sx q[2];
rz(2.807775) q[2];
rz(-0.24091992) q[3];
sx q[3];
rz(-1.4333457) q[3];
sx q[3];
rz(-2.3960927) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
