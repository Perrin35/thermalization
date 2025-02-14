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
rz(-0.24124423) q[0];
sx q[0];
rz(-0.21114199) q[0];
sx q[0];
rz(0.40786064) q[0];
rz(-1.7653699) q[1];
sx q[1];
rz(5.0566109) q[1];
sx q[1];
rz(15.772718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415675) q[0];
sx q[0];
rz(-1.5655101) q[0];
sx q[0];
rz(1.4553817) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32081713) q[2];
sx q[2];
rz(-2.1324369) q[2];
sx q[2];
rz(0.80036847) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3169466) q[1];
sx q[1];
rz(-1.4295409) q[1];
sx q[1];
rz(1.5596175) q[1];
rz(2.519386) q[3];
sx q[3];
rz(-1.4861306) q[3];
sx q[3];
rz(2.0537927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9593418) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(1.135929) q[2];
rz(-0.7302537) q[3];
sx q[3];
rz(-1.3948995) q[3];
sx q[3];
rz(1.6212538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.092904329) q[0];
sx q[0];
rz(-0.95656675) q[0];
sx q[0];
rz(-0.055305716) q[0];
rz(2.7703908) q[1];
sx q[1];
rz(-2.7822045) q[1];
sx q[1];
rz(-0.41195437) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13875554) q[0];
sx q[0];
rz(-1.635972) q[0];
sx q[0];
rz(-2.6145593) q[0];
rz(-pi) q[1];
rz(-2.2949176) q[2];
sx q[2];
rz(-2.7077423) q[2];
sx q[2];
rz(0.84956079) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4476599) q[1];
sx q[1];
rz(-0.74697633) q[1];
sx q[1];
rz(-2.9332317) q[1];
rz(-pi) q[2];
rz(-2.7633414) q[3];
sx q[3];
rz(-0.95767271) q[3];
sx q[3];
rz(-2.48191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61341316) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(-2.2589653) q[2];
rz(1.708185) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9751137) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(-2.574918) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(-1.3272939) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6782129) q[0];
sx q[0];
rz(-0.87317077) q[0];
sx q[0];
rz(0.11059983) q[0];
rz(0.15285531) q[2];
sx q[2];
rz(-2.1899583) q[2];
sx q[2];
rz(2.4094606) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8019405) q[1];
sx q[1];
rz(-1.2627456) q[1];
sx q[1];
rz(-3.1402474) q[1];
rz(-1.0653082) q[3];
sx q[3];
rz(-0.79366854) q[3];
sx q[3];
rz(-2.6738559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3920445) q[2];
sx q[2];
rz(-2.0508524) q[2];
sx q[2];
rz(-0.90847477) q[2];
rz(-0.79489094) q[3];
sx q[3];
rz(-1.1029693) q[3];
sx q[3];
rz(-0.046549646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817151) q[0];
sx q[0];
rz(-0.91232038) q[0];
sx q[0];
rz(-0.6657486) q[0];
rz(-0.84716973) q[1];
sx q[1];
rz(-2.8916292) q[1];
sx q[1];
rz(0.4096823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100135) q[0];
sx q[0];
rz(-2.0635475) q[0];
sx q[0];
rz(-1.3973622) q[0];
x q[1];
rz(-1.9471859) q[2];
sx q[2];
rz(-0.79845286) q[2];
sx q[2];
rz(1.2898766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53654251) q[1];
sx q[1];
rz(-2.4808421) q[1];
sx q[1];
rz(1.9416757) q[1];
x q[2];
rz(2.014854) q[3];
sx q[3];
rz(-1.5448625) q[3];
sx q[3];
rz(-1.3249205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72982558) q[2];
sx q[2];
rz(-0.7080141) q[2];
sx q[2];
rz(1.7626308) q[2];
rz(2.0046558) q[3];
sx q[3];
rz(-1.3552908) q[3];
sx q[3];
rz(2.304145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1788504) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(2.7233997) q[0];
rz(1.7182982) q[1];
sx q[1];
rz(-1.9530674) q[1];
sx q[1];
rz(1.6421912) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20883501) q[0];
sx q[0];
rz(-1.8153738) q[0];
sx q[0];
rz(-0.45758684) q[0];
rz(-pi) q[1];
rz(-1.2659433) q[2];
sx q[2];
rz(-1.8902167) q[2];
sx q[2];
rz(1.4679421) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.114257) q[1];
sx q[1];
rz(-1.789195) q[1];
sx q[1];
rz(-0.059192358) q[1];
x q[2];
rz(1.70056) q[3];
sx q[3];
rz(-1.6356182) q[3];
sx q[3];
rz(0.45634064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8666009) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(-0.17821136) q[2];
rz(0.37402672) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(-1.1427574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571534) q[0];
sx q[0];
rz(-1.3244119) q[0];
sx q[0];
rz(1.4763747) q[0];
rz(0.58552512) q[1];
sx q[1];
rz(-2.245677) q[1];
sx q[1];
rz(2.4077328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7134345) q[0];
sx q[0];
rz(-2.1913986) q[0];
sx q[0];
rz(1.0125068) q[0];
x q[1];
rz(2.3247277) q[2];
sx q[2];
rz(-0.43359633) q[2];
sx q[2];
rz(-2.9006633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7712931) q[1];
sx q[1];
rz(-1.7669189) q[1];
sx q[1];
rz(1.0443772) q[1];
x q[2];
rz(-1.1595585) q[3];
sx q[3];
rz(-2.6542668) q[3];
sx q[3];
rz(-0.05035487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0840941) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(-1.499184) q[2];
rz(-1.5205787) q[3];
sx q[3];
rz(-2.4937544) q[3];
sx q[3];
rz(-0.15629855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-2.2235625) q[0];
sx q[0];
rz(1.3936438) q[0];
rz(-2.5539894) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(-2.844152) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2327691) q[0];
sx q[0];
rz(-0.89134514) q[0];
sx q[0];
rz(2.6772236) q[0];
rz(-pi) q[1];
rz(-3.1226899) q[2];
sx q[2];
rz(-1.3002965) q[2];
sx q[2];
rz(-2.8342385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0175727) q[1];
sx q[1];
rz(-0.67970733) q[1];
sx q[1];
rz(-1.5862096) q[1];
rz(-pi) q[2];
rz(2.8964306) q[3];
sx q[3];
rz(-0.94221598) q[3];
sx q[3];
rz(2.6338284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3191159) q[2];
sx q[2];
rz(-1.7812704) q[2];
sx q[2];
rz(0.3817257) q[2];
rz(2.5413359) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(2.8762347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(1.3375244) q[0];
rz(3.1321101) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(1.7705852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776896) q[0];
sx q[0];
rz(-2.7743917) q[0];
sx q[0];
rz(-2.7504671) q[0];
x q[1];
rz(-0.96641936) q[2];
sx q[2];
rz(-0.76497173) q[2];
sx q[2];
rz(-1.4422095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0095894) q[1];
sx q[1];
rz(-1.6958964) q[1];
sx q[1];
rz(2.3522207) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.480732) q[3];
sx q[3];
rz(-1.9833296) q[3];
sx q[3];
rz(-1.6270669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3235772) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(-2.9991007) q[2];
rz(-2.2583708) q[3];
sx q[3];
rz(-1.764069) q[3];
sx q[3];
rz(1.6787329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0505117) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.8930513) q[0];
rz(2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(1.1936845) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79289397) q[0];
sx q[0];
rz(-1.3137709) q[0];
sx q[0];
rz(-1.2299163) q[0];
rz(-2.2676554) q[2];
sx q[2];
rz(-0.94183445) q[2];
sx q[2];
rz(0.12128092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15577573) q[1];
sx q[1];
rz(-1.078842) q[1];
sx q[1];
rz(-1.1883177) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5709236) q[3];
sx q[3];
rz(-2.8817892) q[3];
sx q[3];
rz(0.75614959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(-1.1129145) q[2];
rz(-0.042996081) q[3];
sx q[3];
rz(-1.6578511) q[3];
sx q[3];
rz(-0.23785166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766895) q[0];
sx q[0];
rz(-1.7273128) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(-1.9211357) q[1];
sx q[1];
rz(-2.7595322) q[1];
sx q[1];
rz(-0.97741309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138151) q[0];
sx q[0];
rz(-1.312698) q[0];
sx q[0];
rz(-1.4679704) q[0];
x q[1];
rz(-1.270411) q[2];
sx q[2];
rz(-1.4255878) q[2];
sx q[2];
rz(-1.6329488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6466684) q[1];
sx q[1];
rz(-2.1416602) q[1];
sx q[1];
rz(2.5350556) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39267003) q[3];
sx q[3];
rz(-1.304198) q[3];
sx q[3];
rz(2.6368898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0105373) q[2];
sx q[2];
rz(-2.8752893) q[2];
sx q[2];
rz(-1.4492501) q[2];
rz(-2.2629755) q[3];
sx q[3];
rz(-2.0802458) q[3];
sx q[3];
rz(-1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286248) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(-1.3799008) q[1];
sx q[1];
rz(-1.5329783) q[1];
sx q[1];
rz(3.0402532) q[1];
rz(0.18086441) q[2];
sx q[2];
rz(-1.7500063) q[2];
sx q[2];
rz(1.0371006) q[2];
rz(2.8158549) q[3];
sx q[3];
rz(-0.67254638) q[3];
sx q[3];
rz(-1.2253958) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
