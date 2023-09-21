OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0599521) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(-0.089846213) q[0];
rz(-pi) q[1];
rz(-0.88959496) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(1.7878469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-0.70322137) q[1];
sx q[1];
rz(-2.1232848) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6775042) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(-0.97809631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-2.1957943) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5058274) q[0];
sx q[0];
rz(-2.8553537) q[0];
sx q[0];
rz(-0.92606996) q[0];
rz(-pi) q[1];
rz(-0.63983812) q[2];
sx q[2];
rz(-2.8176762) q[2];
sx q[2];
rz(-1.6537635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82169689) q[1];
sx q[1];
rz(-0.79274717) q[1];
sx q[1];
rz(-2.4577623) q[1];
rz(-pi) q[2];
rz(-1.3974959) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(-1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2382425) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(-1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.7984614) q[0];
rz(-0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(2.6599191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96268481) q[0];
sx q[0];
rz(-2.509153) q[0];
sx q[0];
rz(0.90895598) q[0];
x q[1];
rz(-0.9573612) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(-2.7797109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8183221) q[1];
sx q[1];
rz(-0.96410492) q[1];
sx q[1];
rz(-0.89241772) q[1];
x q[2];
rz(0.84677245) q[3];
sx q[3];
rz(-1.0792884) q[3];
sx q[3];
rz(1.5848974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-0.552185) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.2447371) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075605) q[0];
sx q[0];
rz(-1.2313156) q[0];
sx q[0];
rz(1.5725122) q[0];
rz(-pi) q[1];
rz(-2.142698) q[2];
sx q[2];
rz(-2.9036387) q[2];
sx q[2];
rz(1.1513125) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(1.0210387) q[1];
x q[2];
rz(-2.5707385) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3969288) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(2.1555107) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0585074) q[2];
sx q[2];
rz(-1.6787046) q[2];
sx q[2];
rz(1.2667058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7170143) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(-1.2574408) q[1];
x q[2];
rz(2.6322707) q[3];
sx q[3];
rz(-1.29042) q[3];
sx q[3];
rz(0.59685368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6333255) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(2.6200263) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-0.26671985) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-2.3430603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079433867) q[0];
sx q[0];
rz(-1.845813) q[0];
sx q[0];
rz(2.9955203) q[0];
rz(-1.2217667) q[2];
sx q[2];
rz(-2.0418228) q[2];
sx q[2];
rz(2.7407321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3244464) q[1];
sx q[1];
rz(-1.518461) q[1];
sx q[1];
rz(-2.1423992) q[1];
rz(-0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-2.7748761) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(-2.6775449) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(-0.52857907) q[0];
rz(1.3299994) q[2];
sx q[2];
rz(-1.5903928) q[2];
sx q[2];
rz(-0.400825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0276427) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(-0.34160683) q[1];
x q[2];
rz(-0.32902284) q[3];
sx q[3];
rz(-2.892422) q[3];
sx q[3];
rz(2.3014625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15526074) q[0];
sx q[0];
rz(-3.0763456) q[0];
sx q[0];
rz(-0.23373993) q[0];
rz(-pi) q[1];
rz(2.5289815) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(-0.14234662) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0563593) q[1];
sx q[1];
rz(-1.7006405) q[1];
sx q[1];
rz(1.4443881) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0184047) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5817225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(-0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50842677) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(-1.7500061) q[0];
rz(-0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(-1.6800113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1074859) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
x q[2];
rz(2.968623) q[3];
sx q[3];
rz(-1.1555472) q[3];
sx q[3];
rz(-3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(0.40965664) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423994) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(-1.3347795) q[0];
sx q[0];
rz(0.34061265) q[0];
x q[1];
rz(-1.2070933) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(-1.9050913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.532383) q[1];
sx q[1];
rz(-1.7222952) q[1];
sx q[1];
rz(-1.7864368) q[1];
rz(-pi) q[2];
rz(0.16898445) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(-2.3224725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60349764) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-0.60097354) q[2];
sx q[2];
rz(-2.8291611) q[2];
sx q[2];
rz(-0.57401382) q[2];
rz(-0.55660558) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
