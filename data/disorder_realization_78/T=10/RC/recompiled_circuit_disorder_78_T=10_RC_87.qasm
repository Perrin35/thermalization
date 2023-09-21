OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8569782) q[0];
sx q[0];
rz(-2.1224788) q[0];
sx q[0];
rz(-0.3509699) q[0];
rz(2.4721488) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(0.66561156) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6751911) q[1];
sx q[1];
rz(-0.89414222) q[1];
sx q[1];
rz(-2.7387709) q[1];
rz(-0.046269429) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-1.083667) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(-0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(1.1670246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3008227) q[0];
sx q[0];
rz(-1.6292028) q[0];
sx q[0];
rz(-1.8286684) q[0];
rz(-pi) q[1];
rz(-1.0300693) q[2];
sx q[2];
rz(-1.7833372) q[2];
sx q[2];
rz(-2.9589047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0868623) q[1];
sx q[1];
rz(-0.80838258) q[1];
sx q[1];
rz(-3.0562835) q[1];
rz(0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(-3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(2.8607821) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1983564) q[0];
sx q[0];
rz(-0.6456607) q[0];
sx q[0];
rz(-0.093430324) q[0];
rz(-1.1655408) q[2];
sx q[2];
rz(-1.3477256) q[2];
sx q[2];
rz(-2.9802393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44753578) q[1];
sx q[1];
rz(-0.60255614) q[1];
sx q[1];
rz(1.4284929) q[1];
x q[2];
rz(-2.879911) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(-2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(2.5337059) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-0.63582173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83946562) q[0];
sx q[0];
rz(-1.7502022) q[0];
sx q[0];
rz(-2.0749712) q[0];
rz(-2.9767838) q[2];
sx q[2];
rz(-2.2909082) q[2];
sx q[2];
rz(2.1317496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3969877) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(-2.8810487) q[1];
rz(0.28663978) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(0.74433792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9277966) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(2.863046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51019788) q[0];
sx q[0];
rz(-0.83980951) q[0];
sx q[0];
rz(0.98548074) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4611545) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.4102175) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7601732) q[1];
sx q[1];
rz(-2.1790494) q[1];
sx q[1];
rz(-0.45254032) q[1];
x q[2];
rz(0.99028011) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(-3.0596628) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(0.10805282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6022588) q[0];
sx q[0];
rz(-0.56319153) q[0];
sx q[0];
rz(-2.5186033) q[0];
x q[1];
rz(0.33424218) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(0.7451171) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9891226) q[1];
sx q[1];
rz(-1.1835386) q[1];
sx q[1];
rz(0.5247922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.990854) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(-0.95668876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(0.89362842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4511787) q[0];
sx q[0];
rz(-2.4126841) q[0];
sx q[0];
rz(-2.6364987) q[0];
x q[1];
rz(-2.8473179) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(2.0725046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(2.8303353) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10223933) q[3];
sx q[3];
rz(-2.3792301) q[3];
sx q[3];
rz(-1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40076462) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-0.056218938) q[0];
rz(-pi) q[1];
rz(2.5904028) q[2];
sx q[2];
rz(-1.6973719) q[2];
sx q[2];
rz(1.7984496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(2.3803821) q[1];
x q[2];
rz(-0.89973255) q[3];
sx q[3];
rz(-1.2601488) q[3];
sx q[3];
rz(1.336738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026924883) q[0];
sx q[0];
rz(-1.9244734) q[0];
sx q[0];
rz(-2.7915519) q[0];
x q[1];
rz(2.3652472) q[2];
sx q[2];
rz(-1.7661957) q[2];
sx q[2];
rz(-2.5459144) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(-2.7705454) q[1];
x q[2];
rz(-3.0036003) q[3];
sx q[3];
rz(-2.237461) q[3];
sx q[3];
rz(2.7033412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(0.38273746) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-0.66463566) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(2.6616667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9057248) q[0];
sx q[0];
rz(-2.1972482) q[0];
sx q[0];
rz(-0.41525526) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2357113) q[2];
sx q[2];
rz(-1.7582298) q[2];
sx q[2];
rz(0.71042176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2952134) q[1];
sx q[1];
rz(-1.0615674) q[1];
sx q[1];
rz(1.9766115) q[1];
rz(-2.2009833) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(-1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(2.4035113) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823572) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(-2.813415) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(-0.76673037) q[3];
sx q[3];
rz(-0.37692108) q[3];
sx q[3];
rz(2.2212096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];