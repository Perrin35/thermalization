OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(-2.0816198) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19883991) q[0];
sx q[0];
rz(-1.4497541) q[0];
sx q[0];
rz(-0.11902703) q[0];
x q[1];
rz(1.5358309) q[2];
sx q[2];
rz(-0.87858534) q[2];
sx q[2];
rz(-3.0193605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4958475) q[1];
sx q[1];
rz(-1.8105257) q[1];
sx q[1];
rz(-1.7868397) q[1];
x q[2];
rz(-2.5422402) q[3];
sx q[3];
rz(-2.9671892) q[3];
sx q[3];
rz(-0.50267437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9609191) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(-1.4291576) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(-0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83577689) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(-2.3961156) q[0];
rz(1.9396793) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(2.1048996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1069477) q[0];
sx q[0];
rz(-1.6634403) q[0];
sx q[0];
rz(1.9361467) q[0];
rz(-pi) q[1];
rz(-2.3747524) q[2];
sx q[2];
rz(-1.5897703) q[2];
sx q[2];
rz(-2.5350867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17192852) q[1];
sx q[1];
rz(-1.8542037) q[1];
sx q[1];
rz(-0.71782459) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5573331) q[3];
sx q[3];
rz(-1.7528894) q[3];
sx q[3];
rz(0.77891536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0006973) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(2.1689283) q[2];
rz(-2.2843212) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(1.2681786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.641441) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(-0.80297536) q[0];
rz(-0.650644) q[1];
sx q[1];
rz(-0.79901189) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652199) q[0];
sx q[0];
rz(-2.1928761) q[0];
sx q[0];
rz(-2.6724044) q[0];
rz(-0.62816633) q[2];
sx q[2];
rz(-1.7587314) q[2];
sx q[2];
rz(-2.4357901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.928339) q[1];
sx q[1];
rz(-1.6633185) q[1];
sx q[1];
rz(1.0800581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.337255) q[3];
sx q[3];
rz(-1.8982049) q[3];
sx q[3];
rz(1.4191607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(-1.6311084) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.889582) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(2.4776283) q[0];
rz(-0.12531677) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(-2.1024316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17500278) q[0];
sx q[0];
rz(-3.1313854) q[0];
sx q[0];
rz(-2.8673579) q[0];
rz(-1.9661994) q[2];
sx q[2];
rz(-1.8186453) q[2];
sx q[2];
rz(1.6958267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81280983) q[1];
sx q[1];
rz(-2.610958) q[1];
sx q[1];
rz(-1.1483436) q[1];
x q[2];
rz(2.6986946) q[3];
sx q[3];
rz(-1.4848466) q[3];
sx q[3];
rz(1.7922873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30109721) q[2];
sx q[2];
rz(-1.7258464) q[2];
sx q[2];
rz(0.077979716) q[2];
rz(-2.0237427) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(-2.1054721) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9341105) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(1.9059937) q[1];
sx q[1];
rz(-1.3256336) q[1];
sx q[1];
rz(1.8280169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395555) q[0];
sx q[0];
rz(-1.3875204) q[0];
sx q[0];
rz(-2.2048143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2143873) q[2];
sx q[2];
rz(-2.5816133) q[2];
sx q[2];
rz(0.72712979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8660248) q[1];
sx q[1];
rz(-1.1236262) q[1];
sx q[1];
rz(0.42110301) q[1];
rz(-0.4540654) q[3];
sx q[3];
rz(-2.3010074) q[3];
sx q[3];
rz(-2.8968352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95973394) q[2];
sx q[2];
rz(-1.9071969) q[2];
sx q[2];
rz(-2.6971297) q[2];
rz(1.2005165) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(-0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926517) q[0];
sx q[0];
rz(-0.3949202) q[0];
sx q[0];
rz(-0.077202395) q[0];
rz(2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(3.107792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2709406) q[0];
sx q[0];
rz(-2.0217289) q[0];
sx q[0];
rz(1.1876039) q[0];
x q[1];
rz(-1.5438429) q[2];
sx q[2];
rz(-2.1183062) q[2];
sx q[2];
rz(1.6545878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2359743) q[1];
sx q[1];
rz(-1.8840569) q[1];
sx q[1];
rz(-2.226023) q[1];
rz(-pi) q[2];
rz(-1.46557) q[3];
sx q[3];
rz(-0.79685539) q[3];
sx q[3];
rz(1.6355255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(-2.1974468) q[3];
sx q[3];
rz(-2.9278432) q[3];
sx q[3];
rz(1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0807226) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(-0.61332214) q[0];
rz(2.2382286) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(-0.14796743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3053235) q[0];
sx q[0];
rz(-0.9137872) q[0];
sx q[0];
rz(2.5227274) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.067909165) q[2];
sx q[2];
rz(-2.7017044) q[2];
sx q[2];
rz(-1.875017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92098713) q[1];
sx q[1];
rz(-1.9927532) q[1];
sx q[1];
rz(0.76442952) q[1];
rz(-pi) q[2];
rz(-0.023588019) q[3];
sx q[3];
rz(-1.0057431) q[3];
sx q[3];
rz(-0.20701359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.938574) q[2];
sx q[2];
rz(-1.4116986) q[2];
sx q[2];
rz(2.1417248) q[2];
rz(1.3153007) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(-0.6704754) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(-1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1647427) q[0];
sx q[0];
rz(-1.688042) q[0];
sx q[0];
rz(-1.4848723) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0664178) q[2];
sx q[2];
rz(-1.3775423) q[2];
sx q[2];
rz(-2.9444864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2439448) q[1];
sx q[1];
rz(-0.49265322) q[1];
sx q[1];
rz(-1.5642464) q[1];
rz(1.3173728) q[3];
sx q[3];
rz(-0.66605036) q[3];
sx q[3];
rz(-2.1061153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8998731) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(-1.0618173) q[2];
rz(3.0070983) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(0.053248052) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(2.9363976) q[0];
rz(2.9070053) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(2.9451784) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043606) q[0];
sx q[0];
rz(-2.4825486) q[0];
sx q[0];
rz(1.0511257) q[0];
rz(0.5242879) q[2];
sx q[2];
rz(-1.7066188) q[2];
sx q[2];
rz(2.6581531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3141296) q[1];
sx q[1];
rz(-1.2478831) q[1];
sx q[1];
rz(-2.9548378) q[1];
rz(2.6077723) q[3];
sx q[3];
rz(-2.2861423) q[3];
sx q[3];
rz(1.8623587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(-1.9632001) q[2];
rz(-0.34934238) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(-1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0013393764) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(2.4358791) q[0];
rz(0.69752518) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(2.2360905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0374091) q[0];
sx q[0];
rz(-1.5539123) q[0];
sx q[0];
rz(2.707858) q[0];
rz(2.414213) q[2];
sx q[2];
rz(-2.670521) q[2];
sx q[2];
rz(1.7076275) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0841964) q[1];
sx q[1];
rz(-1.2603972) q[1];
sx q[1];
rz(1.5211902) q[1];
x q[2];
rz(0.3304146) q[3];
sx q[3];
rz(-1.2252231) q[3];
sx q[3];
rz(-1.6890077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(0.71851292) q[2];
rz(-2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.318442) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.6406583) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(2.4544302) q[2];
sx q[2];
rz(-2.9570815) q[2];
sx q[2];
rz(2.2312192) q[2];
rz(-0.90632306) q[3];
sx q[3];
rz(-2.255848) q[3];
sx q[3];
rz(-2.0391603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
