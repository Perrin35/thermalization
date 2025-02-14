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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(3.0923831) q[1];
sx q[1];
rz(10.13848) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714234) q[0];
sx q[0];
rz(-1.7180301) q[0];
sx q[0];
rz(1.2933613) q[0];
x q[1];
rz(-2.0242516) q[2];
sx q[2];
rz(-2.9988117) q[2];
sx q[2];
rz(2.873024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9992396) q[1];
sx q[1];
rz(-2.0217784) q[1];
sx q[1];
rz(1.5723448) q[1];
rz(-pi) q[2];
rz(-0.034293745) q[3];
sx q[3];
rz(-2.0418797) q[3];
sx q[3];
rz(-0.33609875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.087622341) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-2.7719882) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358148) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(2.5685487) q[0];
rz(-1.0967968) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(0.47164741) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4147316) q[0];
sx q[0];
rz(-1.4607753) q[0];
sx q[0];
rz(0.96745205) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94309965) q[2];
sx q[2];
rz(-0.97979704) q[2];
sx q[2];
rz(-2.2636388) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.009479) q[1];
sx q[1];
rz(-2.1188564) q[1];
sx q[1];
rz(-2.4365042) q[1];
x q[2];
rz(1.0466659) q[3];
sx q[3];
rz(-1.2132304) q[3];
sx q[3];
rz(-0.78525728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(-2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0786781) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(-0.71006376) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(2.5753218) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1120934) q[0];
sx q[0];
rz(-1.5819966) q[0];
sx q[0];
rz(-2.311934) q[0];
x q[1];
rz(0.71193169) q[2];
sx q[2];
rz(-1.997479) q[2];
sx q[2];
rz(-2.349145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33070811) q[1];
sx q[1];
rz(-0.20527923) q[1];
sx q[1];
rz(2.7122871) q[1];
rz(-pi) q[2];
x q[2];
rz(1.503946) q[3];
sx q[3];
rz(-1.8354356) q[3];
sx q[3];
rz(2.8422249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8038586) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(0.13633063) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(1.2135308) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424054) q[0];
sx q[0];
rz(-2.461705) q[0];
sx q[0];
rz(2.9754382) q[0];
rz(-1.1135788) q[1];
sx q[1];
rz(-0.49003092) q[1];
sx q[1];
rz(-0.16214935) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0320764) q[0];
sx q[0];
rz(-1.4776395) q[0];
sx q[0];
rz(1.3901802) q[0];
rz(-1.5461363) q[2];
sx q[2];
rz(-2.1196105) q[2];
sx q[2];
rz(1.8870712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6238193) q[1];
sx q[1];
rz(-0.86408593) q[1];
sx q[1];
rz(-1.0916187) q[1];
x q[2];
rz(-2.6147551) q[3];
sx q[3];
rz(-0.73606811) q[3];
sx q[3];
rz(-1.3730622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98264155) q[2];
sx q[2];
rz(-2.0621767) q[2];
sx q[2];
rz(0.52784935) q[2];
rz(-2.2900901) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(-2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8714137) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(-0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(-2.1606826) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5194395) q[0];
sx q[0];
rz(-1.3596757) q[0];
sx q[0];
rz(0.4163752) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0641685) q[2];
sx q[2];
rz(-1.2384103) q[2];
sx q[2];
rz(1.6695963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4760433) q[1];
sx q[1];
rz(-3.0308423) q[1];
sx q[1];
rz(1.0057463) q[1];
rz(-2.8994336) q[3];
sx q[3];
rz(-2.1091614) q[3];
sx q[3];
rz(2.4360366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(-0.8832461) q[2];
rz(0.20032459) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(-2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14466318) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-0.99217478) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(-1.7209524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.84467) q[0];
sx q[0];
rz(-1.2529182) q[0];
sx q[0];
rz(-2.4852089) q[0];
rz(-1.0384485) q[2];
sx q[2];
rz(-1.4803998) q[2];
sx q[2];
rz(-0.20648512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7758549) q[1];
sx q[1];
rz(-0.72682646) q[1];
sx q[1];
rz(0.50623843) q[1];
x q[2];
rz(3.104131) q[3];
sx q[3];
rz(-1.473663) q[3];
sx q[3];
rz(2.284632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27986032) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(3.1308657) q[3];
sx q[3];
rz(-0.74672943) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456197) q[0];
sx q[0];
rz(-0.27181044) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(1.3019568) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.9702912) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070967) q[0];
sx q[0];
rz(-0.2992054) q[0];
sx q[0];
rz(-3.0181969) q[0];
rz(-pi) q[1];
rz(1.7843855) q[2];
sx q[2];
rz(-1.7860951) q[2];
sx q[2];
rz(2.4395296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.724546) q[1];
sx q[1];
rz(-1.5486716) q[1];
sx q[1];
rz(1.538365) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62452353) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(-0.7398015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0168212) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(0.67145124) q[2];
rz(2.6050383) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(2.0898537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.8272098) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(-0.89624727) q[0];
rz(0.25262901) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(-2.0910697) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80879656) q[0];
sx q[0];
rz(-1.044863) q[0];
sx q[0];
rz(-1.6288847) q[0];
rz(-2.975198) q[2];
sx q[2];
rz(-0.67182589) q[2];
sx q[2];
rz(0.8478129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.412789) q[1];
sx q[1];
rz(-2.0319684) q[1];
sx q[1];
rz(0.16766602) q[1];
rz(-3.0862634) q[3];
sx q[3];
rz(-2.3829975) q[3];
sx q[3];
rz(0.38842312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5184021) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(1.5790348) q[2];
rz(1.3671499) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(-0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(-2.7516464) q[0];
rz(-0.82603106) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(1.0521851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74653502) q[0];
sx q[0];
rz(-1.2647226) q[0];
sx q[0];
rz(-2.2827882) q[0];
x q[1];
rz(2.7446943) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(-1.6610749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4521976) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(-2.2188975) q[1];
x q[2];
rz(0.94433208) q[3];
sx q[3];
rz(-1.1748702) q[3];
sx q[3];
rz(2.348071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(2.6123987) q[2];
rz(0.75096327) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(1.9770812) q[0];
rz(-1.8544082) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(-2.6374292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8541478) q[0];
sx q[0];
rz(-1.7896255) q[0];
sx q[0];
rz(-1.3059421) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0929448) q[2];
sx q[2];
rz(-2.7644483) q[2];
sx q[2];
rz(-2.6312171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.092445651) q[1];
sx q[1];
rz(-1.7610252) q[1];
sx q[1];
rz(-1.3008429) q[1];
x q[2];
rz(0.88621288) q[3];
sx q[3];
rz(-1.4014114) q[3];
sx q[3];
rz(-2.2365295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9265499) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(3.1070993) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(-2.8300986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9813949) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(-1.002671) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(-1.3748693) q[2];
sx q[2];
rz(-2.0215603) q[2];
sx q[2];
rz(0.10822288) q[2];
rz(2.9583951) q[3];
sx q[3];
rz(-2.8092794) q[3];
sx q[3];
rz(2.487779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
