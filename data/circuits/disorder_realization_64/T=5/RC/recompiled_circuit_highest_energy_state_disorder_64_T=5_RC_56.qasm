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
rz(-0.91941961) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(-0.71370178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8827191) q[0];
sx q[0];
rz(-1.2964411) q[0];
sx q[0];
rz(0.15299847) q[0];
x q[1];
rz(-1.6993148) q[2];
sx q[2];
rz(-1.5084195) q[2];
sx q[2];
rz(-0.85278748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9992396) q[1];
sx q[1];
rz(-1.1198143) q[1];
sx q[1];
rz(-1.5723448) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1072989) q[3];
sx q[3];
rz(-2.0418797) q[3];
sx q[3];
rz(2.8054939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0539703) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358148) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(-2.0447958) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(-2.6699452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7268611) q[0];
sx q[0];
rz(-1.4607753) q[0];
sx q[0];
rz(0.96745205) q[0];
x q[1];
rz(0.94309965) q[2];
sx q[2];
rz(-0.97979704) q[2];
sx q[2];
rz(0.87795382) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1321137) q[1];
sx q[1];
rz(-1.0227362) q[1];
sx q[1];
rz(2.4365042) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0949267) q[3];
sx q[3];
rz(-1.9283623) q[3];
sx q[3];
rz(-2.3563354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4196709) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(-2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(-0.56627083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6951272) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(-1.5873853) q[0];
rz(-0.71193169) q[2];
sx q[2];
rz(-1.1441137) q[2];
sx q[2];
rz(0.79244765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0349087) q[1];
sx q[1];
rz(-1.3843754) q[1];
sx q[1];
rz(-1.6572464) q[1];
rz(1.6376466) q[3];
sx q[3];
rz(-1.306157) q[3];
sx q[3];
rz(-0.29936779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8038586) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424054) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(1.1135788) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(2.9794433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55570468) q[0];
sx q[0];
rz(-1.7506208) q[0];
sx q[0];
rz(0.094688133) q[0];
rz(3.1012898) q[2];
sx q[2];
rz(-0.54931123) q[2];
sx q[2];
rz(1.2072762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.051013324) q[1];
sx q[1];
rz(-0.83003488) q[1];
sx q[1];
rz(0.49511893) q[1];
x q[2];
rz(-0.6643296) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(0.20928247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1589511) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(2.2900901) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.270179) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(-0.98091006) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0410515) q[0];
sx q[0];
rz(-1.9773736) q[0];
sx q[0];
rz(1.340614) q[0];
rz(-0.95246117) q[2];
sx q[2];
rz(-0.59788579) q[2];
sx q[2];
rz(0.43274227) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6740711) q[1];
sx q[1];
rz(-1.5115807) q[1];
sx q[1];
rz(-1.6644415) q[1];
rz(2.122287) q[3];
sx q[3];
rz(-1.3634342) q[3];
sx q[3];
rz(0.99121782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(-2.2583466) q[2];
rz(-2.9412681) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(-2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14466318) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(-0.99217478) q[0];
rz(0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(-1.4206402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.84467) q[0];
sx q[0];
rz(-1.2529182) q[0];
sx q[0];
rz(2.4852089) q[0];
x q[1];
rz(-1.0384485) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(0.20648512) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0040505) q[1];
sx q[1];
rz(-2.1909449) q[1];
sx q[1];
rz(1.1637079) q[1];
rz(-1.4735953) q[3];
sx q[3];
rz(-1.5335113) q[3];
sx q[3];
rz(-2.4241222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-0.53708616) q[2];
rz(3.1308657) q[3];
sx q[3];
rz(-0.74672943) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(-2.8017092) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-0.94568959) q[1];
sx q[1];
rz(-1.1713015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.834496) q[0];
sx q[0];
rz(-2.8423873) q[0];
sx q[0];
rz(3.0181969) q[0];
rz(-pi) q[1];
rz(-0.76979678) q[2];
sx q[2];
rz(-0.30210051) q[2];
sx q[2];
rz(3.0506899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.724546) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.538365) q[1];
rz(-pi) q[2];
rz(-0.24089916) q[3];
sx q[3];
rz(-2.5031075) q[3];
sx q[3];
rz(-2.1158259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-2.4701414) q[2];
rz(2.6050383) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3143828) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(0.89624727) q[0];
rz(2.8889636) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80879656) q[0];
sx q[0];
rz(-1.044863) q[0];
sx q[0];
rz(1.512708) q[0];
x q[1];
rz(-0.66506135) q[2];
sx q[2];
rz(-1.6740693) q[2];
sx q[2];
rz(2.2879083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72880367) q[1];
sx q[1];
rz(-2.0319684) q[1];
sx q[1];
rz(-2.9739266) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3837621) q[3];
sx q[3];
rz(-1.5327454) q[3];
sx q[3];
rz(1.1421957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62319055) q[2];
sx q[2];
rz(-0.16885997) q[2];
sx q[2];
rz(1.5625578) q[2];
rz(1.7744428) q[3];
sx q[3];
rz(-2.0953777) q[3];
sx q[3];
rz(-0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.389483) q[0];
sx q[0];
rz(2.7516464) q[0];
rz(-0.82603106) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(1.0521851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5717333) q[0];
sx q[0];
rz(-2.2433407) q[0];
sx q[0];
rz(0.39542266) q[0];
x q[1];
rz(0.40014123) q[2];
sx q[2];
rz(-0.42707983) q[2];
sx q[2];
rz(-2.6838322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6893951) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(2.2188975) q[1];
x q[2];
rz(0.47635079) q[3];
sx q[3];
rz(-0.99925502) q[3];
sx q[3];
rz(1.0494572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(0.52919394) q[2];
rz(0.75096327) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(0.19415893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.1645114) q[0];
rz(-1.8544082) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(-0.50416344) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7499649) q[0];
sx q[0];
rz(-0.34191439) q[0];
sx q[0];
rz(2.2750399) q[0];
rz(-pi) q[1];
rz(-0.37673925) q[2];
sx q[2];
rz(-1.5887056) q[2];
sx q[2];
rz(-1.1056545) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.092445651) q[1];
sx q[1];
rz(-1.7610252) q[1];
sx q[1];
rz(-1.3008429) q[1];
rz(-pi) q[2];
rz(2.9243117) q[3];
sx q[3];
rz(-2.2437527) q[3];
sx q[3];
rz(2.6126044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(-0.034493383) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9813949) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(-1.002671) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(2.6832081) q[2];
sx q[2];
rz(-1.7469363) q[2];
sx q[2];
rz(1.7652702) q[2];
rz(-1.6335842) q[3];
sx q[3];
rz(-1.2442524) q[3];
sx q[3];
rz(2.2942345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
