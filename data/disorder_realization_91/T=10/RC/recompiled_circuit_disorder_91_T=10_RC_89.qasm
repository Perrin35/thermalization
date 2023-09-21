OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(-0.17106549) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2301807) q[0];
sx q[0];
rz(-0.6517621) q[0];
sx q[0];
rz(0.08641152) q[0];
rz(1.5924442) q[2];
sx q[2];
rz(-1.0550371) q[2];
sx q[2];
rz(1.5390918) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5142071) q[1];
sx q[1];
rz(-1.5537019) q[1];
sx q[1];
rz(-1.3196726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29403789) q[3];
sx q[3];
rz(-2.4774744) q[3];
sx q[3];
rz(0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3657637) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-0.20234385) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(1.6888899) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(2.8413049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18314221) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(-1.4248225) q[0];
rz(-pi) q[1];
rz(2.5710201) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(0.74769339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7784361) q[1];
sx q[1];
rz(-1.2769715) q[1];
sx q[1];
rz(2.9086962) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3817915) q[3];
sx q[3];
rz(-1.731589) q[3];
sx q[3];
rz(-2.8652428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(-0.38811362) q[0];
rz(3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(-2.8252576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91352275) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(-1.8158005) q[0];
rz(-pi) q[1];
rz(-0.6109654) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(1.9463469) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12831941) q[1];
sx q[1];
rz(-2.0187906) q[1];
sx q[1];
rz(-2.0112579) q[1];
rz(-pi) q[2];
rz(2.8510677) q[3];
sx q[3];
rz(-1.2447262) q[3];
sx q[3];
rz(2.3323004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(-1.3446993) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.6279189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095703) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(3.0679697) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4430181) q[2];
sx q[2];
rz(-0.56782349) q[2];
sx q[2];
rz(0.77726269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0878144) q[1];
sx q[1];
rz(-0.44856854) q[1];
sx q[1];
rz(-1.5320369) q[1];
x q[2];
rz(-0.73174814) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(1.5665311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058218) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18729678) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(3.0412578) q[0];
x q[1];
rz(-0.15437834) q[2];
sx q[2];
rz(-1.8533684) q[2];
sx q[2];
rz(0.64235657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3225973) q[1];
sx q[1];
rz(-1.7261191) q[1];
sx q[1];
rz(-2.0218693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13637654) q[3];
sx q[3];
rz(-2.4032421) q[3];
sx q[3];
rz(0.34463681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-1.1365183) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14295386) q[0];
sx q[0];
rz(-0.26252258) q[0];
sx q[0];
rz(-2.4130164) q[0];
rz(1.8385356) q[2];
sx q[2];
rz(-1.3053075) q[2];
sx q[2];
rz(0.43648411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(-1.8610958) q[1];
rz(-pi) q[2];
rz(0.45883026) q[3];
sx q[3];
rz(-0.60744951) q[3];
sx q[3];
rz(-1.6096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.3827682) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(1.680254) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(2.2447341) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-2.9607594) q[1];
sx q[1];
rz(1.1175964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5516978) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(-0.19703534) q[0];
rz(-pi) q[1];
rz(1.3300702) q[2];
sx q[2];
rz(-0.72611085) q[2];
sx q[2];
rz(1.7118529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0301789) q[1];
sx q[1];
rz(-1.7708659) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(0.53447978) q[3];
sx q[3];
rz(-1.8090873) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-2.1462671) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156408) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(0.036235972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72139229) q[2];
sx q[2];
rz(-1.0288236) q[2];
sx q[2];
rz(-1.1162356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36775667) q[1];
sx q[1];
rz(-2.3564597) q[1];
sx q[1];
rz(0.30898856) q[1];
rz(-pi) q[2];
rz(2.213845) q[3];
sx q[3];
rz(-1.6258996) q[3];
sx q[3];
rz(-0.4004713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76072389) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-0.6828126) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(1.9314996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34459201) q[0];
sx q[0];
rz(-1.3636149) q[0];
sx q[0];
rz(-1.9318357) q[0];
x q[1];
rz(-2.8898583) q[2];
sx q[2];
rz(-0.98099698) q[2];
sx q[2];
rz(-0.41433197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4267537) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(-2.0183802) q[1];
rz(-pi) q[2];
rz(-2.6415778) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(-2.6422215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.7738336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.66946) q[0];
sx q[0];
rz(-1.2233943) q[0];
sx q[0];
rz(3.0513289) q[0];
rz(-pi) q[1];
rz(2.6752459) q[2];
sx q[2];
rz(-2.9620167) q[2];
sx q[2];
rz(-1.8293543) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29337063) q[1];
sx q[1];
rz(-1.9105517) q[1];
sx q[1];
rz(-2.1128224) q[1];
rz(-pi) q[2];
rz(-1.6465854) q[3];
sx q[3];
rz(-1.4255376) q[3];
sx q[3];
rz(0.69243542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(-2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0859062) q[0];
sx q[0];
rz(-1.1924556) q[0];
sx q[0];
rz(-2.6299155) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-1.5716626) q[2];
sx q[2];
rz(-1.2348839) q[2];
sx q[2];
rz(0.15765794) q[2];
rz(-1.4064895) q[3];
sx q[3];
rz(-2.3621029) q[3];
sx q[3];
rz(2.7030871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];