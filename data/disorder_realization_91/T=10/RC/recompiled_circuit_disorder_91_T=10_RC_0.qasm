OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(5.7167238) q[0];
sx q[0];
rz(9.5958435) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(-2.7984518) q[1];
sx q[1];
rz(-1.8106102) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5509697) q[0];
sx q[0];
rz(-1.5184214) q[0];
sx q[0];
rz(0.64996029) q[0];
rz(0.51585977) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(3.0992103) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-0.25169262) q[1];
sx q[1];
rz(-1.6394872) q[1];
rz(1.7938697) q[3];
sx q[3];
rz(-2.2017456) q[3];
sx q[3];
rz(2.0032721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77582899) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-2.9392488) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99130327) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4703003) q[0];
sx q[0];
rz(-1.6912582) q[0];
sx q[0];
rz(0.60351535) q[0];
rz(-pi) q[1];
rz(1.2533721) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(0.65371338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.049388) q[1];
sx q[1];
rz(-0.37282473) q[1];
sx q[1];
rz(-2.2224109) q[1];
x q[2];
rz(-2.977936) q[3];
sx q[3];
rz(-1.7573342) q[3];
sx q[3];
rz(-1.877762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(-2.1035813) q[2];
rz(-0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(0.38811362) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-2.8252576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65773327) q[0];
sx q[0];
rz(-1.3257926) q[0];
sx q[0];
rz(3.1396975) q[0];
rz(-pi) q[1];
rz(0.6109654) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(-1.9463469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9005147) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(-pi) q[2];
rz(-0.86795904) q[3];
sx q[3];
rz(-2.7083525) q[3];
sx q[3];
rz(-3.0832574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(0.8316935) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.6279189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8320223) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(-3.0679697) q[0];
rz(3.0604826) q[2];
sx q[2];
rz(-2.1334279) q[2];
sx q[2];
rz(2.5155591) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1308243) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(3.1229449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73174814) q[3];
sx q[3];
rz(-1.8879315) q[3];
sx q[3];
rz(1.5750615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18729678) q[0];
sx q[0];
rz(-1.686839) q[0];
sx q[0];
rz(0.10033484) q[0];
x q[1];
rz(-2.0578458) q[2];
sx q[2];
rz(-2.8205928) q[2];
sx q[2];
rz(-0.1333065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0609378) q[1];
sx q[1];
rz(-0.47532156) q[1];
sx q[1];
rz(1.2259543) q[1];
rz(-pi) q[2];
rz(3.0052161) q[3];
sx q[3];
rz(-0.73835056) q[3];
sx q[3];
rz(0.34463681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(0.59801897) q[2];
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
rz(-pi) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(-2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(0.083267033) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71654728) q[0];
sx q[0];
rz(-1.3971359) q[0];
sx q[0];
rz(2.9437149) q[0];
x q[1];
rz(2.8667738) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(-1.2061662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3377933) q[1];
sx q[1];
rz(-0.38227907) q[1];
sx q[1];
rz(1.2804968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5842651) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(2.717358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(2.8921195) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(2.2447341) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5898949) q[0];
sx q[0];
rz(-0.74001827) q[0];
sx q[0];
rz(-2.9445573) q[0];
rz(1.8115225) q[2];
sx q[2];
rz(-2.4154818) q[2];
sx q[2];
rz(1.7118529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.41457957) q[1];
sx q[1];
rz(-1.0158744) q[1];
sx q[1];
rz(-1.8068061) q[1];
rz(1.8459122) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(-0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0573132) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(-2.4826179) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17651672) q[0];
sx q[0];
rz(-2.342431) q[0];
sx q[0];
rz(-1.6060711) q[0];
rz(-pi) q[1];
rz(2.4023513) q[2];
sx q[2];
rz(-2.2695702) q[2];
sx q[2];
rz(-3.0657363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98098552) q[1];
sx q[1];
rz(-1.3541344) q[1];
sx q[1];
rz(-2.3807081) q[1];
rz(1.4790672) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(-1.0969539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(1.210093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4142128) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(2.1754873) q[2];
sx q[2];
rz(-1.3622869) q[2];
sx q[2];
rz(2.1272121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71483892) q[1];
sx q[1];
rz(-0.73298448) q[1];
sx q[1];
rz(-2.0183802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1990511) q[3];
sx q[3];
rz(-2.1575655) q[3];
sx q[3];
rz(1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.575763) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.7738336) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(1.43191) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-0.79968232) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0737338) q[0];
sx q[0];
rz(-1.4859383) q[0];
sx q[0];
rz(-1.2220864) q[0];
x q[1];
rz(-2.9808447) q[2];
sx q[2];
rz(-1.6511917) q[2];
sx q[2];
rz(-2.423167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3696246) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(-0.97009138) q[1];
rz(0.47761376) q[3];
sx q[3];
rz(-2.9778746) q[3];
sx q[3];
rz(1.9660266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.2790537) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(-0.002481133) q[2];
sx q[2];
rz(-2.8056792) q[2];
sx q[2];
rz(0.16028595) q[2];
rz(2.9813319) q[3];
sx q[3];
rz(-2.3370623) q[3];
sx q[3];
rz(-0.2094895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
