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
rz(-1.3309825) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91141191) q[0];
sx q[0];
rz(-2.4898306) q[0];
sx q[0];
rz(-0.08641152) q[0];
x q[1];
rz(-3.1034307) q[2];
sx q[2];
rz(-2.6254203) q[2];
sx q[2];
rz(-1.6463726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5142071) q[1];
sx q[1];
rz(-1.5537019) q[1];
sx q[1];
rz(-1.8219201) q[1];
x q[2];
rz(0.64294502) q[3];
sx q[3];
rz(-1.391198) q[3];
sx q[3];
rz(-2.5760866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.4527028) q[0];
rz(-2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18314221) q[0];
sx q[0];
rz(-0.97226769) q[0];
sx q[0];
rz(-1.7167702) q[0];
rz(-0.57057256) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(0.74769339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.092204658) q[1];
sx q[1];
rz(-2.7687679) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85864752) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(-1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77256569) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(2.7712908) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(-2.753479) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(2.8252576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91352275) q[0];
sx q[0];
rz(-1.5726349) q[0];
sx q[0];
rz(-1.3257922) q[0];
rz(-2.5691368) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(-0.87871694) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2410779) q[1];
sx q[1];
rz(-1.9652275) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(-pi) q[2];
rz(-1.9100788) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(-2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(0.8316935) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261525) q[0];
sx q[0];
rz(-0.8988131) q[0];
sx q[0];
rz(1.5120904) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1349147) q[2];
sx q[2];
rz(-1.5022105) q[2];
sx q[2];
rz(-0.90142957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0878144) q[1];
sx q[1];
rz(-2.6930241) q[1];
sx q[1];
rz(-1.6095557) q[1];
rz(-pi) q[2];
rz(2.6850011) q[3];
sx q[3];
rz(-0.78568229) q[3];
sx q[3];
rz(2.8031138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7164798) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.8378687) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542959) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(-0.10033484) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15437834) q[2];
sx q[2];
rz(-1.2882243) q[2];
sx q[2];
rz(2.4992361) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67700779) q[1];
sx q[1];
rz(-1.1255463) q[1];
sx q[1];
rz(-0.17226179) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0052161) q[3];
sx q[3];
rz(-2.4032421) q[3];
sx q[3];
rz(-2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(2.5435737) q[2];
rz(-0.55001843) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(-2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-3.0583256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71654728) q[0];
sx q[0];
rz(-1.7444567) q[0];
sx q[0];
rz(-0.19787775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8667738) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(-1.9354265) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0373842) q[1];
sx q[1];
rz(-1.4638149) q[1];
sx q[1];
rz(-1.2030829) q[1];
rz(-pi) q[2];
rz(0.45883026) q[3];
sx q[3];
rz(-2.5341431) q[3];
sx q[3];
rz(1.6096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(1.4613387) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7145342) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-1.1175964) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8539124) q[0];
sx q[0];
rz(-0.84830647) q[0];
sx q[0];
rz(-1.393909) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8115225) q[2];
sx q[2];
rz(-0.72611085) q[2];
sx q[2];
rz(-1.7118529) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1114137) q[1];
sx q[1];
rz(-1.7708659) q[1];
sx q[1];
rz(-2.5740037) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53447978) q[3];
sx q[3];
rz(-1.3325053) q[3];
sx q[3];
rz(-1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156408) q[0];
sx q[0];
rz(-2.3693187) q[0];
sx q[0];
rz(3.1053567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4023513) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(-3.0657363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3498889) q[1];
sx q[1];
rz(-2.3096497) q[1];
sx q[1];
rz(-1.2757343) q[1];
rz(-pi) q[2];
rz(-1.4790672) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(1.0969539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(2.4587801) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.9314996) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34459201) q[0];
sx q[0];
rz(-1.3636149) q[0];
sx q[0];
rz(1.9318357) q[0];
rz(-1.9270883) q[2];
sx q[2];
rz(-0.63535832) q[2];
sx q[2];
rz(-0.84745849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2883816) q[1];
sx q[1];
rz(-0.92331159) q[1];
sx q[1];
rz(2.7700469) q[1];
rz(-pi) q[2];
rz(-2.6415778) q[3];
sx q[3];
rz(-0.68272907) q[3];
sx q[3];
rz(-0.4993712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.575763) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(-1.7738336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.43191) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73197094) q[0];
sx q[0];
rz(-0.35847607) q[0];
sx q[0];
rz(1.8147857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6522371) q[2];
sx q[2];
rz(-1.7310206) q[2];
sx q[2];
rz(-0.83934957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.666115) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(0.39132262) q[1];
rz(-pi) q[2];
rz(-2.9959216) q[3];
sx q[3];
rz(-1.6457857) q[3];
sx q[3];
rz(0.86736995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0556864) q[0];
sx q[0];
rz(-1.1924556) q[0];
sx q[0];
rz(-2.6299155) q[0];
rz(1.2790537) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(-2.8056801) q[2];
sx q[2];
rz(-1.5699785) q[2];
sx q[2];
rz(1.7287398) q[2];
rz(1.7351032) q[3];
sx q[3];
rz(-2.3621029) q[3];
sx q[3];
rz(2.7030871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
