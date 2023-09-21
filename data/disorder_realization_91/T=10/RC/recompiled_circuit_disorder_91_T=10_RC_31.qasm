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
rz(2.9705272) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(-1.3309825) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59062293) q[0];
sx q[0];
rz(-1.6231713) q[0];
sx q[0];
rz(-0.64996029) q[0];
rz(0.51585977) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(-0.042382391) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.0099437873) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(1.5021055) q[1];
rz(1.347723) q[3];
sx q[3];
rz(-2.2017456) q[3];
sx q[3];
rz(1.1383206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502894) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0694885) q[0];
sx q[0];
rz(-2.527643) q[0];
sx q[0];
rz(2.9314562) q[0];
x q[1];
rz(-1.8882206) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(0.65371338) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.049388) q[1];
sx q[1];
rz(-2.7687679) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3817915) q[3];
sx q[3];
rz(-1.731589) q[3];
sx q[3];
rz(2.8652428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(-0.84233061) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(0.38811362) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-2.8252576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4916723) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(1.5632167) q[0];
x q[1];
rz(0.57245589) q[2];
sx q[2];
rz(-1.9334031) q[2];
sx q[2];
rz(-2.2628757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1856826) q[1];
sx q[1];
rz(-0.61756402) q[1];
sx q[1];
rz(0.72567852) q[1];
x q[2];
rz(0.86795904) q[3];
sx q[3];
rz(-0.43324019) q[3];
sx q[3];
rz(-3.0832574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(1.3446993) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31878372) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(0.6728234) q[0];
rz(-0.081110031) q[2];
sx q[2];
rz(-1.0081648) q[2];
sx q[2];
rz(-2.5155591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6935389) q[1];
sx q[1];
rz(-1.5539907) q[1];
sx q[1];
rz(2.0190713) q[1];
x q[2];
rz(1.9862595) q[3];
sx q[3];
rz(-0.88298015) q[3];
sx q[3];
rz(0.26879877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(0.41695693) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6131825) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(2.2806703) q[0];
x q[1];
rz(1.2850045) q[2];
sx q[2];
rz(-1.7190061) q[2];
sx q[2];
rz(0.9718026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4645849) q[1];
sx q[1];
rz(-2.0160463) q[1];
sx q[1];
rz(-2.9693309) q[1];
x q[2];
rz(2.4078835) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(1.1250145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6872528) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(2.5435737) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(2.0050744) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(3.0583256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250454) q[0];
sx q[0];
rz(-1.3971359) q[0];
sx q[0];
rz(2.9437149) q[0];
x q[1];
rz(-0.27481885) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(1.9354265) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80379936) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(1.2804968) q[1];
x q[2];
rz(-1.2721328) q[3];
sx q[3];
rz(-1.0335869) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.3827682) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
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
rz(1.5898949) q[0];
sx q[0];
rz(-0.74001827) q[0];
sx q[0];
rz(-0.19703534) q[0];
rz(2.932981) q[2];
sx q[2];
rz(-2.2715855) q[2];
sx q[2];
rz(2.029062) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0301789) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(0.56758893) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44495961) q[3];
sx q[3];
rz(-2.5611454) q[3];
sx q[3];
rz(0.85948932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(0.21729939) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(-0.6589748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12595183) q[0];
sx q[0];
rz(-2.3693187) q[0];
sx q[0];
rz(-0.036235972) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72139229) q[2];
sx q[2];
rz(-1.0288236) q[2];
sx q[2];
rz(-2.0253571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1606071) q[1];
sx q[1];
rz(-1.7874582) q[1];
sx q[1];
rz(-0.76088455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4790672) q[3];
sx q[3];
rz(-2.496521) q[3];
sx q[3];
rz(-2.0446387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(-0.6828126) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3037198) q[0];
sx q[0];
rz(-1.923773) q[0];
sx q[0];
rz(2.9205802) q[0];
rz(-pi) q[1];
rz(0.96610539) q[2];
sx q[2];
rz(-1.3622869) q[2];
sx q[2];
rz(-2.1272121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.853211) q[1];
sx q[1];
rz(-0.92331159) q[1];
sx q[1];
rz(0.37154571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9425415) q[3];
sx q[3];
rz(-0.98402714) q[3];
sx q[3];
rz(-1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.575763) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.43191) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
rz(pi/2) q[2];
sx q[2];
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
x q[1];
rz(-0.46634679) q[2];
sx q[2];
rz(-0.17957598) q[2];
sx q[2];
rz(-1.3122383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4754776) q[1];
sx q[1];
rz(-1.0628504) q[1];
sx q[1];
rz(-0.39132262) q[1];
rz(-0.47761376) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(1.9660266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2852823) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0859062) q[0];
sx q[0];
rz(-1.1924556) q[0];
sx q[0];
rz(-2.6299155) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(1.5699301) q[2];
sx q[2];
rz(-1.2348839) q[2];
sx q[2];
rz(0.15765794) q[2];
rz(-2.9813319) q[3];
sx q[3];
rz(-0.80453034) q[3];
sx q[3];
rz(2.9321032) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];