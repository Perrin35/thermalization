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
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5509697) q[0];
sx q[0];
rz(-1.6231713) q[0];
sx q[0];
rz(0.64996029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.038161909) q[2];
sx q[2];
rz(-0.51617235) q[2];
sx q[2];
rz(-1.6463726) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0099437873) q[1];
sx q[1];
rz(-0.25169262) q[1];
sx q[1];
rz(1.5021055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29403789) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(-2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3657637) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(2.9392488) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99130327) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(1.4527028) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(2.8413049) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9584504) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(-1.7167702) q[0];
x q[1];
rz(-2.5710201) q[2];
sx q[2];
rz(-1.8405481) q[2];
sx q[2];
rz(0.74769339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2806432) q[1];
sx q[1];
rz(-1.3480535) q[1];
sx q[1];
rz(-1.8722948) q[1];
rz(-pi) q[2];
rz(-0.85864752) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(-1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(-2.753479) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(2.8252576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64992031) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(1.578376) q[0];
rz(0.57245589) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(-0.87871694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9005147) q[1];
sx q[1];
rz(-1.9652275) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(-pi) q[2];
rz(-2.2736336) q[3];
sx q[3];
rz(-2.7083525) q[3];
sx q[3];
rz(-0.058335282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1091653) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9008824) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(0.8316935) q[0];
rz(1.3446993) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(-1.6279189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095703) q[0];
sx q[0];
rz(-0.67414588) q[0];
sx q[0];
rz(0.073622965) q[0];
x q[1];
rz(1.0066779) q[2];
sx q[2];
rz(-1.5022105) q[2];
sx q[2];
rz(-2.2401631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6935389) q[1];
sx q[1];
rz(-1.5539907) q[1];
sx q[1];
rz(1.1225213) q[1];
rz(-pi) q[2];
rz(-2.6850011) q[3];
sx q[3];
rz(-2.3559104) q[3];
sx q[3];
rz(2.8031138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(1.224219) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(-1.8378687) q[0];
rz(0.45571348) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-3.0900893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52841016) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(-0.86092237) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0578458) q[2];
sx q[2];
rz(-0.32099989) q[2];
sx q[2];
rz(-0.1333065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3225973) q[1];
sx q[1];
rz(-1.4154735) q[1];
sx q[1];
rz(2.0218693) q[1];
rz(0.73370917) q[3];
sx q[3];
rz(-1.4791616) q[3];
sx q[3];
rz(-2.0165781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(0.59801897) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14295386) q[0];
sx q[0];
rz(-0.26252258) q[0];
sx q[0];
rz(-2.4130164) q[0];
x q[1];
rz(1.8385356) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(-0.43648411) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-0.38227907) q[1];
sx q[1];
rz(1.8610958) q[1];
rz(0.45883026) q[3];
sx q[3];
rz(-0.60744951) q[3];
sx q[3];
rz(-1.6096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(-1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.9607594) q[1];
sx q[1];
rz(-1.1175964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5898949) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(-2.9445573) q[0];
rz(-pi) q[1];
rz(0.20861161) q[2];
sx q[2];
rz(-2.2715855) q[2];
sx q[2];
rz(-2.029062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1114137) q[1];
sx q[1];
rz(-1.7708659) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(0.53447978) q[3];
sx q[3];
rz(-1.3325053) q[3];
sx q[3];
rz(2.0508545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(0.21952595) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(0.21729939) q[0];
rz(0.030933881) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(2.4826179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17651672) q[0];
sx q[0];
rz(-2.342431) q[0];
sx q[0];
rz(1.6060711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73924139) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(0.07585635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3498889) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(1.2757343) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4790672) q[3];
sx q[3];
rz(-0.64507161) q[3];
sx q[3];
rz(-2.0446387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(-2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.9314996) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7970006) q[0];
sx q[0];
rz(-1.3636149) q[0];
sx q[0];
rz(-1.9318357) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8898583) q[2];
sx q[2];
rz(-0.98099698) q[2];
sx q[2];
rz(-2.7272607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4267537) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(-2.0183802) q[1];
rz(-1.9425415) q[3];
sx q[3];
rz(-2.1575655) q[3];
sx q[3];
rz(1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.575763) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(-3.0513289) q[0];
rz(-pi) q[1];
rz(-2.6752459) q[2];
sx q[2];
rz(-2.9620167) q[2];
sx q[2];
rz(-1.3122383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-1.0628504) q[1];
sx q[1];
rz(-2.75027) q[1];
x q[2];
rz(-2.9959216) q[3];
sx q[3];
rz(-1.6457857) q[3];
sx q[3];
rz(0.86736995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-3.1391115) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(-1.7351032) q[3];
sx q[3];
rz(-0.77948979) q[3];
sx q[3];
rz(-0.43850552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];