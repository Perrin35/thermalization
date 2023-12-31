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
rz(-1.6879727) q[1];
sx q[1];
rz(-2.7984518) q[1];
sx q[1];
rz(-1.8106102) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91141191) q[0];
sx q[0];
rz(-2.4898306) q[0];
sx q[0];
rz(0.08641152) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1034307) q[2];
sx q[2];
rz(-0.51617235) q[2];
sx q[2];
rz(-1.49522) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(1.6394872) q[1];
rz(-pi) q[2];
rz(2.8475548) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(-2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.4527028) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(-2.8413049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6712924) q[0];
sx q[0];
rz(-1.4503345) q[0];
sx q[0];
rz(2.5380773) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2533721) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(-2.4878793) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2806432) q[1];
sx q[1];
rz(-1.3480535) q[1];
sx q[1];
rz(1.2692979) q[1];
x q[2];
rz(2.2829451) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(-1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77256569) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(-2.7712908) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(-0.38811362) q[0];
rz(-3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(2.8252576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5726349) q[0];
sx q[0];
rz(1.3257922) q[0];
rz(-pi) q[1];
rz(-0.6109654) q[2];
sx q[2];
rz(-2.4749711) q[2];
sx q[2];
rz(1.1952458) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95591009) q[1];
sx q[1];
rz(-0.61756402) q[1];
sx q[1];
rz(2.4159141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2315138) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(-1.6279189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8228089) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(2.4687693) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1349147) q[2];
sx q[2];
rz(-1.6393822) q[2];
sx q[2];
rz(2.2401631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1308243) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(-3.1229449) q[1];
x q[2];
rz(2.4098445) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(1.5665311) q[3];
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
rz(0.41695693) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.303724) q[0];
rz(0.45571348) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(3.0900893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52841016) q[0];
sx q[0];
rz(-0.15325704) q[0];
sx q[0];
rz(2.2806703) q[0];
x q[1];
rz(-1.8565882) q[2];
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
rz(-1.0609378) q[1];
sx q[1];
rz(-2.6662711) q[1];
sx q[1];
rz(-1.2259543) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6938985) q[3];
sx q[3];
rz(-0.8408635) q[3];
sx q[3];
rz(-2.6134932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(0.59801897) q[2];
rz(2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(0.083267033) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986388) q[0];
sx q[0];
rz(-0.26252258) q[0];
sx q[0];
rz(-2.4130164) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3699049) q[2];
sx q[2];
rz(-0.37479127) q[2];
sx q[2];
rz(0.37116606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6492918) q[1];
sx q[1];
rz(-1.2052844) q[1];
sx q[1];
rz(-0.11458061) q[1];
rz(1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(-0.99029535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.4613387) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516978) q[0];
sx q[0];
rz(-0.74001827) q[0];
sx q[0];
rz(-0.19703534) q[0];
rz(-pi) q[1];
rz(-0.20861161) q[2];
sx q[2];
rz(-0.87000712) q[2];
sx q[2];
rz(-2.029062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2988105) q[1];
sx q[1];
rz(-2.543445) q[1];
sx q[1];
rz(0.36069718) q[1];
rz(-pi) q[2];
rz(-2.696633) q[3];
sx q[3];
rz(-0.58044725) q[3];
sx q[3];
rz(-2.2821033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-0.6589748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4188822) q[0];
sx q[0];
rz(-1.5960777) q[0];
sx q[0];
rz(-0.77194571) q[0];
rz(-pi) q[1];
rz(-0.89491567) q[2];
sx q[2];
rz(-0.96940982) q[2];
sx q[2];
rz(-2.2611157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3498889) q[1];
sx q[1];
rz(-2.3096497) q[1];
sx q[1];
rz(1.8658584) q[1];
rz(-0.068816618) q[3];
sx q[3];
rz(-2.2127082) q[3];
sx q[3];
rz(-1.9300234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34459201) q[0];
sx q[0];
rz(-1.3636149) q[0];
sx q[0];
rz(1.209757) q[0];
x q[1];
rz(0.25173431) q[2];
sx q[2];
rz(-0.98099698) q[2];
sx q[2];
rz(2.7272607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4267537) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(1.1232125) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5217767) q[3];
sx q[3];
rz(-1.8780939) q[3];
sx q[3];
rz(-0.67051552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-0.79968232) q[1];
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
rz(-0.16074796) q[2];
sx q[2];
rz(-1.4904009) q[2];
sx q[2];
rz(-2.423167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4754776) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(-2.75027) q[1];
rz(-pi) q[2];
rz(0.47761376) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(-1.9660266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-0.69520673) q[2];
rz(-2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.5716626) q[2];
sx q[2];
rz(-1.2348839) q[2];
sx q[2];
rz(0.15765794) q[2];
rz(-0.1602608) q[3];
sx q[3];
rz(-2.3370623) q[3];
sx q[3];
rz(-0.2094895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
