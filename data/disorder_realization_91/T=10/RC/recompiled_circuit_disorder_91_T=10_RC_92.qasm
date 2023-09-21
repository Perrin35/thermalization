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
rz(-1.6879727) q[1];
sx q[1];
rz(3.4847335) q[1];
sx q[1];
rz(7.6141678) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59062293) q[0];
sx q[0];
rz(-1.6231713) q[0];
sx q[0];
rz(0.64996029) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1034307) q[2];
sx q[2];
rz(-2.6254203) q[2];
sx q[2];
rz(1.6463726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0806182) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(0.017647839) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4986476) q[3];
sx q[3];
rz(-1.391198) q[3];
sx q[3];
rz(-0.56550607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.3295056) q[2];
rz(1.7261516) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9584504) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(1.4248225) q[0];
rz(-1.8882206) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(-2.4878793) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.092204658) q[1];
sx q[1];
rz(-2.7687679) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(-1.7598011) q[3];
sx q[3];
rz(-1.4100037) q[3];
sx q[3];
rz(0.27634987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(0.31633502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91352275) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(1.3257922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5306273) q[2];
sx q[2];
rz(-2.4749711) q[2];
sx q[2];
rz(-1.9463469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12831941) q[1];
sx q[1];
rz(-2.0187906) q[1];
sx q[1];
rz(-2.0112579) q[1];
rz(1.2315138) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(-2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(2.9807828) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(-2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31878372) q[0];
sx q[0];
rz(-1.524864) q[0];
sx q[0];
rz(0.6728234) q[0];
rz(-2.1349147) q[2];
sx q[2];
rz(-1.6393822) q[2];
sx q[2];
rz(-0.90142957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1308243) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(3.1229449) q[1];
rz(-pi) q[2];
rz(0.73174814) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(1.5750615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(1.224219) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(3.0900893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542959) q[0];
sx q[0];
rz(-1.686839) q[0];
sx q[0];
rz(-3.0412578) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8565882) q[2];
sx q[2];
rz(-1.4225866) q[2];
sx q[2];
rz(2.1697901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4645849) q[1];
sx q[1];
rz(-1.1255463) q[1];
sx q[1];
rz(2.9693309) q[1];
x q[2];
rz(0.73370917) q[3];
sx q[3];
rz(-1.6624311) q[3];
sx q[3];
rz(-1.1250145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6872528) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705868) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(2.9396074) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71654728) q[0];
sx q[0];
rz(-1.7444567) q[0];
sx q[0];
rz(2.9437149) q[0];
x q[1];
rz(-2.3699049) q[2];
sx q[2];
rz(-0.37479127) q[2];
sx q[2];
rz(-2.7704266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3377933) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(-1.2804968) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(-2.1512973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.7588245) q[2];
rz(2.8921195) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(-1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5898949) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(0.19703534) q[0];
x q[1];
rz(-0.859185) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(2.8189916) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7270131) q[1];
sx q[1];
rz(-2.1257183) q[1];
sx q[1];
rz(1.8068061) q[1];
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
rz(-pi) q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(0.21952595) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156408) q[0];
sx q[0];
rz(-2.3693187) q[0];
sx q[0];
rz(-0.036235972) q[0];
rz(2.246677) q[2];
sx q[2];
rz(-0.96940982) q[2];
sx q[2];
rz(0.88047699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36775667) q[1];
sx q[1];
rz(-0.78513297) q[1];
sx q[1];
rz(-2.8326041) q[1];
rz(-pi) q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-0.64507161) q[3];
sx q[3];
rz(1.0969539) q[3];
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
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3808688) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(-2.4587801) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.9314996) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-2.1754873) q[2];
sx q[2];
rz(-1.3622869) q[2];
sx q[2];
rz(-2.1272121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.71483892) q[1];
sx q[1];
rz(-0.73298448) q[1];
sx q[1];
rz(-2.0183802) q[1];
rz(-pi) q[2];
rz(-1.9425415) q[3];
sx q[3];
rz(-0.98402714) q[3];
sx q[3];
rz(-1.1128807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(3.0513289) q[0];
rz(-pi) q[1];
rz(-0.46634679) q[2];
sx q[2];
rz(-2.9620167) q[2];
sx q[2];
rz(-1.8293543) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3696246) q[1];
sx q[1];
rz(-0.63056417) q[1];
sx q[1];
rz(2.1715013) q[1];
rz(-pi) q[2];
rz(-1.6465854) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(-0.69243542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(-2.4485596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-0.002481133) q[2];
sx q[2];
rz(-2.8056792) q[2];
sx q[2];
rz(0.16028595) q[2];
rz(-2.3435076) q[3];
sx q[3];
rz(-1.4555664) q[3];
sx q[3];
rz(1.2496787) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
