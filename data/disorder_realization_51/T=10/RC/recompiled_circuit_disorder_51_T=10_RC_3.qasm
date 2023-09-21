OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7979413) q[0];
sx q[0];
rz(-2.6290253) q[0];
sx q[0];
rz(2.0128065) q[0];
rz(-pi) q[1];
rz(1.2744781) q[2];
sx q[2];
rz(-0.9689435) q[2];
sx q[2];
rz(-1.559343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3095113) q[1];
sx q[1];
rz(-2.32825) q[1];
sx q[1];
rz(1.3002212) q[1];
x q[2];
rz(1.231108) q[3];
sx q[3];
rz(-1.7237444) q[3];
sx q[3];
rz(-1.4121526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.367761) q[2];
rz(-1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(0.071454123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.6921035) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.5140623) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1986952) q[0];
sx q[0];
rz(-1.9778344) q[0];
sx q[0];
rz(1.0242978) q[0];
rz(-0.83346955) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(1.8881063) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4501614) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(1.3732234) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0239003) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(-2.5642455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65511584) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(-0.97186175) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(3.0774975) q[0];
rz(2.8308716) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.4583189) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4401535) q[0];
sx q[0];
rz(-0.15471409) q[0];
sx q[0];
rz(1.2604777) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1737353) q[2];
sx q[2];
rz(-2.1328916) q[2];
sx q[2];
rz(-2.2130655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2787784) q[1];
sx q[1];
rz(-1.1618975) q[1];
sx q[1];
rz(1.5120974) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5879257) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(-2.0103612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-2.4625051) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963592) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.6569998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35614466) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(-1.5201475) q[0];
rz(-0.19214432) q[2];
sx q[2];
rz(-2.0803183) q[2];
sx q[2];
rz(1.1332878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2034519) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(-0.52312619) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37064151) q[3];
sx q[3];
rz(-2.3834043) q[3];
sx q[3];
rz(-1.7565808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0174039) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(-2.1253288) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(1.0884292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50773412) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.9790861) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50786103) q[0];
sx q[0];
rz(-1.5989688) q[0];
sx q[0];
rz(-0.072576056) q[0];
x q[1];
rz(-1.7816254) q[2];
sx q[2];
rz(-1.1846917) q[2];
sx q[2];
rz(0.93036508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9160737) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(-1.4654935) q[1];
rz(1.3624304) q[3];
sx q[3];
rz(-1.2456129) q[3];
sx q[3];
rz(-1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.48012039) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(0.76812569) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(-0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.1434198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(-2.9907945) q[0];
x q[1];
rz(-1.9009695) q[2];
sx q[2];
rz(-2.1483443) q[2];
sx q[2];
rz(1.6799048) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.032635078) q[1];
sx q[1];
rz(-2.7320478) q[1];
sx q[1];
rz(2.0202548) q[1];
rz(-pi) q[2];
rz(-1.0877803) q[3];
sx q[3];
rz(-0.59637585) q[3];
sx q[3];
rz(0.071810178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4986971) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221508) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(1.1332606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5232552) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(-0.32251127) q[0];
rz(-pi) q[1];
rz(1.557017) q[2];
sx q[2];
rz(-2.2829977) q[2];
sx q[2];
rz(0.53879246) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95722317) q[1];
sx q[1];
rz(-1.2979227) q[1];
sx q[1];
rz(-0.89783122) q[1];
rz(-pi) q[2];
rz(2.6071068) q[3];
sx q[3];
rz(-0.30189461) q[3];
sx q[3];
rz(1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(2.6531632) q[2];
rz(-1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(1.2088998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8457444) q[0];
sx q[0];
rz(-0.83202067) q[0];
sx q[0];
rz(-2.8315298) q[0];
x q[1];
rz(0.9135984) q[2];
sx q[2];
rz(-0.29619869) q[2];
sx q[2];
rz(1.2072472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68122411) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(3.1254915) q[1];
rz(-0.49105673) q[3];
sx q[3];
rz(-2.4591082) q[3];
sx q[3];
rz(-1.7562989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(0.87336826) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1677925) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.5100381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841543) q[0];
sx q[0];
rz(-2.7630002) q[0];
sx q[0];
rz(0.5666825) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(2.0575112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1875302) q[1];
sx q[1];
rz(-2.4937951) q[1];
sx q[1];
rz(-2.8026583) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4771677) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.2053306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(0.81595016) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9394768) q[0];
sx q[0];
rz(-2.3411223) q[0];
sx q[0];
rz(-0.70864422) q[0];
x q[1];
rz(1.9071104) q[2];
sx q[2];
rz(-2.2173777) q[2];
sx q[2];
rz(0.99415776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60724466) q[1];
sx q[1];
rz(-1.4048647) q[1];
sx q[1];
rz(-1.0870618) q[1];
x q[2];
rz(-0.55024054) q[3];
sx q[3];
rz(-2.5329258) q[3];
sx q[3];
rz(-0.88702162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6170549) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.9032003) q[1];
sx q[1];
rz(-1.4307784) q[1];
sx q[1];
rz(1.4684114) q[1];
rz(0.94454371) q[2];
sx q[2];
rz(-1.4251475) q[2];
sx q[2];
rz(0.057694358) q[2];
rz(-1.5068549) q[3];
sx q[3];
rz(-2.1422374) q[3];
sx q[3];
rz(-0.16850785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
